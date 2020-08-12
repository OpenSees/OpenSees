/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DistributedDiagonalSOE.cpp,v $

// Written: fmk 
// Created: 05/05
//
// Description: This file contains the implementation for DistributedDiagonalSOE


#include <DistributedDiagonalSOE.h>
#include <DistributedDiagonalSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <AnalysisModel.h>
#include <DOF_Group.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>


DistributedDiagonalSOE::DistributedDiagonalSOE(DistributedDiagonalSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_DistributedDiagonalSOE),
 size(0), A(0), B(0), X(0), vectX(0), vectB(0), isAfactored(false),
 processID(0), numProcesses(0),
 numChannels(0), theChannels(0), localCol(0), 
 myDOFs(0,32), myDOFsShared(0,16), numShared(0), dataShared(0), vectShared(0), theModel(0)
{
    the_Solver.setLinearSOE(*this);
}

DistributedDiagonalSOE::DistributedDiagonalSOE()
:LinearSOE(LinSOE_TAGS_DistributedDiagonalSOE),
 size(0), A(0), B(0), X(0), vectX(0), vectB(0), isAfactored(false),
 processID(0), numProcesses(0),
 numChannels(0), theChannels(0), localCol(0), 
 myDOFs(0,32), myDOFsShared(0,16), numShared(0), dataShared(0), vectShared(0), theModel(0)
{

}


DistributedDiagonalSOE::~DistributedDiagonalSOE()
{
  if (A != 0) delete [] A;
  if (B != 0) delete [] B;
  if (X != 0) delete [] X;
  if (dataShared != 0) delete [] dataShared;
  if (vectX != 0) delete vectX;    
  if (vectB != 0) delete vectB;    
  if (vectShared != 0) delete vectShared;    

  if (theChannels != 0)
    delete [] theChannels;
}


int 
DistributedDiagonalSOE::getNumEqn(void) const
{
  return size;
}


int 
DistributedDiagonalSOE::setSize(Graph &theGraph)
{
  int result = 0;
  size = theGraph.getNumVertex();

  //
  // first we build an ID containing all local DOFs
  //

  myDOFs.resize(size);

  int count = 0;
  Vertex *theVertex;
  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    int vertexTag = theVertex->getTag();
    myDOFs(count) = vertexTag;
    count++;
  }

  static ID otherSize(1);
  ID otherDOFS(0, size/10);

  if (processID != 0) {

    //
    // each process send it's local and receives all other processes IDs (remote IDs).
    // from the remote IDs the local process figures out which DOF's it shares with others.
    // it then sends the shared dofs' out.
    //

    Channel *theChannel = theChannels[0];
    theChannel->recvID(0, 0, otherSize);
    numProcesses = otherSize(0);

    // send local
    otherSize(0) = size;
    theChannel->sendID(0, 0, otherSize);
    if (size != 0)
      theChannel->sendID(0, 0, myDOFs);

    for (int i=0; i<numProcesses-1; i++) {
      // receive remote & check for shared DOFs
      theChannel->recvID(0, 0, otherSize);	
      int numOther = otherSize(0);

      if (numOther != 0) {
	otherDOFS.resize(numOther);
	theChannel->recvID(0, 0, otherDOFS);
	
	for (int j=0; j<numOther; j++) {
	  int otherTag = otherDOFS(j);
	  if (myDOFs.getLocation(otherTag) != -1 && myDOFsShared.getLocation(otherTag) == -1) 
	    myDOFsShared[numShared++] = otherTag;
	}
      }
    }

    // send my shared DOFs
    otherSize(0) = myDOFsShared.Size();
    theChannel->sendID(0, 0, otherSize);
    if (size != 0)
      theChannel->sendID(0, 0, myDOFsShared);

    // recv all shared DOFs
    theChannel->recvID(0, 0, otherSize);	
    numShared = otherSize(0);

    if (numShared != 0) {
      myDOFsShared.resize(numShared);
      theChannel->recvID(0, 0, myDOFsShared);
    } 

  } else { // in absence of a broadcast we need the following

    //
    // each process send it's local and receives all other processes IDs (remote IDs).
    // from the remote IDs the local process figures out which DOF's it shares with others.
    // it then sends the shared dofs' out. lastly receives from P0 list of all shared.
    // NOTE:  if all sent to P0 their ID's; P0 could figure out shared; but scalability issues
    //

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      otherSize(0) = numChannels+1;
      theChannel->sendID(0, 0, otherSize);

      otherSize(0) = size;
      theChannel->sendID(0, 0, otherSize);
      if (size != 0)
	theChannel->sendID(0, 0, myDOFs);
    }

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];

      // receive remote & check for shared DOFs

      theChannel->recvID(0, 0, otherSize);	
      int numOther = otherSize(0);
      if (numOther != 0) {
	otherDOFS.resize(numOther);
	theChannel->recvID(0, 0, otherDOFS);
      }

      // need to send to all others
      for (int k=0; k<numChannels; k++) {
	Channel *theChannel = theChannels[k];
	if (k != j) {
	  theChannel->sendID(0, 0, otherSize);	
	  if (numOther != 0) {
	    theChannel->sendID(0, 0, otherDOFS);	    
	  }
	}
      }
      
      // need to merge with mine 
      for (int l=0; l<numOther; l++) {
	int otherTag = otherDOFS(l);
	if (myDOFs.getLocation(otherTag) != -1 && myDOFsShared.getLocation(otherTag) == -1) { 
	  myDOFsShared[numShared++] = otherTag;
	}
      }
    }

    // now recv each of the shared DOFs & merge with mine to form the primary list
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->recvID(0, 0, otherSize);	
      int numOther = otherSize(0);
      if (numOther != 0) {
	otherDOFS.resize(numOther);
	theChannel->recvID(0, 0, otherDOFS);
      } 
      // need to merge with mine 
      for (int k=0; k<numOther; k++) {
	int otherTag = otherDOFS(k);
	if (myDOFsShared.getLocation(otherTag) == -1) 
	  myDOFsShared(numShared++) = otherTag;
      }
    }

    // now send the shared ID
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      otherSize(0) = numShared;
      theChannel->sendID(0, 0, otherSize);	
      theChannel->sendID(0, 0, myDOFsShared);
    }
  }


  if (A != 0) delete [] A; A = 0;
  if (B != 0) delete [] B; B = 0;
  if (X != 0) delete [] X; X = 0;
  if (dataShared != 0) delete [] dataShared; dataShared = 0;
  if (vectX != 0) delete vectX; vectX = 0;
  if (vectB != 0) delete vectB; vectB = 0;

  A = new double[size];
  B = new double[size];
  X = new double[size];	
  if (X != 0) 
    vectX = new Vector(X,size);
  if (B != 0) 
    vectB = new Vector(B,size);
  
  dataShared = new double[2*numShared];	 // 2 times for A & B
  if (dataShared != 0)
    vectShared = new Vector(dataShared, 2*numShared);


  if (A == 0 || B == 0 || X == 0 || vectX == 0 || vectB == 0 || dataShared == 0 || vectShared == 0) {
    opserr << "ERROR DistributedDiagonalSOE::setSize() - ";
    opserr << " ran out of memory for size: " << size << endln;
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (dataShared != 0) delete [] dataShared;
    size = 0;
    return -1;
  }

  // zero the vectors
  for (int l=0; l<size; l++) {
    A[l] = 0;
    B[l] = 0;
    X[l] = 0;
  }

  //
  // now let's redo the mapping of the dof's .. locally numbered 0 through size
  //

  
  if (theModel == 0) {
    opserr << "WARNING DistributedDiagonalSOE::setSize - no AnalysisModel\n";
  } 
  else {
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    
    while ((dofPtr = theDOFs()) != 0) {
	const ID &theID = dofPtr->getID();
	for (int i=0; i<theID.Size(); i++) {
	  int dof = theID(i);
	  if (dof >= 0) {
	    int newDOF = myDOFs.getLocation(dof);
	    dofPtr->setID(i, newDOF);
	  }
	}   
    }
    
    // iterate through the FE_Element getting them to set their IDs
    FE_EleIter &theEle = theModel->getFEs();
    FE_Element *elePtr;
    while ((elePtr = theEle()) != 0)
	elePtr->setID(); 
  }  

  // invoke setSize() on the Solver
  DistributedDiagonalSolver *the_Solver = (DistributedDiagonalSolver *)this->getSolver();
  the_Solver->setLinearSOE(*this);
  int solverOK = the_Solver->setSize();
  
  if (solverOK < 0) {
    opserr << "WARNING DistributedDiagonalSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    

  return result;
}

int 
DistributedDiagonalSOE::addA(const Matrix &m, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
#ifdef _G3DEBUG
  // check that m and id are of similar size
  int idSize = id.Size();    
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "FullGenLinSOE::addA()	- Matrix and ID not of similar sizes\n";
    return -1;
  }
#endif

  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0) {
	A[pos] += m(i,i);
      }
    }
  } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0) {
	A[pos] -= m(i,i);
      }
    }
  } else {
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0) {
	A[pos] += m(i,i) * fact;
      }
    }
  }	

  return 0;
}
 
    
int 
DistributedDiagonalSOE::addB(const Vector &v, const ID &id, double fact)
{
    
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
#ifdef _G3DEBUG
  // check that m and id are of similar size
  int idSize = id.Size();        
  if (idSize != v.Size() ) {
    opserr << "DistributedDiagonalSOE::addB() -";
    opserr << " Vector and ID not of similar sizes\n";
    return -1;
  }    
#endif
  
  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	B[pos] += v(i);
    }
  } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	B[pos] -= v(i);
    }
  } else {
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	B[pos] += v(i) * fact;
    }
  }	
  return 0;
}


int
DistributedDiagonalSOE::setB(const Vector &v, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
  if (v.Size() != size) {
    opserr << "WARNING DistributedDiagonalSOE::setB() -";
    opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
    return -1;
  }
  
  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<size; i++) {
      B[i] = v(i);
    }
  } else if (fact == -1.0) {
    for (int i=0; i<size; i++) {
      B[i] = -v(i);
    }
  } else {
    for (int i=0; i<size; i++) {
      B[i] = v(i) * fact;
    }
  }	
  return 0;
}

void 
DistributedDiagonalSOE::zeroA(void)
{
  double *Aptr = A;
  for (int i=0; i<size; i++)
    *Aptr++ = 0;
  
  isAfactored = false;
}

void 
DistributedDiagonalSOE::zeroB(void)
{
  double *Bptr = B;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}


void 
DistributedDiagonalSOE::setX(int loc, double value)
{
  if (loc < size && loc >=0)
    X[loc] = value;
}

void 
DistributedDiagonalSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
DistributedDiagonalSOE::getX(void)
{
  if (vectX == 0) {
    opserr << "FATAL DistributedDiagonalSOE::getX - vectX == 0";
    exit(-1);
  }    
  return *vectX;
}

const Vector &
DistributedDiagonalSOE::getB(void)
{
  if (vectB == 0) {
    opserr << "FATAL DistributedDiagonalSOE::getB - vectB == 0";
    exit(-1);
  }        
  return *vectB;
}

double 
DistributedDiagonalSOE::normRHS(void)
{
  double norm =0.0;
  for (int i=0; i<size; i++) {
    double Yi = B[i];
    norm += Yi*Yi;
  }
  return sqrt(norm);
  
}    


int
DistributedDiagonalSOE::setDiagonalSolver(DistributedDiagonalSolver &newSolver)
{
  newSolver.setLinearSOE(*this);
  
  if (size != 0) {
    int solverOK = newSolver.setSize();
    if (solverOK < 0) {
      opserr << "WARNING:DistributedDiagonalSOE::setSolver :";
      opserr << "the new solver could not setSeize() - staying with old\n";
      return -1;
    }
  }
  
  return this->setSolver(newSolver);
}


int 
DistributedDiagonalSOE::sendSelf(int cTag, Channel &theChannel)
{
  processID = 0;
  int sendID =0;
  
  // if P0 check if already sent. If already sent use old processID; if not allocate a new process 
  // id for remote part of object, enlarge channel * to hold a channel * for this remote object.
  
  // if not P0, send current processID
  
  if (processID == 0) {

    // check if already using this object
    bool found = false;
    for (int i=0; i<numChannels; i++)
      if (theChannels[i] == &theChannel) {
	sendID = i+1;
	found = true;
      }

    // if new object, enlarge Channel pointers to hold new channel * & allocate new ID
    if (found == false) {
      int nextNumChannels = numChannels + 1;
      Channel **nextChannels = new Channel *[nextNumChannels];
      if (nextNumChannels == 0) {
	opserr << "DistributedBandGenLinSOE::sendSelf() - failed to allocate channel array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	nextChannels[i] = theChannels[i];
      nextChannels[numChannels] = &theChannel;
      numChannels = nextNumChannels;
      
      if (theChannels != 0)
	delete [] theChannels;
      
      theChannels = nextChannels;
      
  if (localCol != 0)
	delete [] localCol;
      localCol = new ID *[numChannels];
      if (localCol == 0) {
	opserr << "DistributedBandGenLinSOE::sendSelf() - failed to allocate id array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	localCol[i] = 0;    

      // allocate new processID for remote object
      sendID = numChannels;
    }

  } else 
    sendID = processID;

  // send remotes processID
  ID idData(1);
  idData(0) = sendID;
  
  int res = theChannel.sendID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedDiagonalSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  LinearSOESolver *theSoeSolver = this->getSolver();
  if (theSoeSolver != 0) {
    if (theSoeSolver->sendSelf(cTag, theChannel) < 0) {
      opserr <<"WARNING DistributedDiagonalSOE::sendSelf() - failed to send solver\n";
      return -1;
    } 
  } else {
    opserr <<"WARNING DistributedDiagonalSOE::sendSelf() - no solver to send!\n";
    return -1;
  }

  return 0;
}


int 
DistributedDiagonalSOE::recvSelf(int cTag, Channel &theChannel, 
		   FEM_ObjectBroker &theBroker)
{
  ID idData(1);
  int res = theChannel.recvID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedProfileSPDLinSOE::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  DistributedDiagonalSolver *theDistributedDiagonalSolver = new DistributedDiagonalSolver();
  if (theDistributedDiagonalSolver->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr <<"WARNING DistributedBandgenLinSOE::sendSelf() - failed to recv solver\n";
    return -1;
  }
  theDistributedDiagonalSolver->setLinearSOE(*this);
  this->setSolver(*theDistributedDiagonalSolver);
  return 0;
}

int
DistributedDiagonalSOE::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];


  localCol = new ID *[nChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  return 0;
}

int 
DistributedDiagonalSOE::setAnalysisModel(AnalysisModel &theAnalysisModel)
{
  theModel = &theAnalysisModel;
  return 0;
}
