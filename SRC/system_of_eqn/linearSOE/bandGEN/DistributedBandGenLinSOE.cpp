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
                                                                        
// $Revision: 1.4 $
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandGEN/DistributedBandGenLinSOE.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the implementation for BandGenLinSOE


#include <DistributedBandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <BandGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

DistributedBandGenLinSOE::DistributedBandGenLinSOE(BandGenLinSolver &theSolvr)
  :BandGenLinSOE(LinSOE_TAGS_DistributedBandGenLinSOE), 
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0), sizeWork(0), myB(0), myVectB(0)
{
	this->setSolver(theSolvr);
    theSolvr.setLinearSOE(*this);
}

DistributedBandGenLinSOE::DistributedBandGenLinSOE()
  :BandGenLinSOE(LinSOE_TAGS_DistributedBandGenLinSOE), 
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0), sizeWork(0), myB(0), myVectB(0)
{

}


DistributedBandGenLinSOE::~DistributedBandGenLinSOE()
{
  if (theChannels != 0)
    delete [] theChannels;

  if (localCol != 0)
    for (int i=0; i<numChannels; i++)
      if (localCol[i] != 0)
	delete localCol[i];
  delete [] localCol;

  if (workArea != 0)
    delete [] workArea;

  if (myVectB != 0)
    delete myVectB;

  if (myB == 0)
    delete [] myB;
}


int 
DistributedBandGenLinSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  int maxNumSubVertex = 0;

  // if subprocess, collect graph, send it off, 
  // vector back containing size of system, etc.

  if (processID != 0) {
    Channel *theChannel = theChannels[0];
    theGraph.sendSelf(0, *theChannel);
    
    static ID data(3);
    theChannel->recvID(0, 0, data);
    size = data(0);
    numSubD = data(1);
    numSuperD = data(2);

    ID *subMap = new ID(theGraph.getNumVertex());
    localCol[0] = subMap;
    Vertex *vertex;
    VertexIter &theSubVertices = theGraph.getVertices();
    int cnt = 0;
    while((vertex = theSubVertices()) != 0) 
      (*subMap)(cnt++) = vertex->getTag();

    theChannel->sendID(0, 0, *subMap);
  } 
  
  // if main domain, collect graphs from all subdomains,
  // merge into 1, number this one, send to subdomains the
  // id containing dof tags & start id's.
  else {

    // from each distributed soe recv it's graph
    // and merge them into primary graph

    FEM_ObjectBroker theBroker;
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      Graph theSubGraph;
      theSubGraph.recvSelf(0, *theChannel, theBroker);
      theGraph.merge(theSubGraph);

      int numSubVertex = theSubGraph.getNumVertex();
      ID *subMap = new ID(numSubVertex);
      localCol[j] = subMap;
      if (numSubVertex > maxNumSubVertex)
	maxNumSubVertex = numSubVertex;
    }

    size = theGraph.getNumVertex();

    //opserr << "\n\n\n**********************************************\n"; 
    //theGraph.Print(opserr); 
    //opserr << "\n**********************************************\n\n\n"; 
    //
    // determine the number of superdiagonals and subdiagonals
    //

    numSubD = 0;
    numSuperD = 0;

    Vertex *vertexPtr;
    VertexIter &theVertices = theGraph.getVertices();
    
    while ((vertexPtr = theVertices()) != 0) {
      int vertexNum = vertexPtr->getTag();
      const ID &theAdjacency = vertexPtr->getAdjacency();
      for (int i=0; i<theAdjacency.Size(); i++) {
	int otherNum = theAdjacency(i);
	int diff = vertexNum - otherNum;
	if (diff > 0) {
	  if (diff > numSuperD)
	    numSuperD = diff;
	} else 
	  if (diff < numSubD)
	    numSubD = diff;
      }
    }
    numSubD *= -1;
    
    static ID data(3);
    data(0) = size;
    data(1) = numSubD;
    data(2) = numSuperD;

    // to each distributed soe send the size data
    // and merge them into primary graph

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendID(0, 0, data);

      ID *subMap = localCol[j];
      theChannel->recvID(0, 0, *subMap);
    }    
  }

  int numCols = theGraph.getNumVertex();

  if (processID != 0) {
    ID &globalMap = *(localCol[0]);
    ID *localMap = new ID(size);
  
    localMap->Zero();
    for (int k=0; k< globalMap.Size(); k++)
      (*localMap)(globalMap(k)) = k; 
    delete localCol[0];
    localCol[0] = localMap;
  }


   int newSize;
   if (processID != 0)
   	newSize = numCols * (2*numSubD + numSuperD + 1);
   else
   	newSize = size * (2*numSubD + numSuperD + 1);

   // newSize = size * (2*numSubD + numSuperD + 1);

  if (newSize != Asize) { // we have to get another space for A

    if (processID == 0) {
      if (workArea != 0)
	delete [] workArea;

      workArea = new double [newSize];
      sizeWork = newSize;
    }
    sizeWork = numCols * (2*numSubD + numSuperD +1);

    if (A != 0) 
      delete [] A;

    A = new double[newSize];
    
    if (A == 0) {
      opserr << "WARNING DistributedBandGenLinSOE::DistributedBandGenLinSOE :";
      opserr << " ran out of memory for A (size,super,sub) (";
      opserr << size <<", " << numSuperD << ", " << numSubD << ") \n";
      Asize = 0; size = 0; numSubD = 0; numSuperD = 0;
      result= -1;
    }
    else  
      Asize = newSize;
  }
  
  // zero the matrix
  for (int i=0; i<Asize; i++)
    A[i] = 0;
  
  factored = false;


  
  if (size > Bsize) { // we have to get space for the vectors
    
    // delete the old	
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (myB != 0) delete [] myB;
    
    // create the new
    B = new double[size];
    X = new double[size];
    myB = new double[size];
    
    if (B == 0 || X == 0 || myB == 0) {
      opserr << "WARNING DistributedBandGenLinSOE::DistributedBandGenLinSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      Bsize = 0; size = 0; numSubD = 0; numSuperD = 0;
      result = -1;
    }
    else 
      Bsize = size;
  }
  
  // zero the vectors
  for (int j=0; j<size; j++) {
    B[j] = 0;
    X[j] = 0;
    myB[j] = 0;
  }
  
  // get new Vector objects if size has changes
  if (oldSize != size) {
    if (vectX != 0) 
      delete vectX;
    
    if (vectB != 0) 
      delete vectB;
    
    if (myVectB != 0)
      delete myVectB;

    vectX = new Vector(X,size);
    vectB = new Vector(B,size);
    myVectB = new Vector(myB, size);
  }

    // invoke setSize() on the Solver
  LinearSOESolver *theSolvr = this->getSolver();
  if (theSolvr == 0) {
	opserr << "WARNING:DistributedBandGenLinSOE::setSize :";
    opserr << " no solver set\n";
    return -2;
  }    

  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:DistributedBandGenLinSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    
  
  return result;    
}


int 
DistributedBandGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
  // check that m and id are of similar size
  int idSize = id.Size();    
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "BandGenLinSOE::addA()	- Matrix and ID not of similar sizes\n";
    return -1;
  }
  
  int ldA = 2*numSubD + numSuperD + 1;
  ID *theMap = 0;
  if (numChannels > 0)
    theMap = localCol[0];

  if (fact == 1.0) { // do not need to multiply 
    for (int i=0; i<idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
	// double *coliiPtr = A + col*ldA + numSubD + numSuperD;
	
	double *coliiPtr;
	if (processID == 0)
	  coliiPtr = A + col*ldA + numSubD + numSuperD;
	else
	  coliiPtr = A + ((*theMap)(col))*ldA + numSubD + numSuperD;
	
	for (int j=0; j<idSize; j++) {
	  int row = id(j);
	  if (row <size && row >= 0) {		    
	    int diff = col - row;
	    if (diff > 0) {
	      if (diff <= numSuperD) {
		double *APtr = coliiPtr - diff;
		*APtr += m(j,i);
	      }			
	      
	    } else {
	      diff *= -1;
	      if (diff <= numSubD) {
		double *APtr = coliiPtr + diff;
		*APtr += m(j,i);
	      }
	    }
	  }
	}  // for j
      } 
    }  // for i
  } else {
    for (int i=0; i<idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
	// double *coliiPtr = A + col*ldA + numSubD + numSuperD;
	
	double *coliiPtr;
	if (processID == 0)
	  coliiPtr = A + col*ldA + numSubD + numSuperD;
	else
	  coliiPtr = A + (*theMap)(col)*ldA + numSubD + numSuperD;

	for (int j=0; j<idSize; j++) {
	  int row = id(j);
	  if (row <size && row >= 0) {		    
	    int diff = col - row;
	    if (diff > 0) {
	      if (diff <= numSuperD) {
		double *APtr = coliiPtr - diff;
		*APtr += m(j,i) *fact;
	      }
	    } else {
	      diff *= -1;
	      if (diff <= numSubD) {
		double *APtr = coliiPtr + diff;
		*APtr += m(j,i) *fact;
	      }
	    }
	  }
	}  // for j
      } 
    }  // for i
  }    

  return 0;
}

int 
DistributedBandGenLinSOE::solve(void)
{
  static ID result(1);

  //
  // if subprocess send B and A and receive back result X, B & result
  //

  if (processID != 0) {
    Channel *theChannel = theChannels[0];

    // send B
    theChannel->sendVector(0, 0, *myVectB);

    // send A in packets placed in vector X
    //    Vector vectA(A, Asize);    
    if (factored == false) {
      Vector vectA(A, sizeWork);        
      theChannel->sendVector(0, 0, vectA);
    }

    // receive X,B and result
    theChannel->recvVector(0, 0, *vectX);
    theChannel->recvVector(0, 0, *vectB);
    theChannel->recvID(0, 0, result);
    
    factored = true;
  } 

  //
  // if main process, recv B & A from all, solve and send back X, B & result
  //

  else {

    // add P0 contribution to B
    *vectB = *myVectB;

    // receive X and A contribution from subprocess & add them in
    for (int j=0; j<numChannels; j++) {

      // get X & add
      Channel *theChannel = theChannels[j];
      theChannel->recvVector(0, 0, *vectX);
      *vectB += *vectX;

      // get A & add using local map
      if (factored == false) {
	const ID &localMap = *(localCol[j]);
	int ldA = 2*numSubD + numSuperD + 1;    
	int localSize = localMap.Size() * ldA;
	Vector vectA(workArea, localSize);    
	theChannel->recvVector(0, 0, vectA);

	int loc = 0;
	for (int i=0; i<localMap.Size(); i++) {
	  int pos = localMap(i)*ldA;
	  for (int k=0; k<ldA; k++) 
	    A[pos++] += workArea[loc++];
	}    
      }
    }
    
    // solve
    result(0) = this->LinearSOE::solve();

    // send results back
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendVector(0, 0, *vectX);
      theChannel->sendVector(0, 0, *vectB);
      theChannel->sendID(0, 0, result);      
    }
  } 

  return result(0);
}	


int 
DistributedBandGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    
    // check for a quick return 
    if (fact == 0.0)  return 0;

    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "BandSPDLinSOE::addB() - Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		myB[pos] += v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		myB[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		myB[pos] += v(i) * fact;
	}
    }	
    return 0;
}


int
DistributedBandGenLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	opserr << "WARNING DistributedBandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<size; i++) {
	    myB[i] = v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<size; i++) {
	    myB[i] = -v(i);
	}
    } else {
	for (int i=0; i<size; i++) {
	    myB[i] = v(i) * fact;
	}
    }	
    return 0;
}

void 
DistributedBandGenLinSOE::zeroB(void)
{
  double *Bptr = myB;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}


const Vector &
DistributedBandGenLinSOE::getB(void)
{
  if (processID != 0) {
    Channel *theChannel = theChannels[0];

    // send B & recv merged B
    theChannel->sendVector(0, 0, *myVectB);
    theChannel->recvVector(0, 0, *vectB);
  } 

  //
  // if main process, recv B & A from all, solve and send back X, B & result
  //

  else {

    *vectB = *myVectB;

    Vector remoteB(workArea, size);    
    // receive X and A contribution from subprocess & add them in

    for (int j=0; j<numChannels; j++) {

      Channel *theChannel = theChannels[j];
      theChannel->recvVector(0, 0, remoteB);
      *vectB += remoteB;
    }
  
    // send results back
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendVector(0, 0, *vectB);
    }
  } 

  return *vectB;
}	

  
int 
DistributedBandGenLinSOE::sendSelf(int commitTag, Channel &theChannel)
{
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
  
  int res = theChannel.sendID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedBandGenLinSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  LinearSOESolver *theSoeSolver = this->getSolver();
  if (theSoeSolver != 0) {
    if (theSoeSolver->sendSelf(commitTag, theChannel) < 0) {
      opserr <<"WARNING DistributedBandGenLinSOE::sendSelf() - failed to send solver\n";
      return -1;
    } 
  } else {
    opserr <<"WARNING DistributedBandGenLinSOE::sendSelf() - no solver to send!\n";
    return -1;
  }

  return 0;
}


int 
DistributedBandGenLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  ID idData(1);
  int res = theChannel.recvID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedBandGenLinSOE::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  BandGenLinSolver *theBandGenSolver = new BandGenLinLapackSolver();
  if (theBandGenSolver->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr <<"WARNING DistributedBandgenLinSOE::sendSelf() - failed to recv solver\n";
    return -1;
  }
  theBandGenSolver->setLinearSOE(*this);
  this->setSolver(*theBandGenSolver);

  return 0;
}


int
DistributedBandGenLinSOE::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
DistributedBandGenLinSOE::setChannels(int nChannels, Channel **theC)
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
