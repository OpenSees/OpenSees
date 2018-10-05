/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
* University of California, Berkeley, is strictly prohibited.  See   **
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
// $Date: 2009-05-11 20:56:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation for MumpsParallelSOE

#include <MumpsParallelSOE.h>
#include <MumpsParallelSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

MumpsParallelSOE::MumpsParallelSOE(MumpsParallelSolver &theSolvr, int matType)
  :MumpsSOE(theSolvr, LinSOE_TAGS_MumpsParallelSOE, matType), 
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0), 
   sizeWork(0), myB(0), myVectB(0)
{
    theSolvr.setLinearSOE(*this);
}


MumpsParallelSOE::MumpsParallelSOE()
  :MumpsSOE(LinSOE_TAGS_MumpsParallelSOE), 
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0), 
   sizeWork(0), myB(0), myVectB(0)
{

}


MumpsParallelSOE::~MumpsParallelSOE()
{
  if (theChannels != 0)
    delete [] theChannels;

  if (localCol != 0)
    for (int i=0; i<numChannels; i++)
      if (localCol[i] != 0)
	delete localCol[i];
  delete [] localCol;

  if (myB != 0)
    delete [] myB;

  if (myVectB != 0)
    delete myVectB;
}

int 
MumpsParallelSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  int maxNumSubVertex = 0;
  
  // fist itearte through the vertices of the graph to get nnzLoc and n
  int maxVertexTag = -1;
  Vertex *theVertex;
  int newNNZ = 0;
  size = theGraph.getNumVertex();
  int mySize = size;
  //opserr << "MumpsParallelSOE: size : " << size << endln;

  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    int vertexTag = theVertex->getTag();
    if (vertexTag > maxVertexTag)
      maxVertexTag = vertexTag;
    const ID &theAdjacency = theVertex->getAdjacency();
    newNNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
  }

  if (matType !=  0) {

    // symmetric - allows us to reduce nnz by almost half
    newNNZ -= size;
    newNNZ /= 2;
    newNNZ += size;
  }

  nnz = newNNZ;

  if (processID != 0) {

    //
    // if subprocess, send local max vertexTag (n)
    // recv ax n from P0
    //
    static ID data(1);

    data(0) = maxVertexTag;
    Channel *theChannel = theChannels[0];
    theChannel->sendID(0, 0, data);
    theChannel->recvID(0, 0, data);
    
    size = data(0);

  } else {

    //
    // from each distributed soe recv it's max n and compare; return max n to all
    //

    static ID data(1);
    FEM_ObjectBroker theBroker;
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->recvID(0, 0, data);
      if (data(0) > maxVertexTag)
	maxVertexTag = data(0);
    }

    data(0) = maxVertexTag;

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendID(0, 0, data);
    }
    size = maxVertexTag;
  }

  size+=1; // vertices numbered 0 through n-1

  if (nnz > Asize) { // we have to get more space for A and rowA and colA

    if (A != 0) delete [] A;
    if (rowA != 0) delete [] rowA;
    if (colA != 0) delete [] colA;

    A = new double[nnz];    
    rowA = new int[nnz];
    colA = new int[nnz];
      
    for (int i=0; i<nnz; i++) {
      A[i]=0;
      rowA[i]=0;
      colA[i]=0;
    }

    if (rowA == 0 || A == 0 || colA == 0) {
      opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
      opserr << " ran out of memory for A and rowA with nnz = ";
      opserr << nnz << " \n";
      size = 0; Asize = 0; nnz = 0;
	result =  -1;
    } 
    Asize = nnz;
  }

  if (size > Bsize) { // we have to get space for the vectors

    if (B != 0) delete [] B;
    if (X != 0) delete [] X;    
    if (myB != 0) delete [] myB;
    if (workArea != 0) delete [] workArea;
    if (colStartA != 0)  delete [] colStartA;

    // create the new
    B = new double[size];
    X = new double[size];
    myB = new double[size];
    workArea = new double[size];
    colStartA = new int[size+1]; 
    
    if (B == 0 || X == 0 || colStartA == 0 || workArea == 0 || myB == 0) {
      opserr << "WARNING MumpsSOE::MumpsSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; Bsize = 0;
      result =  -1;
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

  // create new Vectors objects
  if (size != oldSize) {
    if (vectX != 0) delete vectX;
    if (vectB != 0) delete vectB;
    if (myVectB != 0) delete myVectB;
    
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);	
    myVectB = new Vector(myB, size);
  }

  // fill in colStartA and rowA
  if (size != 0) {
    colStartA[0] = 0;
    int startLoc = 0;
    int lastLoc = 0;
    for (int a=0; a<size; a++) {
      
      theVertex = theGraph.getVertexPtr(a);
      if (theVertex != 0) {
	
	int vertexTag = theVertex->getTag();
	rowA[lastLoc++] = vertexTag; // place diag in first
	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	
	// now we have to place the entries in the ID into order in rowA
	
	if (matType != 0) {
	  
	  // symmetric
	  for (int i=0; i<idSize; i++) {
	    int row = theAdjacency(i);
	    if (row > vertexTag) {
	      bool foundPlace = false;
	      // find a place in rowA for current col
	      for (int j=startLoc; j<lastLoc; j++)
		if (rowA[j] > row) { 
		  // move the entries already there one further on
		  // and place col in current location
		  for (int k=lastLoc; k>j; k--)
		    rowA[k] = rowA[k-1];
		  rowA[j] = row;
		  foundPlace = true;
		  j = lastLoc;
		}
	      
	      if (foundPlace == false) // put in at the end
		rowA[lastLoc] = row;
	      lastLoc++;
	    }
	  }
	  
	} else {

	  // unsymmetric	  
	  for (int i=0; i<idSize; i++) {
	    int row = theAdjacency(i);
	    bool foundPlace = false;
	    // find a place in rowA for current col
	    for (int j=startLoc; j<lastLoc; j++)
	      if (rowA[j] > row) { 
		// move the entries already there one further on
		// and place col in current location
		for (int k=lastLoc; k>j; k--)
		  rowA[k] = rowA[k-1];
		rowA[j] = row;
		foundPlace = true;
		j = lastLoc;
	      }
	    if (foundPlace == false) // put in at the end
	      rowA[lastLoc] = row;
	    
	    lastLoc++;
	  }
	}
      }
      colStartA[a+1] = lastLoc;
      startLoc = lastLoc;
    }
  }

  // fill in colA
  int count = 0;
  for (int i=0; i<size; i++) {
    for (int k=colStartA[i]; k<colStartA[i+1]; k++)
      colA[count++] = i;
  }

  LinearSOESolver *theSolvr = this->getSolver();

  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:MumpsParallelSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    

  return result;    
}


int 
MumpsParallelSOE::solve(void)
{
  int resSolver = 0;

  //
  // if subprocess send B, solve and recv back X and B
  //

  if (processID != 0) {

    // send B
    Channel *theChannel = theChannels[0];
    theChannel->sendVector(0, 0, *myVectB);

    resSolver =  this->LinearSOE::solve();

    if (resSolver == 0) {
      // receive X,B and result
      theChannel->recvVector(0, 0, *vectX);
      theChannel->recvVector(0, 0, *vectB);
      factored = true;
    }
  } 

  //
  // if main process, recv B & A from all, solve and send back X, B & result
  //
  
  else {
    
    // add P0 contribution to B
    *vectB = *myVectB;
    
    // receive B 
    for (int j=0; j<numChannels; j++) {
      // get X & add
      Channel *theChannel = theChannels[j];
      theChannel->recvVector(0, 0, *vectX);
      *vectB += *vectX;
    }

    // solve
    resSolver = this->LinearSOE::solve();

    // send results back
    if (resSolver == 0) {
      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->sendVector(0, 0, *vectX);
	theChannel->sendVector(0, 0, *vectB);
      }
    }
  } 
  
  return resSolver;
}	



int 
MumpsParallelSOE::addB(const Vector &v, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;

  int idSize = id.Size();    
  // check that m and id are of similar size
  if (idSize != v.Size() ) {
    opserr << "SparseGenColLinSOE::addB() ";
    opserr << " - Vector and ID not of similar sizes\n";
    return -1;
  }    

  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	myB[pos] += v(i);
    }
  } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
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
MumpsParallelSOE::setB(const Vector &v, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;

  //opserr << "MumpsParallelSOE::setB() - start()\n";
  //opserr << v;

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

   //opserr << "MumpsParallelSOE::setB() - end()\n";
  return 0;
}

void 
MumpsParallelSOE::zeroB(void)
{
  double *Bptr = myB;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}


const Vector &
MumpsParallelSOE::getB(void)
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
MumpsParallelSOE::sendSelf(int commitTag, Channel &theChannel)
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
	opserr << "MumpsParallelSOE::sendSelf() - failed to allocate channel array of size: " << 
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
	opserr << "MumpsParallelSOE::sendSelf() - failed to allocate id array of size: " << 
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
  ID idData(2);
  idData(0) = sendID;
  idData(1) = matType;
  
  int res = theChannel.sendID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING MumpsParallelSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  LinearSOESolver *theSoeSolver = this->getSolver();
  if (theSoeSolver != 0) {
    if (theSoeSolver->sendSelf(commitTag, theChannel) < 0) {
      opserr <<"WARNING MumpsParallelSOE::sendSelf() - failed to send solver\n";
      return -1;
    } 
  } else {
    opserr <<"WARNING MumpsParallelSOE::sendSelf() - no solver to send!\n";
    return -1;
  }
    

  return 0;
}


int 
MumpsParallelSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  ID idData(2);
  int res = theChannel.recvID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING MumpsParallelSOE::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);
  matType = idData(1);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;


  MumpsParallelSolver *theSolvr = new MumpsParallelSolver();
  if (theSolvr->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr <<"WARNING MumpsParallelSOE::sendSelf() - failed to recv solver\n";
    return -1;
  }
  
  theSolvr->setLinearSOE(*this);
  this->setSolver(*theSolvr);

  return 0;
}


int
MumpsParallelSOE::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
MumpsParallelSOE::setChannels(int nChannels, Channel **theC)
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
