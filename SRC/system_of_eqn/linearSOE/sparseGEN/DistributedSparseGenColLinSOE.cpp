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
                                                                        
// $Revision: 1.3 $
// $Date: 2009-05-11 20:57:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenColLinSOE.cpp,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation for SparseGenColLinSOE


#include <DistributedSparseGenColLinSOE.h>
#include <DistributedSuperLU.h>
#include <SparseGenColLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

DistributedSparseGenColLinSOE::DistributedSparseGenColLinSOE(SparseGenColLinSolver &theSolvr)
  :SparseGenColLinSOE(theSolvr, LinSOE_TAGS_DistributedSparseGenColLinSOE), 
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0), sizeWork(0), myB(0),
   myVectB(0)
  
{
    theSolvr.setLinearSOE(*this);
}

DistributedSparseGenColLinSOE::DistributedSparseGenColLinSOE()
  :SparseGenColLinSOE(LinSOE_TAGS_DistributedSparseGenColLinSOE), 
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0), sizeWork(0), myB(0),
   myVectB(0)
  
{

}


DistributedSparseGenColLinSOE::~DistributedSparseGenColLinSOE()
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
DistributedSparseGenColLinSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  int maxNumSubVertex = 0;

  // if subprocess, collect graph, send it off, 
  // vector back containing size of system, etc.
  if (processID != 0) {
    Channel *theChannel = theChannels[0];
    theGraph.sendSelf(0, *theChannel);
    
    static ID data(2);
    theChannel->recvID(0, 0, data);
    size = data(0);
    nnz = data(1);

    ID *subMap = new ID(theGraph.getNumVertex());
    localCol[0] = subMap;
    Vertex *vertex;
    VertexIter &theSubVertices = theGraph.getVertices();
    int cnt = 0;
    while((vertex = theSubVertices()) != 0) 
      (*subMap)(cnt++) = vertex->getTag();

    theChannel->sendID(0, 0, *subMap);

    if (nnz > Asize) { // we have to get more space for A and rowA
      if (rowA != 0)
	delete [] rowA;
      
      rowA = new int[nnz];
	
      if (rowA == 0) {
	opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
	opserr << " ran out of memory for A and rowA with nnz = ";
	opserr << nnz << " \n";
	size = 0; Asize = 0; nnz = 0;
	result =  -1;
      } 
    }

    if (size > Bsize) { // we have to get space for the vectors
	
      if (colStartA != 0) 
	delete [] colStartA;
      colStartA = new int[size+1]; 
      
      if (colStartA == 0) {
	opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
	opserr << " ran out of memory for vectors (size) (";
	opserr << size << ") \n";
	size = 0; Bsize = 0;
	result =  -1;
      }    
    }

    ID rowAdata(rowA, nnz);
    ID colStartAdata(colStartA, size+1);
    theChannel->recvID(0, 0, rowAdata);
    theChannel->recvID(0, 0, colStartAdata);
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
  
    //
    // determine the number of non-zeros
    //

    Vertex *theVertex;
    VertexIter &theVertices = theGraph.getVertices();
    nnz = 0;
    while ((theVertex = theVertices()) != 0) {
	const ID &theAdjacency = theVertex->getAdjacency();
	nnz += theAdjacency.Size() +1; // the +1 is for the diag entry
    }

    static ID data(2);
    data(0) = size;
    data(1) = nnz;

    // to each distributed soe send the size data
    // and merge them into primary graph

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendID(0, 0, data);

      ID *subMap = localCol[j];
      theChannel->recvID(0, 0, *subMap);
    }    

    if (nnz > Asize) { // we have to get more space for A and rowA
      if (rowA != 0)
	delete [] rowA;

      if (workArea != 0) 
	delete [] workArea;

      rowA = new int[nnz];
      workArea = new double[nnz];

      sizeWork = nnz;
	
      if (rowA == 0 || workArea == 0) {
	opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
	opserr << " ran out of memory for A and rowA with nnz = ";
	opserr << nnz << " \n";
	size = 0; Asize = 0; nnz = 0;
	result =  -1;
      } 
    }

    if (size > Bsize) { // we have to get space for the vectors
	
      if (colStartA != 0) 
	delete [] colStartA;
      colStartA = new int[size+1]; 
      
      if (colStartA == 0) {
	opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
	opserr << " ran out of memory for vectors (size) (";
	opserr << size << ") \n";
	size = 0; Bsize = 0;
	result =  -1;
      }    
    }

    // fill in colStartA and rowA
    if (size != 0) {
      colStartA[0] = 0;
      int startLoc = 0;
      int lastLoc = 0;
      for (int a=0; a<size; a++) {
	
	theVertex = theGraph.getVertexPtr(a);
	if (theVertex == 0) {
	  opserr << "WARNING:SparseGenColLinSOE::setSize :";
	  opserr << " vertex " << a << " not in graph! - size set to 0\n";
	  size = 0;
	  return -1;
	}
	
	rowA[lastLoc++] = theVertex->getTag(); // place diag in first
	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	
	// now we have to place the entries in the ID into order in rowA
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
	colStartA[a+1] = lastLoc;;	    
	startLoc = lastLoc;
      }
    }


    ID rowAdata(rowA, nnz);
    ID colStartAdata(colStartA, size+1);

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendID(0, 0, rowAdata);
      theChannel->sendID(0, 0, colStartAdata);
    }
  }

  if (nnz > Asize) { // we have to get more space for A and rowA
    if (A != 0) 
      delete [] A;
    
    A = new double[nnz];
    
    if (A == 0 || rowA == 0) {
      opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
      opserr << " ran out of memory for A and rowA with nnz = ";
      opserr << nnz << " \n";
      size = 0; Asize = 0; nnz = 0;
      result =  -1;
    } 
    
    Asize = nnz;
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
    
    if (B == 0 || X == 0 || colStartA == 0 || myB == 0) {
      opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
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
    if (vectX != 0)
      delete vectX;
    
    if (vectB != 0)
      delete vectB;
    
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);	
    myVectB = new Vector(myB, size);
  }

  LinearSOESolver *theSolvr = this->getSolver();
  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:DistributedSparseGenColLinSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    

  return result;    
}


int 
DistributedSparseGenColLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  
	return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "SparseGenColLinSOE::addA() ";
	opserr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply 
      for (int i=0; i<idSize; i++) {
	int col = id(i);
	if (col < size && col >= 0) {
	  int startColLoc = colStartA[col];
	  int endColLoc = colStartA[col+1];
	  for (int j=0; j<idSize; j++) {
	    int row = id(j);
	    if (row <size && row >= 0) {
	      // find place in A using rowA
	      for (int k=startColLoc; k<endColLoc; k++)
		if (rowA[k] == row) {
		  A[k] += m(j,i);
		  k = endColLoc;
		}
	    }
	  }  // for j		
	} 
      }  // for i
    } else {
      for (int i=0; i<idSize; i++) {
	int col = id(i);
	if (col < size && col >= 0) {
	  int startColLoc = colStartA[col];
	  int endColLoc = colStartA[col+1];
	  for (int j=0; j<idSize; j++) {
	    int row = id(j);
	    if (row <size && row >= 0) {
	      // find place in A using rowA
	      for (int k=startColLoc; k<endColLoc; k++)
		if (rowA[k] == row) {
		  A[k] += fact * m(j,i);
		  k = endColLoc;
		}
	    }
	  }  // for j		
	} 
      }  // for i
    }
    return 0;
}

int 
DistributedSparseGenColLinSOE::addB(const Vector &v, const ID &id, double fact)
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
DistributedSparseGenColLinSOE::setB(const Vector &v, double fact)
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
DistributedSparseGenColLinSOE::zeroB(void)
{
  double *Bptr = myB;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}

int 
DistributedSparseGenColLinSOE::solve(void)
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
    if (factored == false) {
      Vector vectA(A, nnz);    
      theChannel->sendVector(0, 0, vectA);
    }

    LinearSOESolver *theSoeSolver = this->getSolver();
    if (theSoeSolver->getClassTag() == SOLVER_TAGS_DistributedSuperLU)
      this->LinearSOE::solve();

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

      if (factored == false) {
	Vector vectA(workArea, nnz);
	theChannel->recvVector(0, 0, vectA);
	for (int i=0; i<nnz; i++)
	  A[i] += workArea[i];	
      }	

      /*
      // get A & add using local map
      const ID &localMap = *(localCol[j]);
      int localSize = localMap.Size() * half_band;
      Vector vectA(workArea, localSize);    
      theChannel->recvVector(0, 0, vectA);

      int loc = 0;
      for (int i=0; i<localMap.Size(); i++) {
	int pos = localMap(i)*half_band;
	for (int k=0; k<half_band; k++) 
	  A[pos++] += workArea[loc++];
      }    
      */
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


const Vector &
DistributedSparseGenColLinSOE::getB(void)
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
DistributedSparseGenColLinSOE::sendSelf(int commitTag, Channel &theChannel)
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
	opserr << "DistributedSparseGenColLinSOE::sendSelf() - failed to allocate channel array of size: " << 
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
	opserr << "DistributedSparseGenColLinSOE::sendSelf() - failed to allocate id array of size: " << 
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
    opserr <<"WARNING DistributedSparseGenColLinSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  return 0;
}


int 
DistributedSparseGenColLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  ID idData(1);
  int res = theChannel.recvID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedSparseGenColLinSOE::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;


  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;
  opserr << "DistributedSparseGenColLinSOE::recvSelf() - DONE\n";

  DistributedSuperLU *theSolvr = new DistributedSuperLU();
  theSolvr->setLinearSOE(*this);
  this->setSolver(*theSolvr);

  return 0;
}


int
DistributedSparseGenColLinSOE::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
DistributedSparseGenColLinSOE::setChannels(int nChannels, Channel **theC)
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
