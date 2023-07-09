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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/DistributedProfileSPDLinSOE.cpp,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation for DistributedProfileSPDLinSOE

#include <DistributedProfileSPDLinSOE.h>
#include <ProfileSPDLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

DistributedProfileSPDLinSOE::DistributedProfileSPDLinSOE(ProfileSPDLinSolver &theSolvr)
  :ProfileSPDLinSOE(theSolvr, LinSOE_TAGS_DistributedProfileSPDLinSOE), 
   processID(0), numChannels(0), theChannels(0), 
   localCol(0), sizeLocal(0), workArea(0), sizeWork(0), myVectB(0), myB(0)
{
    theSolvr.setLinearSOE(*this);
}


DistributedProfileSPDLinSOE::DistributedProfileSPDLinSOE()
  :ProfileSPDLinSOE(LinSOE_TAGS_DistributedProfileSPDLinSOE), 
   processID(0), numChannels(0), theChannels(0), 
   localCol(0), sizeLocal(0), workArea(0), sizeWork(0), myVectB(0), myB(0)
{

}


DistributedProfileSPDLinSOE::~DistributedProfileSPDLinSOE()
{
  if (theChannels != 0)
    delete [] theChannels;

  if (localCol != 0)
    for (int i=0; i<numChannels; i++)
      if (localCol[i] != 0)
	delete localCol[i];
  delete [] localCol;

  if (sizeLocal != 0)
    delete sizeLocal;

  if (workArea != 0)
    delete [] workArea;

  if (myVectB != 0)
    delete myVectB;

  if (myB != 0)
    delete [] myB;
}


int 
DistributedProfileSPDLinSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  int maxNumSubVertex = 0;
  int myProfileSize =0;

  // if subprocess, collect graph, send it off, 
  // vector back containing size of system, etc.
  if (processID != 0) {
    Channel *theChannel = theChannels[0];

    theGraph.sendSelf(0, *theChannel);
    
    static ID data(1);
    theChannel->recvID(0, 0, data);
    size = data(0);

    if (size > Bsize) { 
	if (iDiagLoc != 0) delete [] iDiagLoc;
	iDiagLoc = new int[size];
    }
    
    // receive my iDiagLoad
    ID iLoc(size);    
    theChannel->recvID(0, 0, iLoc);

    // determine iDiagLoc, profileSize & local mapping
    ID *subMap = new ID(theGraph.getNumVertex());
    localCol[0] = subMap;
    Vertex *vertex;
    VertexIter &theSubVertices = theGraph.getVertices();
    int cnt = 0; int loc = 0;  
    while((vertex = theSubVertices()) != 0) {
      int vertexTag = vertex->getTag();
      (*subMap)(cnt) = vertexTag;
      int colHeight =0;
      if (vertexTag == 0) {
	myProfileSize++;
	colHeight= 1;
      } else {
	colHeight = iLoc[vertexTag] - iLoc[vertexTag-1];	    
	myProfileSize += colHeight;
      }
      loc += colHeight;
      iDiagLoc[cnt++] = loc;	      
    }
    (*sizeLocal)(0) = myProfileSize;

    // send local mapping & profile size to P0
    theChannel->sendID(0, 0, *subMap);
    theChannel->sendID(0, 0, *sizeLocal);
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
      
      if (theSubGraph.recvSelf(0, *theChannel, theBroker) < 0) {
	opserr << "WARNING DistributedProfileSPDLinSOE::setSize() : ";
	opserr << " - failed to recv graph\n";
      }

      theGraph.merge(theSubGraph);

      int numSubVertex = theSubGraph.getNumVertex();
      ID *subMap = new ID(numSubVertex);
      localCol[j] = subMap;
      if (numSubVertex > maxNumSubVertex)
	maxNumSubVertex = numSubVertex;
    }

    size = theGraph.getNumVertex();
  
    // check we have enough space in iDiagLoc and iLastCol
    // if not delete old and create new
    if (size != Bsize) { 
	if (iDiagLoc != 0) delete [] iDiagLoc;
	iDiagLoc = new int[size];

	if (iDiagLoc == 0) {
	    opserr << "WARNING DistributedProfileSPDLinSOE::setSize() : ";
	    opserr << " - ran out of memory for iDiagLoc\n";
	    size = 0; Asize = 0;
	    result = -1;
	}
    }

    // zero out iDiagLoc 
    for (int i=0; i<size; i++) {
	iDiagLoc[i] = 0;
    }

    // now we go through the vertices to find the height of each col and
    // width of each row from the connectivity information.
    
    Vertex *vertexPtr;
    VertexIter &theVertices = theGraph.getVertices();

    while ((vertexPtr = theVertices()) != 0) {
	int vertexNum = vertexPtr->getTag();
	const ID &theAdjacency = vertexPtr->getAdjacency();
	int iiDiagLoc = iDiagLoc[vertexNum];
	int *iiDiagLocPtr = &(iDiagLoc[vertexNum]);

	for (int i=0; i<theAdjacency.Size(); i++) {
	    int otherNum = theAdjacency(i);
	    int diff = vertexNum-otherNum;
	    if (diff > 0) {
		if (iiDiagLoc < diff) {
		    iiDiagLoc = diff;
		    *iiDiagLocPtr = diff;
		}
	    } 
	}
    }


    // now go through iDiagLoc, adding 1 for the diagonal element
    // and then adding previous entry to give current location.
    if (iDiagLoc != 0)
	iDiagLoc[0] = 1; // NOTE FORTRAN ARRAY LOCATION - 1 of solvers uses library

    for (int j=1; j<size; j++)
	iDiagLoc[j] = iDiagLoc[j] + 1 + iDiagLoc[j-1];

    static ID data(1);
    data(0) = size;

    ID iLoc(iDiagLoc, size);

    // to each distributed soe send the size data
    // and merge them into primary graph

    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendID(0, 0, data);
      theChannel->sendID(0, 0, iLoc);


      static ID remoteLocalSize(1);
      ID *subMap = localCol[j];
      theChannel->recvID(0, 0, *subMap);
      theChannel->recvID(0, 0, remoteLocalSize);
      (*sizeLocal)(j) = remoteLocalSize(0);
    }    

    myProfileSize = iDiagLoc[size-1];
  }

  profileSize = myProfileSize;
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


  // check if we need more space to hold A
  // if so then go get it
  if (profileSize != Asize) { 

    // create some space for storing data sent from processes
    if (processID == 0) {
      if (workArea != 0)
	delete [] workArea;

      workArea = new double [profileSize];
      sizeWork = profileSize;
    }
    
    // delete old space
    if (A != 0)
      delete [] A;
    
    // get new space
    A = new double[profileSize];
    
    if (A == 0) {
      opserr << "DistributedProfileSPDLinSOE::ProfileSPDLinSOE :";
      opserr << " ran out of memory for A (size,Profile) (";
      opserr << size <<", " << profileSize << ") \n";
      size = 0;  Asize = 0;  profileSize = 0;
      result = -1;
    }
    else 
      Asize = profileSize;
  }
  
  // zero the matrix
  for (int k=0; k<profileSize; k++)
    A[k] = 0;
  
  isAfactored = false;
  isAcondensed = false;    
  
  if (size > Bsize) { // we have to get another space for A
    
    // delete the old
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (myB != 0) delete [] myB;
    
    // create the new	
    B = new double[size];
    X = new double[size];
    myB = new double[size];

    if (B == 0 || X == 0 || myB == 0) {
      opserr << "DistributedProfileSPDLinSOE::ProfileSPDLinSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; Bsize = 0;
      result = -1;
    }
  }
  
  // zero the vectors
  for (int l=0; l<size; l++) {
    B[l] = 0;
    X[l] = 0;
    myB[l] = 0;

  }
  
  if (size != oldSize) {
    
    if (vectX != 0)
      delete vectX;

    if (vectB != 0)
      delete vectB;

    if (myVectB != 0)
      delete myVectB;
    
    if (myVectB != 0)
      delete myVectB;

    vectX = new Vector(X,size);
    vectB = new Vector(B,size);
    myVectB = new Vector(myB, size);
    
    if (size > Bsize)
      Bsize = size;
  }
  
  // invoke setSize() on the Solver
  LinearSOESolver *the_Solver = this->getSolver();
  int solverOK = the_Solver->setSize();
  if (solverOK < 0) {
    opserr << "WARNING DistributedProfileSPDLinSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    

  return result;    
}


int 
DistributedProfileSPDLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
  // check that m and id are of similar size
  int idSize = id.Size();    
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "DistributedProfileSPDLinSOE::addA()	- Matrix and ID not of similar sizes\n";
    return -1;
  }

  ID *theMap = 0;
  if (numChannels > 0) 
    theMap = localCol[0];
  
  if (fact == 1.0) { // do not need to multiply 
    for (int i=0; i<idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
	//	double *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
	double *coliiPtr;
	if (processID == 0)
	  coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
	else
	  coliiPtr = &A[iDiagLoc[(*theMap)(col)] -1];

	for (int j=0; j<idSize; j++) {
	  int row = id(j);
	  if (row <size && row >= 0 && 
	      row <= col) { 
	    
	    // we only add upper and inside profile
	    double *APtr = coliiPtr + (row-col);
	    *APtr += m(j,i);
	  }
	}  // for j
      } 
    }  // for i
  } else {
    for (int i=0; i<idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
	//  ydouble *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 		
	double *coliiPtr;
	if (processID == 0)
	  coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
	else
	  coliiPtr = &A[iDiagLoc[(*theMap)(col)] -1];
			
	for (int j=0; j<idSize; j++) {
	  int row = id(j);
	  if (row <size && row >= 0 && 
	      row <= col) { 
	    
	    // we only add upper and inside profile
	    double *APtr = coliiPtr + (row-col);
	    *APtr += m(j,i) * fact;
	  }
	}  // for j
      } 
    }  // for i
    
  }
  return 0;
}

int 
DistributedProfileSPDLinSOE::solve(void)
{
  static ID result(1);

  //
  // if subprocess send B and A and receive back result X, B & result
  //

  if (processID != 0) {
    Channel *theChannel = theChannels[0];

    // send B
    theChannel->sendVector(0, 0, *myVectB);

    if (isAfactored == false) {
      // send A in packets placed in vector X
      Vector vectA(A, (*sizeLocal)(0));    
      theChannel->sendVector(0, 0, vectA);
    }
    // receive X,B and result
    theChannel->recvVector(0, 0, *vectX);
    theChannel->recvVector(0, 0, *vectB);
    theChannel->recvID(0, 0, result);
    isAfactored = true;
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
      if (isAfactored == false) {
	const ID &localMap = *(localCol[j]);
	int localSize = (*sizeLocal)(j);
	Vector vectA(workArea, localSize);    
	theChannel->recvVector(0, 0, vectA);
	
	int loc = 0;
	for (int i=0; i<localMap.Size(); i++) {
	  int col = localMap(i);
	  int colSize, pos;
	  
	  if (col == 0) {
	    colSize = 1;
	    pos = 0;
	  }
	  else {
	    pos = iDiagLoc[col-1];
	    colSize = iDiagLoc[col] - iDiagLoc[col-1];
	  }
	  for (int k=0; k<colSize; k++) {
	    A[pos++] += workArea[loc++];
	  }
	}
      }    
    }

    // solve
    result(0) = this->LinearSOE::solve();

    //    opserr << *vectX;

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
DistributedProfileSPDLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    
    // check for a quick return 
    if (fact == 0.0)  return 0;

    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "DistributedProfileSPDLinSOE::addB() -";
	opserr << " Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		myB[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		myB[pos] -= v(i);
	}
    } else {
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		myB[pos] += v(i) * fact;
	}
    }	
    return 0;
}

int
DistributedProfileSPDLinSOE::setB(const Vector &v, double fact)
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
DistributedProfileSPDLinSOE::zeroB(void)
{
  double *Bptr = myB;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}


const Vector &
DistributedProfileSPDLinSOE::getB(void)
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

    // receive X and A contribution from subprocess & add them in

    *vectB = *myVectB;
    Vector remoteB(workArea, size);    

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
DistributedProfileSPDLinSOE::sendSelf(int commitTag, Channel &theChannel)
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
	opserr << "DistributedProfileSPDLinSOE::sendSelf() - failed to allocate channel array of size: " << 
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
	opserr << "DistributedProfileSPDLinSOE::sendSelf() - failed to allocate id array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	localCol[i] = 0;    

      if (sizeLocal != 0)
	delete sizeLocal;
      
      sizeLocal = new ID(numChannels);

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
    opserr <<"WARNING DistributedProfileSPDLinSOE::sendSelf() - failed to send data\n";
    return -1;
  }
  
  return 0;
}


int 
DistributedProfileSPDLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

  ID idData(1);
  int res = theChannel.recvID(0, commitTag, idData);
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

  if (sizeLocal != 0)
    delete sizeLocal;

  sizeLocal = new ID(numChannels);

  return 0;
}


int
DistributedProfileSPDLinSOE::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
DistributedProfileSPDLinSOE::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;


  if (sizeLocal != 0)
    delete sizeLocal;

  sizeLocal = new ID(numChannels);

  return 0;
}
