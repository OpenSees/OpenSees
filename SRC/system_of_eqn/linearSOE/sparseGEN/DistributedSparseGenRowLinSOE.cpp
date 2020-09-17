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
                                                                        
// $Revision: 1.2 $
// $Date: 2006-10-02 20:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenRowLinSOE.cpp,v $
                                                                        
// Written: fmk 
// Created: 04/05
//
// Description: This file contains the implementation for DistributedSparseGenRowLinSOE

#include <DistributedSparseGenRowLinSOE.h>
#include <DistributedSparseGenRowLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

DistributedSparseGenRowLinSOE::DistributedSparseGenRowLinSOE(DistributedSparseGenRowLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_DistributedSparseGenRowLinSOE),
 n(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
 vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false),
 numP(0), rank(-1), startRow(-1), endRow(-1), numRows(0),
 rows(0), otherProcessesRows(0), otherProcessesRowStart(0)
{
    the_Solver.setLinearSOE(*this);
}


DistributedSparseGenRowLinSOE::~DistributedSparseGenRowLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (rowStartA != 0) delete [] rowStartA;
    if (colA != 0) delete []colA;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;        
    if (otherProcessesRows != 0) delete [] otherProcessesRows;
    if (otherProcessesRowStart != 0) delete [] otherProcessesRowStart;
}


int
DistributedSparseGenRowLinSOE::getNumEqn(void) const
{
    return n;
}

int 
DistributedSparseGenRowLinSOE::setSize(Graph &theGraph)
{
  int result = 0;
  /*
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
    int numShared = otherSize(0);
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

  //
  // now primary list has been created we recv from each process the edges corresponding
  // to shared vertices.
  // 



  // 
  // determine system size
  //
  
  // first determine the largest in graph on local
  int myN = 0;
  VertexIter &theVertices = theGraph.getVertices();
  Vertex *theVertex;
  while ((theVertex = theVertices()) != 0) {
    if (theVertex->getTag() > myN)
      myN = theVertex->getTag();
  }
    
  // each process sends local, n is max of them all
  n = myN;
  for (int i=0; i<numP; i++) {
    static ID rankSendRecv(1);
    if (rank == i) {
      rankSendRecv = myN;
      theChannel->sendID(0,0,rankSendRecv);
    } else {
      theChannel->recvID(0,0,rankSendRecv);
      if (rankSendRecv(0) > n)
	n = rankSendRecv(0);
    }
  }
  
  //
  // now determine rows assigned to each process
  //
  
  numRows = n/numP;
  startRow = numRows*rank;
    endRow = numRows*(rank+1);
  } else {
    startRow = numRows*(numP-1);
    endRow = n;
  }
  
  numRows = endRow-startRow;

    //
    // now create a graph of all rows assigned to processor & any rows
    // that are included from local graph (theGraph passed in method call)
    //

    //
    // graph must obviously contain the original plus any new Vertices from startRow
    // through endRow not already in Graph; rows ID must also contain any vertex in Graph
    // but not in current startRow through endRow.
    //

    Graph myGraph(theGraph);
    rows.resize(numRows);

    for (int i=0; i<numRows; i++) {
      int row = startRow + i;
      rows(i) = row;
      if ((theVertex = myGraph.getVertexPtr(row)) == 0) {
	theVertex = new Vertex(row, row);
	if (myGraph.addVertex(theVertex, false) == false) {
	  opserr << "WARNING DistributedDistributedSparseGenRowLinSOE::addVertex - error adding vertex\n";
	}
      }
    }
      
    int numAllRows = numRows;
    VertexIter &theVertices2 = myGraph.getVertices();
    while ((theVertex = theVertices2()) != 0) {      
      int row = theVertex->getTag();
      int rowLoc = rows.getLocation(row);
      if (rowLoc == -1) {
	rows[numAllRows] = row;
	numAllRows++;
      }
    }

    rowStartA = new int[numAllRows+1]; 

    // graph must also contain all edges that exist in other graphs for vertices
    // in numRows; we will only add edges & not vertices; graph will end up undirected
    // which is why we can't use addEdge Graph method but must add edges to Vertices

    for (int i=0; i<numP; i++) {
      if (rank == i) {

	//
	// of course we only need to send the vertices of the graph that we are not responsible for!
	//

	Graph graphToSend;
	Vertex *theVertex;
	VertexIter &theVertices = theGraph.getVertices();	
	while ((theVertex = theVertices()) != 0) {
	  int vertexTag = theVertex->getTag();
	  if (vertexTag < startRow && vertexTag >= endRow) {
	    Vertex *copyVertex = new Vertex(*theVertex);
	    graphToSend.addVertex(copyVertex, false);
	  }
	}

	// now send this new graph
	graphToSend.sendSelf(0, *theChannel);
      } else {

	// receive the graph, for each vertex in recvd graph, if vertex is in myGraph
	// add any edges not already present in newly recvd graph to vertex in myGraph.

	FEM_ObjectBroker theBroker;
	Graph otherGraph;
	otherGraph.recvSelf(0, *theChannel, theBroker);
	VertexIter &theVertices = otherGraph.getVertices();
	Vertex *otherVertex;
	while ((otherVertex = theVertices()) != 0) {
	  int otherTag = otherVertex->getTag();
	  int otherTagLoc = rows.getLocation(otherTag);
	  if (otherTagLoc >= 0) {
	    theVertex = myGraph.getVertexPtr(otherTag);
	    const ID &otherAdjacency = otherVertex->getAdjacency();
	    for (int j=0; j<otherAdjacency.Size(); j++) {
	      theVertex->addEdge(otherAdjacency(j));
	    }
	  }
	}
      }
    }

    int newNNZ = 0;
    VertexIter &theVertices3 = myGraph.getVertices();
    while ((theVertex = theVertices3()) != 0) {
	const ID &theAdjacency = theVertex->getAdjacency();
	newNNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
    }
    nnz = newNNZ;

    if (newNNZ > Asize) { // we have to get more space for A and colA
	if (A != 0) 
	    delete [] A;
	if (colA != 0)
	    delete [] colA;
	
	A = new double[newNNZ];
	colA = new int[newNNZ];
	
        if (A == 0 || colA == 0) {
            opserr << "WARNING DistributedSparseGenRowLinSOE::DistributedSparseGenRowLinSOE :";
	    opserr << " ran out of memory for A and colA with nnz = ";
	    opserr << newNNZ << " \n";
	    n = 0; Asize = 0; nnz = 0;
	    result =  -1;
        } 
	
	Asize = newNNZ;
    }

    // zero the matrix
    for (int i=0; i<Asize; i++)
	A[i] = 0;
	
    factored = false;
    
    if (n > Bsize) { // we have to get space for the vectors
	
	// delete the old	
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;
	if (rowStartA != 0) delete [] rowStartA;

	// create the new
	B = new double[n];
	X = new double[n];
	
        if (B == 0 || X == 0) {
            opserr << "WARNING DistributedSparseGenRowLinSOE::DistributedSparseGenRowLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << n << ") \n";
	    n = 0; Bsize = 0;
	    result =  -1;
        }
	else
	    Bsize = n;
    }

    // zero the vectors
    for (int j=0; j<n; j++) {
	B[j] = 0;
	X[j] = 0;
    }
    
    // create new Vectors objects
    if (n != oldSize) {
	if (vectX != 0)
	    delete vectX;

	if (vectB != 0)
	    delete vectB;
	
	vectX = new Vector(X, n);
	vectB = new Vector(B, n);	
    }

    // fill in rowStartA and colA
    if (n != 0) {
      rowStartA[0] = 0;
      int startLoc = 0;
      int lastLoc = 0;
      for (int a=0; a<numAllRows; a++) {

	theVertex = theGraph.getVertexPtr(rows(a));
	if (theVertex == 0) {
	  opserr << "WARNING:DistributedSparseGenRowLinSOE::setSize :";
	  opserr << " vertex " << a << " not in graph! - n set to 0\n";
	  n = 0;
	  return -1;
	}

	colA[lastLoc++] = theVertex->getTag(); // place diag in first
	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	
	// now we have to place the entries in the ID into order in colA
	for (int i=0; i<idSize; i++) {

	  int row = theAdjacency(i);
	  bool foundPlace = false;
	  // find a place in colA for current col
	  for (int j=startLoc; j<lastLoc; j++)
	    if (colA[j] > row) { 
	      // move the entries already there one further on
	      // and place col in current location
	      for (int k=lastLoc; k>j; k--)
		
		colA[k] = colA[k-1];
	      colA[j] = row;
	      foundPlace = true;
	      j = lastLoc;
	    }
	  if (foundPlace == false) // put in at the end
	    colA[lastLoc] = row;

	  lastLoc++;
	}
	rowStartA[a+1] = lastLoc;;	    
	startLoc = lastLoc;
      }
    }
    
    // invoke setSize() on the Solver    
    LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:DistributedSparseGenRowLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    
  */
    return result;
}

int 
DistributedSparseGenRowLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  
	return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar n
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "DistributedSparseGenRowLinSOE::addA() ";
	opserr << " - Matrix and ID not of similar ns\n";
	return -1;
    }
    
    for (int i=0; i<idSize; i++) {
      int row = id(i);
      if (row < n && row >= 0) {
	int localRow = rows.getLocation(row);
	if (localRow == -1) 
	  opserr << "DistributedSparseGenRowLinSOE::addA() - you goofed\n";
	int startRowLoc = rowStartA[localRow];
	int endRowLoc = rowStartA[localRow+1];
	for (int j=0; j<idSize; j++) {
	  int col = id(j);
	  if ((col < n) && (col >= 0)) {
	    // find place in A using colA
	    for (int k=startRowLoc; k<endRowLoc; k++)
	      if (colA[k] == col) {
		A[k] += m(i,j) * fact;
		k = endRowLoc;
	      }
	  }
	}  // for j		
      } 
    }  // for i

    return 0;
}

    
int 
DistributedSparseGenRowLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar n
    if (idSize != v.Size() ) {
	opserr << "DistributedSparseGenRowLinSOE::addB() ";
	opserr << " - Vector and ID not of similar ns\n";
	return -1;
    }    

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos < n && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <n && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <n && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	

    return 0;
}


int
DistributedSparseGenRowLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != n) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible ns " << n << " and " << v.Size() << endln;
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<n; i++) {
	    B[i] = v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<n; i++) {
	    B[i] = -v(i);
	}
    } else {
	for (int i=0; i<n; i++) {
	    B[i] = v(i) * fact;
	}
    }	
    return 0;
}

void 
DistributedSparseGenRowLinSOE::zeroA(void)
{
    double *Aptr = A;
    for (int i=0; i<Asize; i++)
	*Aptr++ = 0;

    factored = false;
}
	
void 
DistributedSparseGenRowLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<n; i++)
	*Bptr++ = 0;
}

void 
DistributedSparseGenRowLinSOE::setX(int loc, double value)
{
    if (loc < n && loc >=0)
	X[loc] = value;
}

void 
DistributedSparseGenRowLinSOE::setX(const Vector &x)
{
  if (x.Size() == n && vectX != 0)
    *vectX = x;
}

const Vector &
DistributedSparseGenRowLinSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL DistributedSparseGenRowLinSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
DistributedSparseGenRowLinSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL DistributedSparseGenRowLinSOE::getB - vectB == 0";
	exit(-1);
    }        
    return *vectB;
}

double 
DistributedSparseGenRowLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<n; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    


int
DistributedSparseGenRowLinSOE::setDistributedSparseGenRowSolver(DistributedSparseGenRowLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (n != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:DistributedSparseGenRowLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
DistributedSparseGenRowLinSOE::sendSelf(int cTag, Channel &theChannel)
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

  return 0;
}


int 
DistributedSparseGenRowLinSOE::recvSelf(int cTag, Channel &theChannel, 
		   FEM_ObjectBroker &theBroker)
{
  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  return 0;
}

int
DistributedSparseGenRowLinSOE::setChannels(int nChannels, Channel **theC)
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



