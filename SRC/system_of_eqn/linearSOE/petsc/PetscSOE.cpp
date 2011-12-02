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
// $Date: 2006-01-13 00:00:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/PetscSOE.cpp,v $
                                                                        
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the implementation for PetscSOE


#include "PetscSOE.h"
#include "PetscSolver.h"
#include <petscvec.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

PetscSOE::PetscSOE(PetscSolver &theSOESolver, int bs)
:LinearSOE(theSOESolver, LinSOE_TAGS_PetscSOE),
 isFactored(0),size(0), processID(0), numProcesses(0), B(0), X(0), 
 indices(0), vectX(0), vectB(0), A(0), x(0), b(0), blockSize(bs),
 numChannels(0), theChannels(0), localCol(0)
{
  theSOESolver.setLinearSOE(*this);
}


int
PetscSOE::getNumEqn(void) const
{
  return size;
}
    
PetscSOE::~PetscSOE()
{
  if (theChannels != 0)
    delete [] theChannels;

  if (localCol != 0)
    for (int i=0; i<numChannels; i++)
      if (localCol[i] != 0)
	delete localCol[i];
  delete [] localCol;
  
  if (vectX != 0) delete vectX;  
  if (vectB != 0) delete vectB;  
  if (B != 0) delete [] B;
  if (X != 0) delete [] X;
  
  // invoke the petsc destructors
  if (A != 0) MatDestroy(A);
  if (b != 0) VecDestroy(b);
  if (x != 0) VecDestroy(x);
}


int 
PetscSOE::setSize(Graph &theGraph)
{
  PetscInitialize(0, PETSC_NULL, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(PETSC_COMM_WORLD, &processID);
  
  int result = 0;
  int ierr = 0;

    // 
    // first determine system size
    //

    if (numProcesses == 1) {
      
      // if single process, the system size is size of graph
      size = theGraph.getNumVertex();
      isFactored = 0;

    } else {

      // first determine local max
      size = 0;
      isFactored = 0;
      VertexIter &theVertices = theGraph.getVertices();
      Vertex *theVertex;
      while ((theVertex = theVertices()) != 0) {
	int vertexTag = theVertex->getTag();
	if (vertexTag > size) 
	  size = vertexTag;
      }

      static ID data(1);

      // all local max's sent to P0 which determines the max
      // and informs all others

      if (processID != 0) {
	Channel *theChannel = theChannels[0];
	
	data(0) = size;
	theChannel->sendID(0, 0, data);
	theChannel->recvID(0, 0, data);
	
	size = data(0);
      } else {

	for (int j=0; j<numChannels; j++) {
	  Channel *theChannel = theChannels[j];
	  theChannel->recvID(0, 0, data);
	  if (data(0) > size)
	    size = data(0);
	}
	data(0) = size;
	for (int j=0; j<numChannels; j++) {
	  Channel *theChannel = theChannels[j];
	  theChannel->sendID(0, 0, data);
	}
      }
      size = size+1; // vertices numbered 0 through n-1
    }

    // invoke the petsc destructors
    if (A != 0) MatDestroy(A);
    if (b != 0) VecDestroy(b);
    if (x != 0) VecDestroy(x);
    
    //
    // now we create the opensees vector objects
    //

    // delete the old vectors
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    
    // create the new
    B = new double[size];
    X = new double[size];
    
    if (B == 0 || X == 0) {
      opserr << "WARNING PetscSOE::PetscSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; 
      result = -1;
    }
    
    // zero the vectors
    for (int j=0; j<size; j++) {
      B[j] = 0;
      X[j] = 0;
    }
    
    if (vectX != 0) 
      delete vectX;
    
    if (vectB != 0) 
      delete vectB;
    
    vectX = new Vector(X, size);
    vectB = new Vector(B, size);

    // 
    // now create petsc matrix and vector objects
    //

    if (numProcesses == 1) {
    
      // we determine the number of non-zeros & number of nonzero's
      // in each row of A
      int *rowA = new int[size];  // will contain no of non-zero entries
      // in each row
      
      int NNZ = 0;
      for (int a=0; a<size; a++) {
	Vertex *theVertex = theGraph.getVertexPtr(a);
	if (theVertex == 0) {
	  opserr << "WARNING:PetscSOE::setSize :";
	  opserr << " vertex " << a << " not in graph! - size set to 0\n";
	  size = 0;
	  return -1;
	}
	
	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	
	NNZ += idSize +1;
	rowA[a] = idSize +1;  // +1 for the diagonal entry
      }

      // 
      // Call Petsc VecCreate & MatCreate; NOTE: using previously allocated storage for vectors
      //      

      //      ierr = PetscOptionsGetInt(PETSC_NULL, "-n", &size, &flg); CHKERRQ(ierr);
      
      if (blockSize == 1) {
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 0, rowA, &A); CHKERRQ(ierr);
      } else {
	ierr = MatCreateSeqBAIJ(PETSC_COMM_SELF, blockSize, size,size, 0, rowA, &A); CHKERRQ(ierr);
      }
      
      ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, size, X, &x); CHKERRQ(ierr); 
      ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, size, B, &b); CHKERRQ(ierr); 

      // invoke setSize() on the Solver
      LinearSOESolver *tSolver = this->getSolver();
      int solverOK = tSolver->setSize();
      if (solverOK < 0) {
	opserr << "WARNING:PetscSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
      }    
      
      // clear up the memory we used
      delete [] rowA;


    } else {

      // 
      // Call Petsc VecCreate & MatCreate; NOTE: using previously allocated storage for vectors
      // 
      //

      ierr = PetscOptionsGetInt(PETSC_NULL, "-n", &size, &flg); CHKERRQ(ierr);
      ierr = MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,size, size, &A); CHKERRQ(ierr);
      ierr = MatSetFromOptions(A);CHKERRQ(ierr);
      ierr = MatGetOwnershipRange(A, &startRow, &endRow);CHKERRQ(ierr);

      ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, endRow-startRow, size, &X[startRow],  &x); CHKERRQ(ierr);
      ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, endRow-startRow, size, &B[startRow],  &b); CHKERRQ(ierr);

      // invoke setSize() on the Solver
      LinearSOESolver *tSolver = this->getSolver();
      int solverOK = tSolver->setSize();
      if (solverOK < 0) {
	opserr << "WARNING:PetscSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
      }    
    }

    return result;    
}


int 
PetscSOE::addA(const Matrix &m, const ID &id, double fact)
{
  isFactored = 0;

    // check for a quick return 
    if (fact == 0.0)  return 0;

    
    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "PetscSOE::addA() - Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    int n = id.Size();
    int row;
    int col;
    double value;
    for (int i=0; i<n; i++) {
      row = id(i);
      if (row >= 0) {
	for (int j=0; j<n; j++) {
	  col = id(j);
	  if (col >= 0) {
	    value = m(i,j)*fact;
	    int ierr = MatSetValues(A,1,&row,1,&col,&value,ADD_VALUES); 
	    if (ierr) opserr << processID << " " << row << " " << col << endln; 
	    CHKERRQ(ierr); 
	  }
	}
      }
    }

    return 0;
}

    
int 
PetscSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "BandGenLinSOE::addB()	- Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	
    return 0;
}


int
PetscSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incomptable sizes " << size << " and " << v.Size() << endln;
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
PetscSOE::zeroA(void)
{
  isFactored = 0;
  MatZeroEntries(A);
}
	
void 
PetscSOE::zeroB(void)
{
  double *Bptr = B;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}


const Vector &
PetscSOE::getX(void)
{
  if (vectX == 0) {
    opserr << "FATAL PetscSOE::getX - vectX == 0!";
    exit(-1);
  }  
  return *vectX;
}


const Vector &
PetscSOE::getB(void)
{
  if (vectB == 0) {
    opserr << "FATAL PetscSOE::getB - vectB == 0!";
    exit(-1);
  }    
  return *vectB;
}


double 
PetscSOE::normRHS(void)
{
  this->getB();
  double norm =0.0;
  double *Bptr = B;
  for (int i=0; i<size; i++) {
    double Yi = *Bptr++;
    norm += Yi*Yi;
  }
  return sqrt(norm);
}    


void 
PetscSOE::setX(int loc, double value)
{
  if (loc < size && loc >= 0)
	X[loc] = value;
}

void 
PetscSOE::setX(const Vector &xData)
{
  if (xData.Size() == size && vectX != 0)
    *vectX = xData;
}

int
PetscSOE::setSolver(PetscSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:PetscSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return solverOK;
	}
    }	
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
PetscSOE::sendSelf(int cTag, Channel &theChannel)
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
PetscSOE::recvSelf(int cTag, Channel &theChannel, 
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
PetscSOE::setChannels(int nChannels, Channel **theC)
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
