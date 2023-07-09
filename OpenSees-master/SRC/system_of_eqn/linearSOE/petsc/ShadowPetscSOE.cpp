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
// $Date: 2003-02-14 23:02:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/ShadowPetscSOE.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/petsc/ShadowPetscSOE.C
//
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the implementation for BandGenLinSOE

#include <string.h>
#include "ShadowPetscSOE.h"
#include "PetscSolver.h"
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

ShadowPetscSOE::ShadowPetscSOE(PetscSolver &theSOESolver, int bs)
:LinearSOE(theSOESolver, LinSOE_TAGS_ShadowPetscSOE),
 theSOE(theSOESolver, bs), theSolver(&theSOESolver), myRank(0),
 sendBuffer(0), blockSize(bs)
{

  MPI_Comm_rank(PETSC_COMM_WORLD, &myRank);
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcessors);
  // MPI_Comm_dup(PETSC_COMM_WORLD, &theComm);
  if (myRank != 0) {
    opserr << " ShadowPetscSOE::ShadowPetscSOE - must be rank 0\n";
  }
  sendBuffer = (void *)(&sendData[0]);
  MPI_Barrier(PETSC_COMM_WORLD);

/************
  // send data on the PetscSolver
  int petscMethod;
  int petscPre;

  KSPType method = theSOESolver.method;
  PCType preconditioner = theSOESolver.preconditioner;

  if (strcmp(method,KSPCG) == 0)
    petscMethod = 1;
  else if (strcmp(method,KSPGMRES) == 0)
    petscMethod = 2;
  else if (strcmp(method,KSPBCGS) == 0)
    petscMethod = 3;
  else if (strcmp(method,KSPCGS) == 0)
    petscMethod = 4;
  else if (strcmp(method,KSPTFQMR) == 0)
    petscMethod = 5;
  else if (strcmp(method,KSPTCQMR) == 0)
    petscMethod = 6;
  else if (strcmp(method,KSPCR) == 0)
    petscMethod = 7;
  else if (strcmp(method, KSPLSQR) == 0)
    petscMethod = 8;
  else if (strcmp(method, KSPRICHARDSON) == 0)
    petscMethod = 9;
  else if (strcmp(method,KSPCHEBYCHEV) == 0)
    petscMethod = 10;
  else if (strcmp(method, KSPPREONLY) == 0)
    petscMethod = 11;
  else {
    opserr << "ShadowPetscSOE::ShadowPetscSOE - unknown KSP method\n";
    petscMethod = 12;
  }

  if (strcmp(preconditioner,PCNONE) == 0)
    petscPre = 1;
  else if (strcmp(preconditioner,PCJACOBI) == 0)
    petscPre = 2;
  else if (strcmp(preconditioner,PCSOR) == 0)
    petscPre = 3;
  else if (strcmp(preconditioner,PCEISENSTAT) == 0)
    petscPre = 4;
  else if (strcmp(preconditioner,PCBJACOBI) == 0)
    petscPre = 5;
  else if (strcmp(preconditioner,PCASM) == 0)
    petscPre = 6;
  else if (strcmp(preconditioner,PCILU) == 0)
    petscPre = 7;
  else if (strcmp(preconditioner,PCICC) == 0)
    petscPre = 8;  
  else if (strcmp(preconditioner,PCBGS) == 0)
    petscPre = 9;
  else if (strcmp(preconditioner,PCMG) == 0)
    petscPre = 10;
  else if (strcmp(preconditioner,PCSHELL) == 0)
    petscPre = 11;
  else if (strcmp(preconditioner,PCLU) == 0)
    petscPre = 12;
  else {
    opserr << "ActorPetscSOE::ActorPetscSOE - unknown PC method\n";
    petscPre = 13;
  }

  // now send the data to the remote actor objects
  sendData[0] = petscMethod;
  sendData[1] = petscPre;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);    
**********/

}


int
ShadowPetscSOE::getNumEqn(void) const
{
    return theSOE.size;
}
    
ShadowPetscSOE::~ShadowPetscSOE()
{
  sendData[0] = 0;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
}


int 
ShadowPetscSOE::solve(void)
{
  sendData[0] = 1;
  sendData[1] = theSOE.isFactored;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);
  return theSOE.solve();
}

int 
ShadowPetscSOE::setSize(Graph &theGraph)
{
  int n = theGraph.getNumVertex();
  int size = n;
  int N = n;

  // fist itearte through the vertices of the graph to get nnz,
  // the number of non-zeros.
  Vertex *theVertex;
  int NNZ = 0;
  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    const ID &theAdjacency = theVertex->getAdjacency();
    NNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
  }

  // we determine the number of non-zeros & number of nonzero's
  // in each row of A
  // create two integer arrays. One containing the column indices for each
  // entry in A stored in order by row and column number. 
  // The other a pointer into this array for each row.
  int *rowStartA = new int[size];  // start locations of rows in colA
  int *colA = new int[NNZ]; // column locations, stored row-wise

  // fill in rowStartA and colA
  if (size != 0) {
    rowStartA[0] = 0;
    int startLoc = 0;
    int lastLoc = 0;
    for (int a=0; a<size; a++) {

      theVertex = theGraph.getVertexPtr(a);
      if (theVertex == 0) {
	opserr << "WARNING:SparseGenLinSOE::setSize :";
	opserr << " vertex " << a << " not in graph! - size set to 0\n";
	size = 0;
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

  // now for each of the SOE's we determine how to invoke setSizeParallel()
  sendData[0] = 2;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);

  int dnz = 0;
  int onz = 0;
  int numRowsTyp = N/blockSize/numProcessors;
  numRowsTyp *= blockSize;
  int numRowsLast = numRowsTyp + (N - numProcessors*numRowsTyp); // lastProc extra rows
  int *dnnz = new int[numRowsLast];
  int *onnz = new int[numRowsLast];

  // first determine start and end rows of diag block
  int numRows = numRowsLast;
  int endRow = N-1;
  int startRow = endRow-numRows +1; 
  
  int result = 0;
  // we have to go last to first
  for (int i=numProcessors-1; i>=0; i--) {

    // for each processor determine onnz and dnnz info in colA[]
    for (int j=0; j<numRows; j++) {
      dnnz[j] = 0; 
      onnz[j] = 0;
      int rowStart = rowStartA[startRow+j];
      int nextRowStart = rowStartA[startRow+j+1];
      for (int k=rowStart; k<nextRowStart; k++) {
	int col = colA[k];
	if (col < startRow || col > endRow)
	  onnz[j] += 1;
	else
	  dnnz[j] += 1;
      }
    }

    // now send the data
    if (i != 0) {  // remote object
      sendData[0] = 2;
      sendData[1] = numRows;
      sendData[2] = n;
      
      int tag = 100;
      MPI_Send(sendBuffer, 3, MPI_INT, i, tag, PETSC_COMM_WORLD);

      tag = 101;
      void *buffer = (void *)dnnz;
      MPI_Send(buffer, numRows, MPI_INT, i, tag, PETSC_COMM_WORLD);

      tag = 102;
      buffer = (void *)onnz;
      MPI_Send(buffer, numRows, MPI_INT, i, tag, PETSC_COMM_WORLD);
    }      

    // increment startRow, endRow and numRows if necessary
    endRow = startRow-1;
    numRows = numRowsTyp;
    startRow = endRow-numRows +1; 
  }

  // we broadcast again before we start setSizeParallel()
  // this is because Petsc all processes need to call setup at same
  // time .. if don't we hang
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);  
  result = theSOE.setSizeParallel(numRows, n, dnz, dnnz, onz, onnz);

  // free up the memory used
  /*
  delete [] colA;
  delete [] rowStartA;
  delete [] onnz;
  delete [] dnnz;
  */

  return result;
}



int 
ShadowPetscSOE::addA(const Matrix &m, const ID &id, double fact)
{
  return theSOE.addA(m, id, fact);
}

    
int 
ShadowPetscSOE::addB(const Vector &v, const ID &id, double fact)
{
  return theSOE.addB(v, id, fact);
}


int
ShadowPetscSOE::setB(const Vector &v, double fact)
{
  return theSOE.setB(v, fact);
}


void 
ShadowPetscSOE::zeroA(void)
{
  sendData[0] = 3;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);
  theSOE.zeroA();
}
	
void 
ShadowPetscSOE::zeroB(void)
{
  sendData[0] = 4;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);
  theSOE.zeroB();
}


const Vector &
ShadowPetscSOE::getX(void)
{
  sendData[0] = 5;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);


  double *theData =0;
  double *X = theSOE.X;

  // put the local data into the X array
  int high, low;
  int ierr = VecGetOwnershipRange(theSOE.x, &low, &high);CHKERRA(ierr);       

  ierr = VecGetArray(theSOE.x, &theData); CHKERRA(ierr);       
  opserr << "ShadowPetscSOE::getX - low:high " << low << " " << high << endln;
  for (int i=low; i<=high; i++)
    X[i] = theData[i-low];
  ierr = VecRestoreArray(theSOE.x, &theData); CHKERRA(ierr);             

  // now from each actor object get it's local copy
  // and place in the array
  int tag = 99;
  for (i=1; i<numProcessors; i++) {
    MPI_Status status;
    MPI_Recv(sendBuffer, 3, MPI_INT, i, tag, PETSC_COMM_WORLD, &status);
    low = sendData[0];
    high = sendData[1];
  opserr << "ShadowPetscSOE::getX - low:high " << low << " " << high << endln;
    double *dataBuf = &X[low];
    void *buffer = (void *)dataBuf;
    MPI_Recv(buffer, high-low, MPI_DOUBLE, i, tag, PETSC_COMM_WORLD, &status);
  }

  return *(theSOE.vectX);
}


const Vector &
ShadowPetscSOE::getB(void)
{
  sendData[0] = 6;
  MPI_Bcast(sendBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);

  // STOP ****** some more work here
  return theSOE.getB();
}


double 
ShadowPetscSOE::normRHS(void)
{
  theSOE.getB();
  double norm =0.0;
  int size = theSOE.size;
  double *Bptr = theSOE.B;
  for (int i=0; i<size; i++) {
    double Yi = *Bptr++;
    norm += Yi*Yi;
  }
  return sqrt(norm);
}    


void 
ShadowPetscSOE::setX(int loc, double value)
{
  theSOE.setX(loc, value);
}

int
ShadowPetscSOE::setSolver(PetscSolver &newSolver)
{
  opserr << "ShadowPetscSOE::setSolver - not yet working\n";
  return -1;
}


int 
ShadowPetscSOE::sendSelf(int cTag, Channel &theChannel)
{
  opserr << "WARNING ShadowPetscSOE::sendSelf - does not send itself YET\n";
  return 0;
}


int 
ShadowPetscSOE::recvSelf(int cTag, 
			 Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  opserr << "WARNING ShadowPetscSOE::sendSelf - does not receive itself YET\n";
  return 0;
}





