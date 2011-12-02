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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/itpack/ItpackLinSOE.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/umfGEN/ItpackLinSOE.h
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for 
// ItpackLinSolver. It solves the ItpackLinSOEobject by calling
// UMFPACK2.2.1 routines.
//
// What: "@(#) ItpackLinSolver.h, revA"

#include <fstream.h>

#include <ItpackLinSOE.h>
#include <ItpackLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

ItpackLinSOE::ItpackLinSOE(ItpackLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_ItpackLinSOE),
 size(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
 vectX(0), vectB(0), Asize(0), Bsize(0), Aformed(false)
{
  the_Solver.setLinearSOE(*this);
}

ItpackLinSOE::~ItpackLinSOE()
{
  if (A != 0) delete [] A;
  if (B != 0) delete [] B;
  if (X != 0) delete [] X;
  if (rowStartA != 0) delete [] rowStartA;
  if (colA != 0) delete [] colA;
  if (vectX != 0) delete vectX;    
  if (vectB != 0) delete vectB;        
}

int
ItpackLinSOE::getNumEqn(void) const
{
  return size;
}

int 
ItpackLinSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  size = theGraph.getNumVertex();
  
  // fist itearte through the vertices of the graph to get nnz
  Vertex *theVertex;
  int newNNZ = 0;
  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    const ID &theAdjacency = theVertex->getAdjacency();
    newNNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
  }
  nnz = newNNZ;
  
  opserr << "ItpackLinSOE::setSize - n " << size << " nnz " << nnz << endln;
  
  if (nnz > Asize) { // we have to get more space for A and colA
    
    if (A != 0) 
      delete [] A;
    if (colA != 0)
      delete [] colA;
    
    A = new double[nnz];
    colA = new int[nnz];
    
    if (A == 0 || colA == 0) {
      opserr << "WARNING ItpackLinSOE::ItpackLinSOE :";
      opserr << " ran out of memory for A and colA with nnz = ";
      opserr << newNNZ << " \n";
      size = 0; Asize = 0; nnz = 0;
      result =  -1;
    } 
    
    Asize = nnz;
  }
  
  // zero the matrix
  for (int i=0; i<Asize; i++)
    A[i] = 0;
  
  if (size > Bsize) { // we have to get space for the vectors
    
    // delete the old	
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (rowStartA != 0) delete [] rowStartA;
    
    // create the new
    B = new double[size];
    X = new double[size];
    rowStartA = new int[size+1]; 
    
    if (B == 0 || X == 0 || rowStartA == 0) {
      opserr << "WARNING ItpackLinSOE::ItpackLinSOE :";
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
  }
  
  // create new Vectors objects
  if (size != oldSize) {
    if (vectX != 0)
      delete vectX;
    
    if (vectB != 0)
      delete vectB;
    
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);	
  }
  
  // fill in rowStartA and colA
  if (size != 0) {
    rowStartA[0] = 0;
    int startLoc = 0;
    int lastLoc = 0;
    for (int a=0; a<size; a++) {
      
      theVertex = theGraph.getVertexPtr(a);
      if (theVertex == 0) {
	opserr << "WARNING:ItpackLinSOE::setSize :";
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
      rowStartA[a+1] = lastLoc;	    
      startLoc = lastLoc;
    }
  }
  
  /*
  // fill out index
  int *indexRowPtr = &index[0];
  int *indexColPtr = &index[nnz];
  for (int ii=0; ii<size; ii++) {
    int rowBegin = rowStartA[ii];
    int rowEnd = rowStartA[ii+1] -1;
    for (int j=rowBegin; j<=rowEnd; j++) {
      int row = ii+1;  // + 1 for fortarn indexing
      int col = colA[j] +1; // +1 for fortran indexing
      *indexRowPtr++ = row;
      *indexColPtr++ = col;
    }
  }
  */

  Aformed = false;

  // invoke setSize() on the Solver    
  LinearSOESolver *the_Solver = this->getSolver();
  int solverOK = the_Solver->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:ItpackLinSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    
  return result;
}

int 
ItpackLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  
    return 0;
  
  int idSize = id.Size();
  
  // check that m and id are of similar size
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "ItpackLinSOE::addA() ";
    opserr << " - Matrix and ID not of similar sizes\n";
    return -1;
  }
  
  if (fact == 1.0) { // do not need to multiply 
    for (int i=0; i<idSize; i++) {
      int row = id(i);
      if (row < size && row >= 0) {
	int startRowLoc = rowStartA[row];
	int endRowLoc = rowStartA[row+1];
	for (int j=0; j<idSize; j++) {
	  int col = id(j);
	  if (col <size && col >= 0) {
	    // find place in A using colA
	    for (int k=startRowLoc; k<endRowLoc; k++)
	      if (colA[k] == col) {
		A[k] += m(i,j);
		k = endRowLoc;
	      }
	  }
	}  // for j		
      } 
    }  // for i
  }
  else {
    for (int i=0; i<idSize; i++) {
      int row = id(i);
      if (row < size && row >= 0) {
	int startRowLoc = rowStartA[row];
	int endRowLoc = rowStartA[row+1];
	for (int j=0; j<idSize; j++) {
	  int col = id(j);
	  if (col <size && col >= 0) {
	    // find place in A using colA
	    for (int k=startRowLoc; k<endRowLoc; k++)
	      if (colA[k] == col) {
		A[k] += fact * m(i,j);
		k = endRowLoc;
	      }
	  }
	}  // for j		
      } 
    }  // for i
  }
  
  return 0;
}

int
ItpackLinSOE::addB(const Vector &v, const ID &id, double fact)
{

  // check for a quick return 
  if (fact == 0.0)
    return 0;
  
  int idSize = id.Size();    
  // check that m and id are of similar size
  if (idSize != v.Size() ) {
    opserr << "ItpackLinSOE::addB() ";
    opserr << " - Vector and ID not of similar sizes\n";
    return -1;
  }    
  
  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	B[pos] += v(i);
    }
  }
  else if (fact == -1.0) { // do not need to multiply if fact == -1.0
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	B[pos] -= v(i) * fact;
    }
  }
  else {
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	B[pos] += v(i) * fact;
    }
  }	
  
  return 0;
}

int
ItpackLinSOE::setB(const Vector &v, double fact)
{
  if (v.Size() != size) {
    opserr << "WARNING ItpackLinSOE::setB() -";
    opserr << " incomptable sizes " << size << " and " << v.Size() << endln;
    return -1;
  }
  
  if (fact == 0.0) {
    for (int i=0; i<size; i++)
      B[i] = 0.0;
  }
  else if (fact == 1.0) {
    for (int i=0; i<size; i++)
      B[i] = v(i);
  }
  else if (fact == -1.0) {
    for (int i=0; i<size; i++)
      B[i] = -v(i);
  }
  else {
    for (int i=0; i<size; i++)
      B[i] = v(i) * fact;
  }	
  return 0;
}

void 
ItpackLinSOE::zeroA(void)
{
  double *Aptr = A;
  for (int i=0; i<Asize; i++)
    *Aptr++ = 0;
  Aformed = false;
}

void 
ItpackLinSOE::zeroB(void)
{
  double *Bptr = B;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}

void 
ItpackLinSOE::setX(int loc, double value)
{
  if (loc < size && loc >=0)
    X[loc] = value;
}

void
ItpackLinSOE::setX(const Vector &x)
{
  if (x.Size() != size) {
    opserr << "WARNING ItpackLinSOE::setX() -";
    opserr << " incomptable sizes " << size << " and " << x.Size() << endln;
  }
  else {
    for (int i = 0; i < size; i++)
      X[i] = x(i);
  }
}

const Vector &
ItpackLinSOE::getX(void)
{
  if (vectX == 0) {
    opserr << "FATAL ItpackLinSOE::getX - vectX == 0";
    exit(-1);
  }
  return *vectX;
}

const Vector &
ItpackLinSOE::getB(void)
{
  if (vectB == 0) {
    opserr << "FATAL ItpackLinSOE::getB - vectB == 0";
    exit(-1);
  }        
  return *vectB;
}

double 
ItpackLinSOE::normRHS(void)
{
  double norm =0.0;
  for (int i=0; i<size; i++) {
    double Yi = B[i];
    norm += Yi*Yi;
  }
  return sqrt(norm);
}    

int
ItpackLinSOE::setItpackLinSolver(ItpackLinSolver &newSolver)
{
  newSolver.setLinearSOE(*this);
  
  if (size != 0) {
    int solverOK = newSolver.setSize();
    if (solverOK < 0) {
      opserr << "WARNING:ItpackLinSOE::setSolver :";
      opserr << "the new solver could not setSeize() - staying with old\n";
      return -1;
    }
  }
  
  return this->LinearSOE::setSolver(newSolver);
}

int 
ItpackLinSOE::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
ItpackLinSOE::recvSelf(int cTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  return -1;
}
