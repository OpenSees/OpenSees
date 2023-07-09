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
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DiagonalSOE.cpp,v $

// Written: fmk 
// Created: February 1997
// Revision: A
//
// Description: This file contains the implementation for DiagonalSOE


#include <DiagonalSOE.h>
#include <DiagonalSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

DiagonalSOE::DiagonalSOE(DiagonalSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_DiagonalSOE),
 size(0), A(0), B(0), X(0), vectX(0), vectB(0), isAfactored(false)
{
    the_Solver.setLinearSOE(*this);
}


DiagonalSOE::DiagonalSOE(int N, DiagonalSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_DiagonalSOE),
 size(0), A(0), B(0), X(0), vectX(0), vectB(0), isAfactored(false)
{
  if (size > 0) {
    size = N;
    A = new double[size];
    B = new double[size];
    X = new double[size];	
    
    if (A == 0 || B == 0 || X == 0) {
      opserr << "ERROR DiagonalSOE::DiagonalSOE :";
      opserr << " ran out of memory for size: " << size << endln;
      if (A != 0) delete [] A;
      if (B != 0) delete [] B;
      if (X != 0) delete [] X;
      size = 0;
    }
    
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);
    
    if (vectB == 0 || vectX == 0) {
      opserr << "ERROR DiagonalSOE::DiagonalSOE :";
      opserr << " ran out of memory for size: " << size << endln;
      if (A != 0) delete [] A;
      if (B != 0) delete [] B;
      if (X != 0) delete [] X;
      size = 0;
    }
  }
  the_Solver.setLinearSOE(*this);        
  
  int solverOK = the_Solver.setSize();
  if (solverOK < 0) {
    opserr << "WARNING DiagonalSOE::DiagonalSOE :";
    opserr << " solver failed setSize() in constructor\n";
  }
}

    
DiagonalSOE::~DiagonalSOE()
{
  if (A != 0) delete [] A;
  if (B != 0) delete [] B;
  if (X != 0) delete [] X;
  if (vectX != 0) delete vectX;    
  if (vectB != 0) delete vectB;    
}


int 
DiagonalSOE::getNumEqn(void) const
{
  return size;
}

int 
DiagonalSOE::setSize(Graph &theGraph)
{
  int oldSize = size;
  int result = 0;
  size = theGraph.getNumVertex();
  
  // check we have enough space in iDiagLoc and iLastCol
  // if not delete old and create new
  if (size > oldSize) { 
    
    if (A != 0) delete [] A; A = 0;
    if (B != 0) delete [] B; B = 0;
    if (X != 0) delete [] X; X = 0;
    A = new double[size];
    B = new double[size];
    X = new double[size];	
    
    if (A == 0 || B == 0 || X == 0) {
      opserr << "ERROR DiagonalSOE::setSize() - ";
      opserr << " ran out of memory for size: " << size << endln;
      if (A != 0) delete [] A;
      if (B != 0) delete [] B;
      if (X != 0) delete [] X;
      size = 0;
      return -1;
    }
  }

  if (size != oldSize && size != 0) {
    if (vectX != 0) delete vectX; vectX = 0;
    if (vectB != 0) delete vectB; vectB = 0;
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);
    
    if (vectB == 0 || vectX == 0) {
      opserr << "ERROR DiagonalSOE::setSize() - ";
      opserr << " ran out of memory for size: " << size << endln;
      if (A != 0) delete [] A;
      if (B != 0) delete [] B;
      if (X != 0) delete [] X;
      size = 0;
      return -1;
    }
  }    
    
  // zero the vectors
  for (int l=0; l<size; l++) {
    A[l] = 0;
    B[l] = 0;
    X[l] = 0;
  }
    
  // invoke setSize() on the Solver
  LinearSOESolver *the_Solver = this->getSolver();
  int solverOK = the_Solver->setSize();
  if (solverOK < 0) {
    opserr << "WARNING DiagonalSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    
  
  return result;
}

int 
DiagonalSOE::addA(const Matrix &m, const ID &id, double fact)
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
      if (pos <size && pos >= 0)
	A[pos] += m(i,i);
    }
  } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	A[pos] -= m(i,i);
    }
  } else {
    for (int i=0; i<id.Size(); i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	A[pos] += m(i,i) * fact;
    }
  }	

  return 0;
}
 
    
int 
DiagonalSOE::addB(const Vector &v, const ID &id, double fact)
{
    
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
#ifdef _G3DEBUG
  // check that m and id are of similar size
  int idSize = id.Size();        
  if (idSize != v.Size() ) {
    opserr << "DiagonalSOE::addB() -";
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
DiagonalSOE::setB(const Vector &v, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
  if (v.Size() != size) {
    opserr << "WARNING DiagonalSOE::setB() -";
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
DiagonalSOE::zeroA(void)
{
  double *Aptr = A;
  for (int i=0; i<size; i++)
    *Aptr++ = 0;
  
  isAfactored = false;
}

void 
DiagonalSOE::zeroB(void)
{
  double *Bptr = B;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}

int
DiagonalSOE::formAp(const Vector &p, Vector &Ap)
{
  double *Aptr = A;
  for (int i = 0; i < size; i++, Aptr++)
    Ap(i) = (*Aptr) * p(i);

  return 0;
}

void 
DiagonalSOE::setX(int loc, double value)
{
  if (loc < size && loc >=0)
    X[loc] = value;
}

void 
DiagonalSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
DiagonalSOE::getX(void)
{
  if (vectX == 0) {
    opserr << "FATAL DiagonalSOE::getX - vectX == 0";
    exit(-1);
  }    
  return *vectX;
}

const Vector &
DiagonalSOE::getB(void)
{
  if (vectB == 0) {
    opserr << "FATAL DiagonalSOE::getB - vectB == 0";
    exit(-1);
  }        
  return *vectB;
}

double 
DiagonalSOE::normRHS(void)
{
  double norm =0.0;
  for (int i=0; i<size; i++) {
    double Yi = B[i];
    norm += Yi*Yi;
  }
  return sqrt(norm);
  
}    


int
DiagonalSOE::setDiagonalSolver(DiagonalSolver &newSolver)
{
  newSolver.setLinearSOE(*this);
  
  if (size != 0) {
    int solverOK = newSolver.setSize();
    if (solverOK < 0) {
      opserr << "WARNING:DiagonalSOE::setSolver :";
      opserr << "the new solver could not setSeize() - staying with old\n";
      return -1;
    }
  }
  
  return this->setSolver(newSolver);
}


int 
DiagonalSOE::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}


int 
DiagonalSOE::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}
