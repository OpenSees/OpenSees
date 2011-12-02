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
// $Date: 2004-10-05 00:17:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/SymBandEigenSOE.cpp,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for
// SymBandEigenSOE, which stores a symmetric banded matrix, A, for
// standard eigenvalue computations.

#include <SymBandEigenSOE.h>
#include <SymBandEigenSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>

SymBandEigenSOE::SymBandEigenSOE(SymBandEigenSolver &theSolvr,  
				 AnalysisModel &aModel)
:EigenSOE(theSolvr, EigenSOE_TAGS_SymBandEigenSOE),
 size(0), numSuperD(0), A(0), 
 Asize(0), factored(false), theModel(&aModel), M(0), Msize(0)
{
  theSolvr.setEigenSOE(*this);
}

int
SymBandEigenSOE::getNumEqn(void) const
{
  return size;
}
    
SymBandEigenSOE::~SymBandEigenSOE()
{
  if (A != 0)
    delete [] A;
  if (M != 0)
    delete [] M;
}

int 
SymBandEigenSOE::setSize(Graph &theGraph)
{
  int result = 0;
  size = theGraph.getNumVertex();
  
  // determine the number of superdiagonals and subdiagonals
  
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
	if (diff < -numSuperD)
	  numSuperD = -diff;
    }
  }
  
  int newSize = size*(numSuperD+1);
  if (newSize > Asize) { // we have to get another space for A
    
    if (A != 0) 
      delete [] A;
    
    A = new double[newSize];
    
    if (A == 0) {
      opserr << "SymBandEigenSOE::setSize() -- ran out of memory for A, size = " <<
	size <<  " and numSuperD = : " << numSuperD << endln;
			      
      Asize = 0; size = 0; numSuperD = 0;
      result= -1;
    }
    else  
      Asize = newSize;
  }
  
  // zero the matrix
  for (int i = 0; i < Asize; i++)
    A[i] = 0.0;
  
  factored = false;
  
  // invoke setSize() on the Solver
  EigenSolver *theSolvr = this->getSolver();
  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "SymBandEigenSOE::setSize() -- solver failed in setSize()\n";
			    
    return solverOK;
  } 
  
  return result;    
}

int 
SymBandEigenSOE::addA(const Matrix &m, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)
    return 0;

  // check that m and id are of similar size
  int idSize = id.Size();    
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "SymBandEigenSOE::addA() -- Matrix and ID not of similar sizes,\n";
    return -1;
  }

  if (fact == 1.0) { // do not need to multiply 
    for (int i = 0; i < idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
	double *coliiPtr = A +(col+1)*(numSuperD+1) - 1;
	int minColRow = col - (numSuperD+1) + 1;
	for (int j = 0; j < idSize; j++) {
	  int row = id(j);
	  if (row <size && row >= 0 && 
	      row <= col && row >= minColRow) { // only add upper
	    double *APtr = coliiPtr + (row-col);
	    *APtr += m(j,i);
	  }
	}  // for j
      } 
    }  // for i
  } else {
    for (int i = 0; i < idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
	double *coliiPtr = A +(col+1)*(numSuperD+1) - 1;
	int minColRow = col - (numSuperD+1) +1;
	for (int j = 0; j < idSize; j++) {
	  int row = id(j);
	  if (row < size && row >= 0 && 
	      row <= col && row >= minColRow) { // only add upper
	    double *APtr = coliiPtr + (row-col);
	    *APtr += m(j,i)*fact;
	  }
	}  // for j
      } 
    }  // for i
  }

  return 0;
}

void 
SymBandEigenSOE::zeroA(void)
{
  double *Aptr = A;
  
  for (int i = 0; i < Asize; i++)
    *Aptr++ = 0;
  
  factored = false;
}

int 
SymBandEigenSOE::addM(const Matrix &m, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)
    return 0;

  // check memory for M has been allocated
  if (M == 0 || Msize != size) {
    if (M != 0)
      delete [] M;
    M = new double[size];
    Msize = size;
    if (M == 0) {
      opserr << "WARNING SymBandEigenSOE::addM() - out of memory creating M for size: " << size << endln;
      return -1;
    }
    for (int i=0; i<size; i++)
      M[i] = 0.0;
  }

  // check that m and id are of similar size
  int idSize = id.Size();    
  if (idSize != m.noRows() && idSize != m.noCols()) {
    opserr << "WARING: SymBandEigenSOE::addM() -- Matrix and ID not of similar sizes!!\n";
    return -1;
  }

  for (int i = 0; i <idSize; i++) {
    int loc = id(i);
    if (loc >= 0)
      M[loc] += fact * m(i,i);
  }


  bool issueWarning = false;
  for (int ii=0; ii<idSize; ii++) 
    for (int jj=0; jj<idSize; jj++)
      if (ii!=jj)
	if (m(ii,jj) != 0.0)
	  issueWarning = true;
  if (issueWarning == true) {
      opserr << "WARNING SymBandEigenSOE::addM() - m passed was not diagonal, only diagonal entries added\n";
  }
  
  return 0;
}

void 
SymBandEigenSOE::zeroM(void)
{
  if (M != 0)
    for (int i=0; i<size; i++)
      M[i] = 0.0;

  return;
}

int 
SymBandEigenSOE::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}
    
int 
SymBandEigenSOE::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  return 0;
}
