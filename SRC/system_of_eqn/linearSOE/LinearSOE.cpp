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
// $Date: 2009-05-14 22:45:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/LinearSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
//
// Description: This file contains the implementation of LinearSOE.
//
// What: "@(#) LinearSOE.C, revA"

#include<LinearSOE.h>
#include<LinearSOESolver.h>
#include<Matrix.h>

LinearSOE::LinearSOE(LinearSOESolver &theLinearSOESolver, int classtag)
    :MovableObject(classtag), theModel(0), theSolver(&theLinearSOESolver)
{

}

LinearSOE::LinearSOE(int classtag)
:MovableObject(classtag), theModel(0), theSolver(0)
{

}


LinearSOE::~LinearSOE()
{
  if (theSolver != 0)
    delete theSolver;
}

int 
LinearSOE::solve(void)
{
  if (theSolver != 0)
    return (theSolver->solve());
  else 
    return -1;
}

int
LinearSOE::formAp(const Vector &p, Vector &Ap)
{
  return 0;
}

double
LinearSOE::getDeterminant(void)
{
  if (theSolver != 0)
    return theSolver->getDeterminant();
  else 
    return 0;
}



int 
LinearSOE::setSolver(LinearSOESolver &newSolver)
{
    theSolver = &newSolver;
    return 0;
}

LinearSOESolver *
LinearSOE::getSolver(void)
{
    return theSolver;
}

int 
LinearSOE::setLinks(AnalysisModel &theModel)
{
    this->theModel = &theModel;
    return 0;
}


int
LinearSOE::addA(const Matrix &) {
  return -1;
}

int
LinearSOE::addColA(const Vector &col, int colIndex, double fact) {
  return -1;
}

int LinearSOE::saveSparseA(OPS_Stream& output, int baseIndex) {
  Matrix* A = const_cast<Matrix*>(this->getA());
  
  if (A == nullptr) {
    return -1;
  }
  
  int rows = A->noRows();
  int cols = A->noCols();
  int nnz = rows * cols;
  
  // Assume the header is already written to output stream
  
  output << rows << " " << cols << " " << nnz << "\n";
  
  // Write all elements with base index
  int nnz_written = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double val = (*A)(i,j);
      output << i + baseIndex << " " << j + baseIndex << " " << val << "\n";
      nnz_written++;
    }
  }
  if (nnz_written != nnz) {
    return -1;
  }
  return 0;
}