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
#include<Vector.h>
#include<ID.h>
#include<OPS_Globals.h>
#include<fstream>
#include<stdint.h>
#include<stdio.h>

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
  const Matrix* A = this->getA();
  
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

int LinearSOE::getSparseA(ID& rowIndices, ID& colIndices, Vector& values, int baseIndex) {
  const Matrix* A = this->getA();
  
  if (A == nullptr) {
    return -1;
  }
  
  int rows = A->noRows();
  int cols = A->noCols();
  int nnz = rows * cols;
  
  // Resize vectors to hold all elements
  rowIndices.resize(nnz);
  colIndices.resize(nnz);
  values.resize(nnz);
  
  // Fill vectors with all elements
  int idx = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      rowIndices(idx) = i + baseIndex;
      colIndices(idx) = j + baseIndex;
      values(idx) = (*A)(i,j);
      idx++;
    }
  }
  
  return 0;
}

int LinearSOE::getSparseA(std::vector<int>& rowIndices, std::vector<int>& colIndices, std::vector<double>& values, int baseIndex) {
  const Matrix* A = this->getA();
  
  if (A == nullptr) {
    return -1;
  }
  
  int rows = A->noRows();
  int cols = A->noCols();
  int nnz = rows * cols;
  
  // Resize vectors to hold all elements
  rowIndices.resize(nnz);
  colIndices.resize(nnz);
  values.resize(nnz);
  
  // Fill vectors with all elements
  int idx = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      rowIndices[idx] = i + baseIndex;
      colIndices[idx] = j + baseIndex;
      values[idx] = (*A)(i,j);
      idx++;
    }
  }
  
  return 0;
}

// Writes the entries of A that getSparseA() returns to a binary file, in the byte order of
// the machine:
//   int64   number of rows, number of columns, number of entries (nnz), baseIndex
//   int32   row index of each entry     (nnz of them)
//   int32   column index of each entry  (nnz of them)
//   double  value of each entry         (nnz of them)
// The values are written exactly, and each array is written in one piece, so a large matrix
// is written in about the time the disk takes, where saveSparseA() formats every entry as text.
int LinearSOE::saveSparseABinary(const char *fileName, int baseIndex) {
  static_assert(sizeof(int) == 4, "printA -binary writes the indices as 32-bit integers");

  std::vector<int> rowIndices, colIndices;
  std::vector<double> values;
  if (this->getSparseA(rowIndices, colIndices, values, baseIndex) != 0)
    return -1;

  const size_t nnz = values.size();
  if (rowIndices.size() != nnz || colIndices.size() != nnz)
    return -1;

  std::ofstream file(fileName, std::ios::out | std::ios::binary | std::ios::trunc);
  if (!file.is_open()) {
    opserr << "WARNING: LinearSOE::saveSparseABinary() - failed to open file: " << fileName << endln;
    return -1;
  }

  const int64_t size = this->getNumEqn();
  const int64_t header[4] = {size, size, (int64_t)nnz, baseIndex};
  file.write((const char *)header, sizeof(header));
  if (nnz > 0) {
    file.write((const char *)rowIndices.data(), nnz * sizeof(int));
    file.write((const char *)colIndices.data(), nnz * sizeof(int));
    file.write((const char *)values.data(), nnz * sizeof(double));
  }
  file.close();
  if (file.fail()) {
    opserr << "WARNING: LinearSOE::saveSparseABinary() - failed to write file: " << fileName << endln;
    remove(fileName);  // not leave a file shorter than its header says
    return -1;
  }

  return 0;
}