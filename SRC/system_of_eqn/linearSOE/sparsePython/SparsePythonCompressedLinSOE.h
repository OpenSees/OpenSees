#ifndef SparsePythonCompressedLinSOE_h
#define SparsePythonCompressedLinSOE_h

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

// Written: Gustavo A. Araujo R. (garaujor@stanford.edu, gaaraujor@gmail.com)
// Created: 11/25
//
// Description: This file contains the class definition for
// SparsePythonCompressedLinSOE. The system stores a general sparse matrix in
// either Compressed Sparse Row (CSR) or Compressed Sparse Column (CSC) layout
// and delegates the linear solve to a Python callable provided to the
// associated SparsePythonCompressedLinSolver instance.

#include <LinearSOE.h>
#include <Vector.h>
#include <vector>

#include "SparsePythonCommon.h"

class Channel;
class FEM_ObjectBroker;
class Graph;
class ID;
class Matrix;
class SparsePythonCompressedLinSolver;

class SparsePythonCompressedLinSOE : public LinearSOE
{
  public:
    using StorageScheme = SparsePythonStorageScheme;

    using MatrixStatus = SparsePythonMatrixStatus;

    SparsePythonCompressedLinSOE(SparsePythonCompressedLinSolver &theSolver,
                                 StorageScheme scheme = StorageScheme::CSR);
    explicit SparsePythonCompressedLinSOE(StorageScheme scheme = StorageScheme::CSR);
    ~SparsePythonCompressedLinSOE() override;

    StorageScheme getStorageScheme() const { return storageScheme; }

    int getNumEqn(void) const override;
    int setSize(Graph &theGraph) override;

    int addA(const Matrix &, const ID &, double fact = 1.0) override;
    int addB(const Vector &, const ID &, double fact = 1.0) override;
    int setB(const Vector &, double fact = 1.0) override;

    void zeroA(void) override;
    void zeroB(void) override;

    const Vector &getX(void) override;
    const Vector &getB(void) override;
    double normRHS(void) override;

    void setX(int loc, double value) override;
    void setX(const Vector &X) override;
    int solve(void) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    int setPythonSolver(SparsePythonCompressedLinSolver &newSolver);

    MatrixStatus getMatrixStatus(void) const;

    const std::vector<int> &getIndexPtr(void) const { return indexPtr; }
    const std::vector<int> &getIndices(void) const { return indices; }
    const std::vector<double> &getValues(void) const { return values; }
    std::vector<double> &getValues(void) { return values; }
    std::vector<double> &getRHS(void) { return rhs; }
    std::vector<double> &getSolution(void) { return sol; }
    const std::vector<double> &getRHS(void) const { return rhs; }
    const std::vector<double> &getSolution(void) const { return sol; }

    std::size_t getNNZ(void) const { return values.size(); }

  protected:
    void updateVectorViews();
    void buildCSRPattern(Graph &theGraph);
    void buildCSCPattern(Graph &theGraph);

    std::vector<int> &accessIndexPtr() { return indexPtr; }
    std::vector<int> &accessIndices() { return indices; }

    const std::vector<int> &getRowPtr() const { return indexPtr; } // compatibility alias
    const std::vector<int> &getColInd() const { return indices; }  // compatibility alias

  private:
    std::vector<int> indexPtr;
    std::vector<int> indices;
    std::vector<double> values;
    std::vector<double> rhs;
    std::vector<double> sol;

    Vector xVector;
    Vector bVector;

    MatrixStatus matrixStatus;
    StorageScheme storageScheme;

    friend class SparsePythonCompressedLinSolver;
};

#endif /* SparsePythonCompressedLinSOE_h */


