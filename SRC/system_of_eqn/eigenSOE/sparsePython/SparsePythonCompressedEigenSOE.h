#ifndef SparsePythonCompressedEigenSOE_h
#define SparsePythonCompressedEigenSOE_h

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
// SparsePythonCompressedEigenSOE. The system stores both stiffness (K) and mass (M)
// matrices in compressed sparse format (CSR or CSC) with shared sparsity pattern
// and delegates the eigen solve to a Python callable provided to the associated
// SparsePythonCompressedEigenSolver instance.

#include <EigenSOE.h>
#include <Vector.h>
#include <vector>

#include "SparsePythonEigenCommon.h"

class Channel;
class FEM_ObjectBroker;
class Graph;
class ID;
class Matrix;
class SparsePythonCompressedEigenSolver;

class SparsePythonCompressedEigenSOE : public EigenSOE
{
  public:
    using StorageScheme = SparsePythonEigenStorageScheme;
    using MatrixStatus = SparsePythonEigenMatrixStatus;

    SparsePythonCompressedEigenSOE(SparsePythonCompressedEigenSolver &theSolver,
                                    StorageScheme scheme = StorageScheme::CSR);
    explicit SparsePythonCompressedEigenSOE(StorageScheme scheme = StorageScheme::CSR);
    ~SparsePythonCompressedEigenSOE() override;

    StorageScheme getStorageScheme() const { return storageScheme; }

    int getNumEqn(void) const;
    int setSize(Graph &theGraph) override;

    int addA(const Matrix &, const ID &, double fact = 1.0) override;
    int addM(const Matrix &, const ID &, double fact = 1.0) override;

    void zeroA(void) override;
    void zeroM(void) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    int setPythonSolver(SparsePythonCompressedEigenSolver &newSolver);

    MatrixStatus getMatrixStatus(void) const;

    // Accessors for shared sparsity pattern
    const std::vector<int> &getIndexPtr(void) const { return indexPtr; }
    const std::vector<int> &getIndices(void) const { return indices; }

    // Accessors for K matrix (stiffness) values
    const std::vector<double> &getKValues(void) const { return kValues; }
    std::vector<double> &getKValues(void) { return kValues; }

    // Accessors for M matrix (mass) values
    const std::vector<double> &getMValues(void) const { return mValues; }
    std::vector<double> &getMValues(void) { return mValues; }

    std::size_t getNNZ(void) const { return kValues.size(); }

  protected:
    void updateMatrixStatus();
    void buildCSRPattern(Graph &theGraph);
    void buildCSCPattern(Graph &theGraph);

    std::vector<int> &accessIndexPtr() { return indexPtr; }
    std::vector<int> &accessIndices() { return indices; }

  private:
    // Shared sparsity pattern (same for K and M)
    std::vector<int> indexPtr;
    std::vector<int> indices;

    // K matrix (stiffness) values
    std::vector<double> kValues;

    // M matrix (mass) values
    std::vector<double> mValues;

    MatrixStatus matrixStatus;
    StorageScheme storageScheme;

    friend class SparsePythonCompressedEigenSolver;
};

#endif /* SparsePythonCompressedEigenSOE_h */

