#ifndef SparsePythonCompressedEigenSolver_h
#define SparsePythonCompressedEigenSolver_h

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
// SparsePythonCompressedEigenSolver. The solver delegates the eigen solve to a Python
// callable that operates directly on the sparse buffers owned by
// SparsePythonCompressedEigenSOE.

#include "SparsePythonCompressedEigenSOE.h"
#include "SparsePythonEigenCommon.h"

#include <EigenSolver.h>
#include <Vector.h>
#include <string>

struct _object;
typedef _object PyObject;

class Channel;
class FEM_ObjectBroker;
class SparsePythonCompressedEigenSolver : public EigenSolver
{
  public:
    SparsePythonCompressedEigenSolver();
    explicit SparsePythonCompressedEigenSolver(PyObject *callable,
                                                const char *methodName = "solve");
    ~SparsePythonCompressedEigenSolver() override;

    int solve(int numModes, bool generalized, bool findSmallest = true) override;
    int setSize(void) override;

    virtual const Vector &getEigenvector(int mode) override;
    virtual double getEigenvalue(int mode) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    int setEigenSOE(SparsePythonCompressedEigenSOE &theSOE);
    int setPythonCallable(PyObject *callable, const char *methodName = "solve");

    PyObject *getPythonCallable() const { return solverObject; }
    const std::string &getMethodName() const { return methodName; }

  private:
    int callPythonSolver(int numModes, bool generalized, bool findSmallest);

    SparsePythonCompressedEigenSOE *theSOE;
    PyObject *solverObject;
    std::string methodName;

    // Storage for eigenvalues and eigenvectors
    // Flat buffers for zero-copy transfer from Python
    std::vector<double> eigenvaluesBuffer;      // Flat buffer for eigenvalues
    std::vector<double> eigenvectorsBuffer;     // Flat buffer: numModes * numEqn (row-major: mode-major)
    std::vector<Vector> eigenvectorWrappers;    // Vector wrappers pointing to slices of eigenvectorsBuffer
    int numComputedModes;
};

#endif /* SparsePythonCompressedEigenSolver_h */

