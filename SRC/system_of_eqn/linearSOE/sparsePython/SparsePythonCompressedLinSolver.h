#ifndef SparsePythonCompressedLinSolver_h
#define SparsePythonCompressedLinSolver_h

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
// SparsePythonCompressedLinSolver. The solver delegates the factorization /
// back-substitution to a Python callable that operates directly on the sparse
// buffers owned by SparsePythonCompressedLinSOE.

#include "SparsePythonCompressedLinSOE.h"
#include "SparsePythonCommon.h"

#include <LinearSOESolver.h>
#include <string>

struct _object;
typedef _object PyObject;

class Channel;
class FEM_ObjectBroker;
class SparsePythonCompressedLinSolver : public LinearSOESolver
{
  public:
    SparsePythonCompressedLinSolver();
    explicit SparsePythonCompressedLinSolver(PyObject *callable,
                                             const char *methodName = "solve");
    ~SparsePythonCompressedLinSolver() override;

    int solve(void) override;
    int setSize(void) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    int setLinearSOE(SparsePythonCompressedLinSOE &theSOE);
    int setPythonCallable(PyObject *callable, const char *methodName = "solve");
    void setWritableFlags(const SparsePythonWritableFlags &flags) { writableFlags = flags; }

    PyObject *getPythonCallable() const { return solverObject; }
    const std::string &getMethodName() const { return methodName; }
    const SparsePythonWritableFlags &getWritableFlags() const { return writableFlags; }

  private:
    int callPythonSolver();

    SparsePythonCompressedLinSOE *theSOE;
    PyObject *solverObject;
    std::string methodName;
    SparsePythonWritableFlags writableFlags;
};
#endif /* SparsePythonCompressedLinSolver_h */