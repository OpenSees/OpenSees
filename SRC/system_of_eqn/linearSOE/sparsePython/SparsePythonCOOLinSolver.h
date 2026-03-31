#ifndef SparsePythonCOOLinSolver_h
#define SparsePythonCOOLinSolver_h

#include "SparsePythonCOOLinSOE.h"
#include "SparsePythonCommon.h"

#include <LinearSOESolver.h>
#include <string>

struct _object;
typedef _object PyObject;

class Channel;
class FEM_ObjectBroker;

class SparsePythonCOOLinSolver : public LinearSOESolver
{
  public:
    SparsePythonCOOLinSolver();
    explicit SparsePythonCOOLinSolver(PyObject *callable, const char *methodName = "solve");
    ~SparsePythonCOOLinSolver() override;

    int solve(void) override;
    int setSize(void) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    int setLinearSOE(SparsePythonCOOLinSOE &theSOE);
    int setPythonCallable(PyObject *callable, const char *methodName = "solve");
    void setWritableFlags(const SparsePythonWritableFlags &flags) { writableFlags = flags; }

    PyObject *getPythonCallable() const { return solverObject; }
    const std::string &getMethodName() const { return methodName; }
    const SparsePythonWritableFlags &getWritableFlags() const { return writableFlags; }

  private:
    int callPythonSolver();

    SparsePythonCOOLinSOE *theSOE;
    PyObject *solverObject;
    std::string methodName;
    SparsePythonWritableFlags writableFlags;
};

#endif /* SparsePythonCOOLinSolver_h */
