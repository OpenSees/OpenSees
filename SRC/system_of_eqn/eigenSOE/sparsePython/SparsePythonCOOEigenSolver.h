#ifndef SparsePythonCOOEigenSolver_h
#define SparsePythonCOOEigenSolver_h

#include <EigenSolver.h>
#include <Vector.h>

#include <string>
#include <vector>

class SparsePythonCOOEigenSOE;
struct _object;
typedef _object PyObject;

class Channel;
class FEM_ObjectBroker;

class SparsePythonCOOEigenSolver : public EigenSolver
{
  public:
    SparsePythonCOOEigenSolver();
    explicit SparsePythonCOOEigenSolver(PyObject *callable, const char *methodName = "solve");
    ~SparsePythonCOOEigenSolver() override;

    int solve(int numModes, bool generalized, bool findSmallest = true) override;
    int setSize(void) override;

    const Vector &getEigenvector(int mode) override;
    double getEigenvalue(int mode) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    int setEigenSOE(SparsePythonCOOEigenSOE &theSOE);
    int setPythonCallable(PyObject *callable, const char *methodName = "solve");

  private:
    int callPythonSolver(int numModes, bool generalized, bool findSmallest);

    SparsePythonCOOEigenSOE *theSOE;
    PyObject *solverObject;
    std::string methodName;

    std::vector<double> eigenvalues;
    std::vector<Vector> eigenvectors;
    int numComputedModes;
};

#endif /* SparsePythonCOOEigenSolver_h */


