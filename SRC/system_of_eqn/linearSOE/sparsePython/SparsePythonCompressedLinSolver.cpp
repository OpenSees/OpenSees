#include "SparsePythonCompressedLinSolver.h"

#include "SparsePythonCompressedLinSOE.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <cstring>
#include <elementAPI.h>
#include <Python.h>
#include "SparsePythonCommon.h"

SparsePythonCompressedLinSolver::SparsePythonCompressedLinSolver()
    : LinearSOESolver(SOLVER_TAGS_SparsePythonCompressedLinSolver),
      theSOE(nullptr),
      solverObject(nullptr),
      methodName("solve")
{
}

SparsePythonCompressedLinSolver::SparsePythonCompressedLinSolver(PyObject *callable,
                                                                 const char *method)
    : LinearSOESolver(SOLVER_TAGS_SparsePythonCompressedLinSolver),
      theSOE(nullptr),
      solverObject(nullptr),
      methodName(method != nullptr ? method : "solve")
{
    setPythonCallable(callable, methodName.c_str());
}

SparsePythonCompressedLinSolver::~SparsePythonCompressedLinSolver()
{
    // NOTE: This is generally unsafe to call Python C-API from destructors,
    // as destructors may run after Py_Finalize(). This will be addressed later.
    // For now, we assume the user will call ops.wipe() at the end of their script
    // to ensure resources are released before exiting Python, making it safe to
    // decrement refcounts in destructors.
    if (solverObject != nullptr && Py_IsInitialized()) {
        PyGILState_STATE gilState = PyGILState_Ensure();
        Py_XDECREF(solverObject);
        PyGILState_Release(gilState);
    }
}

int
SparsePythonCompressedLinSolver::setLinearSOE(SparsePythonCompressedLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    return 0;
}

int
SparsePythonCompressedLinSolver::setPythonCallable(PyObject *callable, const char *method)
{
    if (method != nullptr) {
        methodName = method;
    }

    // Acquire GIL before any Python API calls
    PyGILState_STATE gilState = PyGILState_Ensure();

    if (callable == nullptr) {
        Py_XDECREF(solverObject);
        solverObject = nullptr;
        PyGILState_Release(gilState);
        return 0;
    }

    Py_XINCREF(callable);
    Py_XDECREF(solverObject);
    solverObject = callable;

    PyGILState_Release(gilState);
    return 0;
}


int
SparsePythonCompressedLinSolver::solve(void)
{
    if (theSOE == nullptr) {
        opserr << "SparsePythonCompressedLinSolver::solve - no associated SOE\n";
        return -1;
    }
    if (solverObject == nullptr) {
        opserr << "SparsePythonCompressedLinSolver::solve - no Python callable bound\n";
        return -2;
    }
    if (theSOE->getNumEqn() == 0) {
        return 0;
    }

    return callPythonSolver();
}

int
SparsePythonCompressedLinSolver::setSize(void)
{
    if (theSOE == nullptr) {
        opserr << "SparsePythonCompressedLinSolver::setSize - no associated SOE\n";
        return -1;
    }
    return 0;
}

int
SparsePythonCompressedLinSolver::sendSelf(int, Channel &)
{
    return 0;
}

int
SparsePythonCompressedLinSolver::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

int
SparsePythonCompressedLinSolver::callPythonSolver()
{
    // Pull the SOE storage; indexPtr/indices are immutable views, values/rhs/sol are writable.
    const std::vector<int> &indexPtr = theSOE->getIndexPtr();
    const std::vector<int> &indices = theSOE->getIndices();
    std::vector<double> &values = theSOE->getValues();
    std::vector<double> &rhs = theSOE->getRHS();
    std::vector<double> &sol = theSOE->getSolution();

    
    if (indexPtr.size() < 2) {
        opserr << "SparsePythonCompressedLinSolver::callPythonSolver - invalid sparse structure\n";
        return -3;
    }

    // Ensure we hold the GIL while touching Python APIs.
    PyGILState_STATE gilState = PyGILState_Ensure();

    // Lookup the bound solve method; bail if missing or not callable.
    PyObjectHolder method(PyObject_GetAttrString(solverObject, methodName.c_str()));
    if (method.get() == nullptr || !PyCallable_Check(method.get())) {
        opserr << "SparsePythonCompressedLinSolver::callPythonSolver - attribute '"
               << methodName.c_str() << "' not callable\n";
        PyGILState_Release(gilState);
        return -4;
    }

    // Build memoryviews that hand the raw SOE buffers to Python without copying.
    const Py_ssize_t ptrBytes = static_cast<Py_ssize_t>(indexPtr.size() * sizeof(int));
    const Py_ssize_t indBytes = static_cast<Py_ssize_t>(indices.size() * sizeof(int));
    const Py_ssize_t valuesBytes = static_cast<Py_ssize_t>(values.size() * sizeof(double));
    const Py_ssize_t rhsBytes = static_cast<Py_ssize_t>(rhs.size() * sizeof(double));
    const Py_ssize_t solBytes = static_cast<Py_ssize_t>(sol.size() * sizeof(double));

    PyObjectHolder pyPtr(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(const_cast<int *>(indexPtr.data())),
        ptrBytes,
        PyBUF_READ));

    PyObjectHolder pyInd(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(const_cast<int *>(indices.data())),
        indBytes,
        PyBUF_READ));

    PyObjectHolder pyValues(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(values.data()),
        valuesBytes,
        writableFlags.values ? PyBUF_WRITE : PyBUF_READ));

    PyObjectHolder pyRhs(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(rhs.data()),
        rhsBytes,
        writableFlags.rhs ? PyBUF_WRITE : PyBUF_READ));

    PyObjectHolder pySol(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(sol.data()),
        solBytes,
        PyBUF_WRITE));

    if (pyPtr.get() == nullptr || pyInd.get() == nullptr ||
        pyValues.get() == nullptr || pyRhs.get() == nullptr || pySol.get() == nullptr) {
        opserr << "SparsePythonCompressedLinSolver::callPythonSolver - failed to create memory views\n";
        PyGILState_Release(gilState);
        return -5;
    }

    // Scalar metadata handed to the Python solver (as all-caps strings).
    PyObjectHolder matrixStatus(
        PyUnicode_FromString(theSOE->getMatrixStatus().to_str()));
    PyObjectHolder storageScheme(
        PyUnicode_FromString(theSOE->getStorageScheme().to_str()));

    PyObjectHolder kwargs(PyDict_New());
    const char *ptrKeyword = "index_ptr";
    const char *indicesKeyword = "indices";
    const char *valuesKeyword = "values";
    const char *rhsKeyword = "rhs";
    const char *xKeyword = "x";
    const char *matrixStatusKeyword = "matrix_status";
    const char *storageSchemeKeyword = "storage_scheme";
    const char *numEqnKeyword = "num_eqn";
    const char *nnzKeyword = "nnz";

    if (kwargs.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), ptrKeyword, pyPtr.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), indicesKeyword, pyInd.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), valuesKeyword, pyValues.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), rhsKeyword, pyRhs.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), xKeyword, pySol.get()) != 0 ||
        matrixStatus.get() == nullptr ||
        storageScheme.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), matrixStatusKeyword, matrixStatus.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), storageSchemeKeyword, storageScheme.get()) != 0) {
        opserr << "SparsePythonCompressedLinSolver::callPythonSolver - failed to populate kwargs\n";
        PyGILState_Release(gilState);
        return -6;
    }

    // Add additional metadata expected by the solver (num equations, number of nonzeros).
    PyObjectHolder numEqn(PyLong_FromLong(theSOE->getNumEqn()));
    PyObjectHolder nnz(PyLong_FromSize_t(theSOE->getNNZ()));

    if (numEqn.get() == nullptr || nnz.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), numEqnKeyword, numEqn.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), nnzKeyword, nnz.get()) != 0) {
        opserr << "SparsePythonCompressedLinSolver::callPythonSolver - failed to set metadata kwargs\n";
        PyGILState_Release(gilState);
        return -7;
    }

    // The Python callable expects no positional args; we pass everything via kwargs.
    PyObject *args = PyTuple_New(0);
    if (args == nullptr) {
        opserr << "SparsePythonCompressedLinSolver::callPythonSolver - failed to create args tuple\n";
        PyGILState_Release(gilState);
        return -8;
    }
    PyObjectHolder argsHolder(args);

    // PyObject_Call must be made with GIL held. CuPy/NumPy will release the GIL
    // internally during their native GPU/CPU computations, so we don't need to
    // do anything special here.
    PyObjectHolder result(PyObject_Call(method.get(), argsHolder.get(), kwargs.get()));
    if (result.get() == nullptr) {
        opserr << "SparsePythonCompressedLinSolver::callPythonSolver - Python callable raised an exception\n";
        PyErr_Print();
        PyGILState_Release(gilState);
        return -9;
    }

    // Interpret the return value: None -> success, numeric -> explicit status code.
    int status = 0;
    if (result.get() != Py_None) {
        if (PyLong_Check(result.get())) {
            status = static_cast<int>(PyLong_AsLong(result.get()));
        } else if (PyFloat_Check(result.get())) {
            status = static_cast<int>(PyFloat_AsDouble(result.get()));
        } else {
            status = 0;
        }
    }

    PyGILState_Release(gilState);
    return status;
}
