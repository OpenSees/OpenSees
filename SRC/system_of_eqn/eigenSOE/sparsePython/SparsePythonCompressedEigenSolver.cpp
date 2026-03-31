#include "SparsePythonCompressedEigenSolver.h"

#include "SparsePythonCompressedEigenSOE.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <Vector.h>
#include <cstring>
#include <elementAPI.h>
#include <Python.h>
#include "SparsePythonEigenCommon.h"

SparsePythonCompressedEigenSolver::SparsePythonCompressedEigenSolver()
    : EigenSolver(EigenSOLVER_TAGS_SparsePythonCompressedEigenSolver),
      theSOE(nullptr),
      solverObject(nullptr),
      methodName("solve"),
      numComputedModes(0)
{
}

SparsePythonCompressedEigenSolver::SparsePythonCompressedEigenSolver(PyObject *callable,
                                                   const char *method)
    : EigenSolver(EigenSOLVER_TAGS_SparsePythonCompressedEigenSolver),
      theSOE(nullptr),
      solverObject(nullptr),
      methodName(method != nullptr ? method : "solve"),
      numComputedModes(0)
{
    setPythonCallable(callable, methodName.c_str());
}

SparsePythonCompressedEigenSolver::~SparsePythonCompressedEigenSolver()
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
SparsePythonCompressedEigenSolver::setEigenSOE(SparsePythonCompressedEigenSOE &theEigenSOE)
{
    theSOE = &theEigenSOE;
    return 0;
}

int
SparsePythonCompressedEigenSolver::setPythonCallable(PyObject *callable, const char *method)
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
SparsePythonCompressedEigenSolver::solve(int numModes, bool generalized, bool findSmallest)
{
    if (theSOE == nullptr) {
        opserr << "SparsePythonCompressedEigenSolver::solve - no associated SOE\n";
        return -1;
    }
    if (solverObject == nullptr) {
        opserr << "SparsePythonCompressedEigenSolver::solve - no Python callable bound\n";
        return -2;
    }
    if (theSOE->getNumEqn() == 0) {
        return 0;
    }

    return callPythonSolver(numModes, generalized, findSmallest);
}

int
SparsePythonCompressedEigenSolver::setSize(void)
{
    if (theSOE == nullptr) {
        opserr << "SparsePythonCompressedEigenSolver::setSize - no associated SOE\n";
        return -1;
    }
    return 0;
}

const Vector &
SparsePythonCompressedEigenSolver::getEigenvector(int mode)
{
    // mode is 1-indexed (matches ArpackSolver and analysis code)
    if (mode <= 0 || mode > numComputedModes) {
        opserr << "SparsePythonCompressedEigenSolver::getEigenvector - mode " << mode
               << " out of range (1 - " << numComputedModes << ")\n";
        static Vector dummy;
        return dummy;
    }
    return eigenvectorWrappers[static_cast<std::size_t>(mode - 1)];
}

double
SparsePythonCompressedEigenSolver::getEigenvalue(int mode)
{
    // mode is 1-indexed (matches ArpackSolver and analysis code)
    if (mode <= 0 || mode > numComputedModes) {
        opserr << "SparsePythonCompressedEigenSolver::getEigenvalue - mode " << mode
               << " out of range (1 - " << numComputedModes << ")\n";
        return 0.0;
    }
    return eigenvaluesBuffer[static_cast<std::size_t>(mode - 1)];
}

int
SparsePythonCompressedEigenSolver::sendSelf(int, Channel &)
{
    return 0;
}

int
SparsePythonCompressedEigenSolver::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

int
SparsePythonCompressedEigenSolver::callPythonSolver(int numModes, bool generalized, bool findSmallest)
{
    // Pull the SOE storage (shared sparsity pattern, separate values for K and M)
    const std::vector<int> &indexPtr = theSOE->getIndexPtr();
    const std::vector<int> &indices = theSOE->getIndices();
    std::vector<double> &kValues = theSOE->getKValues();
    std::vector<double> &mValues = theSOE->getMValues();

    if (indexPtr.size() < 2) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - invalid sparse structure\n";
        return -3;
    }

    // Ensure we hold the GIL while touching Python APIs.
    PyGILState_STATE gilState = PyGILState_Ensure();

    // Lookup the bound solve method; bail if missing or not callable.
    PyObjectHolder method(PyObject_GetAttrString(solverObject, methodName.c_str()));
    if (method.get() == nullptr || !PyCallable_Check(method.get())) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - attribute '"
               << methodName.c_str() << "' not callable\n";
        PyGILState_Release(gilState);
        return -4;
    }

    // Build memoryviews for shared sparsity pattern
    const Py_ssize_t ptrBytes = static_cast<Py_ssize_t>(indexPtr.size() * sizeof(int));
    const Py_ssize_t indBytes = static_cast<Py_ssize_t>(indices.size() * sizeof(int));

    PyObjectHolder pyPtr(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(const_cast<int *>(indexPtr.data())),
        ptrBytes,
        PyBUF_READ));

    PyObjectHolder pyInd(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(const_cast<int *>(indices.data())),
        indBytes,
        PyBUF_READ));

    // Build memoryviews for K and M values
    const Py_ssize_t kValuesBytes = static_cast<Py_ssize_t>(kValues.size() * sizeof(double));
    const Py_ssize_t mValuesBytes = static_cast<Py_ssize_t>(mValues.size() * sizeof(double));

    PyObjectHolder pyKValues(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(kValues.data()),
        kValuesBytes,
        PyBUF_READ));

    PyObjectHolder pyMValues(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(mValues.data()),
        mValuesBytes,
        PyBUF_READ));

    if (pyPtr.get() == nullptr || pyInd.get() == nullptr || 
        pyKValues.get() == nullptr || pyMValues.get() == nullptr) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - failed to create memory views\n";
        PyGILState_Release(gilState);
        return -5;
    }

    // Pre-allocate buffers for zero-copy transfer from Python
    const int numEqn = theSOE->getNumEqn();
    eigenvaluesBuffer.resize(static_cast<std::size_t>(numModes));
    eigenvectorsBuffer.resize(static_cast<std::size_t>(numModes) * static_cast<std::size_t>(numEqn));

    // Create writable memoryviews for eigenvalues and eigenvectors buffers
    const Py_ssize_t eigenvaluesBytes = static_cast<Py_ssize_t>(eigenvaluesBuffer.size() * sizeof(double));
    const Py_ssize_t eigenvectorsBytes = static_cast<Py_ssize_t>(eigenvectorsBuffer.size() * sizeof(double));

    PyObjectHolder pyEigenvaluesBuffer(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(eigenvaluesBuffer.data()),
        eigenvaluesBytes,
        PyBUF_WRITE));

    PyObjectHolder pyEigenvectorsBuffer(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(eigenvectorsBuffer.data()),
        eigenvectorsBytes,
        PyBUF_WRITE));

    if (pyEigenvaluesBuffer.get() == nullptr || pyEigenvectorsBuffer.get() == nullptr) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - failed to create writable memoryviews for results\n";
        PyGILState_Release(gilState);
        return -5;
    }

    // Scalar metadata
    PyObjectHolder matrixStatus(
        PyUnicode_FromString(theSOE->getMatrixStatus().to_str()));
    PyObjectHolder storageScheme(
        PyUnicode_FromString(theSOE->getStorageScheme().to_str()));

    PyObjectHolder kwargs(PyDict_New());
    const char *indexPtrKeyword = "index_ptr";
    const char *indicesKeyword = "indices";
    const char *kValuesKeyword = "k_values";
    const char *mValuesKeyword = "m_values";
    const char *matrixStatusKeyword = "matrix_status";
    const char *storageSchemeKeyword = "storage_scheme";
    const char *numEqnKeyword = "num_eqn";
    const char *nnzKeyword = "nnz";
    const char *numModesKeyword = "num_modes";
    const char *generalizedKeyword = "generalized";
    const char *findSmallestKeyword = "find_smallest";
    const char *eigenvaluesKeyword = "eigenvalues";
    const char *eigenvectorsKeyword = "eigenvectors";

    if (kwargs.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), indexPtrKeyword, pyPtr.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), indicesKeyword, pyInd.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), kValuesKeyword, pyKValues.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), mValuesKeyword, pyMValues.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), eigenvaluesKeyword, pyEigenvaluesBuffer.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), eigenvectorsKeyword, pyEigenvectorsBuffer.get()) != 0 ||
        matrixStatus.get() == nullptr ||
        storageScheme.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), matrixStatusKeyword, matrixStatus.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), storageSchemeKeyword, storageScheme.get()) != 0) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - failed to populate kwargs\n";
        PyGILState_Release(gilState);
        return -6;
    }

    // Add additional metadata
    PyObjectHolder numEqnObj(PyLong_FromLong(numEqn));
    PyObjectHolder nnz(PyLong_FromSize_t(theSOE->getNNZ()));
    PyObjectHolder pyNumModes(PyLong_FromLong(numModes));
    PyObjectHolder pyGeneralized(PyBool_FromLong(generalized ? 1 : 0));
    PyObjectHolder pyFindSmallest(PyBool_FromLong(findSmallest ? 1 : 0));

    if (numEqnObj.get() == nullptr || nnz.get() == nullptr ||
        pyNumModes.get() == nullptr || pyGeneralized.get() == nullptr || pyFindSmallest.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), numEqnKeyword, numEqnObj.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), nnzKeyword, nnz.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), numModesKeyword, pyNumModes.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), generalizedKeyword, pyGeneralized.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), findSmallestKeyword, pyFindSmallest.get()) != 0) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - failed to set metadata kwargs\n";
        PyGILState_Release(gilState);
        return -7;
    }

    // Call Python solver - it should write directly to the writable buffers
    PyObject *args = PyTuple_New(0);
    if (args == nullptr) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - failed to create args tuple\n";
        PyGILState_Release(gilState);
        return -8;
    }
    PyObjectHolder argsHolder(args);

    // PyObject_Call must be made with GIL held. CuPy/NumPy will release the GIL
    // internally during their native GPU/CPU computations, so we don't need to
    // do anything special here.
    PyObjectHolder result(PyObject_Call(method.get(), args, kwargs.get()));

    if (result.get() == nullptr) {
        opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - Python call failed\n";
        PyErr_Print();
        PyGILState_Release(gilState);
        return -9;
    }

    // Check return value - Python should write directly to buffers and return None or status code
    int status = 0;
    if (result.get() != Py_None) {
        if (PyLong_Check(result.get())) {
            status = static_cast<int>(PyLong_AsLong(result.get()));
            if (status < 0) {
                PyGILState_Release(gilState);
                numComputedModes = 0;
                return status;
            }
        } else {
            opserr << "SparsePythonCompressedEigenSolver::callPythonSolver - solver returned unexpected type (expected None or integer status)\n";
            PyGILState_Release(gilState);
            return -10;
        }
    }

    // Create Vector wrappers for eigenvectors pointing to the flat buffer
    // Buffer layout: eigenvectorsBuffer[mode * numEqn + eqn] (mode-major, row-major)
    eigenvectorWrappers.clear();
    eigenvectorWrappers.reserve(static_cast<std::size_t>(numModes));
    for (int mode = 0; mode < numModes; ++mode) {
        double *modeData = eigenvectorsBuffer.data() + mode * numEqn;
        eigenvectorWrappers.emplace_back(Vector(modeData, numEqn));
    }

    PyGILState_Release(gilState);
    numComputedModes = numModes;
    return status;
}

