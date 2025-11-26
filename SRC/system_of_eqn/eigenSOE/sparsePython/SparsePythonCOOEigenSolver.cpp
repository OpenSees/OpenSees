#include "SparsePythonCOOEigenSolver.h"

#include "SparsePythonCOOEigenSOE.h"
#include "SparsePythonEigenCommon.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <Vector.h>
#include <elementAPI.h>
#include <Python.h>

namespace {

struct PyObjectHolder {
    PyObjectHolder() = default;
    explicit PyObjectHolder(PyObject *obj) : ptr(obj) {}
    ~PyObjectHolder() { Py_XDECREF(ptr); }
    PyObject *get() const { return ptr; }
    PyObject *release() {
        PyObject *tmp = ptr;
        ptr = nullptr;
        return tmp;
    }
    void reset(PyObject *obj = nullptr) {
        if (ptr == obj) {
            return;
        }
        Py_XDECREF(ptr);
        ptr = obj;
    }

  private:
    PyObject *ptr{nullptr};
};

} // namespace

SparsePythonCOOEigenSolver::SparsePythonCOOEigenSolver()
    : EigenSolver(EigenSOLVER_TAGS_SparsePythonCOOEigenSolver),
      theSOE(nullptr),
      solverObject(nullptr),
      methodName("solve"),
      numComputedModes(0)
{
}

SparsePythonCOOEigenSolver::SparsePythonCOOEigenSolver(PyObject *callable, const char *method)
    : EigenSolver(EigenSOLVER_TAGS_SparsePythonCOOEigenSolver),
      theSOE(nullptr),
      solverObject(nullptr),
      methodName(method != nullptr ? method : "solve"),
      numComputedModes(0)
{
    setPythonCallable(callable, methodName.c_str());
}

SparsePythonCOOEigenSolver::~SparsePythonCOOEigenSolver()
{
    if (solverObject != nullptr) {
        if (Py_IsInitialized()) {
            // Acquire GIL before decrementing reference count
            // This is critical: Py_XDECREF must be called with GIL held
            PyGILState_STATE gil_state = PyGILState_Ensure();
            Py_XDECREF(solverObject);
            PyGILState_Release(gil_state);
        }
        // If Python is not initialized, we cannot safely call Py_XDECREF
        // The object will be cleaned up when Python interpreter shuts down
        solverObject = nullptr;
    }
}

int
SparsePythonCOOEigenSolver::setEigenSOE(SparsePythonCOOEigenSOE &theEigenSOE)
{
    theSOE = &theEigenSOE;
    return 0;
}

int
SparsePythonCOOEigenSolver::setPythonCallable(PyObject *callable, const char *method)
{
    if (method != nullptr) {
        methodName = method;
    }

    if (callable == nullptr) {
        Py_XDECREF(solverObject);
        solverObject = nullptr;
        return 0;
    }

    Py_XINCREF(callable);
    Py_XDECREF(solverObject);
    solverObject = callable;
    return 0;
}

int
SparsePythonCOOEigenSolver::solve(int numModes, bool generalized, bool findSmallest)
{
    if (theSOE == nullptr) {
        opserr << "SparsePythonCOOEigenSolver::solve - no associated SOE\n";
        return -1;
    }
    if (solverObject == nullptr) {
        opserr << "SparsePythonCOOEigenSolver::solve - no Python callable bound\n";
        return -2;
    }
    if (theSOE->getNumEqn() == 0) {
        return 0;
    }

    return callPythonSolver(numModes, generalized, findSmallest);
}

int
SparsePythonCOOEigenSolver::setSize(void)
{
    if (theSOE == nullptr) {
        opserr << "SparsePythonCOOEigenSolver::setSize - no associated SOE\n";
        return -1;
    }
    return 0;
}

const Vector &
SparsePythonCOOEigenSolver::getEigenvector(int mode)
{
    if (mode <= 0 || mode > numComputedModes) {
        opserr << "SparsePythonCOOEigenSolver::getEigenvector - mode " << mode
               << " out of range (1 - " << numComputedModes << ")\n";
        static Vector dummy;
        return dummy;
    }
    return eigenvectors[static_cast<std::size_t>(mode - 1)];
}

double
SparsePythonCOOEigenSolver::getEigenvalue(int mode)
{
    if (mode <= 0 || mode > numComputedModes) {
        opserr << "SparsePythonCOOEigenSolver::getEigenvalue - mode " << mode
               << " out of range (1 - " << numComputedModes << ")\n";
        return 0.0;
    }
    return eigenvalues[static_cast<std::size_t>(mode - 1)];
}

int
SparsePythonCOOEigenSolver::sendSelf(int, Channel &)
{
    return 0;
}

int
SparsePythonCOOEigenSolver::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

int
SparsePythonCOOEigenSolver::callPythonSolver(int numModes, bool generalized, bool findSmallest)
{
    const std::vector<int> &rowIdx = theSOE->getRowIndices();
    const std::vector<int> &colIdx = theSOE->getColIndices();
    std::vector<double> &kValues = theSOE->getKValues();
    std::vector<double> &mValues = theSOE->getMValues();

    if (rowIdx.size() != colIdx.size() || rowIdx.size() != kValues.size() || rowIdx.size() != mValues.size()) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - inconsistent COO data\n";
        return -3;
    }

    PyGILState_STATE gilState = PyGILState_Ensure();

    PyObjectHolder method(PyObject_GetAttrString(solverObject, methodName.c_str()));
    if (method.get() == nullptr || !PyCallable_Check(method.get())) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - attribute '"
               << methodName.c_str() << "' not callable\n";
        PyGILState_Release(gilState);
        return -4;
    }

    const Py_ssize_t nnzBytes = static_cast<Py_ssize_t>(rowIdx.size() * sizeof(int));
    const Py_ssize_t kValuesBytes = static_cast<Py_ssize_t>(kValues.size() * sizeof(double));
    const Py_ssize_t mValuesBytes = static_cast<Py_ssize_t>(mValues.size() * sizeof(double));

    PyObjectHolder pyRows(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(const_cast<int *>(rowIdx.data())),
        nnzBytes,
        PyBUF_READ));

    PyObjectHolder pyCols(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(const_cast<int *>(colIdx.data())),
        nnzBytes,
        PyBUF_READ));

    PyObjectHolder pyKValues(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(kValues.data()),
        kValuesBytes,
        PyBUF_READ));

    PyObjectHolder pyMValues(PyMemoryView_FromMemory(
        reinterpret_cast<char *>(mValues.data()),
        mValuesBytes,
        PyBUF_READ));

    if (pyRows.get() == nullptr || pyCols.get() == nullptr ||
        pyKValues.get() == nullptr || pyMValues.get() == nullptr) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - failed to create memory views\n";
        PyGILState_Release(gilState);
        return -5;
    }

    PyObjectHolder matrixStatus(PyUnicode_FromString(theSOE->getMatrixStatus().to_str()));
    PyObjectHolder storageScheme(PyUnicode_FromString(theSOE->getStorageScheme().to_str()));

    PyObjectHolder kwargs(PyDict_New());
    const char *rowIndicesKeyword = "row_indices";
    const char *colIndicesKeyword = "col_indices";
    const char *kValuesKeyword = "k_values";
    const char *mValuesKeyword = "m_values";
    const char *matrixStatusKeyword = "matrix_status";
    const char *storageSchemeKeyword = "storage_scheme";
    const char *numEqnKeyword = "num_eqn";
    const char *nnzKeyword = "nnz";
    const char *numModesKeyword = "num_modes";
    const char *generalizedKeyword = "generalized";
    const char *findSmallestKeyword = "find_smallest";

    if (kwargs.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), rowIndicesKeyword, pyRows.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), colIndicesKeyword, pyCols.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), kValuesKeyword, pyKValues.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), mValuesKeyword, pyMValues.get()) != 0 ||
        matrixStatus.get() == nullptr ||
        storageScheme.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), matrixStatusKeyword, matrixStatus.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), storageSchemeKeyword, storageScheme.get()) != 0) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - failed to populate kwargs\n";
        PyGILState_Release(gilState);
        return -6;
    }

    PyObjectHolder numEqn(PyLong_FromLong(theSOE->getNumEqn()));
    PyObjectHolder nnz(PyLong_FromSize_t(theSOE->getNNZ()));
    PyObjectHolder pyNumModes(PyLong_FromLong(numModes));
    PyObjectHolder pyGeneralized(PyBool_FromLong(generalized ? 1 : 0));
    PyObjectHolder pyFindSmallest(PyBool_FromLong(findSmallest ? 1 : 0));

    if (numEqn.get() == nullptr || nnz.get() == nullptr ||
        pyNumModes.get() == nullptr || pyGeneralized.get() == nullptr || pyFindSmallest.get() == nullptr ||
        PyDict_SetItemString(kwargs.get(), numEqnKeyword, numEqn.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), nnzKeyword, nnz.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), numModesKeyword, pyNumModes.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), generalizedKeyword, pyGeneralized.get()) != 0 ||
        PyDict_SetItemString(kwargs.get(), findSmallestKeyword, pyFindSmallest.get()) != 0) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - failed to set metadata kwargs\n";
        PyGILState_Release(gilState);
        return -7;
    }

    PyObject *args = PyTuple_New(0);
    if (args == nullptr) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - failed to create args tuple\n";
        PyGILState_Release(gilState);
        return -8;
    }
    PyObjectHolder argsHolder(args);

    PyObjectHolder result(PyObject_Call(method.get(), args, kwargs.get()));
    PyGILState_Release(gilState);

    if (result.get() == nullptr) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - Python call failed\n";
        PyErr_Print();
        return -9;
    }

    if (!PyTuple_Check(result.get()) || PyTuple_Size(result.get()) != 2) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - result must be tuple (eigenvalues, eigenvectors)\n";
        return -10;
    }

    PyObject *pyEigenvalues = PyTuple_GetItem(result.get(), 0);
    PyObject *pyEigenvectors = PyTuple_GetItem(result.get(), 1);

    if (pyEigenvalues == nullptr || pyEigenvectors == nullptr) {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - failed to extract result components\n";
        return -11;
    }

    Py_ssize_t numEigen = 0;
    if (PyList_Check(pyEigenvalues)) {
        numEigen = PyList_Size(pyEigenvalues);
        eigenvalues.resize(static_cast<std::size_t>(numEigen));
        for (Py_ssize_t i = 0; i < numEigen; ++i) {
            PyObject *item = PyList_GetItem(pyEigenvalues, i);
            if (item == nullptr || !PyFloat_Check(item)) {
                opserr << "SparsePythonCOOEigenSolver::callPythonSolver - invalid eigenvalue at index " << i << "\n";
                return -12;
            }
            eigenvalues[static_cast<std::size_t>(i)] = PyFloat_AsDouble(item);
        }
    } else {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - eigenvalues must be a list\n";
        return -13;
    }

    if (PyList_Check(pyEigenvectors)) {
        Py_ssize_t numVecs = PyList_Size(pyEigenvectors);
        if (numVecs != numEigen) {
            opserr << "SparsePythonCOOEigenSolver::callPythonSolver - number of eigenvectors ("
                   << numVecs << ") != number of eigenvalues (" << numEigen << ")\n";
            return -14;
        }
        eigenvectors.clear();
        eigenvectors.reserve(static_cast<std::size_t>(numVecs));
        const int size = theSOE->getNumEqn();
        for (Py_ssize_t i = 0; i < numVecs; ++i) {
            PyObject *vec = PyList_GetItem(pyEigenvectors, i);
            if (vec == nullptr || !PyList_Check(vec)) {
                opserr << "SparsePythonCOOEigenSolver::callPythonSolver - invalid eigenvector at index " << i << "\n";
                return -15;
            }
            Py_ssize_t vecSize = PyList_Size(vec);
            if (vecSize != size) {
                opserr << "SparsePythonCOOEigenSolver::callPythonSolver - eigenvector " << i
                       << " size (" << vecSize << ") != system size (" << size << ")\n";
                return -16;
            }
            Vector eigenvec(size);
            for (Py_ssize_t j = 0; j < vecSize; ++j) {
                PyObject *item = PyList_GetItem(vec, j);
                if (item == nullptr || !PyFloat_Check(item)) {
                    opserr << "SparsePythonCOOEigenSolver::callPythonSolver - invalid eigenvector element at ("
                           << i << ", " << j << ")\n";
                    return -17;
                }
                eigenvec(static_cast<int>(j)) = PyFloat_AsDouble(item);
            }
            eigenvectors.push_back(eigenvec);
        }
    } else {
        opserr << "SparsePythonCOOEigenSolver::callPythonSolver - eigenvectors must be a list of lists\n";
        return -18;
    }

    numComputedModes = static_cast<int>(numEigen);
    return 0;
}


