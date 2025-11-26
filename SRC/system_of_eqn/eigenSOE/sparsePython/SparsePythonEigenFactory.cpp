#include "SparsePythonEigenFactory.h"

#include "SparsePythonCompressedEigenSOE.h"
#include "SparsePythonCompressedEigenSolver.h"
#include "SparsePythonCOOEigenSOE.h"
#include "SparsePythonCOOEigenSolver.h"
#include "SparsePythonEigenCommon.h"

#include <OPS_Globals.h>
#include <elementAPI.h>

#include <Python.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <string>

namespace {

bool
IsSparsePythonEigenType(const char *type)
{
    return strcmp(type, "PythonSparse") == 0 ||
           strcmp(type, "PythonCompressedSparseEigen") == 0;
}

std::string
ToLower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool
ParseSchemeToken(const char *token, SparsePythonEigenStorageScheme &scheme)
{
    if (token == nullptr) {
        return false;
    }
    std::string normalized = ToLower(token);
    if (normalized == "csr" || normalized == "row") {
        scheme = SparsePythonEigenStorageScheme::CSR;
        return true;
    }
    if (normalized == "csc" || normalized == "column" || normalized == "col") {
        scheme = SparsePythonEigenStorageScheme::CSC;
        return true;
    }
    if (normalized == "coo") {
        scheme = SparsePythonEigenStorageScheme::COO;
        return true;
    }
    return false;
}

} // namespace

void *
OPS_SparsePythonEigenSolver()
{
    const char *expectedSyntax = "eigen 'PythonSparse' numModes {'solver': SolverObject, 'scheme': 'CSR'|'CSC'|'COO'}";

    // Note: This function is called AFTER numModes has been read in PythonAnalysisBuilder.cpp
    // So we don't need to reset or read the type string again - just get the dict argument

    if (!Py_IsInitialized()) {
        Py_Initialize();
    }

    // Get the dictionary argument (should be the next argument after numModes)
    void *dictPtr = OPS_GetVoidPtr();
    if (dictPtr == nullptr) {
        opserr << "WARNING: eigen PythonSparse - requires a dictionary argument as third parameter" << endln;
        opserr << "Expected syntax: " << expectedSyntax << endln;
        return nullptr;
    }

    PyObject *dict = static_cast<PyObject *>(dictPtr);
    if (!PyDict_Check(dict)) {
        opserr << "WARNING: eigen PythonSparse - third argument must be a dictionary" << endln;
        opserr << "Expected syntax: " << expectedSyntax << endln;
        return nullptr;
    }

    // Extract 'solver' (required)
    PyObject *solverObj = PyDict_GetItemString(dict, "solver");
    if (solverObj == nullptr) {
        opserr << "WARNING: eigen PythonSparse - dictionary must contain 'solver' key" << endln;
        opserr << "Expected syntax: " << expectedSyntax << endln;
        return nullptr;
    }

    // Extract 'scheme' (optional, defaults to CSR)
    SparsePythonEigenStorageScheme scheme = SparsePythonEigenStorageScheme::CSR;
    PyObject *schemeObj = PyDict_GetItemString(dict, "scheme");
    if (schemeObj != nullptr) {
        if (!PyUnicode_Check(schemeObj)) {
            opserr << "WARNING: eigen PythonSparse - 'scheme' must be a string" << endln;
            return nullptr;
        }
        const char *schemeStr = PyUnicode_AsUTF8(schemeObj);
        if (schemeStr == nullptr || !ParseSchemeToken(schemeStr, scheme)) {
            opserr << "WARNING: eigen PythonSparse - unknown storage scheme '"
                   << (schemeStr != nullptr ? schemeStr : "null") << "' (expected CSR, CSC, or COO)" << endln;
            return nullptr;
        }
    }

    if (scheme == SparsePythonEigenStorageScheme::COO) {
        SparsePythonCOOEigenSolver *solver = new SparsePythonCOOEigenSolver();
        if (solver->setPythonCallable(solverObj, "solve") != 0) {
            opserr << "WARNING: eigen PythonSparse - failed to set Python callable" << endln;
            delete solver;
            return nullptr;
        }
        SparsePythonCOOEigenSOE *soe = new SparsePythonCOOEigenSOE(*solver);
        return static_cast<void *>(soe);
    }

    SparsePythonCompressedEigenSolver *solver = new SparsePythonCompressedEigenSolver();
    if (solver->setPythonCallable(solverObj, "solve") != 0) {
        opserr << "WARNING: eigen PythonSparse - failed to set Python callable" << endln;
        delete solver;
        return nullptr;
    }

    SparsePythonCompressedEigenSOE *soe = new SparsePythonCompressedEigenSOE(*solver, scheme);
    return static_cast<void *>(soe);
}


