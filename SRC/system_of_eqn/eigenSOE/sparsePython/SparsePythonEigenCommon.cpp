#include "SparsePythonEigenCommon.h"

// Define static const members for SparsePythonEigenStorageScheme
const SparsePythonEigenStorageScheme SparsePythonEigenStorageScheme::CSR(SparsePythonEigenStorageScheme::Value::CSR);
const SparsePythonEigenStorageScheme SparsePythonEigenStorageScheme::CSC(SparsePythonEigenStorageScheme::Value::CSC);
const SparsePythonEigenStorageScheme SparsePythonEigenStorageScheme::COO(SparsePythonEigenStorageScheme::Value::COO);

const char* SparsePythonEigenStorageScheme::to_str() const {
    switch (value_) {
        case Value::CSR:
            return "CSR";
        case Value::CSC:
            return "CSC";
        case Value::COO:
            return "COO";
        default:
            return "CSR";  // fallback
    }
}

// Define static const members for SparsePythonEigenMatrixStatus
const SparsePythonEigenMatrixStatus SparsePythonEigenMatrixStatus::UNCHANGED(SparsePythonEigenMatrixStatus::Value::UNCHANGED);
const SparsePythonEigenMatrixStatus SparsePythonEigenMatrixStatus::COEFFICIENTS_CHANGED(SparsePythonEigenMatrixStatus::Value::COEFFICIENTS_CHANGED);
const SparsePythonEigenMatrixStatus SparsePythonEigenMatrixStatus::STRUCTURE_CHANGED(SparsePythonEigenMatrixStatus::Value::STRUCTURE_CHANGED);

const char* SparsePythonEigenMatrixStatus::to_str() const {
    switch (value_) {
        case Value::UNCHANGED:
            return "UNCHANGED";
        case Value::COEFFICIENTS_CHANGED:
            return "COEFFICIENTS_CHANGED";
        case Value::STRUCTURE_CHANGED:
            return "STRUCTURE_CHANGED";
        default:
            return "UNCHANGED";  // fallback
    }
}

