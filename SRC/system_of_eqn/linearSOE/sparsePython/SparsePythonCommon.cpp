#include "SparsePythonCommon.h"

// Define static const members for SparsePythonStorageScheme
const SparsePythonStorageScheme SparsePythonStorageScheme::CSR(SparsePythonStorageScheme::Value::CSR);
const SparsePythonStorageScheme SparsePythonStorageScheme::CSC(SparsePythonStorageScheme::Value::CSC);
const SparsePythonStorageScheme SparsePythonStorageScheme::COO(SparsePythonStorageScheme::Value::COO);

const char* SparsePythonStorageScheme::to_str() const {
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

// Define static const members for SparsePythonMatrixStatus
const SparsePythonMatrixStatus SparsePythonMatrixStatus::UNCHANGED(SparsePythonMatrixStatus::Value::UNCHANGED);
const SparsePythonMatrixStatus SparsePythonMatrixStatus::COEFFICIENTS_CHANGED(SparsePythonMatrixStatus::Value::COEFFICIENTS_CHANGED);
const SparsePythonMatrixStatus SparsePythonMatrixStatus::STRUCTURE_CHANGED(SparsePythonMatrixStatus::Value::STRUCTURE_CHANGED);

const char* SparsePythonMatrixStatus::to_str() const {
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

