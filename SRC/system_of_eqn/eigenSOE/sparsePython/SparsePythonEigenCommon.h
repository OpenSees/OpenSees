#ifndef SparsePythonEigenCommon_h
#define SparsePythonEigenCommon_h

/**
 * Storage scheme class: acts like an enum with string conversion.
 */
class SparsePythonEigenStorageScheme {
private:
    enum class Value : int {
        CSR = 0,
        CSC = 1,
        COO = 2
    };
    Value value_;

public:
    // Static const members for easy access (defined in .cpp)
    static const SparsePythonEigenStorageScheme CSR;
    static const SparsePythonEigenStorageScheme CSC;
    static const SparsePythonEigenStorageScheme COO;

    // Constructors
    SparsePythonEigenStorageScheme() : value_(Value::CSR) {}
    explicit SparsePythonEigenStorageScheme(Value v) : value_(v) {}
    SparsePythonEigenStorageScheme(int v) : value_(static_cast<Value>(v)) {}

    // Implicit conversion operators for comparisons and assignments
    operator Value() const { return value_; }
    operator int() const { return static_cast<int>(value_); }

    // Comparison operators
    bool operator==(const SparsePythonEigenStorageScheme& other) const {
        return value_ == other.value_;
    }
    bool operator!=(const SparsePythonEigenStorageScheme& other) const {
        return value_ != other.value_;
    }

    // String conversion for Python
    const char* to_str() const;
};

/**
 * Matrix status class: tracks changes for potential factorization reuse.
 * Useful for shift-and-invert methods (e.g., ARPACK mode 3) that factorize
 * (K - sigma*M) and can reuse symbolic factorization if only coefficients change.
 */
class SparsePythonEigenMatrixStatus {
private:
    enum class Value : int {
        UNCHANGED = 0,
        COEFFICIENTS_CHANGED = 1,
        STRUCTURE_CHANGED = 2
    };
    Value value_;

public:
    // Static const members for easy access (defined in .cpp)
    static const SparsePythonEigenMatrixStatus UNCHANGED;
    static const SparsePythonEigenMatrixStatus COEFFICIENTS_CHANGED;
    static const SparsePythonEigenMatrixStatus STRUCTURE_CHANGED;

    // Constructors
    SparsePythonEigenMatrixStatus() : value_(Value::UNCHANGED) {}
    explicit SparsePythonEigenMatrixStatus(Value v) : value_(v) {}
    SparsePythonEigenMatrixStatus(int v) : value_(static_cast<Value>(v)) {}

    // Implicit conversion operators for comparisons and assignments
    operator Value() const { return value_; }
    operator int() const { return static_cast<int>(value_); }

    // Comparison operators
    bool operator==(const SparsePythonEigenMatrixStatus& other) const {
        return value_ == other.value_;
    }
    bool operator!=(const SparsePythonEigenMatrixStatus& other) const {
        return value_ != other.value_;
    }

    // String conversion for Python
    const char* to_str() const;
};

#endif /* SparsePythonEigenCommon_h */

