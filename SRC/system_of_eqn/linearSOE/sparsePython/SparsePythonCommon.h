#ifndef SparsePythonCommon_h
#define SparsePythonCommon_h

#include <Python.h>

/**
 * Storage scheme class: acts like an enum with string conversion.
 */
class SparsePythonStorageScheme {
private:
    enum class Value : int {
        CSR = 0,
        CSC = 1,
        COO = 2
    };
    Value value_;

public:
    // Static const members for easy access (defined in .cpp)
    static const SparsePythonStorageScheme CSR;
    static const SparsePythonStorageScheme CSC;
    static const SparsePythonStorageScheme COO;

    // Constructors
    SparsePythonStorageScheme() : value_(Value::CSR) {}
    explicit SparsePythonStorageScheme(Value v) : value_(v) {}
    SparsePythonStorageScheme(int v) : value_(static_cast<Value>(v)) {}

    // Implicit conversion operators for comparisons and assignments
    operator Value() const { return value_; }
    operator int() const { return static_cast<int>(value_); }

    // Comparison operators
    bool operator==(const SparsePythonStorageScheme& other) const {
        return value_ == other.value_;
    }
    bool operator!=(const SparsePythonStorageScheme& other) const {
        return value_ != other.value_;
    }

    // String conversion for Python
    const char* to_str() const;
};

/**
 * Matrix status class: acts like an enum with string conversion.
 */
class SparsePythonMatrixStatus {
private:
    enum class Value : int {
        UNCHANGED = 0,
        COEFFICIENTS_CHANGED = 1,
        STRUCTURE_CHANGED = 2
    };
    Value value_;

public:
    // Static const members for easy access (defined in .cpp)
    static const SparsePythonMatrixStatus UNCHANGED;
    static const SparsePythonMatrixStatus COEFFICIENTS_CHANGED;
    static const SparsePythonMatrixStatus STRUCTURE_CHANGED;

    // Constructors
    SparsePythonMatrixStatus() : value_(Value::UNCHANGED) {}
    explicit SparsePythonMatrixStatus(Value v) : value_(v) {}
    SparsePythonMatrixStatus(int v) : value_(static_cast<Value>(v)) {}

    // Implicit conversion operators for comparisons and assignments
    operator Value() const { return value_; }
    operator int() const { return static_cast<int>(value_); }

    // Comparison operators
    bool operator==(const SparsePythonMatrixStatus& other) const {
        return value_ == other.value_;
    }
    bool operator!=(const SparsePythonMatrixStatus& other) const {
        return value_ != other.value_;
    }

    // String conversion for Python
    const char* to_str() const;
};

/**
 * Flags controlling which buffers are exported as writable to Python solvers.
 * By default, all buffers except 'x' (solution) are read-only for safety.
 * Advanced users can enable writable buffers for in-place modifications.
 */
struct SparsePythonWritableFlags {
    bool values = false;  // Matrix coefficients (values array)
    bool rhs = false;     // Right-hand side vector
    // Note: index_ptr, indices, row, col are always read-only (structural)
    // Note: x (solution) is always writable (output buffer)

    SparsePythonWritableFlags() = default;
    SparsePythonWritableFlags(bool v, bool r) : values(v), rhs(r) {}
};

/**
 * RAII helper that manages reference counting for transient PyObject* instances.
 * Acquires the GIL before calling Py_XDECREF to ensure thread safety.
 */
struct PyObjectHolder {
    PyObjectHolder() = default;
    explicit PyObjectHolder(PyObject *obj) : ptr(obj) {}
    ~PyObjectHolder() {
        // Acquire GIL before calling Py_XDECREF to ensure thread safety
        // Note: We don't set ptr to nullptr here because the object is being destroyed
        if (ptr != nullptr && Py_IsInitialized()) {
            PyGILState_STATE gilState = PyGILState_Ensure();
            Py_XDECREF(ptr);
            PyGILState_Release(gilState);
        }
    }
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
        // reset() should only be called while GIL is held
        Py_XDECREF(ptr);
        ptr = obj;
    }
  private:
    PyObject *ptr{nullptr};
};

#endif /* SparsePythonCommon_h */
