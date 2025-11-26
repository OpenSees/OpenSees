#ifndef SparsePythonFactory_h
#define SparsePythonFactory_h

/**
 * Parse the `system PythonSparse` command and construct the corresponding SOE.
 *
 * Expected syntax:
 *   `system 'PythonSparse' {'solver': SolverObject, 'scheme': 'CSR'|'CSC'|'COO', 'writable': 'values'|'rhs'|'values,rhs'|'all'|'none'}`
 *
 * The supplied Python object must expose a `solve` method. When invoked the
 * solver receives memoryviews pointing directly to the SOE storage:
 *   - CSR/CSC schemes share the `index_ptr` (int32) and `indices` (int32) arrays
 *   - COO exposes `row` and `col` (int32) arrays instead
 *   - All schemes provide `values`, `rhs`, `x` (float64) along with metadata
 *     (`num_eqn`, `nnz`, `matrix_status`, `storage_scheme`).
 *
 * The callable must operate in place on these buffers and return either `None`
 * for success or an integer status/error code. DO NOT rebind the keywords;
 * wrap the memoryviews with NumPy/JAX/CuPy arrays so updates propagate directly
 * to the underlying C++ storage.
 *
 * The returned pointer owns a Python-backed `LinearSOE` instance and should be
 * treated as a `LinearSOE*` by the caller.
 */
void *OPS_SparsePythonSolver();

#endif /* SparsePythonFactory_h */
