#ifndef SparsePythonEigenFactory_h
#define SparsePythonEigenFactory_h

/**
 * Parse the `eigen PythonSparse` command and construct the corresponding EigenSOE.
 *
 * Expected syntax:
 *   `eigen PythonSparse numModes {'solver': SolverObject, 'scheme': 'CSR'|'CSC'|'COO'}`
 *
 * The supplied Python object must expose a `solve` method. When invoked the
 * solver receives memoryviews pointing directly to the EigenSOE storage:
 *   - Shared sparsity pattern: `index_ptr` (int32) and `indices` (int32) arrays
 *   - Separate value arrays: `k_values` (float64) for stiffness, `m_values` (float64) for mass
 *   - Metadata: `num_eqn`, `nnz`, `matrix_status`, `storage_scheme`, `num_modes`, `generalized`, `find_smallest`
 *
 * The solver must return a tuple (eigenvalues, eigenvectors) where:
 *   - eigenvalues: list of floats
 *   - eigenvectors: list of lists (each inner list is an eigenvector)
 *
 * The returned pointer owns a Python-backed `EigenSOE` instance and should be
 * treated as an `EigenSOE*` by the caller.
 */
void *OPS_SparsePythonEigenSolver();

#endif /* SparsePythonEigenFactory_h */

