"""
Simple benchmark comparing SciPy eigen solver with native OpenSees eigen solvers.

Model geometry from Michael H. Scott's blog post:
https://portwooddigital.com/2021/12/19/three-dimensional-meshing/

"""

from __future__ import annotations

import csv
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_linalg

# Import OpenSeesPy
SCRIPT_PATH = Path(__file__).resolve()
FILENAME = SCRIPT_PATH.name
SCRIPT_DIR = SCRIPT_PATH.parent
REPO_ROOT = SCRIPT_PATH.parents[2]
OPENSEESPY_BUILD = REPO_ROOT / "build"
print(f"[{FILENAME}] Importing OpenSeesPy from: {OPENSEESPY_BUILD}")
sys.path.append(str(OPENSEESPY_BUILD))
import opensees as ops


# -----------------------------------------------------------------------------
# SciPy Eigen Solver Definition
# -----------------------------------------------------------------------------


class SciPyGeneralizedEigenSolver:
    """SciPy-based generalized eigenvalue solver."""
    
    def __init__(self, maxiter=None, tol=0.0):
        self.maxiter = maxiter
        self.tol = tol
        self._k_matrix = None  # Cache stiffness matrix
        self._m_matrix = None  # Cache mass matrix
    
    def solve(self, **kwargs):
        # Extract buffers from kwargs
        index_ptr: memoryview = kwargs['index_ptr']  # int32, read-only
        indices: memoryview = kwargs['indices']  # int32, read-only
        k_values: memoryview = kwargs['k_values']  # float64, read-only
        m_values: memoryview = kwargs['m_values']  # float64, read-only
        eigenvalues: memoryview = kwargs['eigenvalues']  # float64, writable
        eigenvectors: memoryview = kwargs['eigenvectors']  # float64, writable (flat buffer)
        num_eqn: int = kwargs['num_eqn']
        nnz: int = kwargs['nnz']
        matrix_status: str = kwargs['matrix_status']
        num_modes: int = kwargs['num_modes']
        find_smallest: bool = kwargs['find_smallest']
        
        # Wrap memoryviews in NumPy arrays
        indptr = np.frombuffer(index_ptr, dtype=np.int32, count=num_eqn + 1)
        idx = np.frombuffer(indices, dtype=np.int32, count=nnz)
        k_data = np.frombuffer(k_values, dtype=np.float64, count=nnz)
        m_data = np.frombuffer(m_values, dtype=np.float64, count=nnz)
        
        # Rebuild matrices if structure changed (assume CSR storage scheme)
        if matrix_status == 'STRUCTURE_CHANGED' or self._k_matrix is None:
            self._k_matrix = sp.csr_matrix((k_data, idx, indptr), shape=(num_eqn, num_eqn))
            self._m_matrix = sp.csr_matrix((m_data, idx, indptr), shape=(num_eqn, num_eqn))
        elif matrix_status == 'COEFFICIENTS_CHANGED':
            # Update only the values
            self._k_matrix.data[:] = k_data
            self._m_matrix.data[:] = m_data
        
        # Solve generalized eigenvalue problem Kx = λMx
        eigsh_kwargs = {
            'k': num_modes,
            'M': self._m_matrix,
            'which': 'LM',
        }
        if find_smallest:
            eigsh_kwargs['sigma'] = 0.0  # shift-invert method for smallest eigenvalues
        if self.maxiter is not None:
            eigsh_kwargs['maxiter'] = self.maxiter
        if self.tol > 0.0:
            eigsh_kwargs['tol'] = self.tol
        
        eigenvalues_result, eigenvectors_result = sp_linalg.eigsh(self._k_matrix, **eigsh_kwargs)
        
        # Write eigenvalues directly to the writable buffer
        eigenvalues_buf = np.frombuffer(eigenvalues, dtype=np.float64, count=num_modes)
        eigenvalues_buf[:] = eigenvalues_result[:num_modes]
        
        # Write eigenvectors directly to the flat buffer (mode-major layout)
        # Buffer layout: eigenvectors[mode * num_eqn + eqn]
        # Transpose and flatten for efficient single assignment (eigenvectors_result is (num_eqn, num_modes))
        eigenvectors_buf = np.frombuffer(eigenvectors, dtype=np.float64, count=num_modes * num_eqn)
        eigenvectors_buf[:] = eigenvectors_result.T.flatten()
        
        return None  # Success


# -----------------------------------------------------------------------------
# Model Geometry and Material Properties
# -----------------------------------------------------------------------------

BAR_LENGTH = 10.0  # inches
BAR_HEIGHT = 2.0  # inches
BAR_THICKNESS = 1.0  # inches

ELASTIC_MODULUS = 29_000.0  # kip / in^2
POISSON_RATIO = 0.3
STEEL_DENSITY = 0.284e-3 / 386.4  # kip s^2 / in^4

# -----------------------------------------------------------------------------
# Model Building Functions
# -----------------------------------------------------------------------------


def build_solid_bar_model(nx: int, ny: int, nz: int) -> Optional[int]:
    """Create a structured hexahedral mesh using block3D.
    
    Returns the node tag at the far corner (x=BAR_LENGTH, y=BAR_THICKNESS/2, z=BAR_HEIGHT/2)
    for eigenvector reporting.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, ELASTIC_MODULUS, POISSON_RATIO, STEEL_DENSITY)

    eleType = "stdBrick"
    eleArgs = 1
    ops.block3D(
        nx, ny, nz, 1, 1, eleType, eleArgs,
        1, 0.0, -BAR_THICKNESS / 2.0, -BAR_HEIGHT / 2.0,
        2, BAR_LENGTH, -BAR_THICKNESS / 2.0, -BAR_HEIGHT / 2.0,
        3, BAR_LENGTH, BAR_THICKNESS / 2.0, -BAR_HEIGHT / 2.0,
        4, 0.0, BAR_THICKNESS / 2.0, -BAR_HEIGHT / 2.0,
        5, 0.0, -BAR_THICKNESS / 2.0, BAR_HEIGHT / 2.0,
        6, BAR_LENGTH, -BAR_THICKNESS / 2.0, BAR_HEIGHT / 2.0,
        7, BAR_LENGTH, BAR_THICKNESS / 2.0, BAR_HEIGHT / 2.0,
        8, 0.0, BAR_THICKNESS / 2.0, BAR_HEIGHT / 2.0,
    )

    ops.fixX(0.0, 1, 1, 1)
    
    # Find the node at the far corner for eigenvector reporting
    far_corner_node = None
    for node in ops.getNodeTags():
        x = ops.nodeCoord(node, 1)
        y = ops.nodeCoord(node, 2)
        z = ops.nodeCoord(node, 3)
        if np.isclose(
            [x, y, z],
            [BAR_LENGTH, BAR_THICKNESS / 2.0, BAR_HEIGHT / 2.0],
            atol=1e-9,
        ).all():
            far_corner_node = node
            break
    
    return far_corner_node


# -----------------------------------------------------------------------------
# Solver Run Functions
# -----------------------------------------------------------------------------


def run_gen_band_arpack(num_modes: int) -> Optional[List[float]]:
    """Run genBandArpack eigen solver."""
    return ops.eigen("genBandArpack", num_modes)


def run_full_gen_lapack(num_modes: int) -> Optional[List[float]]:
    """Run fullGenLapack eigen solver."""
    return ops.eigen("fullGenLapack", num_modes)


def run_scipy_eigsh(num_modes: int) -> Optional[List[float]]:
    """Run SciPy eigsh eigen solver."""
    solver = SciPyGeneralizedEigenSolver(maxiter=None, tol=0.0)
    return ops.eigen("PythonSparse", num_modes, {"solver": solver, "scheme": "CSR"})


# -----------------------------------------------------------------------------
# Benchmark Data Structures
# -----------------------------------------------------------------------------


@dataclass
class BenchmarkRow:
    solver_name: str
    mesh_factor: float
    mesh_size: float
    num_elements: int
    num_nodes: int
    num_equations: int
    status: int
    time_seconds: float
    eigenvalue_mode1: Optional[float] = None
    eigenvalue_mode5: Optional[float] = None
    eigenvec_mode1_x: Optional[float] = None
    eigenvec_mode1_y: Optional[float] = None
    eigenvec_mode1_z: Optional[float] = None
    eigenvec_mode5_x: Optional[float] = None
    eigenvec_mode5_y: Optional[float] = None
    eigenvec_mode5_z: Optional[float] = None
    eigenvalues_all: Optional[List[float]] = None  # First 5 eigenvalues for console output


CSV_HEADER = (
    "solver",
    "mesh_factor",
    "mesh_c",
    "num_elements",
    "num_nodes",
    "num_equations",
    "status",
    "time_seconds",
    "eigenvalue_mode1",
    "eigenvalue_mode5",
    "eigenvec_mode1_x",
    "eigenvec_mode1_y",
    "eigenvec_mode1_z",
    "eigenvec_mode5_x",
    "eigenvec_mode5_y",
    "eigenvec_mode5_z",
)


# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------


def counts_from_mesh_factor(factor: float) -> Tuple[float, Tuple[int, int, int]]:
    """Return mesh size and brick counts for a given refinement factor."""
    mesh_size = BAR_THICKNESS / factor

    def count(dim: float) -> int:
        return max(1, int(math.ceil(dim / mesh_size)))

    return mesh_size, (count(BAR_LENGTH), count(BAR_THICKNESS), count(BAR_HEIGHT))


# -----------------------------------------------------------------------------
# Benchmark Execution
# -----------------------------------------------------------------------------


def run_benchmark(
    solver_name: str,
    solver_run: Callable[[int], Optional[List[float]]],
    mesh_factor: float,
    mesh_size: float,
    counts: Tuple[int, int, int],
    num_modes: int = 5,
) -> BenchmarkRow:
    """Run benchmark for a single eigen solver and mesh configuration."""
    nx, ny, nz = counts
    far_corner_node = build_solid_bar_model(nx, ny, nz)
    
    # Initialize to get system size
    try:
        ops.system("Diagonal")  # to obtain the number of equations
        ops.analysis("Static", "-noWarnings")
        ops.analyze(1)
        neq = ops.systemSize()
    except AttributeError:
        neq = -1
    
    # Run eigen analysis
    start_time = time.perf_counter()
    eigenvalues = None
    status = -1
    
    try:
        eigenvalues = solver_run(num_modes)
        status = 0 if eigenvalues is not None and len(eigenvalues) > 0 else -1
    except Exception as e:
        print(f"Error in eigen analysis: {e}", file=sys.stderr)
        status = -1
    
    time_seconds = time.perf_counter() - start_time
    
    # Get system information
    num_elements = len(ops.getEleTags())
    num_nodes = len(ops.getNodeTags())
    
    # Extract eigenvalues for modes 1 and 5
    eigenvalue_mode1 = None
    eigenvalue_mode5 = None
    eigenvalues_all = None
    
    if eigenvalues is not None and len(eigenvalues) >= 5:
        eigenvalue_mode1 = eigenvalues[0]
        eigenvalue_mode5 = eigenvalues[4]
        eigenvalues_all = eigenvalues[:5]
    elif eigenvalues is not None and len(eigenvalues) > 0:
        eigenvalue_mode1 = eigenvalues[0]
        eigenvalues_all = list(eigenvalues) + [None] * (5 - len(eigenvalues))
    
    # Extract eigenvectors for modes 1 and 5 at the far corner node
    eigenvec_mode1 = None
    eigenvec_mode5 = None
    
    if status == 0 and far_corner_node is not None:
        try:
            # Mode 1 (1-based indexing in OpenSees)
            if eigenvalues is not None and len(eigenvalues) >= 1:
                eigenvec_mode1 = ops.nodeEigenvector(far_corner_node, 1)
            # Mode 5 (1-based indexing in OpenSees)
            if eigenvalues is not None and len(eigenvalues) >= 5:
                eigenvec_mode5 = ops.nodeEigenvector(far_corner_node, 5)
        except Exception:
            # If we can't get eigenvectors, mark as failed
            status = -1
    
    if eigenvec_mode1 is not None:
        eigenvec_mode1_x, eigenvec_mode1_y, eigenvec_mode1_z = eigenvec_mode1
    else:
        eigenvec_mode1_x = eigenvec_mode1_y = eigenvec_mode1_z = None
    
    if eigenvec_mode5 is not None:
        eigenvec_mode5_x, eigenvec_mode5_y, eigenvec_mode5_z = eigenvec_mode5
    else:
        eigenvec_mode5_x = eigenvec_mode5_y = eigenvec_mode5_z = None

    return BenchmarkRow(
        solver_name=solver_name,
        mesh_factor=mesh_factor,
        mesh_size=mesh_size,
        num_elements=num_elements,
        num_nodes=num_nodes,
        num_equations=neq,
        status=status,
        time_seconds=time_seconds,
        eigenvalue_mode1=eigenvalue_mode1,
        eigenvalue_mode5=eigenvalue_mode5,
        eigenvec_mode1_x=eigenvec_mode1_x,
        eigenvec_mode1_y=eigenvec_mode1_y,
        eigenvec_mode1_z=eigenvec_mode1_z,
        eigenvec_mode5_x=eigenvec_mode5_x,
        eigenvec_mode5_y=eigenvec_mode5_y,
        eigenvec_mode5_z=eigenvec_mode5_z,
        eigenvalues_all=eigenvalues_all,
    )


# -----------------------------------------------------------------------------
# Main Function
# -----------------------------------------------------------------------------


def main():
    """Run benchmarks and output CSV."""
    if len(sys.argv) < 2:
        print(f"Usage: python {FILENAME} <output.csv>")
        sys.exit(1)

    # Handle output CSV path: if simple filename, store in script directory
    output_arg = sys.argv[1]
    if '/' not in output_arg and '\\' not in output_arg:
        # Simple filename, store in script directory
        output_csv = SCRIPT_DIR / output_arg
    else:
        # Path provided, use as-is
        output_csv = Path(output_arg)
    
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    # Solver configurations: (name, run_function)
    solvers = [
        ("genBandArpack", run_gen_band_arpack),
        ("fullGenLapack", run_full_gen_lapack),
        ("SciPyEigsh", run_scipy_eigsh),
    ]

    # Known problematic solver/mesh factor thresholds that cause segmentation faults or are too large
    # Format: solver_name -> mesh_factor limit (skip if mesh_factor >= limit)
    SKIP_LIMITS = {
        "fullGenLapack": 4.0,
    }

    # Mesh refinement factors
    mesh_factors = [1.0, 2.0, 3.0, 4.0]

    print("\n=== Python Sparse Eigen Solver Benchmark ===")
    print(f"Comparing: {[s[0] for s in solvers]}")
    print(f"Mesh factors: {mesh_factors}")
    print(f"Results will be written to: {output_csv}\n")

    with output_csv.open("w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(CSV_HEADER)

        for mesh_factor in mesh_factors:
            mesh_size, counts = counts_from_mesh_factor(mesh_factor)
            mesh_info_printed = False

            for solver_name, solver_run in solvers:
                # Skip if mesh_factor exceeds the limit for this solver
                limit = SKIP_LIMITS.get(solver_name, float('inf'))
                if mesh_factor >= limit:
                    # Write skipped row
                    writer.writerow((
                        solver_name,
                        mesh_factor,
                        mesh_size,
                        -1,  # num_elements
                        -1,  # num_nodes
                        -1,  # num_equations
                        -999,  # status
                        float("nan"),
                        "",  # eigenvalue_mode1
                        "",  # eigenvalue_mode5
                        "",  # eigenvec_mode1_x
                        "",  # eigenvec_mode1_y
                        "",  # eigenvec_mode1_z
                        "",  # eigenvec_mode5_x
                        "",  # eigenvec_mode5_y
                        "",  # eigenvec_mode5_z
                    ))
                    csvfile.flush()
                    
                    status_str = "⊘"
                    # Print mesh info after first solver completes, then its status on one line
                    if not mesh_info_printed:
                        print(f"\nMesh factor {mesh_factor:.1f}: (skipped)")
                        print(f"  Running {solver_name}... {status_str} (skipped - either known segfault or too large)")
                        mesh_info_printed = True
                    else:
                        print(f"  Running {solver_name}... {status_str} (skipped - either known segfault or too large)")
                    continue
                
                try:
                    # For first solver, don't print "Running..." yet - we'll print it with mesh info
                    if mesh_info_printed:
                        print(f"  Running {solver_name}...", end=" ", flush=True)
                    
                    row = run_benchmark(
                        solver_name=solver_name,
                        solver_run=solver_run,
                        mesh_factor=mesh_factor,
                        mesh_size=mesh_size,
                        counts=counts,
                    )
                    
                    writer.writerow((
                        row.solver_name,
                        row.mesh_factor,
                        row.mesh_size,
                        row.num_elements,
                        row.num_nodes,
                        row.num_equations,
                        row.status,
                        row.time_seconds,
                        row.eigenvalue_mode1 if row.eigenvalue_mode1 is not None else "",
                        row.eigenvalue_mode5 if row.eigenvalue_mode5 is not None else "",
                        row.eigenvec_mode1_x if row.eigenvec_mode1_x is not None else "",
                        row.eigenvec_mode1_y if row.eigenvec_mode1_y is not None else "",
                        row.eigenvec_mode1_z if row.eigenvec_mode1_z is not None else "",
                        row.eigenvec_mode5_x if row.eigenvec_mode5_x is not None else "",
                        row.eigenvec_mode5_y if row.eigenvec_mode5_y is not None else "",
                        row.eigenvec_mode5_z if row.eigenvec_mode5_z is not None else "",
                    ))
                    csvfile.flush()
                    
                    status_str = "✓" if row.status == 0 else "✗"
                    # Print mesh info after first solver completes, then its status on one line
                    if not mesh_info_printed:
                        print(f"\nMesh factor {mesh_factor:.1f}: {row.num_elements} elements, {row.num_nodes} nodes, {row.num_equations} equations")
                        # Format first 5 eigenvalues for display
                        if row.eigenvalues_all:
                            eig_str = ", ".join(
                                f"{e:.6e}" if e is not None else "N/A"
                                for e in row.eigenvalues_all
                            )
                        else:
                            eig_str = "N/A"
                        print(f"  Running {solver_name}... {status_str} ({row.time_seconds:.3f}s)  λ=[{eig_str}]")
                        mesh_info_printed = True
                    else:
                        # Format first 5 eigenvalues for display
                        if row.eigenvalues_all:
                            eig_str = ", ".join(
                                f"{e:.6e}" if e is not None else "N/A"
                                for e in row.eigenvalues_all
                            )
                        else:
                            eig_str = "N/A"
                        print(f"{status_str} ({row.time_seconds:.3f}s)  λ=[{eig_str}]")
                except Exception as e:
                    print(f"✗ Error: {e}")
                    # Write error row
                    writer.writerow((
                        solver_name,
                        mesh_factor,
                        mesh_size,
                        -1,  # num_elements
                        -1,  # num_nodes
                        -1,  # num_equations
                        -999,  # status
                        float("nan"),
                        "",  # eigenvalue_mode1
                        "",  # eigenvalue_mode5
                        "",  # eigenvec_mode1_x
                        "",  # eigenvec_mode1_y
                        "",  # eigenvec_mode1_z
                        "",  # eigenvec_mode5_x
                        "",  # eigenvec_mode5_y
                        "",  # eigenvec_mode5_z
                    ))
                    csvfile.flush()

    print(f"\n✓ Results written to {output_csv}")


if __name__ == "__main__":
    main()

