"""
Simple benchmark comparing CuPyCG solver with native OpenSees solvers.

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

# Import OpenSeesPy
SCRIPT_PATH = Path(__file__).resolve()
FILENAME = SCRIPT_PATH.name
SCRIPT_DIR = SCRIPT_PATH.parent
REPO_ROOT = SCRIPT_PATH.parents[2]
OPENSEESPY_BUILD = REPO_ROOT / "build"
print(f"[{FILENAME}] Importing OpenSeesPy from: {OPENSEESPY_BUILD}")
sys.path.append(str(OPENSEESPY_BUILD))
import opensees as ops

# Import CuPy for GPU solvers
import cupy as cp
import cupyx.scipy.sparse.linalg


# -----------------------------------------------------------------------------
# CuPy CG Solver Definition
# -----------------------------------------------------------------------------


class CuPyCGSolver:
    """CuPy-based conjugate gradient solver for linear systems."""
    
    def __init__(self, rtol=1e-5, atol=1e-12, maxiter=None):
        self.rtol = rtol
        self.atol = atol
        self.maxiter = maxiter
        self.A = None  # Cache the sparse matrix
    
    def solve(self, **kwargs):
        # Extract buffers from kwargs
        index_ptr: memoryview = kwargs['index_ptr']  # int32, read-only memoryview
        indices: memoryview = kwargs['indices']  # int32, read-only memoryview
        values: memoryview = kwargs['values']  # float64, read-only memoryview
        rhs: memoryview = kwargs['rhs']  # float64, read-only memoryview
        x: memoryview = kwargs['x']  # float64, writeable memoryview
        num_eqn: int = kwargs['num_eqn']
        nnz: int = kwargs['nnz']
        matrix_status: str = kwargs['matrix_status']  # UNCHANGED, STRUCTURE_CHANGED, COEFFICIENTS_CHANGED
        
        # Wrap memoryviews using zero-copy numpy views
        indptr = np.frombuffer(index_ptr, dtype=np.int32, count=num_eqn + 1)
        idx = np.frombuffer(indices, dtype=np.int32, count=nnz)
        vals = np.frombuffer(values, dtype=np.float64, count=nnz)
        
        # Rebuild matrix if structure changed, update values if coefficients changed
        if matrix_status == 'STRUCTURE_CHANGED' or self.A is None:
            # Copy the entire CSR matrix to the GPU
            values_gpu = cp.asarray(vals)
            indices_gpu = cp.asarray(idx)
            index_ptr_gpu = cp.asarray(indptr)
            self.A = cp.sparse.csr_matrix((values_gpu, indices_gpu, index_ptr_gpu), shape=(num_eqn, num_eqn))
        elif matrix_status == 'COEFFICIENTS_CHANGED':
            # Update the values of the CSR matrix on the GPU
            values_gpu = cp.asarray(vals)
            self.A.data[:] = values_gpu  # in-place update
        else:
            # If UNCHANGED, do nothing
            pass
        
        # Wrap RHS for solving
        rhs_buf = np.frombuffer(rhs, dtype=np.float64, count=num_eqn)
        rhs_gpu = cp.asarray(rhs_buf)
        
        # Solve using conjugate gradient (without preconditioning) on the GPU
        x_gpu, info = cupyx.scipy.sparse.linalg.cg(self.A, rhs_gpu, tol=self.rtol, atol=self.atol, maxiter=self.maxiter)
        
        # Copy result back to CPU buffer
        x_buf = np.frombuffer(x, dtype=np.float64, count=num_eqn)
        x_buf[:] = cp.asnumpy(x_gpu)  # in-place update
        return -int(info)  # Return the info from the solver


# -----------------------------------------------------------------------------
# Model Geometry and Material Properties
# -----------------------------------------------------------------------------

BAR_LENGTH = 10.0  # inches
BAR_HEIGHT = 2.0  # inches
BAR_THICKNESS = 1.0  # inches

ELASTIC_MODULUS = 29_000.0  # kip / in^2
POISSON_RATIO = 0.3
STEEL_DENSITY = 0.284e-3 / 386.4  # kip s^2 / in^4
YIELD_STRESS = 50.0  # kip / in^2

# -----------------------------------------------------------------------------
# Model Building Functions
# -----------------------------------------------------------------------------


def build_solid_bar_model(nx: int, ny: int, nz: int) -> None:
    """Create a structured hexahedral mesh using block3D."""
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


def apply_load_pattern(total_load: float) -> None:
    """Apply a vertical load to all nodes on the free face."""
    tStart, tEnd, period = 0.0, 6.0, 4.0
    ops.timeSeries('Trig', 1, tStart, tEnd, period, '-factor', 1.0)
    ops.pattern("Plain", 1, 1)
    far_x_nodes = [
        node for node in ops.getNodeTags()
        if math.isclose(ops.nodeCoord(node, 1), BAR_LENGTH, abs_tol=1e-9)
    ]
    if not far_x_nodes:
        raise ValueError("No nodes found on the far face (x = BAR_LENGTH)")
    load_per_node = total_load / len(far_x_nodes)
    for node in far_x_nodes:
        ops.load(node, 0.0, 0.0, -load_per_node)


# -----------------------------------------------------------------------------
# Analysis Configuration
# -----------------------------------------------------------------------------


def configure_static_analysis(solver_setup: Callable[[], None], numberer: str, num_steps: int, tol: float, max_iter: int) -> None:
    """Configure the analysis for a given solver."""
    ops.constraints("Plain")
    ops.numberer(numberer)
    solver_setup()
    ops.integrator("LoadControl", 1.0 / num_steps)
    ops.test("NormUnbalance", tol, max_iter)
    ops.algorithm("ModifiedNewton", "-FactorOnce")
    ops.analysis("Static")


# -----------------------------------------------------------------------------
# Solver Setup Functions
# -----------------------------------------------------------------------------


def setup_bandspd() -> None:
    """Setup BandSPD solver."""
    ops.system("BandSPD")


def setup_umfpack() -> None:
    """Setup UmfPack solver."""
    ops.system("UmfPack")


def setup_cupy_cg() -> None:
    """Setup CuPy CG solver."""
    solver = CuPyCGSolver(rtol=1.0e-7, atol=1.0e-12, maxiter=None)
    ops.system("PythonSparse", {"solver": solver, "scheme": "CSR"})


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
    displacement_x: Optional[float] = None
    displacement_y: Optional[float] = None
    displacement_z: Optional[float] = None


CSV_HEADER = (
    "solver",
    "mesh_factor",
    "mesh_c",
    "num_elements",
    "num_nodes",
    "num_equations",
    "status",
    "time_seconds",
    "displacement_x",
    "displacement_y",
    "displacement_z",
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
    solver_setup: Callable[[], None],
    numberer: str,
    mesh_factor: float,
    mesh_size: float,
    counts: Tuple[int, int, int],
    num_steps: int = 5,
    tol: float = 1.0e-7,
    max_iter: int = 50,
) -> BenchmarkRow:
    """Run benchmark for a single solver and mesh configuration."""
    nx, ny, nz = counts
    build_solid_bar_model(nx, ny, nz)

    # Calculate load
    total_load = 1.25 * YIELD_STRESS * (BAR_THICKNESS * BAR_HEIGHT**2) / (6 * BAR_LENGTH)
    apply_load_pattern(total_load)

    # Static analysis
    configure_static_analysis(solver_setup, numberer, num_steps, tol, max_iter)
    start_time = time.perf_counter()
    status = ops.analyze(num_steps)
    time_seconds = time.perf_counter() - start_time

    # Get system information
    num_elements = len(ops.getEleTags())
    num_nodes = len(ops.getNodeTags())
    try:
        neq = ops.systemSize()
    except AttributeError:
        neq = -1

    # Get displacement at far corner
    displacement = None
    if status == 0:
        for node in ops.getNodeTags():
            x = ops.nodeCoord(node, 1)
            y = ops.nodeCoord(node, 2)
            z = ops.nodeCoord(node, 3)
            if np.isclose(
                [x, y, z],
                [BAR_LENGTH, BAR_THICKNESS / 2.0, BAR_HEIGHT / 2.0],
                atol=1e-9,
            ).all():
                dx = ops.nodeDisp(node, 1)
                dy = ops.nodeDisp(node, 2)
                dz = ops.nodeDisp(node, 3)
                displacement = (dx, dy, dz)
                break

    if displacement is not None:
        dx, dy, dz = displacement
    else:
        dx = dy = dz = None

    return BenchmarkRow(
        solver_name=solver_name,
        mesh_factor=mesh_factor,
        mesh_size=mesh_size,
        num_elements=num_elements,
        num_nodes=num_nodes,
        num_equations=neq,
        status=status,
        time_seconds=time_seconds,
        displacement_x=dx,
        displacement_y=dy,
        displacement_z=dz,
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

    # Solver configurations: (name, setup_function, numberer)
    solvers = [
        ("BandSPD", setup_bandspd, "RCM"),
        ("UmfPack", setup_umfpack, "Plain"),
        ("CuPyCG", setup_cupy_cg, "RCM"),
    ]

    # Known problematic solver/mesh factor thresholds that cause segmentation faults or are too large
    # Format: solver_name -> mesh_factor limit (skip if mesh_factor >= limit)
    SKIP_LIMITS = {
        "UmfPack": 14.0,
        "BandSPD": 20.0,
    }

    # Mesh refinement factors
    mesh_factors = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0]

    print("\n=== Python Sparse Solver Benchmark ===")
    print(f"Comparing: {[s[0] for s in solvers]}")
    print(f"Mesh factors: {mesh_factors}")
    print(f"Results will be written to: {output_csv}\n")

    with output_csv.open("w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(CSV_HEADER)

        for mesh_factor in mesh_factors:
            mesh_size, counts = counts_from_mesh_factor(mesh_factor)
            mesh_info_printed = False

            for solver_name, solver_setup, numberer in solvers:
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
                        "",
                        "",
                        "",
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
                        solver_setup=solver_setup,
                        numberer=numberer,
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
                        row.displacement_x if row.displacement_x is not None else "",
                        row.displacement_y if row.displacement_y is not None else "",
                        row.displacement_z if row.displacement_z is not None else "",
                    ))
                    csvfile.flush()
                    
                    status_str = "✓" if row.status == 0 else "✗"
                    # Print mesh info after first solver completes, then its status on one line
                    if not mesh_info_printed:
                        print(f"\nMesh factor {mesh_factor:.1f}: {row.num_elements} elements, {row.num_nodes} nodes, {row.num_equations} equations")
                        print(f"  Running {solver_name}... {status_str} ({row.time_seconds:.3f}s)")
                        mesh_info_printed = True
                    else:
                        print(f"{status_str} ({row.time_seconds:.3f}s)")
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
                        "",
                        "",
                        "",
                    ))
                    csvfile.flush()

    print(f"\n✓ Results written to {output_csv}")


if __name__ == "__main__":
    main()

