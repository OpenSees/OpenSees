# MINPACK Solver Module for OpenSees

## Overview
This module embeds the classic **MINPACK** library into OpenSees, providing battle-tested solvers for **nonlinear equations** and **least-squares problems** in pure C.

## Solvers Included
| Function | Purpose | Key Feature |
| :--- | :--- | :--- |
| **`HYBRD`** | Solve systems of nonlinear equations `F(x) = 0` | No derivative required (uses internal finite-difference) |
| **`HYBRJ`** | Solve systems of nonlinear equations `F(x) = 0` | **Requires Jacobian** (faster/more robust if provided) |
| **`LMDIF`** | Solve nonlinear least-squares problems `min \|F(x)\|²` | No derivative required |
| **`LMDER`** | Solve nonlinear least-squares problems `min \|F(x)\|²` | **Requires Jacobian** |
| **`CHKDER`** | Check user-provided Jacobian correctness | Essential for debugging `HYBRJ`/`LMDER` |

## Why Use This?
1.  **Reliability**: Algorithms are the reference implementation from **MINPACK** (Argonne National Laboratory), used for decades in scientific computing.
2.  **Zero Dependencies**: Pure C code. No external BLAS/LAPACK or other libraries needed.
3.  **License Compatibility**: Licensed under the **BSD-3-Clause** license, identical to OpenSees core.
4.  **New Capabilities**: Enables **parameter inversion**, **model calibration**, and complex equilibrium searches directly within OpenSees.

## How to Integrate into Your Build
1.  Add `OTHER/minpack/minpack.c` to your compiler's source file list.
2.  Ensure the directory containing `OTHER/minpack/` is in your compiler's include path.
3.  Include the header in your C/C++ code: `#include "minpack.h"`.
4.  Link the math library (`-lm` on Unix-like systems).

## Quick Example
See the full example in [`example.cpp`](./example.cpp) which solves the equilibrium equations of a two-bar truss using `HYBRD`.

## Reference
- Original MINPACK documentation: [https://www.netlib.org/minpack/](https://www.netlib.org/minpack/)
- SciPy's `fsolve` and `leastsq` are high-level wrappers of these same routines.