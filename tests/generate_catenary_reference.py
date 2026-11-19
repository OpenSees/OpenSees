#!/usr/bin/env python
"""
Generate reference results for the old CatenaryCable element (3-DOF per node)
using OpenSeesPy from the ``~/engineering-venv`` virtual environment.

The script builds the 3-node, 2-element verification example from
Salehi Ahmad Abad et al. (2013) / OpenSeesWiki, runs a static analysis,
and writes the node displacements and element end-forces to
``catenary_reference.json``.

Usage (from the engineering-venv):
    source ~/engineering-venv/bin/activate
    python tests/generate_catenary_reference.py

Reference:
    Salehi Ahmad Abad, M., Shooshtari, A., Esmaeili, V., & Naghavi Riabi, A. (2013).
    Nonlinear analysis of cable structures under general loadings.
    Finite Elements in Analysis and Design, 73, 11–19.
    https://doi.org/10.1016/j.finel.2013.05.002
"""

import json
import math
import pathlib
import sys

# ---------------------------------------------------------------------------
#  Use OLD OpenSeesPy (3-DOF per node catenary)
# ---------------------------------------------------------------------------
try:
    import openseespy.opensees as ops
except ModuleNotFoundError:
    print("ERROR: OpenSeesPy not found. Activate the ~/engineering-venv first.")
    print("    source ~/engineering-venv/bin/activate")
    sys.exit(1)

# ---------------------------------------------------------------------------
#  Output path (same directory as this script)
# ---------------------------------------------------------------------------
_OUTPUT = pathlib.Path(__file__).resolve().parent / "catenary_reference.json"


def build_and_solve():
    """
    Build the 3-node, 2-element catenary model from the OpenSeesWiki example
    (x=30 case) and run the static analysis.

    Returns
    -------
    dict
        A dictionary with node displacements and element forces.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    # -----------------------------------------------------------------------
    # Geometry  (x=30 case from the wiki)
    # -----------------------------------------------------------------------
    x = 30.0
    ops.node(1, 0.0, 0.0, 90.0)
    ops.node(2, x / 2, 0.0, 40.0)
    ops.node(3, x, 60.0, 30.0)

    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)
    ops.fix(3, 1, 1, 1)

    # -----------------------------------------------------------------------
    # Cable parameters (from the wiki example)
    # -----------------------------------------------------------------------
    w3 = -0.00001          # weight per unit length
    E = 3.0e7              # Young's modulus
    A = 1.0                # cross-sectional area
    L0_total = 100.0       # total unstretched cable length
    alpha = 6.5e-6         # coefficient of thermal expansion
    temperature_change = 100.0  # temperature change
    rho = w3 / 9.81        # mass per unit length
    error_tol = 1e-6       # within-element convergence tolerance
    N_substeps = 20        # within-element sub-steps
    mass_type = 0          # lumped mass

    # Two cable elements, each L0_total/2 long
    ops.element(
        "CatenaryCable", 1, 1, 2,
        w3, E, A, L0_total / 2, alpha, temperature_change, rho,
        error_tol, N_substeps, mass_type
    )
    ops.element(
        "CatenaryCable", 2, 2, 3,
        w3, E, A, L0_total / 2, alpha, temperature_change, rho,
        error_tol, N_substeps, mass_type
    )

    # -----------------------------------------------------------------------
    # Gravity load
    # -----------------------------------------------------------------------
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 2, 1)
    ops.eleLoad("-ele", 1, 2, "-type", "-beamUniform", 0.0, 0.0, -1.0)

    # -----------------------------------------------------------------------
    # Static analysis
    # -----------------------------------------------------------------------
    ops.system("FullGeneral")
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.test("NormDispIncr", 1.0e-5, 100, 1)
    ops.integrator("LoadControl", 1.0 / 10)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.analyze(10)

    # -----------------------------------------------------------------------
    # Collect results
    # -----------------------------------------------------------------------
    data = {
        "description": (
            "Reference results for CatenaryCable element verification. "
            "Model: 3-node, 2-element catenary from Salehi Ahmad Abad (2013) / "
            "OpenSeesWiki example (x=30). "
            "Generated with old OpenSeesPy (3-DOF per node)."
        ),
        "parameters": {
            "x_span": x,
            "L0_total": L0_total,
            "weight_per_unit_length": w3,
            "E": E,
            "A": A,
            "alpha": alpha,
            "temperature_change": temperature_change,
            "rho": rho,
            "error_tol": error_tol,
            "N_substeps": N_substeps,
            "mass_type": mass_type,
        },
        "node_displacements": {
            str(tag): [ops.nodeDisp(tag, i + 1) for i in range(3)]
            for tag in (1, 2, 3)
        },
        "element_forces": {
            str(tag): [float(v) for v in ops.eleForce(tag)]
            for tag in (1, 2)
        },
    }

    return data


def main():
    print("Building old CatenaryCable model (3-DOF per node)...")
    data = build_and_solve()

    # Check the known analytical value as a sanity check
    d2 = data["node_displacements"]["2"]
    expected = (8.58693, 0.0, 2.82578)
    print(f"Node 2 displacement: {d2}")
    for got, exp in zip(d2, expected):
        if not math.isclose(got, exp, abs_tol=1e-4):
            print(
                f"WARNING: Node-2 displacement mismatch: "
                f"got {got}, expected {exp}"
            )

    with open(_OUTPUT, "w") as fh:
        json.dump(data, fh, indent=2)

    print(f"Reference results written to {_OUTPUT}")


if __name__ == "__main__":
    main()