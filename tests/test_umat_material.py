"""Regression test for the generic `nDMaterial UMAT` wrapper.

The wrapper (SRC/material/nD/UmatMaterial.cpp) runs an Abaqus-format UMAT
subroutine as an OpenSees nDMaterial. No constitutive model ships with it,
so this test brings its own: tests/umat_elastic.f90, a stateless isotropic
linear-elastic UMAT. The test compiles that file into a shared library with
gfortran, binds it with `-library`, and checks that a single stdBrick driven
by the UMAT reproduces the built-in ElasticIsotropic material to round-off.

Covered:
  * symbol resolution through `-library` (and the gfortran `_` suffix)
  * the Voigt reorder between OpenSees [11 22 33 12 23 31] and Abaqus
    [11 22 33 12 13 23] on strain in, stress out and the tangent
  * trial/commit handling over several Newton load steps
  * the `-check_tangent` diagnostic running without disturbing the result

Skipped when gfortran is not on PATH, or on Windows (the Windows symbol
binding path of the wrapper is compile-guarded only and not yet exercised).
"""

import os
import shutil
import subprocess
import sys
import tempfile

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E_MOD = 30000.0
NU = 0.3
N_STEPS = 5

HERE = os.path.dirname(os.path.abspath(__file__))
UMAT_SOURCE = os.path.join(HERE, "umat_elastic.f90")

pytestmark = pytest.mark.skipif(
    sys.platform.startswith("win"),
    reason="UMAT wrapper symbol binding is not exercised on Windows yet",
)


@pytest.fixture(scope="module")
def umat_library():
    """Compile umat_elastic.f90 into a shared library, return its path."""
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran not found on PATH")
    build_dir = tempfile.mkdtemp(prefix="umat_elastic_")
    library = os.path.join(build_dir, "libumat_elastic.so")
    subprocess.run(
        [gfortran, "-shared", "-fPIC", "-O2", "-o", library, UMAT_SOURCE],
        check=True,
    )
    yield library
    shutil.rmtree(build_dir, ignore_errors=True)


def _build_brick(define_material):
    """Unit cube stdBrick, base fixed, loaded on the top face.

    `define_material` is called once the model is open and must define
    nDMaterial tag 1. The load has vertical, lateral and shear components
    so that every stress component and every entry of the tangent takes
    part.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    define_material()
    coords = {
        1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0),
        3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
        5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0),
        7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
    }
    for tag, xyz in coords.items():
        ops.node(tag, *xyz)
    for tag in (1, 2, 3, 4):
        ops.fix(tag, 1, 1, 1)

    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in (5, 6, 7, 8):
        ops.load(tag, 3.0, -2.0, -10.0)

    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / N_STEPS)
    ops.analysis("Static")


def _run_and_collect():
    for _ in range(N_STEPS):
        assert ops.analyze(1) == 0
    disp = [ops.nodeDisp(tag, dof) for tag in (5, 6, 7, 8) for dof in (1, 2, 3)]
    stress = ops.eleResponse(1, "stresses")
    assert len(stress) == 8 * 6
    return disp, stress


def _run_reference():
    _build_brick(lambda: ops.nDMaterial("ElasticIsotropic", 1, E_MOD, NU))
    try:
        return _run_and_collect()
    finally:
        ops.wipe()


def _run_umat(library, extra_args=()):
    def define_umat():
        ops.nDMaterial(
            "UMAT", 1,
            "-symbol", "umat_elastic",
            "-nprops", 2, "-props", E_MOD, NU,
            "-nstatv", 0,
            "-library", library,
            *extra_args,
        )

    _build_brick(define_umat)
    try:
        return _run_and_collect()
    finally:
        ops.wipe()


def _assert_close(got, want, rel=1.0e-9):
    scale = max(abs(v) for v in want)
    assert scale > 0.0
    for g, w in zip(got, want):
        assert abs(g - w) <= rel * scale, (g, w)


def test_umat_elastic_matches_elastic_isotropic(umat_library):
    disp_ref, stress_ref = _run_reference()
    disp_umat, stress_umat = _run_umat(umat_library)
    _assert_close(disp_umat, disp_ref)
    _assert_close(stress_umat, stress_ref)


def test_umat_check_tangent_does_not_change_response(umat_library):
    disp_ref, stress_ref = _run_reference()
    disp_umat, stress_umat = _run_umat(umat_library, ("-check_tangent", 1.0e-7))
    _assert_close(disp_umat, disp_ref)
    _assert_close(stress_umat, stress_ref)


def test_umat_unknown_symbol_is_rejected(umat_library):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    try:
        with pytest.raises(Exception):
            ops.nDMaterial(
                "UMAT", 1,
                "-symbol", "no_such_umat",
                "-nprops", 2, "-props", E_MOD, NU,
                "-nstatv", 0,
                "-library", umat_library,
            )
    finally:
        ops.wipe()
