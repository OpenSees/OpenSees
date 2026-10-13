import math

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _build_static_system(material_modulus=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, material_modulus)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")


def test_bandgeneral_rejects_nonfinite_rhs_without_crashing():
    try:
        _build_static_system()
        ops.load(2, math.nan)
        assert ops.analyze(1) != 0
    finally:
        ops.wipe()


def test_bandgeneral_rejects_nonfinite_matrix_without_crashing():
    try:
        _build_static_system(math.nan)
        ops.load(2, 1.0)
        assert ops.analyze(1) != 0
    finally:
        ops.wipe()
