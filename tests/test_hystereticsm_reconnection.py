try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


POSITIVE_ENVELOPE = (
    170.0908879721198, 0.0008256150522648083,
    806.5664389032247, 0.3299875444832722,
    68.03635518884792, 0.33779606261927597,
    0.0, 0.33794004766395014,
)
NEGATIVE_ENVELOPE = (
    -170.0908879721198, -0.0008256150522648083,
    -806.5664389032247, -0.3299875444832722,
    -68.03635518884792, -0.33779606261927597,
    0.0, -0.33794004766395014,
)
RETAINED_MAXIMUM = 0.0009391377421725002
COMMITTED_DEFORMATIONS = (
    RETAINED_MAXIMUM,
    0.0008689933708300856,
    -0.0004230444166321932,
)
EPSILON = 1.0e-12


def advance(increment):
    ops.integrator("DisplacementControl", 2, 2, increment)
    assert ops.analyze(1) == 0


def reconnection_residual(sign):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 0, 1)
    ops.uniaxialMaterial(
        "HystereticSM", 1,
        "-posEnv", *POSITIVE_ENVELOPE,
        "-negEnv", *NEGATIVE_ENVELOPE,
        "-pinch", 0.08, 0.98,
        "-damage", 0.0, 0.0,
        "-beta", 0.15,
        "-degEnv", 0.0, 0.0,
    )
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 2)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 1.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-8, 100)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 2, EPSILON)
    ops.analysis("Static")

    current = 0.0
    for deformation in COMMITTED_DEFORMATIONS:
        target = sign * deformation
        advance(target - current)
        current = target

    retained_extreme = sign * RETAINED_MAXIMUM
    inside = retained_extreme - sign * EPSILON
    advance(inside - current)
    force_inside = ops.eleResponse(1, "material", 1, "stress")[0]
    tangent_inside = ops.eleResponse(1, "material", 1, "tangent")[0]

    advance(retained_extreme - inside)
    force_at_extreme = ops.eleResponse(1, "material", 1, "stress")[0]
    ops.wipe()

    return (
        force_at_extreme
        - force_inside
        - tangent_inside * (retained_extreme - inside)
    )


def test_positive_reload_reconnects_continuously():
    assert abs(reconnection_residual(1.0)) < 1.0e-9


def test_negative_reload_reconnects_continuously():
    assert abs(reconnection_residual(-1.0)) < 1.0e-9
