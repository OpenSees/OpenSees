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
YIELD_ROTATION = POSITIVE_ENVELOPE[1]
PARTIAL_RELOAD_STEP = YIELD_ROTATION * 1.0e-11
TANGENT_PROBE = YIELD_ROTATION * 1.0e-7
ZERO_CROSSING_TARGET = 0.00055


def advance(increment):
    ops.integrator("DisplacementControl", 2, 2, increment)
    assert ops.analyze(1) == 0


def build_model():
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


def material_response():
    stress = ops.eleResponse(1, "material", 1, "stress")[0]
    tangent = ops.eleResponse(1, "material", 1, "tangent")[0]
    return stress, tangent


def run_deformation_history(history):
    build_model()
    current = 0.0
    responses = []
    for target in history:
        advance(target - current)
        current = target
        responses.append(material_response())
    ops.wipe()
    return responses


def reconnection_residual(sign):
    build_model()
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


def partial_reload_residual(sign, reload_step):
    peak = sign * 2.0 * YIELD_ROTATION
    reversal = sign * 1.9 * YIELD_ROTATION
    increment = sign * reload_step
    before, after = run_deformation_history(
        (peak, reversal, reversal + increment)
    )[-2:]

    assert sign * before[0] > 0.0
    return after[0] - before[0] - before[1] * increment, before, after


def zero_force_reload_rotation(sign):
    prefix = tuple(sign * deformation for deformation in COMMITTED_DEFORMATIONS)
    direction = sign
    probe_target = prefix[-1] + direction * EPSILON
    stress, tangent = run_deformation_history((*prefix, probe_target))[-1]
    zero_force_rotation = probe_target - stress / tangent

    target = sign * ZERO_CROSSING_TARGET
    assert direction * (zero_force_rotation - prefix[-1]) > 0.0
    assert direction * (target - zero_force_rotation) > 0.0
    return zero_force_rotation


def zero_crossing_responses(sign):
    prefix = tuple(sign * deformation for deformation in COMMITTED_DEFORMATIONS)
    target = sign * ZERO_CROSSING_TARGET
    direction = sign
    zero = zero_force_reload_rotation(sign)
    offset = abs(target - prefix[-1]) * 1.0e-6
    intermediate_targets = (
        (),
        (zero,),
        (zero - direction * offset,),
        (zero + direction * offset,),
    )

    return [
        run_deformation_history((*prefix, *intermediate, target))[-1]
        for intermediate in intermediate_targets
    ]


def test_positive_reload_reconnects_continuously():
    assert abs(reconnection_residual(1.0)) < 1.0e-9


def test_negative_reload_reconnects_continuously():
    assert abs(reconnection_residual(-1.0)) < 1.0e-9


def test_positive_partial_unload_tiny_reload_is_continuous():
    residual, _, _ = partial_reload_residual(1.0, PARTIAL_RELOAD_STEP)
    assert abs(residual) < 1.0e-9

    _, before, after = partial_reload_residual(1.0, TANGENT_PROBE)
    finite_difference = (after[0] - before[0]) / TANGENT_PROBE
    assert abs(finite_difference - after[1]) < 1.0e-3


def test_negative_partial_unload_tiny_reload_is_continuous():
    residual, _, _ = partial_reload_residual(-1.0, PARTIAL_RELOAD_STEP)
    assert abs(residual) < 1.0e-9

    _, before, after = partial_reload_residual(-1.0, TANGENT_PROBE)
    finite_difference = (after[0] - before[0]) / -TANGENT_PROBE
    assert abs(finite_difference - after[1]) < 1.0e-3


def test_positive_full_reversal_zero_crossing_is_subdivision_independent():
    responses = zero_crossing_responses(1.0)
    stresses = [stress for stress, _ in responses]
    tangents = [tangent for _, tangent in responses]
    assert max(stresses) - min(stresses) < 1.0e-9
    assert max(tangents) - min(tangents) < 1.0e-7


def test_negative_full_reversal_zero_crossing_is_subdivision_independent():
    responses = zero_crossing_responses(-1.0)
    stresses = [stress for stress, _ in responses]
    tangents = [tangent for _, tangent in responses]
    assert max(stresses) - min(stresses) < 1.0e-9
    assert max(tangents) - min(tangents) < 1.0e-7
