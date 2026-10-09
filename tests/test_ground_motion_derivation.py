import math

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _configure_transient_analysis():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def _build_dummy_system():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.fix(1, 1, 1, 0)
    ops.fix(2, 0, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    ops.element("truss", 1, 1, 2, 1.0, 1)


def test_multiple_support_velocity_does_not_create_acceleration():
    _build_dummy_system()
    pressure = 361.008
    ops.timeSeries("Constant", 1, "-factor", pressure)
    ops.pattern("MultipleSupport", 1)
    ops.groundMotion(1, "Plain", "-vel", 1)
    ops.imposedMotion(1, 3, 1)
    _configure_transient_analysis()

    assert ops.analyze(1, 0.1) == 0
    assert math.isclose(ops.nodeVel(1, 3), pressure)
    assert math.isclose(ops.nodeAccel(1, 3), 0.0, abs_tol=1.0e-12)


def _uniform_excitation_response(series_option, values, num_steps=None):
    ops.wipe()
    if num_steps is None:
        num_steps = len(values) - 1
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.mass(2, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 10.0)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Path", 1, "-dt", 0.1, "-values", *values)
    ops.pattern("UniformExcitation", 1, 1, series_option, 1)
    _configure_transient_analysis()

    response = []
    for _ in range(num_steps):
        assert ops.analyze(1, 0.1) == 0
        response.append(ops.nodeDisp(2, 1))
    return response


def test_uniform_excitation_derives_acceleration_from_velocity():
    acceleration = [0.0, 1.0, -1.0, 0.0]
    velocity = [0.0, 0.05, 0.05, 0.0]

    from_acceleration = _uniform_excitation_response("-accel", acceleration)
    from_velocity = _uniform_excitation_response("-vel", velocity)

    assert len(from_acceleration) == len(from_velocity)
    for expected, actual in zip(from_acceleration, from_velocity):
        assert math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)


def test_uniform_excitation_derives_acceleration_from_displacement():
    acceleration = [0.0, 1.0, -1.0, 0.0, 0.0]
    displacement = [0.0, 0.0025, 0.0075, 0.01, 0.01]

    from_acceleration = _uniform_excitation_response("-accel", acceleration, 3)
    from_displacement = _uniform_excitation_response("-disp", displacement, 3)

    assert len(from_acceleration) == len(from_displacement)
    for expected, actual in zip(from_acceleration, from_displacement):
        assert math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)


def test_explicit_acceleration_is_preserved_for_multiple_support():
    _build_dummy_system()
    acceleration = 3.25
    ops.timeSeries("Constant", 1, "-factor", acceleration)
    ops.pattern("MultipleSupport", 1)
    ops.groundMotion(1, "Plain", "-accel", 1)
    ops.imposedMotion(1, 3, 1)
    _configure_transient_analysis()

    assert ops.analyze(1, 0.1) == 0
    assert math.isclose(ops.nodeAccel(1, 3), acceleration)
