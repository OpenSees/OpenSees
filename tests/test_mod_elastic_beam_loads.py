"""Load equilibrium and fixed-end checks for the remaining gap in issue #1612."""

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


@pytest.fixture(autouse=True)
def clean_domain():
    ops.wipe()
    yield
    ops.wipe()


def _model(kind="ModElasticBeam3d", modifiers=(4., 4., 2.), fixed_bending=False):
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0., 0., 0.)
    ops.node(2, 10., 0., 0.)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    if fixed_bending:
        # Keep the axial equation free to avoid a zero-equation solver.
        ops.fix(2, 0, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0., 0., 1.)
    extra = (*modifiers, *modifiers) if kind == "ModElasticBeam3d" else ()
    ops.element(kind, 1, 1, 2, .01, 2.e11, 8.e10, 1.e-4, 1.e-4, 1.e-4, *extra, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)


def _load(axis, start, end, a, b):
    first, last = [0.] * 3, [0.] * 3
    first[axis], last[axis] = start, end
    # Command component order is local y, z, x.
    ops.eleLoad("-ele", 1, "-type", "-beamUniform",
                first[1], first[2], first[0], a, b, last[1], last[2], last[0])


def _analyze(factor=1.):
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", factor)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()


def _integral(coefficients, start, end, a, b):
    # Exact polynomial antiderivative in x/L; independent of element quadrature.
    if a == b:
        return 0.
    slope = (end - start) / (b - a)
    intercept = start - slope * a
    return 10. * sum(c * (intercept * (b**(i+1) - a**(i+1)) / (i+1)
                         + slope * (b**(i+2) - a**(i+2)) / (i+2))
                     for i, c in enumerate(coefficients))


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("start,end", [(-2., -2.), (-2., 0.), (0., -2.), (-2., 3.)])
@pytest.mark.parametrize("a,b", [(0., 1.), (.2, .8), (.4, .4)])
@pytest.mark.parametrize("modifiers", [(4., 4., 2.), (3.6, 3.6, 2.4)])
def test_cantilever_force_and_moment_balance(axis, start, end, a, b, modifiers):
    _model(modifiers=modifiers)
    _load(axis, start, end, a, b)
    _analyze(.5)
    expected = [0.] * 6
    expected[axis] = -.5 * _integral([1.], start, end, a, b)
    if axis != 0:
        moment = .5 * _integral([0., 10.], start, end, a, b)
        expected[5 if axis == 1 else 4] = -moment if axis == 1 else moment
    assert ops.nodeReaction(1) == pytest.approx(expected, abs=1.e-9)


@pytest.mark.parametrize("axis", [1, 2])
@pytest.mark.parametrize("start,end", [(-2., 0.), (0., -2.), (-2., 3.), (-2., -2.)])
@pytest.mark.parametrize("a,b", [(0., 1.), (.2, .8)])
def test_fixed_end_forces_in_ordinary_beam_limit(axis, start, end, a, b):
    _model(fixed_bending=True)
    _load(axis, start, end, a, b)
    _analyze()
    shears = ([1., 0., -3., 2.], [0., 0., 3., -2.])
    moments = ([0., 10., -20., 10.], [0., 0., -10., 10.])
    for index, tag in enumerate((1, 2)):
        expected = [0.] * 6
        expected[axis] = -_integral(shears[index], start, end, a, b)
        moment = _integral(moments[index], start, end, a, b)
        expected[5 if axis == 1 else 4] = -moment if axis == 1 else moment
        assert ops.nodeReaction(tag) == pytest.approx(expected, abs=1.e-9)


def test_ordinary_beam_limit_displacements_and_load_superposition():
    results = []
    for kind in ("elasticBeamColumn", "ModElasticBeam3d"):
        ops.wipe()
        _model(kind=kind)
        _load(1, -2., 3., .2, .8)
        _load(2, 0., -4., 0., 1.)
        _load(0, 5., -1., .1, .9)
        _analyze()
        results.append((ops.nodeDisp(2), ops.eleForce(1)))
    for actual, expected in zip(results[1], results[0]):
        assert actual == pytest.approx(expected, rel=1.e-10, abs=1.e-10)


def test_partial_uniform_matches_existing_full_uniform_load():
    results = []
    for partial in (False, True):
        ops.wipe()
        _model(modifiers=(3.6, 3.6, 2.4))
        if partial:
            for axis, value in enumerate((3., -2., 4.)):
                _load(axis, value, value, 0., 1.)
        else:
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", -2., 4., 3.)
        _analyze()
        results.append((ops.nodeDisp(2), ops.eleForce(1)))
    for actual, expected in zip(results[1], results[0]):
        assert actual == pytest.approx(expected, rel=1.e-10, abs=1.e-10)
