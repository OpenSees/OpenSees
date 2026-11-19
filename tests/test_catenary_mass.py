"""Clough-style cable mass must conserve mass and be node-order independent."""

import math

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


def _model(reverse=False, height=30., rho=1., mass_type=2, free_end=False):
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0., 0., 0.)
    ops.node(2, 100., 0., height)
    ops.fix(1, 1, 1, 1)
    if not free_end:
        ops.fix(2, 1, 1, 1)
    ends = (2, 1) if reverse else (1, 2)
    ops.element("CatenaryCable", 1, *ends, -10., 2.e11, 1.e-3,
                110., 0., 0., rho, 1.e-9, 10, mass_type)
    return ends


def _nodal_masses(acceleration, reverse=False, height=30., rho=1., mass_type=2):
    ops.wipe()
    _model(reverse, height, rho, mass_type)
    ops.reactions()
    static = [ops.nodeReaction(tag) for tag in (1, 2)]
    for tag in (1, 2):
        for dof, value in enumerate(acceleration, start=1):
            ops.setNodeAccel(tag, dof, value, "-commit")
    ops.reactions("-dynamic")
    return [[(ops.nodeReaction(tag, dof+1) - static[tag-1][dof]) / value
             for dof, value in enumerate(acceleration)] for tag in (1, 2)]


@pytest.mark.parametrize("height", [-30., 0., 30.])
@pytest.mark.parametrize("rho", [0., 1., 2.5])
@pytest.mark.parametrize("mass_type", [0, 2, 3])
def test_rigid_acceleration_conserves_mass_and_node_order(height, rho, mass_type):
    acceleration = (1., -2., 3.)
    forward = _nodal_masses(acceleration, height=height, rho=rho, mass_type=mass_type)
    backward = _nodal_masses(acceleration, reverse=True, height=height,
                              rho=rho, mass_type=mass_type)
    for dof in range(3):
        assert forward[0][dof] + forward[1][dof] == pytest.approx(110. * rho, abs=1.e-8)
    for actual, expected in zip(backward, forward):
        assert actual == pytest.approx(expected, abs=1.e-8)


@pytest.mark.parametrize("reverse", [False, True])
def test_each_end_mass_uses_its_own_tension(reverse):
    ends = _model(reverse=reverse)
    forces = ops.eleForce(1)
    tensions = [math.sqrt(sum(f*f for f in forces[i:i+3])) for i in (0, 3)]
    expected = {tag: 110. * tension / sum(tensions) for tag, tension in zip(ends, tensions)}
    masses = _nodal_masses((1., 2., -3.), reverse=reverse)
    for tag in (1, 2):
        assert masses[tag-1] == pytest.approx([expected[tag]] * 3, abs=1.e-8)


@pytest.mark.parametrize("query_force_first", [False, True])
def test_clough_mass_assembly(query_force_first):
    _model(free_end=True)
    if query_force_first:
        ops.eleForce(1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("GimmeMCK", 1., 0., 0.)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    assert ops.analyze(1, 0.) == 0
    mass = ops.printA("-ret")
    force = ops.eleForce(1)
    tensions = [math.sqrt(sum(f*f for f in force[i:i+3])) for i in (0, 3)]
    nodal_mass = 110. * tensions[1] / sum(tensions)
    assert mass == pytest.approx([nodal_mass if i == j else 0.
                                 for i in range(3) for j in range(3)], abs=1.e-8)
