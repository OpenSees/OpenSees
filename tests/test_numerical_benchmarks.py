"""Small analytical benchmarks shared by every supported CI platform."""

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


def _static():
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()


@pytest.mark.parametrize("load", [-1000., 1000.])
def test_cantilever_displacement_reaction_and_strain_energy(load):
    length, modulus, inertia = 3., 2.e11, 8.e-6
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0., 0.)
    ops.node(2, length, 0.)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.element("elasticBeamColumn", 1, 1, 2, .01, modulus, inertia, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0., load, 0.)
    _static()
    displacement = ops.nodeDisp(2, 2)
    assert displacement == pytest.approx(load * length**3 / (3. * modulus * inertia))
    assert ops.nodeDisp(2, 3) == pytest.approx(load * length**2 / (2. * modulus * inertia))
    assert ops.nodeReaction(1) == pytest.approx([0., -load, -load * length], abs=1.e-9)
    # Independently integrated M(x)^2/(2EI), compared with nodal work.
    energy = load**2 * length**3 / (6. * modulus * inertia)
    assert .5 * load * displacement == pytest.approx(energy)


@pytest.mark.parametrize("consistent", [False, True])
def test_truss_mass_assembly_and_rigid_translation(consistent):
    length, density = 2., 3.
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.)
    ops.node(2, length)
    ops.uniaxialMaterial("Elastic", 1, 1000.)
    ops.element("truss", 1, 1, 2, 1., 1, "-rho", density, "-cMass", int(consistent))
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("GimmeMCK", 1., 0., 0.)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    assert ops.analyze(1, 0.) == 0
    mass = ops.printA("-ret")
    expected = [2., 1., 1., 2.] if consistent else [3., 0., 0., 3.]
    assert mass == pytest.approx(expected)
    assert sum(mass) == pytest.approx(density * length)


def test_single_dof_frequency_and_undamped_energy():
    stiffness, mass, displacement = 100., 2., .01
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.)
    ops.node(2, 0.)
    ops.fix(1, 1)
    ops.mass(2, mass)
    ops.uniaxialMaterial("Elastic", 1, stiffness)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, stiffness * displacement)
    _static()
    assert ops.eigen("-fullGenLapack", 1)[0] == pytest.approx(stiffness / mass)
    ops.remove("loadPattern", 1)
    ops.wipeAnalysis()
    ops.setTime(0.)
    # Set equilibrium initial acceleration for unforced vibration.
    ops.setNodeAccel(2, 1, -stiffness * displacement / mass, "-commit")
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", .5, .25)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    energy = .5 * stiffness * displacement**2
    dt = 2. * math.pi / math.sqrt(stiffness / mass) / 200.
    for _ in range(200):
        assert ops.analyze(1, dt) == 0
        kinetic = .5 * mass * ops.nodeVel(2, 1)**2
        strain = .5 * stiffness * ops.nodeDisp(2, 1)**2
        assert kinetic + strain == pytest.approx(energy, rel=1.e-8)
