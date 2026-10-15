"""Ground-motion inertia load of the elastic beams with element mass.

Under UniformExcitation, addInertiaLoadToUnbalance stores the inertia load
-M*R*ag in the element load vector. getResistingForce() subtracts it when
rho != 0, and getResistingForceIncInertia() calls getResistingForce(), so the
load vector must not be subtracted a second time there. ElasticBeam3d is the
control: it always subtracted the load once.

Reference values (ag = 2.0, rho = 0.8 per unit length):
  rigid-body probe  a beam free only in the excitation direction has no
                    stiffness in rigid translation, so every node's relative
                    acceleration is exactly -ag = -2.0 (double count: -4.0).
  base reaction     all nodes fixed in the excitation direction, lumped mass:
                    the dynamic reactions in that direction sum to
                    rho*L*ag = 0.8*3.0*2.0 = 4.8 (double count: 9.6).
  differential      a cantilever with element -mass (lumped) reproduces, step
                    for step, the same cantilever carrying the identical lumped
                    masses as nodal masses (double count: exactly 2x).
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


AG, E, G, A, I, J, RHO = 2.0, 2.0e5, 8.0e4, 0.02, 1.5e-4, 3.0e-4, 0.8
AV = 0.8 * A
NEL, L_BEAM = 2, 3.0

ELEMENTS = [("elasticBeamColumn", 2), ("elasticBeamColumn", 3),
            ("ElasticTimoshenkoBeam", 2), ("ElasticTimoshenkoBeam", 3)]
IDS = ["ElasticBeam2d", "ElasticBeam3d(control)",
       "ElasticTimoshenkoBeam2d", "ElasticTimoshenkoBeam3d"]


@pytest.fixture(autouse=True)
def clean_domain():
    ops.wipe()
    yield
    ops.wipe()


def _model(dim):
    ops.wipe()
    if dim == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.geomTransf("Linear", 1)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)


def _section(name, dim):
    if name == "elasticBeamColumn":
        return (A, E, I) if dim == 2 else (A, E, G, J, I, I)
    # ElasticTimoshenkoBeam: 2D E G A Iz Avy ; 3D E G A J Iz Iy Avy Avz
    return (E, G, A, I, AV) if dim == 2 else (E, G, A, J, I, I, AV, AV)


def _element(name, dim, tag, i, j, element_mass, cmass=False):
    args = [name, tag, i, j, *_section(name, dim), 1]
    if element_mass:
        args += ["-mass", RHO]
    if cmass:
        args.append("-cMass")
    ops.element(*args)


def _nodes(dim, nel, length):
    h = length / nel
    for k in range(nel + 1):
        ops.node(k + 1, *((k * h, 0.0) if dim == 2 else (k * h, 0.0, 0.0)))


def _ground_motion_y():
    ops.timeSeries("Constant", 1, "-factor", AG)
    ops.pattern("UniformExcitation", 1, 2, "-accel", 1)


def _transient():
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.analysis("Transient")


@pytest.mark.parametrize("cmass", [False, True], ids=["lumped", "consistent"])
@pytest.mark.parametrize("name,dim", ELEMENTS, ids=IDS)
def test_rigid_body_relative_accel_is_minus_ag(name, dim, cmass):
    _model(dim)
    _nodes(dim, NEL, L_BEAM)
    for k in range(NEL + 1):
        ops.fix(k + 1, *((1, 0, 1) if dim == 2 else (1, 0, 1, 1, 1, 1)))
    for e in range(1, NEL + 1):
        _element(name, dim, e, e, e + 1, element_mass=True, cmass=cmass)
    _ground_motion_y()
    _transient()
    assert ops.analyze(1, 0.01) == 0
    acc = [ops.nodeAccel(k + 1, 2) for k in range(NEL + 1)]
    worst = max(abs(a + AG) for a in acc)
    assert worst <= 1e-9 * AG, (
        f"relative y-acceleration {min(acc):+.6f}..{max(acc):+.6f}, expected {-AG}")


@pytest.mark.parametrize("name,dim", ELEMENTS, ids=IDS)
def test_dynamic_base_reaction_is_total_mass_times_ag(name, dim):
    _model(dim)
    _nodes(dim, NEL, L_BEAM)
    ops.fix(1, *([1] * (3 if dim == 2 else 6)))
    for k in range(2, NEL + 2):
        ops.fix(k, *((0, 1, 0) if dim == 2 else (0, 1, 0, 0, 0, 0)))
    for e in range(1, NEL + 1):
        _element(name, dim, e, e, e + 1, element_mass=True)
    _ground_motion_y()
    _transient()
    assert ops.analyze(1, 0.01) == 0
    ops.reactions("-dynamic")
    total = sum(ops.nodeReaction(k + 1, 2) for k in range(NEL + 1))
    expected = RHO * L_BEAM * AG
    assert abs(total - expected) <= 1e-9 * expected, (
        f"sum of dynamic y-reactions {total:.6f}, expected rho*L*ag = {expected:.6f}")


def _cantilever(name, dim, element_mass, nel=4):
    _model(dim)
    _nodes(dim, nel, L_BEAM)
    h = L_BEAM / nel
    ops.fix(1, *([1] * (3 if dim == 2 else 6)))
    for e in range(1, nel + 1):
        _element(name, dim, e, e, e + 1, element_mass=element_mass)
    if not element_mass:
        for k in range(2, nel + 2):
            m = RHO * h * (0.5 if k == nel + 1 else 1.0)
            ops.mass(k, *([m, m, 0.0] if dim == 2 else [m, m, m, 0.0, 0.0, 0.0]))
    _ground_motion_y()
    _transient()
    w1 = math.sqrt(ops.eigen("-fullGenLapack", 1)[0])
    nstep = 100
    dt = 1.5 * (2.0 * math.pi / w1) / nstep
    out = []
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0
        out.append(ops.nodeDisp(nel + 1, 2))
    return out


@pytest.mark.parametrize("name,dim", ELEMENTS, ids=IDS)
def test_element_mass_equals_nodal_mass_under_ground_motion(name, dim):
    elem = _cantilever(name, dim, element_mass=True)
    nodal = _cantilever(name, dim, element_mass=False)
    scale = max(abs(u) for u in nodal)
    assert scale > 0.0
    worst = max(abs(a - b) for a, b in zip(elem, nodal))
    assert worst <= 1e-9 * scale, (
        f"element -mass departs from the same masses as nodal masses by "
        f"{worst / scale:.3e} of peak")
