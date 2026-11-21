"""DDM sensitivity of elasticBeamColumn (ElasticBeam2d, ElasticBeam3d) checked against central finite differences.

The finite differences perturb the same element parameters with updateParameter: static cases rerun the analysis
from the initial state, transient cases rebuild the model and set the parameter before rayleigh, because the
committed stiffness used by betaKcomm is stored by the element at commit and is not reset by ops.reset(). Covered: static frames with nodal, point,
uniform, partial uniform and (2D) temperature loads, all moment releases, with and without shear deformation
(3D: also with a shear area in one bending plane only, and a cantilever checked against the closed form);
transient analyses with consistent mass and the Rayleigh damping variants (betaK, betaKinit, betaKcomm),
including the mass density rho."""

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


# ---------------------------------------------------------------- helpers

def _params(spec):
    """spec: list of (element tag, parameter name); returns the parameter tags and current values."""
    tags, vals = [], []
    for k, (ele, name) in enumerate(spec):
        tag = 100 + k
        ops.parameter(tag, "element", ele, name)
        tags.append(tag)
        vals.append(ops.getParamValue(tag))
    return tags, vals


def _compare(g, fd, ele, name, tol, scale_all):
    """Error of the normalized DDM derivative g = v du/dv against the finite differences fd, relative to the
    larger of its own size and 1e-3 of the largest normalized derivative in the model (scale_all): finite
    differences cannot resolve the derivative of a parameter with a very weak effect to a small relative error.
    A parameter on which the response does not depend (e.g. G of a beam released at both ends) must give a
    zero DDM derivative."""
    scale = max(abs(x) for x in fd)
    if scale < 1e-9 * scale_all:
        assert max(abs(x) for x in g) < 1e-9 * scale_all, f"element {ele}, parameter {name}: nonzero DDM derivative"
        return
    err = max(abs(a - b) for a, b in zip(g, fd)) / max(scale, 1e-3 * scale_all)
    assert err < tol, f"element {ele}, parameter {name}: relative error {err:.2e}"


def _static_response(nodes, ndf):
    ops.reset()
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(n, d) for n in nodes for d in range(1, ndf + 1)]


def _check_static(nodes, ndf, spec, rel=1e-4, tol=1e-6):
    tags, vals = _params(spec)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    ops.sensitivityAlgorithm("-computeByCommand")
    _static_response(nodes, ndf)
    ops.computeGradients()
    ddm = [[ops.sensNodeDisp(n, d, t) for n in nodes for d in range(1, ndf + 1)] for t in tags]
    scale_all = max(abs(v * x) for t, v, row in zip(tags, vals, ddm) for x in row)   # largest v * du/dv
    for t, v, g, (ele, name) in zip(tags, vals, ddm, spec):
        h = rel * abs(v)
        ops.updateParameter(t, v + h)
        up = _static_response(nodes, ndf)
        ops.updateParameter(t, v - h)
        dn = _static_response(nodes, ndf)
        ops.updateParameter(t, v)
        fd = [(a - b) / (2 * h) for a, b in zip(up, dn)]
        _compare([v * x for x in g], [v * x for x in fd], ele, name, tol, scale_all)


def _transient_run(build, rayleigh, nodes, dof, spec, nsteps, dt, set_value=None, ddm=False):
    """Builds the model, optionally sets one parameter (index, value) before rayleigh, and runs the analysis.
    Returns the displacements and accelerations of the nodes in direction dof at each step and, with ddm=True,
    their DDM derivatives."""
    ops.wipe()
    build()
    tags, vals = _params(spec)
    if set_value is not None:
        ops.updateParameter(tags[set_value[0]], set_value[1])
    ops.rayleigh(*rayleigh)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    if ddm:
        ops.sensitivityAlgorithm("-computeByCommand")
    resp, sens = [], []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        resp.append([ops.nodeDisp(n, dof) for n in nodes] + [ops.nodeAccel(n, dof) for n in nodes])
        if ddm:
            ops.computeGradients()
            sens.append([[ops.sensNodeDisp(n, dof, t) for n in nodes] + [ops.sensNodeAccel(n, dof, t) for n in nodes]
                         for t in tags])
    return resp, sens, vals


def _check_transient(build, rayleigh, nodes, dof, spec, nsteps=40, dt=0.002, rel=1e-5, tol=1e-5):
    _, ddm, vals = _transient_run(build, rayleigh, nodes, dof, spec, nsteps, dt, ddm=True)
    n = len(nodes)
    for k, (v, (ele, name)) in enumerate(zip(vals, spec)):
        h = rel * abs(v)
        up, _, _ = _transient_run(build, rayleigh, nodes, dof, spec, nsteps, dt, set_value=(k, v + h))
        dn, _, _ = _transient_run(build, rayleigh, nodes, dof, spec, nsteps, dt, set_value=(k, v - h))
        fd = [[(a - b) / (2 * h) for a, b in zip(u, d)] for u, d in zip(up, dn)]
        g = [step[k] for step in ddm]
        for part in (slice(0, n), slice(n, 2 * n)):          # displacements, accelerations
            scale = max(abs(x) for row in fd for x in row[part])
            assert scale > 0, f"response does not depend on {name} of element {ele}"
            err = max(abs(a - b) for rg, rf in zip(g, fd) for a, b in zip(rg[part], rf[part])) / scale
            assert err < tol, f"element {ele}, parameter {name}: relative error {err:.2e}"


# ---------------------------------------------------------------- 2D

def _frame2d(release=0, shear=False):
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for n, (x, y) in enumerate([(0., 0.), (0., 3.), (4., 3.5), (4., 0.)], 1):
        ops.node(n, x, y)
    ops.fix(1, 1, 1, 1)
    ops.fix(4, 1, 1, 0)
    ops.geomTransf("Linear", 1)
    E, G = 2.0e8, 8.0e7
    sec = [(0.012, 1.2e-4, 0.006), (0.010, 0.9e-4, 0.005), (0.011, 1.0e-4, 0.0055)]
    conn = [(1, 2), (2, 3), (3, 4)]
    for e, ((i, j), (A, I, Av)) in enumerate(zip(conn, sec), 1):
        props = [A, E, I, G, Av] if shear else [A, E, I]
        opts = ["-alpha", 1.2e-5, "-depth", 0.3]
        if e == 2:
            opts += ["-release", release]
        ops.element("elasticBeamColumn", e, i, j, *props, 1, *opts)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10., -5., 2.)
    ops.eleLoad("-ele", 1, "-type", "-beamPoint", 3., 0.4, 1.)           # Py, xL, Px
    ops.eleLoad("-ele", 2, "-type", "-beamUniform", -6., 1.)            # Wy, Wx
    ops.eleLoad("-ele", 3, "-type", "-beamUniform", 2., 0.5, 0.2, 0.7)  # partial: Wy, Wx, aL, bL
    ops.eleLoad("-ele", 2, "-type", "-beamTemp", 15., -5.)              # Ttop, Tbottom
    names = ["E", "A", "I"] + (["G", "Av"] if shear else [])
    return [(e, p) for e in (1, 2, 3) for p in names]


@pytest.mark.parametrize("shear", [False, True])
@pytest.mark.parametrize("release", [0, 1, 2, 3])
def test_static_2d(release, shear):
    spec = _frame2d(release, shear)
    _check_static([2, 3, 4], 3, spec)


def _beam2d_dynamic():
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    n_el, L = 4, 8.0
    for n in range(n_el + 1):
        ops.node(n + 1, n * L / n_el, 0.)
    ops.fix(1, 1, 1, 0)
    ops.fix(n_el + 1, 0, 1, 0)
    ops.geomTransf("Linear", 1)
    for e in range(1, n_el + 1):
        ops.element("elasticBeamColumn", e, e, e + 1, 0.01, 2.0e8, 1.0e-4, 1, "-mass", 80., "-cMass")
    ops.timeSeries("Trig", 1, 0., 10., 0.05)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0., -10., 0.)


@pytest.mark.parametrize("rayleigh", [(0., 0., 0., 0.), (0.5, 1e-3, 0., 0.), (0.5, 0., 1e-3, 0.), (0.5, 0., 0., 1e-3)],
                         ids=["none", "betaK", "betaKinit", "betaKcomm"])
def test_transient_2d(rayleigh):
    spec = [(e, p) for e in (1, 2) for p in ("E", "I", "rho")]
    _check_transient(_beam2d_dynamic, rayleigh, [2, 3, 4], 2, spec)


# ---------------------------------------------------------------- 3D

def _frame3d(releasez=0, releasey=0, shear=False):
    """shear: False (no shear deformation), True (Avy and Avz), "y" (only Avy, Avz = 0) or "z" (only Avz)."""
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for n, xyz in enumerate([(0., 0., 0.), (0., 0., 3.), (4., 1., 3.5), (4., 1., 0.)], 1):
        ops.node(n, *xyz)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(4, 1, 1, 1, 1, 0, 0)
    ops.geomTransf("Linear", 1, 1., 0., 0.)    # columns (along z)
    ops.geomTransf("Linear", 2, 0., 0., 1.)    # beam
    E, G = 2.0e8, 8.0e7
    sec = [(0.012, 2.0e-5, 1.0e-4, 1.2e-4, 0.006, 0.005), (0.010, 1.5e-5, 0.8e-4, 0.9e-4, 0.005, 0.004),
           (0.011, 1.8e-5, 0.9e-4, 1.0e-4, 0.0055, 0.0045)]
    conn, transf = [(1, 2), (2, 3), (3, 4)], [1, 2, 1]
    for e, ((i, j), (A, J, Iy, Iz, Avy, Avz), tr) in enumerate(zip(conn, sec, transf), 1):
        if shear == "y":
            Avz = 0.
        elif shear == "z":
            Avy = 0.
        props = [A, E, G, J, Iy, Iz] + ([Avy, Avz] if shear else [])
        opts = ["-releasez", releasez, "-releasey", releasey] if e == 2 else []
        ops.element("elasticBeamColumn", e, i, j, *props, tr, *opts)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10., -4., 3., 1., 2., -1.)
    ops.eleLoad("-ele", 1, "-type", "-beamPoint", 3., -2., 0.4, 1.)            # Py, Pz, xL, Px
    ops.eleLoad("-ele", 2, "-type", "-beamUniform", -6., 2., 1.)               # Wy, Wz, Wx
    ops.eleLoad("-ele", 3, "-type", "-beamUniform", 2., -1., 0.5, 0.2, 0.7)    # partial: Wy, Wz, Wx, aL, bL
    names = ["E", "A", "Iz", "Iy", "G", "J"] + {False: [], True: ["Avy", "Avz"], "y": ["Avy"], "z": ["Avz"]}[shear]
    return [(e, p) for e in (1, 2, 3) for p in names]


@pytest.mark.parametrize("shear", [False, True])
@pytest.mark.parametrize("releasez,releasey", [(0, 0), (1, 0), (2, 0), (0, 3), (3, 1)])
def test_static_3d(releasez, releasey, shear):
    spec = _frame3d(releasez, releasey, shear)
    _check_static([2, 3, 4], 6, spec)


@pytest.mark.parametrize("shear", ["y", "z"])
@pytest.mark.parametrize("releasez,releasey", [(0, 0), (1, 1), (2, 2)])
def test_static_3d_one_shear_area(releasez, releasey, shear):
    """Shear deformation in one bending plane only (Avy > 0 and Avz = 0, or the reverse)."""
    spec = _frame3d(releasez, releasey, shear)
    _check_static([2, 3, 4], 6, spec)


@pytest.mark.parametrize("section", [False, True])
@pytest.mark.parametrize("Avy,Avz", [(0.0, 0.0), (0.006, 0.0), (0.0, 0.005), (0.006, 0.005)])
def test_cantilever_3d_shear_areas(Avy, Avz, section):
    """Tip deflections of a cantilever, P L^3/(3 E I) + P L/(G Av) in each bending plane, with Newton iterations,
    so that getResistingForce() and getTangentStiff() must agree; Av = 0 means no shear deformation in that plane.
    section=True builds the element from an Elastic section (with shear: alphaY = Avy/A, alphaZ = Avz/A), whose
    shear areas the element reads from the section response codes."""
    L, A, E, G, J, Iy, Iz, Py, Pz = 3., 0.01, 2.0e8, 8.0e7, 1.0e-4, 1.0e-4, 2.0e-4, 1., 1.
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0., 0., 0.)
    ops.node(2, L, 0., 0.)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0., 0., 1.)
    if not section:
        ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, J, Iy, Iz, Avy, Avz, 1)
    else:
        shear = [Avy / A, Avz / A] if Avy > 0 or Avz > 0 else []
        ops.section("Elastic", 1, E, A, Iz, Iy, G, J, *shear)
        ops.element("elasticBeamColumn", 1, 1, 2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0., Py, Pz, 0., 0., 0.)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormUnbalance", 1e-10, 5)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    uy = Py*L**3/(3*E*Iz) + (Py*L/(G*Avy) if Avy > 0 else 0.)
    uz = Pz*L**3/(3*E*Iy) + (Pz*L/(G*Avz) if Avz > 0 else 0.)
    assert ops.nodeDisp(2, 2) == pytest.approx(uy, rel=1e-10)
    assert ops.nodeDisp(2, 3) == pytest.approx(uz, rel=1e-10)


def _beam3d_dynamic():
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    n_el, L = 4, 8.0
    for n in range(n_el + 1):
        ops.node(n + 1, n * L / n_el, 0., 0.)
    ops.fix(1, 1, 1, 1, 1, 0, 0)
    ops.fix(n_el + 1, 0, 1, 1, 0, 0, 0)
    ops.geomTransf("Linear", 1, 0., 0., 1.)
    for e in range(1, n_el + 1):
        ops.element("elasticBeamColumn", e, e, e + 1, 0.01, 2.0e8, 8.0e7, 2.0e-5, 0.8e-4, 1.0e-4, 1,
                    "-mass", 80., "-cMass")
    ops.timeSeries("Trig", 1, 0., 10., 0.05)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0., -10., 6., 1., 0., 0.)


PARAMS_3D = {2: ("E", "Iz", "rho"), 3: ("E", "Iy", "rho"), 4: ("G", "J", "A", "rho")}   # parameters per checked dof


@pytest.mark.parametrize("rayleigh", [(0., 0., 0., 0.), (0.5, 1e-3, 0., 0.), (0.5, 0., 1e-3, 0.), (0.5, 0., 0., 1e-3)],
                         ids=["none", "betaK", "betaKinit", "betaKcomm"])
@pytest.mark.parametrize("dof", [2, 3, 4], ids=["uy", "uz", "rx"])
def test_transient_3d(rayleigh, dof):
    spec = [(e, p) for e in (1, 2) for p in PARAMS_3D[dof]]
    _check_transient(_beam3d_dynamic, rayleigh, [2, 3, 4], dof, spec)
