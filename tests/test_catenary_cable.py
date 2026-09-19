import math
import xml.etree.ElementTree as ET

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


def _run_catenary_model(ndf):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", ndf)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 10.0, 0.0, 0.0)
    ops.fix(1, *([1] * ndf))
    ops.fix(2, *([0, 1, 1] + [1] * (ndf - 3)))
    ops.element(
        "CatenaryCable",
        1,
        1,
        2,
        1.0,
        2.0e11,
        0.01,
        10.0,
        0.0,
        0.0,
        0.0,
        1.0e-8,
        10,
        0,
    )
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, *([1000.0, 0.0, 0.0] + [0.0] * (ndf - 3)))
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    try:
        assert ops.analyze(1) == 0
        force = ops.eleResponse(1, "force")
        assert len(force) == 2 * ndf
        if ndf == 6:
            assert force[3:6] == [0.0, 0.0, 0.0]
            assert force[9:12] == [0.0, 0.0, 0.0]
    finally:
        ops.wipe()


def test_catenary_cable_supports_three_dof_nodes():
    _run_catenary_model(3)


def test_catenary_cable_supports_six_dof_nodes():
    _run_catenary_model(6)


def _verification_model(ndf):
    # OpenSeesWiki CatenaryCableElement example (x = 30).
    # https://opensees.berkeley.edu/wiki/index.php/CatenaryCableElement
    # Broader reference/equilibrium coverage suggested by apalazzi in PR #1796.
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", ndf)
    ops.node(1, 0.0, 0.0, 90.0)
    ops.node(2, 15.0, 0.0, 40.0)
    ops.node(3, 30.0, 60.0, 30.0)
    ops.fix(1, *([1] * ndf))
    ops.fix(2, *([0, 1, 0] + [1] * (ndf - 3)))
    ops.fix(3, *([1] * ndf))
    for tag, start, end in ((1, 1, 2), (2, 2, 3)):
        ops.element("CatenaryCable", tag, start, end,
                    -1.0e-5, 3.0e7, 1.0, 50.0, 6.5e-6, 100.0,
                    0.0, 1.0e-6, 20, 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 1, 2, "-type", "-beamUniform", 0.0, 0.0, -1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-5, 100)
    ops.integrator("LoadControl", 0.1)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(10) == 0


@pytest.mark.parametrize("ndf", [3, 6])
def test_documented_verification_displacement(ndf):
    _verification_model(ndf)
    assert ops.nodeDisp(2)[:3] == pytest.approx([8.58693, 0.0, 2.82578], abs=1.0e-4)


def test_three_and_six_dof_results_match():
    _verification_model(3)
    displacements = {tag: ops.nodeDisp(tag) for tag in (1, 2, 3)}
    forces = {tag: ops.eleResponse(tag, "force") for tag in (1, 2)}
    _verification_model(6)
    for tag in (1, 2, 3):
        assert ops.nodeDisp(tag)[:3] == pytest.approx(displacements[tag], abs=1.0e-9)
    for tag in (1, 2):
        force = ops.eleResponse(tag, "force")
        assert len(force) == 12
        assert force[:3] + force[6:9] == pytest.approx(forces[tag], rel=1.0e-9, abs=1.0e-8)
        assert force[3:6] + force[9:12] == [0.0] * 6


def _symmetric_cable(ndf, offset=0, weight=-10.0, rho=0.0, mass_type=0):
    ops.model("basic", "-ndm", 3, "-ndf", ndf)
    ops.node(offset + 1, 0.0, float(offset), 0.0)
    ops.node(offset + 2, 100.0, float(offset), 0.0)
    ops.fix(offset + 1, *([1] * ndf))
    ops.fix(offset + 2, *([1] * ndf))
    ops.element("CatenaryCable", offset + 1, offset + 1, offset + 2,
                weight, 2.0e11, 1.0e-3, 102.0, 0.0, 0.0, rho, 1.0e-9, 10, mass_type)


@pytest.mark.parametrize("ndf", [3, 6])
def test_self_weight_equilibrium_and_elastic_catenary_tension(ndf):
    _symmetric_cable(ndf)
    # Domain::addElement solves the cable's internal equilibrium. Query the
    # fully restrained model without invoking a zero-equation linear solver.
    ops.reactions()
    reactions = [ops.nodeReaction(tag) for tag in (1, 2)]
    assert reactions[0][2] == pytest.approx(510.0, abs=1.0e-8)
    assert reactions[1][2] == pytest.approx(510.0, abs=1.0e-8)
    assert reactions[0][0] + reactions[1][0] == pytest.approx(0.0, abs=1.0e-8)

    # Exact span relation for a symmetric elastic cable with uniform load per
    # unstretched length: x = 2 H/w asinh(w L0/(2 H)) + H L0/(EA).
    low, high = 1.0, 2.0e8
    for _ in range(100):
        tension = (low + high) / 2.0
        span = 2.0 * tension / 10.0 * math.asinh(510.0 / tension)
        span += tension * 102.0 / 2.0e8
        if span < 100.0:
            low = tension
        else:
            high = tension
    assert abs(reactions[0][0]) == pytest.approx(tension, rel=1.0e-6)
    if ndf == 6:
        assert reactions[0][3:] + reactions[1][3:] == [0.0] * 6


@pytest.mark.parametrize("order", [(3, 6), (6, 3)])
def test_mixed_dof_elements_keep_separate_force_storage(order):
    for index, ndf in enumerate(order):
        _symmetric_cable(ndf, offset=10 * index, weight=-10.0 * (index + 1))
    expected = {tag: ops.eleResponse(tag, "force") for tag in (1, 11)}
    for tag in (11, 1, 11, 1):
        ndf = order[0 if tag == 1 else 1]
        force = ops.eleResponse(tag, "force")
        assert len(force) == 2 * ndf
        assert force == pytest.approx(expected[tag], abs=1.0e-9)
        assert force[2] + force[ndf + 2] == pytest.approx(1020.0 * (1 if tag == 1 else 2))
        if ndf == 6:
            assert force[3:6] + force[9:12] == [0.0] * 6


@pytest.mark.parametrize("ndf", [3, 6])
@pytest.mark.parametrize("mass_type, nodal_mass", [(0, 51.0), (1, 51.0), (3, 34.0)])
def test_translational_mass_assembly(ndf, mass_type, nodal_mass):
    ops.model("basic", "-ndm", 3, "-ndf", ndf)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 100.0, 0.0, 0.0)
    ops.fix(1, *([1] * ndf))
    # Added nodal inertias make the rotational equations solvable; the cable
    # must contribute zero to these entries and to all off-diagonal entries.
    if ndf == 6:
        ops.mass(2, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0)
    ops.element("CatenaryCable", 1, 1, 2,
                -10.0, 2.0e11, 1.0e-3, 102.0, 0.0, 0.0, 1.0, 1.0e-9, 10, mass_type)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("GimmeMCK", 1.0, 0.0, 0.0)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    assert ops.analyze(1, 0.0) == 0
    diagonal = [nodal_mass] * 3 + ([2.0, 3.0, 4.0] if ndf == 6 else [])
    expected = [diagonal[i] if i == j else 0.0 for i in range(ndf) for j in range(ndf)]
    assert ops.printA("-ret") == pytest.approx(expected, abs=1.0e-10)


@pytest.mark.parametrize("ndf", [3, 6])
@pytest.mark.parametrize("mass_type", [0, 3])
def test_inertia_force_mapping(ndf, mass_type):
    _symmetric_cable(ndf, rho=1.0, mass_type=mass_type)
    ops.reactions()
    static = [ops.nodeReaction(tag) for tag in (1, 2)]
    acceleration = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    for tag in (1, 2):
        values = acceleration[tag - 1] + ([100.0, 200.0, 300.0] if ndf == 6 else [])
        for dof, value in enumerate(values, start=1):
            ops.setNodeAccel(tag, dof, value, "-commit")
    ops.reactions("-dynamic")
    for index, tag in enumerate((1, 2)):
        reaction = ops.nodeReaction(tag)
        if mass_type == 0:
            expected = [51.0 * value for value in acceleration[index]]
        else:
            expected = [17.0 * (2.0 * a + b) for a, b in
                        zip(acceleration[index], acceleration[1 - index])]
        inertia = [reaction[i] - static[index][i] for i in range(3)]
        assert inertia == pytest.approx(expected, abs=1.0e-9)
        if ndf == 6:
            assert reaction[3:] == [0.0] * 3


@pytest.mark.parametrize("ndf", [3, 6])
def test_force_recorder_labels_and_values(ndf, tmp_path):
    _symmetric_cable(ndf)
    force = ops.eleResponse(1, "force")
    output = tmp_path / "cable.xml"
    ops.recorder("Element", "-xml", str(output), "-precision", 16, "-ele", 1, "force")
    ops.record()
    ops.remove("recorders")
    root = ET.parse(output).getroot()
    assert len(root.findall("ElementOutput")) == 1
    assert [item.text for item in root.iter("ResponseType")] == [
        f"f{i}" for i in range(1, 2 * ndf + 1)
    ]
    data = [float(value) for value in root.find("Data").text.split()]
    assert data == pytest.approx(force, rel=1.0e-12, abs=1.0e-12)


def test_mixed_dof_recorder_elements_are_siblings(tmp_path):
    _symmetric_cable(3)
    _symmetric_cable(6, offset=10, weight=-20.0)
    force = ops.eleResponse(1, "force") + ops.eleResponse(11, "force")
    output = tmp_path / "mixed.xml"
    ops.recorder("Element", "-xml", str(output), "-precision", 16,
                 "-ele", 1, 11, "force")
    ops.record()
    ops.remove("recorders")
    root = ET.parse(output).getroot()
    elements = root.findall("ElementOutput")
    assert [element.get("eleTag") for element in elements] == ["1", "11"]
    assert [len(element.findall("ResponseType")) for element in elements] == [6, 12]
    data = [float(value) for value in root.find("Data").text.split()]
    assert data == pytest.approx(force, rel=1.0e-12, abs=1.0e-12)


def test_cable_with_beam_allows_nonzero_node_rotations():
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 100.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 2, 1, 2, 0.01, 2.0e11, 8.0e10,
                1.0e-4, 1.0e-4, 1.0e-4, 1)
    ops.element("CatenaryCable", 1, 1, 2,
                -10.0, 2.0e11, 1.0e-3, 102.0, 0.0, 0.0, 0.0, 1.0e-9, 10, 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 100.0, 0.0, 50.0, 0.0, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-8, 50)
    ops.integrator("LoadControl", 0.1)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(10) == 0
    assert abs(ops.nodeDisp(2, 4)) > 1.0e-6
    force = ops.eleResponse(1, "force")
    assert len(force) == 12
    assert all(math.isfinite(value) for value in force)
    assert force[3:6] + force[9:12] == [0.0] * 6
    ops.reactions()
    assert ops.nodeReaction(1, 4) == pytest.approx(-50.0, abs=1.0e-6)
