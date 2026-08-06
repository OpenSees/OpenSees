try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


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
