import os
import sys
TEST_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"
INTERPRETER_PATH = TEST_DIR + "../interpreter/"
sys.path.append(INTERPRETER_PATH)

import opensees as opy


def test_recorder_time_step_can_handle_fp_precision():
    import tempfile
    opy.model('basic', '-ndm', 2, '-ndf', 3)
    opy.node(1, 0.0, 0.0)
    opy.node(2, 0.0, 5.0)
    opy.fix(2, 0, 1, 0)
    opy.fix(1, 1, 1, 1)
    opy.equalDOF(2, 1, 2)
    opy.mass(2, 1.0, 0.0, 0.0)
    opy.geomTransf('Linear', 1, '-jntOffset')
    opy.element('elasticBeamColumn', 1, 1, 2, 1.0, 1e+06, 0.00164493, 1)
    opy.timeSeries('Path', 1, '-dt', 0.1, '-values', 0.0, -0.001, 0.001, -0.015, 0.033, 0.105, 0.18)
    opy.pattern('UniformExcitation', 1, 1, '-accel', 1)
    opy.rayleigh(0.0, 0.0159155, 0.0, 0.0)
    opy.wipeAnalysis()
    opy.algorithm('Newton')
    opy.system('SparseSYM')
    opy.numberer('RCM')
    opy.constraints('Transformation')
    opy.integrator('Newmark', 0.5, 0.25)
    opy.analysis('Transient')
    opy.test('EnergyIncr', 1e-07, 10, 0, 2)
    node_rec_ffp = tempfile.NamedTemporaryFile(delete=False).name
    ele_rec_ffp = tempfile.NamedTemporaryFile(delete=False).name
    rdt = 0.01
    adt = 0.001
    opy.recorder('Node', '-file', node_rec_ffp, '-precision', 16, '-dT', rdt, '-rTolDt', 0.00001, '-time', '-node', 1, '-dof', 1, 'accel')
    opy.recorder('Element', '-file', ele_rec_ffp, '-precision', 16, '-dT', rdt, '-rTolDt', 0.00001, '-time', '-ele', 1, 'force')

    opy.record()
    for i in range(1100):
        opy.analyze(1, adt)
        opy.getTime()
    opy.wipe()

    a = open(node_rec_ffp).read().splitlines()
    for i in range(len(a) - 1):
        dt = float(a[i + 1].split()[0]) - float(a[i].split()[0])
        assert abs(dt - 0.01) < adt * 0.1, (i, dt)
    a = open(ele_rec_ffp).read().splitlines()
    for i in range(len(a) - 1):
        dt = float(a[i + 1].split()[0]) - float(a[i].split()[0])
        assert abs(dt - 0.01) < adt * 0.1, (i, dt)


if __name__ == '__main__':
    test_recorder_time_step_can_handle_fp_precision()
