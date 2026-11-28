# Bathe3D.py
#
# OpenSeesPy version of the Bathe & Bolourchi example in the
# CorotationalBattini documentation (see Bathe3D.tcl for the full Tcl test,
# which also compares the original Corotational transformation).
#
# A cantilever curved in the x-z plane (radius 100, 45 degrees, 8 elements,
# unit square section, E = 1e7) is loaded at the tip by a force of 600 normal
# to its plane, applied in a single load step.
#
# Reference (CorotationalBattini / xara Corotational03, Euler-Bernoulli):
#   tip (ux, uy, uz) = (-23.7823, 53.6463, -13.6962), 12 Newton iterations
#
# Source: Bathe, K.-J. and Bolourchi, S. (1979), IJNME 14(7):961-986; used
# (Table 8) by Perez, C.M. and Filippou, F.C. (2024) "On nonlinear geometric
# transformations of finite elements", IJNME 125(17), e7506.
#
# usage: python Bathe3D.py

import openseespy.opensees as ops
from math import pi, sin, cos

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)
ne, R = 8, 100.0
for n in range(1, ne + 2):
    a = (n - 1) * pi / 4.0 / ne
    ops.node(n, R * sin(a), 0.0, R * (1 - cos(a)))
ops.fix(1, 1, 1, 1, 1, 1, 1)
E, G, A = 1.0e7, 0.5e7, 1.0
I = 1.0 / 12.0
J = 2 * I * 0.8436
ops.geomTransf('CorotationalBattini', 1, 0, 1, 0)
for e in range(1, ne + 1):
    ops.element('elasticBeamColumn', e, e, e + 1, A, E, G, J, I, I, 1)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(ne + 1, 0.0, 600.0, 0.0, 0.0, 0.0, 0.0)
ops.system('FullGeneral'); ops.numberer('Plain'); ops.constraints('Plain')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0)
ops.analysis('Static')
ops.analyze(1)
print(ops.nodeDisp(ne + 1))   # [-23.78, 53.65, -13.70, ...]
