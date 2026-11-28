# Rollup2D.py
#
# OpenSeesPy version of the 2d roll-up example in the Corotational
# documentation (see Rollup2D.tcl for the full Tcl test).
#
# A cantilever of length 10 (EI = 100, 10 elements) is loaded by an end
# moment M = 2*pi*EI/L. The exact solution is a full circle with the tip back
# at the support: tip displacement (-10, 0) and tip rotation 2*pi.
#
# Source: 2d version of the cantilever of Section 4.2.1 of Perez, C.M. and
# Filippou, F.C. (2024) "On nonlinear geometric transformations of finite
# elements", IJNME 125(17), e7506.
#
# usage: python Rollup2D.py

import openseespy.opensees as ops
from math import pi

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
L, E, A, I, ne = 10.0, 1.0e4, 1.0, 1.0e-2, 10
for n in range(1, ne + 2):
    ops.node(n, (n - 1) * L / ne, 0.0)
ops.fix(1, 1, 1, 1)
ops.geomTransf('Corotational', 1)
for e in range(1, ne + 1):
    ops.element('elasticBeamColumn', e, e, e + 1, A, E, I, 1)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(ne + 1, 0.0, 0.0, 2 * pi * E * I / L)
ops.system('FullGeneral'); ops.numberer('Plain'); ops.constraints('Plain')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 3.0)
ops.analysis('Static')
ops.analyze(3)
print(ops.nodeDisp(ne + 1))   # [-10.0, 0.0, 6.283185...]
