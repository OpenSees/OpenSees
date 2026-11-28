# PDeltaExample2D.py
#
# Example of the PDelta documentation page (geomTransf/PDelta.rst):
# moment amplification and implied buckling load of a cantilever column,
# 1, 2, 4 and 8 elements. Expected output:
#   exact amplification = 1.986
#   elements 1 / 2 / 4 / 8: amplification 1.698 / 1.893 / 1.961 / 1.980,
#   implied buckling load / Pcr 1.216 / 1.060 / 1.020 / 1.010
#
# cantilever column under axial load P = 0.5 Pcr plus a small lateral tip load H,
# divided into 1, 2, 4 and 8 elements
import math
import openseespy.opensees as ops

L, E, A, I = 10.0, 1.0e4, 100.0, 1.0
Pcr = math.pi**2*E*I/(4*L**2)            # Euler load of a cantilever
P = 0.5*Pcr
H = 1.0e-3

# exact amplification of the lateral tip displacement
u = math.pi/2*math.sqrt(P/Pcr)
print(f"exact amplification = {3*(math.tan(u) - u)/u**3:.3f}\n")
print("elements   amplification   implied buckling load / Pcr")

for ne in [1, 2, 4, 8]:
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    for i in range(ne + 1):
        ops.node(i + 1, 0.0, i*L/ne)
    ops.fix(1, 1, 1, 1)

    ops.geomTransf('PDelta', 1)
    for i in range(1, ne + 1):
        ops.element('elasticBeamColumn', i, i, i + 1, A, E, I, 1)
    tip = ne + 1

    # 1) axial load, applied in 10 steps and then held constant
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(tip, 0.0, -P, 0.0)
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.test('NormDispIncr', 1.0e-12, 20)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)
    ops.analysis('Static')
    ops.analyze(10)
    ops.loadConst('-time', 0.0)

    # 2) small lateral load at the tip
    ops.pattern('Plain', 2, 1)
    ops.load(tip, H, 0.0, 0.0)
    ops.integrator('LoadControl', 1.0)
    ops.analyze(1)

    # amplification over the first-order displacement H L^3 / 3EI
    amp = ops.nodeDisp(tip, 1)/(H*L**3/(3*E*I))
    print(f"{ne:8d}   {amp:13.3f}   {(P/Pcr)/(1.0 - 1.0/amp):27.3f}")
