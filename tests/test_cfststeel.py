import math

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

PARAMETERS = [206000., 305., .02, -342.07, 26.58, -74.87, 3.16,
              4.07, 6.38, 3.12, 2.36, 5.82, 1., 1., 1., 1., 1.,
              116.39, 11.81, 1, 4467.68, 38.39]

def test_cfststeel_elastic_response():
    ops.wipe()
    try:
        ops.uniaxialMaterial('CFSTsteel', 1, *PARAMETERS)
        ops.testUniaxialMaterial(1)
        ops.setStrain(1.e-6)
        assert math.isclose(ops.getStress(), .206, rel_tol=1.e-10)
        assert math.isclose(ops.getTangent(), 206000., rel_tol=1.e-10)
    finally:
        ops.wipe()

def test_cfststeel_cyclic_response_is_finite():
    ops.wipe()
    try:
        ops.uniaxialMaterial('CFSTsteel', 1, *PARAMETERS)
        ops.testUniaxialMaterial(1)
        last = 0.
        for peak in [.003, -.004, .01, -.02, .015, -.03, .02, -.04]:
            for j in range(1, 101):
                ops.setStrain(last + (peak-last)*j/100.)
                assert math.isfinite(ops.getStress())
                assert math.isfinite(ops.getTangent())
            last = peak
    finally:
        ops.wipe()
