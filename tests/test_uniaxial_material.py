import math

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def test_elastic_uniaxial_material_response():
    modulus = 1000.0

    ops.wipe()
    try:
        ops.uniaxialMaterial("Elastic", 1, modulus)
        ops.testUniaxialMaterial(1)

        for strain in (-0.002, 0.0, 0.003):
            ops.setStrain(strain)
            assert math.isclose(ops.getStrain(), strain, rel_tol=0.0, abs_tol=1.0e-15)
            assert math.isclose(ops.getStress(), modulus * strain, rel_tol=0.0, abs_tol=1.0e-12)
            assert math.isclose(ops.getTangent(), modulus, rel_tol=0.0, abs_tol=1.0e-12)
    finally:
        ops.wipe()
