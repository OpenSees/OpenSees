import opensees as ops

ops.wipe()

ops.uniaxialMaterial("Elastic", 1, 1000.);
ops.testUniaxialMaterial(1);

for strain in [0.01, 0.02, 0.03, 0.04, 0.05]:
    ops.setStrain(strain);
    print("strain: ", str(ops.getStrain()), " stress: ", str(ops.getStress()), " tangent: ", str(ops.getTangent()));

ops.uniaxialMaterial("Elastic", 2, 1000.);
ops.uniaxialMaterial("Parallel", 3, 1, 2);
ops.testUniaxialMaterial(3);

for strain in [0.01, 0.02, 0.03, 0.04, 0.05]:
    ops.setStrain(strain);
    print("strain: ", str(ops.getStrain()), " stress: ", str(ops.getStress()), " tangent: ", str(ops.getTangent()));
