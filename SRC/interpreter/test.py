import opensees

opensees.wipe()

opensees.uniaxialMaterial("Elastic", 1, 1000.);
opensees.testUniaxialMaterial(1);

for strain in [0.01, 0.02, 0.03, 0.04, 0.05]:
    opensees.setStrain(strain);
    print "strain: " + str(opensees.getStrain()) + " stress: " + str(opensees.getStress()) + " tangent: " + str(opensees.getTangent());

opensees.uniaxialMaterial("Elastic", 2, 1000.);
opensees.uniaxialMaterial("Parallel", 3, 1, 2);
opensees.testUniaxialMaterial(3);

for strain in [0.01, 0.02, 0.03, 0.04, 0.05]:
    opensees.setStrain(strain);
    print "strain: " + str(opensees.getStrain()) + " stress: " + str(opensees.getStress()) + " tangent: " + str(opensees.getTangent());
