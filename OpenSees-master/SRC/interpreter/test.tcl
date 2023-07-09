
uniaxialMaterial Elastic 1 1000.
testUniaxialMaterial 1

foreach strain {1.0e-3 1.1e-2 1.0} {
    setStrain $strain
    puts "strain: [getStrain] stress: [getStress] tangent: [getTangent]"
}

uniaxialMaterial Elastic 2 1000.
uniaxialMaterial Parallel 3 2 1 1000.
testUniaxialMaterial 3

foreach strain {1.0e-3 1.1e-2 1.0} {
    setStrain $strain
    puts "strain: [getStrain] stress: [getStress] tangent: [getTangent]"
}
