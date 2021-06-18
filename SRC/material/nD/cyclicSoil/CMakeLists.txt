target_sources(material
    PRIVATE
	MultiaxialCyclicPlasticity.cpp
	MultiaxialCyclicPlasticityPlaneStrain.cpp
	MultiaxialCyclicPlasticity3D.cpp
	MultiaxialCyclicPlasticityAxiSymm.cpp
    PUBLIC
	MultiaxialCyclicPlasticity.h
	MultiaxialCyclicPlasticityPlaneStrain.h
	MultiaxialCyclicPlasticity3D.h
	MultiaxialCyclicPlasticityAxiSymm.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))
