#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
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

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
