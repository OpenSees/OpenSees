#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
    ConcreteL01.cpp
    ConcreteZ01.cpp
    TendonL01.cpp
    SteelZ01.cpp
    ReinforcedConcretePlaneStress.cpp
    FAReinforcedConcretePlaneStress.cpp
    FAFourSteelRCPlaneStress.cpp
    RAFourSteelRCPlaneStress.cpp
    PrestressedConcretePlaneStress.cpp
    FAPrestressedConcretePlaneStress.cpp
    FAFourSteelPCPlaneStress.cpp
    RAFourSteelPCPlaneStress.cpp
    PUBLIC
    ConcreteL01.h
    ConcreteZ01.h
    TendonL01.h
    SteelZ01.h
    ReinforcedConcretePlaneStress.h
    FAReinforcedConcretePlaneStress.h
    FAFourSteelRCPlaneStress.h
    RAFourSteelRCPlaneStress.h
    PrestressedConcretePlaneStress.h
    FAPrestressedConcretePlaneStress.h
    FAFourSteelPCPlaneStress.h
    RAFourSteelPCPlaneStress.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
