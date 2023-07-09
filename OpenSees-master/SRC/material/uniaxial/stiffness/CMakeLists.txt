#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
        ConstantStiffnessDegradation.cpp
        DuctilityStiffnessDegradation.cpp
        EnergyStiffnessDegradation.cpp
        PincheiraStiffnessDegradation.cpp
        StiffnessDegradation.cpp
    PUBLIC
        ConstantStiffnessDegradation.h
        DuctilityStiffnessDegradation.h
        EnergyStiffnessDegradation.h
        PincheiraStiffnessDegradation.h
        StiffnessDegradation.h
)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

