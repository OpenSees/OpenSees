#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
    StrengthDegradation.cpp
    ACIStrengthDegradation.cpp 
    ConstantStrengthDegradation.cpp
    DuctilityStrengthDegradation.cpp 
    EnergyStrengthDegradation.cpp
    PetrangeliStrengthDegradation.cpp 
    SectionStrengthDegradation.cpp
    PUBLIC
    StrengthDegradation.h
    ACIStrengthDegradation.h 
    ConstantStrengthDegradation.h
    DuctilityStrengthDegradation.h 
    EnergyStrengthDegradation.h
    PetrangeliStrengthDegradation.h 
    SectionStrengthDegradation.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})


