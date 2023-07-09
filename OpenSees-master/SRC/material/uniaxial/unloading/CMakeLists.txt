#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
    PRIVATE
    UnloadingRule.cpp
    ConstantUnloadingRule.cpp 
    EnergyUnloadingRule.cpp
    KarsanUnloadingRule.cpp 
    TakedaUnloadingRule.cpp
    PUBLIC
    UnloadingRule.h
    ConstantUnloadingRule.h 
    EnergyUnloadingRule.h
    KarsanUnloadingRule.h 
    TakedaUnloadingRule.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})


