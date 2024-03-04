#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
    ReinforcedConcreteLayeredMembraneSection.cpp
    LayeredMembraneSection.cpp
    PUBLIC
    ReinforcedConcreteLayeredMembraneSection.h
    LayeredMembraneSection.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

