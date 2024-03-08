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
    ElasticMembraneSection.cpp
    PUBLIC
    ReinforcedConcreteLayeredMembraneSection.h
    LayeredMembraneSection.h
    ElasticMembraneSection.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

