#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
    ReinforcedConcreteLayerMembraneSection01.cpp
    ReinforcedConcreteLayerMembraneSection02.cpp
    PUBLIC
    ReinforcedConcreteLayerMembraneSection01.h
    ReinforcedConcreteLayerMembraneSection02.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

