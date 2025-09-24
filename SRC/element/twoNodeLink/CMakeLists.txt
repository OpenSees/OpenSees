#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
    Inerter.cpp
    LinearElasticSpring.cpp
    TwoNodeLink.cpp
    TwoNodeLinkSection.cpp
    PUBLIC
    Inerter.h
    LinearElasticSpring.h
    TwoNodeLink.h
    TwoNodeLinkSection.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
