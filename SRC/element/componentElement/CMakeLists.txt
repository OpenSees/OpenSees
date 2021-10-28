#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Element_componentElement OBJECT)

target_sources(OPS_Element_componentElement
    PRIVATE
        ComponentElement2d.cpp
    PUBLIC
        ComponentElement2d.h
)
target_include_directories(OPS_Element_componentElement PUBLIC ${CMAKE_CURRENT_LIST_DIR})

