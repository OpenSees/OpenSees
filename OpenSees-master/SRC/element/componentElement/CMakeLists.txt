#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_componentElement OBJECT)

#target_sources(OPS_Element_componentElement
target_sources(OPS_Element
    PRIVATE
        ComponentElement2d.cpp
        ComponentElement3d.cpp    
    PUBLIC
        ComponentElement2d.h
        ComponentElement3d.h
)
#target_include_directories(OPS_Element_componentElement PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

