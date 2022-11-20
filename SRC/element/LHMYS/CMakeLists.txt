#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Element_LHMYS OBJECT)

target_sources(OPS_Element_LHMYS
    PRIVATE
        BeamColumnwLHNMYS.cpp
        GPYS2d.cpp
    PUBLIC
        BeamColumnwLHNMYS.h
        GPYS2d.h
)
target_include_directories(OPS_Element_LHMYS PUBLIC ${CMAKE_CURRENT_LIST_DIR})

