#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_catenaryCable OBJECT)

target_sources(OPS_Element
    PRIVATE
        CatenaryCable.cpp
    PUBLIC
        CatenaryCable.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

