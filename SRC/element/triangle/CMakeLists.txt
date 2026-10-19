#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
    Tri31.cpp
    Tri31UP.cpp
    PUBLIC
    Tri31.h
    Tri31UP.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
