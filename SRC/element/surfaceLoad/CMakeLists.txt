#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        SurfaceLoad.cpp
        TriSurfaceLoad.cpp
    PUBLIC
        SurfaceLoad.h
        TriSurfaceLoad.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

