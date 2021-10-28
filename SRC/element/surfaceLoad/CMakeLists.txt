#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Element_surfaceLoad OBJECT)

target_sources(OPS_Element_surfaceLoad
    PRIVATE
        SurfaceLoad.cpp
        TriSurfaceLoad.cpp
    PUBLIC
        SurfaceLoad.h
        TriSurfaceLoad.h
)
target_include_directories(OPS_Element_surfaceLoad PUBLIC ${CMAKE_CURRENT_LIST_DIR})

