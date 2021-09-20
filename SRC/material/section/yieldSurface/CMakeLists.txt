#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Section_YieldSurface OBJECT)

target_sources(OPS_Section_YieldSurface
    PRIVATE
        SoilFootingSection2d.cpp
        YieldSurfaceSection2d.cpp
        YS_Section2D01.cpp
        YS_Section2D02.cpp
    PUBLIC
        SoilFootingSection2d.h
        YieldSurfaceSection2d.h
        YS_Section2D01.h
        YS_Section2D02.h
)
target_include_directories(OPS_Section_YieldSurface PUBLIC ${CMAKE_CURRENT_LIST_DIR})

