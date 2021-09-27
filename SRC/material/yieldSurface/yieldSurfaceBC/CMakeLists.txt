#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material_YieldSurface
    PRIVATE
        Attalla2D.cpp
        ElTawil2D.cpp
        ElTawil2DUnSym.cpp
        Hajjar2D.cpp
        NullYS2D.cpp
        Orbison2D.cpp
        YieldSurface_BC2D.cpp
        YieldSurface_BC.cpp
    PUBLIC
        Attalla2D.h
        ElTawil2D.h
        ElTawil2DUnSym.h
        Hajjar2D.h
        NullYS2D.h
        Orbison2D.h
        YieldSurface_BC2D.h
        YieldSurface_BC.h
)
target_include_directories(OPS_Material_YieldSurface PUBLIC ${CMAKE_CURRENT_LIST_DIR})

