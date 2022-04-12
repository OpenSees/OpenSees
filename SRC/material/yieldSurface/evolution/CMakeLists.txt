#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material_YieldSurface
    PRIVATE
        BkStressLimSurface2D.cpp
        BoundingSurface2D.cpp
        CombinedIsoKin2D01.cpp
        CombinedIsoKin2D02.cpp
        Isotropic2D01.cpp
        Kinematic2D01.cpp
        Kinematic2D02.cpp
        NullEvolution.cpp
        PeakOriented2D01.cpp
        PeakOriented2D02.cpp
        PlasticHardening2D.cpp
        YS_Evolution2D.cpp
        YS_Evolution.cpp
    PUBLIC
        BkStressLimSurface2D.h
        BoundingSurface2D.h
        CombinedIsoKin2D01.h
        CombinedIsoKin2D02.h
        Isotropic2D01.h
        Kinematic2D01.h
        Kinematic2D02.h
        NullEvolution.h
        PeakOriented2D01.h
        PeakOriented2D02.h
        PlasticHardening2D.h
        YS_Evolution2D.h
        YS_Evolution.h
)
target_include_directories(OPS_Material_YieldSurface PUBLIC ${CMAKE_CURRENT_LIST_DIR})

