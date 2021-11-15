#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        ForceBeamColumn2d.cpp
        ForceBeamColumn3d.cpp
        ForceBeamColumn2dThermal.cpp
        ElasticForceBeamColumn2d.cpp
        ElasticForceBeamColumn3d.cpp
        ElasticForceBeamColumnWarping2d.cpp
	ForceBeamColumnCBDI2d.cpp
        ForceBeamColumnCBDI3d.cpp	
        ForceBeamColumnWarping2d.cpp
        BeamIntegration.cpp
        LobattoBeamIntegration.cpp
        LegendreBeamIntegration.cpp
        RadauBeamIntegration.cpp
        NewtonCotesBeamIntegration.cpp
        HingeMidpointBeamIntegration.cpp
        HingeEndpointBeamIntegration.cpp
        HingeRadauBeamIntegration.cpp
        HingeRadauTwoBeamIntegration.cpp
        UserDefinedBeamIntegration.cpp
        SimpsonBeamIntegration.cpp
        UserDefinedHingeIntegration.cpp
        DistHingeIntegration.cpp
        RegularizedHingeIntegration.cpp
        TrapezoidalBeamIntegration.cpp
        CompositeSimpsonBeamIntegration.cpp
        FixedLocationBeamIntegration.cpp
        LowOrderBeamIntegration.cpp
        MidDistanceBeamIntegration.cpp
        #GaussQBeamIntegration.cpp
	ChebyshevBeamIntegration.cpp
    PUBLIC
        ForceBeamColumn2d.h
        ForceBeamColumn3d.h
        ForceBeamColumn2dThermal.h
        ElasticForceBeamColumn2d.h
        ElasticForceBeamColumn3d.h
        ElasticForceBeamColumnWarping2d.h
        ForceBeamColumnCBDI2d.h
        ForceBeamColumnCBDI3d.h	
        ForceBeamColumnWarping2d.h
        BeamIntegration.h
        LobattoBeamIntegration.h
        LegendreBeamIntegration.h
        RadauBeamIntegration.h
        NewtonCotesBeamIntegration.h
        HingeMidpointBeamIntegration.h
        HingeEndpointBeamIntegration.h
        HingeRadauBeamIntegration.h
        HingeRadauTwoBeamIntegration.h
        UserDefinedBeamIntegration.h
        SimpsonBeamIntegration.h
        UserDefinedHingeIntegration.h
        DistHingeIntegration.h
        RegularizedHingeIntegration.h
        TrapezoidalBeamIntegration.h
        CompositeSimpsonBeamIntegration.h
        FixedLocationBeamIntegration.h
        LowOrderBeamIntegration.h
        MidDistanceBeamIntegration.h
        #GaussQBeamIntegration.h
	ChebyshevBeamIntegration.h
)

target_sources(OPS_ElementFortran
    PRIVATE
        gaussq.f
        d1mach.f
)

#target_link_libraries(OPS_Element_forceBeamColumn PUBLIC OPS_Element_forceBeamColumn_f)
target_include_directories(OPS_Element PUBLIC $(CMAKE_CURRENT_LIST_DIR))

