#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
    PRIVATE
        ConstantSeries.cpp
        EarthquakePattern.cpp
        LinearSeries.cpp
        LoadPattern.cpp
        LoadPatternIter.cpp
        MultiSupportPattern.cpp
        PathSeries.cpp
        PathTimeSeries.cpp
        PulseSeries.cpp
        RectangularSeries.cpp
        SimpsonTimeSeriesIntegrator.cpp
        TimeSeries.cpp
        TimeSeriesIntegrator.cpp
        TrapezoidalTimeSeriesIntegrator.cpp
        TriangleSeries.cpp
        TrigSeries.cpp
        RampSeries.cpp
	MPAccSeries.cpp
        UniformExcitation.cpp
    PUBLIC
        ConstantSeries.h
        EarthquakePattern.h
        LinearSeries.h
        LoadPattern.h
        LoadPatternIter.h
        MultiSupportPattern.h
        PathSeries.h
        PathTimeSeries.h
        PulseSeries.h
        RectangularSeries.h
        SimpsonTimeSeriesIntegrator.h
        TimeSeries.h
        TimeSeriesIntegrator.h
        TrapezoidalTimeSeriesIntegrator.h
        TriangleSeries.h
        TrigSeries.h
        RampSeries.h
	MPAccSeries.h
        UniformExcitation.h
)

target_sources(OPS_Domain
    PRIVATE
        FireLoadPattern.cpp
        PathTimeSeriesThermal.cpp
        PeerMotion.cpp
        PeerNGAMotion.cpp
    PUBLIC
        FireLoadPattern.h
        PathTimeSeriesThermal.h
        PeerMotion.h
        PeerNGAMotion.h
)

# The following files are currently being added to the OPS_Tcl target
# in the root-level CMakeLists.txt file

#target_sources(OPS_Tcl
#    PRIVATE
#        TclSeriesCommand.cpp
#        TclSeriesIntegratorCommand.cpp
#        TclPatternCommand.cpp
#)



target_sources(OPS_Reliability
    PRIVATE
    SimulatedRandomProcessSeries.cpp
    DiscretizedRandomProcessSeries.cpp    
    PUBLIC
    SimulatedRandomProcessSeries.h
    DiscretizedRandomProcessSeries.h    
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(drm)


