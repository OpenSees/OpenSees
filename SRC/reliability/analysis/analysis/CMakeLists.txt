
add_subdirectory(system)

target_sources(OPS_ReliabilityTcl
    PRIVATE
        DP_RSM_Sim.cpp
        DP_RSM_Sim_TimeVariant.cpp
        MonteCarloResponseAnalysis.cpp
        MultiDimVisPrincPlane.cpp
      PUBLIC
        DP_RSM_Sim.h
        DP_RSM_Sim_TimeVariant.h
        MonteCarloResponseAnalysis.h
        MultiDimVisPrincPlane.h
)
      
target_sources(OPS_Reliability
    PRIVATE
        BivariateDecomposition.cpp
        ExperimentalPointRule1D.cpp
        FORMAnalysis.cpp
        FOSMAnalysis.cpp
	# FragilityAnalysis.cpp
        GFunVisualizationAnalysis.cpp
        GridPlane.cpp
        ImportanceSamplingAnalysis.cpp
	#        MVFOSMAnalysis.cpp
        # OptimizationAnalysis.cpp
        OrthogonalPlaneSamplingAnalysis.cpp
        OutCrossingAnalysis.cpp
        # ParametricReliabilityAnalysis.cpp
        PrincipalAxis.cpp
        ReliabilityAnalysis.cpp
        RespSurfaceSimulation.cpp
        SORMAnalysis.cpp
        SamplingAnalysis.cpp
        SurfaceDesign.cpp
        SystemAnalysis.cpp
        UniformExperimentalPointRule1D.cpp
        UnivariateDecomposition.cpp
    PUBLIC
        BivariateDecomposition.h
        ExperimentalPointRule1D.h
        FORMAnalysis.h
        FOSMAnalysis.h
        # FragilityAnalysis.h
        GFunVisualizationAnalysis.h
        GridPlane.h
        ImportanceSamplingAnalysis.h
        # MVFOSMAnalysis.h
        # OptimizationAnalysis.h
        OrthogonalPlaneSamplingAnalysis.h
        OutCrossingAnalysis.h
        # ParametricReliabilityAnalysis.h
        PrincipalAxis.h
        ReliabilityAnalysis.h
        RespSurfaceSimulation.h
        SORMAnalysis.h
        SamplingAnalysis.h
        SurfaceDesign.h
        SystemAnalysis.h
        UniformExperimentalPointRule1D.h
        UnivariateDecomposition.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

