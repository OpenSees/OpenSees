target_sources(OPS_Reliability
    PRIVATE
        FunctionEvaluator.cpp
	# MatlabEvaluator.cpp
    PUBLIC
        FunctionEvaluator.h
        # MatlabEvaluator.h
)

target_sources(OPS_ReliabilityTcl
    PRIVATE
        TclEvaluator.cpp
    PUBLIC
        TclEvaluator.h
)

target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))
target_include_directories(OPS_ReliabilityTcl PUBLIC $(CMAKE_CURRENT_LIST_DIR))

