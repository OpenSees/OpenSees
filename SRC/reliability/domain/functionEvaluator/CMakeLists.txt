target_sources(OPS_Reliability
    PRIVATE
        FunctionEvaluator.cpp
	# MatlabEvaluator.cpp
    PUBLIC
        FunctionEvaluator.h
        # MatlabEvaluator.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

