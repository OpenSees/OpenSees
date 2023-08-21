
target_sources(OPS_Reliability
    PRIVATE
        OptimalityConditionReliabilityConvergenceCheck.cpp
        ReliabilityConvergenceCheck.cpp
        StandardReliabilityConvergenceCheck.cpp
    PUBLIC
        OptimalityConditionReliabilityConvergenceCheck.h
        ReliabilityConvergenceCheck.h
        StandardReliabilityConvergenceCheck.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

