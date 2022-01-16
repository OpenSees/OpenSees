target_sources(OPS_Reliability
    PRIVATE
        LimitStateFunction.cpp
        LimitStateFunctionIter.cpp
        PerformanceFunction.cpp
        PerformanceFunctionIter.cpp
    PUBLIC
        LimitStateFunction.h
        LimitStateFunctionIter.h
        PerformanceFunction.h
        PerformanceFunctionIter.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

