
target_sources(OPS_Reliability
    PRIVATE
        GradientProjectionSearchDirection.cpp
        HLRFSearchDirection.cpp
        PolakHeSearchDirectionAndMeritFunction.cpp
        SQPsearchDirectionMeritFunctionAndHessian.cpp
        SearchDirection.cpp
    PUBLIC
        GradientProjectionSearchDirection.h
        HLRFSearchDirection.h
        PolakHeSearchDirectionAndMeritFunction.h
        SQPsearchDirectionMeritFunctionAndHessian.h
        SearchDirection.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

