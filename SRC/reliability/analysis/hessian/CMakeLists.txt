
target_sources(OPS_Reliability
    PRIVATE
        FiniteDifferenceHessian.cpp
        HessianEvaluator.cpp
    PUBLIC
        FiniteDifferenceHessian.h
        HessianEvaluator.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

