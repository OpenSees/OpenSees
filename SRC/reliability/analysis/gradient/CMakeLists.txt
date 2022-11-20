
target_sources(OPS_Reliability
    PRIVATE
        FiniteDifferenceGradient.cpp
        GradientEvaluator.cpp
        ImplicitGradient.cpp
    PUBLIC
        FiniteDifferenceGradient.h
        GradientEvaluator.h
        ImplicitGradient.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

