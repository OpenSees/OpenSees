
target_sources(OPS_Reliability
    PRIVATE
        FiniteDifferenceGradGEvaluator.cpp
        GradGEvaluator.cpp
        OpenSeesGradGEvaluator.cpp
    PUBLIC
        FiniteDifferenceGradGEvaluator.h
        GradGEvaluator.h
        OpenSeesGradGEvaluator.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

