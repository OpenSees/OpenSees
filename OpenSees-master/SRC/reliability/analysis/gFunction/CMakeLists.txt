
target_sources(OPS_Reliability
    PRIVATE
        BasicGFunEvaluator.cpp
        GFunEvaluator.cpp
        OpenSeesGFunEvaluator.cpp
    PUBLIC
        BasicGFunEvaluator.h
        GFunEvaluator.h
        OpenSeesGFunEvaluator.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

