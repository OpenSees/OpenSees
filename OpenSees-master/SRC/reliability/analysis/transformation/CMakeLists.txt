
target_sources(OPS_Reliability
    PRIVATE
        AllIndependentTransformation.cpp
        NatafProbabilityTransformation.cpp
        ProbabilityTransformation.cpp
    PUBLIC
        AllIndependentTransformation.h
        NatafProbabilityTransformation.h
        ProbabilityTransformation.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

