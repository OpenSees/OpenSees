
target_sources(OPS_Reliability
    PRIVATE
        FindDesignPointAlgorithm.cpp
        SearchWithStepSizeAndStepDirection.cpp
    PUBLIC
        FindDesignPointAlgorithm.h
        SearchWithStepSizeAndStepDirection.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

