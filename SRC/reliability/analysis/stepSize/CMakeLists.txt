
target_sources(OPS_Reliability
    PRIVATE
        ArmijoStepSizeRule.cpp
        FixedStepSizeRule.cpp
        StepSizeRule.cpp
    PUBLIC
        ArmijoStepSizeRule.h
        FixedStepSizeRule.h
        StepSizeRule.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

