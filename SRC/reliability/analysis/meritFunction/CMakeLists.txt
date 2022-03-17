
target_sources(OPS_Reliability
    PRIVATE
        AdkZhangMeritFunctionCheck.cpp
        CriteriaReductionMeritFunctionCheck.cpp
        MeritFunctionCheck.cpp
    PUBLIC
        AdkZhangMeritFunctionCheck.h
        CriteriaReductionMeritFunctionCheck.h
        MeritFunctionCheck.h
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

