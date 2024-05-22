
target_sources(OPS_Reliability
    PRIVATE
        CStdLibRandGenerator.cpp
        RandomNumberGenerator.cpp
    PUBLIC
        CStdLibRandGenerator.h
        RandomNumberGenerator.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

