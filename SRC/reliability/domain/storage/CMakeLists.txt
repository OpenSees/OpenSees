target_sources(OPS_Reliability
    PRIVATE
        FORMStorage.cpp
        ReliabilityStorage.cpp
    PUBLIC
        FORMStorage.h
        ReliabilityStorage.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

