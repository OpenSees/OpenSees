
target_sources(OPS_Reliability
    PRIVATE
        CorrelatedStandardNormal.cpp
        MatrixOperations.cpp
    PUBLIC
        CorrelatedStandardNormal.h
        MatrixOperations.h
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

