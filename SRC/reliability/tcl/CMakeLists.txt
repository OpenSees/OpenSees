target_sources(OPS_ReliabilityTcl
    PRIVATE
        TclReliabilityBuilder.cpp
    PUBLIC
        TclReliabilityBuilder.h    
)
target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

