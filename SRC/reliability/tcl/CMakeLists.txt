target_sources(OPS_Reliability
    PRIVATE
        TclReliabilityBuilder.cpp
    PUBLIC
        TclReliabilityBuilder.h    
)
target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

