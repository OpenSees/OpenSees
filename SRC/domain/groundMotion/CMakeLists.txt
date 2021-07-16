
target_sources(G3_Domain
    PRIVATE    
        GroundMotion.cpp
    GroundMotionRecord.cpp
    InterpolatedGroundMotion.cpp
    PUBLIC
        GroundMotion.h
    GroundMotionRecord.h
    InterpolatedGroundMotion.h
)

target_include_directories(OPS_Domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
