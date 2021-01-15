
target_sources(domain
    PRIVATE	
    	GroundMotion.cpp
	GroundMotionRecord.cpp
	InterpolatedGroundMotion.cpp
    PUBLIC
    	GroundMotion.h
	GroundMotionRecord.h
	InterpolatedGroundMotion.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
