add_library(utility "")

target_sources(utility
    PRIVATE
	Timer.cpp 
	FileIter.cpp 
	File.cpp 
	SimulationInformation.cpp 
	StringContainer.cpp 
    PUBLIC
	Timer.h 
	FileIter.h 
	File.h 
	SimulationInformation.h 
	StringContainer.h 
)

target_include_directories(utility PUBLIC $(CMAKE_CURRENT_LIST_DIR))
