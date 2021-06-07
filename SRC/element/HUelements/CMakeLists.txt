target_sources(element
    PRIVATE
	MultipleShearSpring.cpp
	MultipleNormalSpring.cpp
	DBESI0.cpp
	DBESI1.cpp
	YamamotoBiaxialHDR.cpp
	KikuchiBearing.cpp
    PUBLIC
	MultipleShearSpring.h
	MultipleNormalSpring.h
	DBESI0.h
	DBESI1.h
	YamamotoBiaxialHDR.h
	KikuchiBearing.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
