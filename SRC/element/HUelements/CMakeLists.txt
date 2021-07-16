target_sources(OPS_Element
    PRIVATE
	MultipleShearSpring.cpp
	MultipleNormalSpring.cpp
	DBESI0.C
	DBESI1.C
	YamamotoBiaxialHDR.cpp
	KikuchiBearing.cpp
    PUBLIC
	MultipleShearSpring.h
	MultipleNormalSpring.h
	#DBESI0.h
	#DBESI1.h
	YamamotoBiaxialHDR.h
	KikuchiBearing.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
