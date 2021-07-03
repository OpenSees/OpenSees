target_sources(element
    PRIVATE
	FourNodeQuad.cpp
	ConstantPressureVolumeQuad.cpp
	EnhancedQuad.cpp
	FourNodeQuad3d.cpp
	NineNodeMixedQuad.cpp
	FourNodeQuadWithSensitivity.cpp
    PUBLIC
	FourNodeQuad.h
	ConstantPressureVolumeQuad.h
	EnhancedQuad.h
	FourNodeQuad3d.h
	NineNodeMixedQuad.h
	FourNodeQuadWithSensitivity.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
