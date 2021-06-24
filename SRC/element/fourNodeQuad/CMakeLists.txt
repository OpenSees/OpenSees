target_sources(OPS_Element
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

target_include_directories(OPS_Element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
