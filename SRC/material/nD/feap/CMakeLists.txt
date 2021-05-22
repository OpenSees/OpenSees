
target_sources(material
    PRIVATE
	feapCommon.f
	FeapMaterial01.cpp
	matl01.f
	FeapMaterial02.cpp
	matl02.f
	FeapMaterial03.cpp
	matl03.f
    PUBLIC
	FeapMaterial01.h
	FeapMaterial02.h
	FeapMaterial03.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))
