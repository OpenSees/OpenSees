
target_sources(material
    PRIVATE
	DrainHardeningMaterial.cpp
	anal00.f
	common00.f
	DrainBilinearMaterial.cpp
	DrainClough1Material.cpp
	DrainClough2Material.cpp
	DrainPinch1Material.cpp
    PUBLIC
	DrainHardeningMaterial.h
	DrainBilinearMaterial.h
	DrainClough1Material.h
	DrainClough2Material.h
	DrainPinch1Material.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))


