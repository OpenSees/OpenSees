
target_sources(material
    PRIVATE
	Pinching.cpp
	Bilinear.cpp
	Clough.cpp
	CloughHenry.cpp
	CloughDamage.cpp
	PinchingDamage.cpp
    PUBLIC
	Pinching.h
	Bilinear.h
	Clough.h
	CloughHenry.h
	CloughDamage.h
	PinchingDamage.h
)

target_include_directories(material PUBLIC $(CMAKE_CURRENT_LIST_DIR))


