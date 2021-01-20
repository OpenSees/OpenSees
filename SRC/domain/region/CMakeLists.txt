
target_sources(domain
    PRIVATE
	MeshRegion.cpp
    PUBLIC
	MeshRegion.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
