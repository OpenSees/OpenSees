target_sources(element
    PRIVATE
	FourNodeTetrahedron.cpp
    PUBLIC
	FourNodeTetrahedron.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
