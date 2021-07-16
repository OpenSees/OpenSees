
target_sources(G3_Domain
    PRIVATE
	MeshRegion.cpp
    PUBLIC
	MeshRegion.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

