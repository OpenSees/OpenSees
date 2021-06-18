target_sources(element
    PRIVATE
	Inerter.cpp
	LinearElasticSpring.cpp
	TwoNodeLink.cpp
    PUBLIC
	Inerter.h
	LinearElasticSpring.h
	TwoNodeLink.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
