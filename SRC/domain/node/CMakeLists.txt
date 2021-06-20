
target_sources(libg3
    PRIVATE
	Node.cpp
	NodalLoad.cpp
    PUBLIC
	Node.h
	NodalLoad.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
