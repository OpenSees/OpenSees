
target_sources(libg3
    PRIVATE
	DomainPartitioner.cpp
    PUBLIC
	DomainPartitioner.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
