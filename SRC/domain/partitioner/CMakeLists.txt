
target_sources(G3_Domain
    PRIVATE
	DomainPartitioner.cpp
    PUBLIC
	DomainPartitioner.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
