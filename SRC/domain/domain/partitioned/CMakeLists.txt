
target_sources(domain
    PRIVATE
	PartitionedDomainEleIter.cpp 
	PartitionedDomainSubIter.cpp
	PartitionedDomain.cpp
    PUBLIC
	PartitionedDomainEleIter.h
	PartitionedDomainSubIter.h
	PartitionedDomain.h
)

target_include_directories(domain PUBLIC $(CMAKE_CURRENT_LIST_DIR))
