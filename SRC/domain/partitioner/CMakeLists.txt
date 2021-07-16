
target_sources(G3_Domain
    PRIVATE
    DomainPartitioner.cpp
    PUBLIC
    DomainPartitioner.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
