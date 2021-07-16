target_sources(G3_Analysis
    PRIVATE
      DomainDecompAlgo.cpp
    PUBLIC
      DomainDecompAlgo.h
)

target_include_directories(G3_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
