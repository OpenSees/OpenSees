target_sources(analysis
    PRIVATE
      DomainDecompAlgo.cpp
    PUBLIC
      DomainDecompAlgo.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
