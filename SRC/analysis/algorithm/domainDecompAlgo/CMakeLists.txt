target_sources(LG3_Analysis
    PRIVATE
      DomainDecompAlgo.cpp
    PUBLIC
      DomainDecompAlgo.h
)

target_include_directories(LG3_Analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))
