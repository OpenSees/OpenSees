
target_sources(LG3_Analysis
    PRIVATE
      SolutionAlgorithm.cpp
    PUBLIC
      SolutionAlgorithm.h
)

#target_include_directories(LG3_Analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(domainDecompAlgo)
add_subdirectory(equiSolnAlgo)
add_subdirectory(eigenAlgo)
