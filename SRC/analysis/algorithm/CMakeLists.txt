
target_sources(analysis
    PRIVATE
      SolutionAlgorithm.cpp
    PUBLIC
      SolutionAlgorithm.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(domainDecompAlgo)
add_subdirectory(equiSolnAlgo)
add_subdirectory(eigenAlgo)
