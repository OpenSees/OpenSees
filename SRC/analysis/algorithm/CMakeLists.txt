
target_sources(G3_Analysis
    PRIVATE
      SolutionAlgorithm.cpp
    PUBLIC
      SolutionAlgorithm.h
)

#target_include_directories(G3_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(domainDecompAlgo)
add_subdirectory(equiSolnAlgo)
add_subdirectory(eigenAlgo)
