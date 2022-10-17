#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Analysis
    PRIVATE
      SolutionAlgorithm.cpp
    PUBLIC
      SolutionAlgorithm.h
)

#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(domainDecompAlgo)
add_subdirectory(equiSolnAlgo)
add_subdirectory(eigenAlgo)
