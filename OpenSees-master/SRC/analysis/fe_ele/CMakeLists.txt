#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      FE_Element.cpp
    PUBLIC
      FE_Element.h
)

#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(penalty)
add_subdirectory(lagrange)
add_subdirectory(transformation)
