#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Section_Repres OBJECT)

target_include_directories(OPS_Section_Repres PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(cell)
add_subdirectory(patch)
add_subdirectory(reinfBar)
add_subdirectory(reinfLayer)

