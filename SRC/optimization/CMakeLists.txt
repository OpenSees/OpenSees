#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Optimization OBJECT)

target_include_directories(OPS_Optimization PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(./domain/)
add_subdirectory(./SNOPT/)
add_subdirectory(./tcl/)
