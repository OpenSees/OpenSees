#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#add_library(OPS_Reliability OBJECT)
add_library(OPS_ReliabilityTcl OBJECT)

target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_ReliabilityTcl PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(./analysis/)
add_subdirectory(./domain/)
add_subdirectory(./FEsensitivity/)
add_subdirectory(./tcl/)
