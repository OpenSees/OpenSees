#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_subdirectory(algorithm)
add_subdirectory(numberer)
add_subdirectory(dof_grp)
add_subdirectory(fe_ele)
add_subdirectory(handler)
add_subdirectory(model)
add_subdirectory(integrator)
add_subdirectory(analysis)

#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})


