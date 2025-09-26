#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_executable(xara_repl main.cpp)
set_target_properties(xara_repl PROPERTIES 
  OUTPUT_NAME "xara"
)
if (WIN32)
  target_link_libraries(xara_repl PRIVATE Netapi32.lib)
  target_link_libraries(xara_repl PRIVATE ${TCL_LIBRARIES} OpenSeesRT)
else()
  target_link_libraries(xara_repl PRIVATE ${TCL_LIBRARY} OpenSeesRT)
endif()

add_executable(xara_static static.cpp)
set_target_properties(xara_static PROPERTIES 
  OUTPUT_NAME "xaras"
)
target_link_libraries(xara_static PUBLIC ${TCL_LIBRARY} OpenSeesRT_Static OPS_Runtime)