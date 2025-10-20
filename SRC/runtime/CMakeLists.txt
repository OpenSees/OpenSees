#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_library(OpenSeesRT SHARED)
add_library(OpenSeesRT_Static STATIC)

target_sources(OpenSeesRT PRIVATE "OpenSeesRT.cpp")
target_sources(OpenSeesRT_Static PRIVATE "OpenSeesRT.cpp")
target_compile_definitions(OpenSeesRT PUBLIC USE_TCL_STUBS _TCL85)
set_property(TARGET OpenSeesRT PROPERTY POSITION_INDEPENDENT_CODE 1)

target_compile_options(OPS_Runtime PRIVATE
  "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-unused-parameter>"
  "$<$<CXX_COMPILER_ID:MSVC>:/W3>"
  "$<$<CXX_COMPILER_ID:MSVC>:/bigobj>"
  "$<$<CXX_COMPILER_ID:MSVC>:/wd4100>"
)

target_include_directories(OPS_Runtime PUBLIC
  ${OPS_SRC_DIR}/runtime/commands
  ${OPS_SRC_DIR}/runtime/parsing/
  ${OPS_SRC_DIR}/runtime/commands/modeling
  ${OPS_SRC_DIR}/runtime/runtime/
)
target_link_libraries(OPS_Runtime PRIVATE G3 OPS_Algorithm)

add_subdirectory(commands)
add_subdirectory(runtime)
add_subdirectory(parsing)


if (DEFINED OPENSEESRT_VERSION)
  set_property(
    SOURCE OpenSeesRT.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS 
           "OPENSEESRT_VERSION=\"${OPENSEESRT_VERSION}\""
  )
endif()


#
# Tcl
#
if (${TCL_INCLUDE_PATH})
  target_include_directories(OpenSeesRT  PUBLIC ${TCL_INCLUDE_PATH})
  target_include_directories(OPS_Runtime PUBLIC ${TCL_INCLUDE_PATH})
endif()

target_link_libraries(OpenSeesRT 
  PRIVATE 
  OPS_Runtime 
  OPS_Renderer 
  OPS_Algorithm 
  ${TCL_STUB_LIBRARY}
)

target_link_libraries(OpenSeesRT_Static
  PUBLIC 
  OPS_Runtime 
  OPS_Renderer 
  OPS_Algorithm
)

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/packages/UCSD_CompGeoMech/CMakeLists.txt)
  add_subdirectory(packages/UCSD_CompGeoMech)
endif ()
