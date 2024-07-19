#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_compile_definitions(_NOGRAPHICS)

add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)

target_sources(OPS_Renderer
  PRIVATE
    PlainMap.cpp
    Renderer.cpp
  PUBLIC
    PlainMap.h
    Renderer.h
)

target_include_directories(OPS_Renderer PUBLIC ${CMAKE_CURRENT_LIST_DIR})
