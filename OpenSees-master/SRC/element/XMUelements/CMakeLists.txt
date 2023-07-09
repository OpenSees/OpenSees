#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
  PRIVATE
      AC3D8HexWithSensitivity.cpp
      ASI3D8QuadWithSensitivity.cpp
      AV3D4QuadWithSensitivity.cpp
      VS3D4QuadWithSensitivity.cpp
  PUBLIC
      AC3D8HexWithSensitivity.h
      ASI3D8QuadWithSensitivity.h
      AV3D4QuadWithSensitivity.h
      VS3D4QuadWithSensitivity.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
