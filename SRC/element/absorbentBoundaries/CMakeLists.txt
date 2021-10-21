#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Element_absorbentBoundaries OBJECT)

target_sources(OPS_Element_absorbentBoundaries
  PRIVATE
    ASDAbsorbingBoundary2D.cpp
    ASDAbsorbingBoundary3D.cpp
    LysmerTriangle.cpp
  PUBLIC
    ASDAbsorbingBoundary2D.h
    ASDAbsorbingBoundary3D.h
    LysmerTriangle.h
)
target_include_directories(OPS_Element_absorbentBoundaries PUBLIC ${CMAKE_CURRENT_LIST_DIR})

