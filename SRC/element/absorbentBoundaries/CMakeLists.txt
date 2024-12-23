#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#add_library(OPS_Element_absorbentBoundaries OBJECT)

#target_sources(OPS_Element_absorbentBoundaries

target_sources(OPS_Element
  PRIVATE
    ASDAbsorbingBoundary2D.cpp
    ASDAbsorbingBoundary3D.cpp
    LysmerTriangle.cpp
    FSIFluidElement2D.cpp
    FSIInterfaceElement2D.cpp
    FSIFluidBoundaryElement2D.cpp
  PUBLIC
    ASDAbsorbingBoundary2D.h
    ASDAbsorbingBoundary3D.h
    LysmerTriangle.h
    FSIFluidElement2D.h
    FSIInterfaceElement2D.h
    FSIFluidBoundaryElement2D.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

