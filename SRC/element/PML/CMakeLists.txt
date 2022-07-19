#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#opensees_elements(OPS_Element_PML#
#  SOURCES

target_sources(OPS_Element
  PRIVATE
    PML2D.cpp
    PML3D.cpp
  PUBLIC    
    PML3D.h
    PML2D.h
)

target_sources(OPS_ElementFortran
  PRIVATE
    pml_3d.f
    pml_2d.f
)

#target_link_libraries(OPS_Element_PML PUBLIC OPS_Element_PML_f)
#target_link_libraries(OPS_Element PUBLIC OPS_Element_PML_f)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

