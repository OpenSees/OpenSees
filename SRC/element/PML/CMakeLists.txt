
opensees_elements(OPS_Element_PML
  SOURCES
    PML2D.cpp
    PML2D.h
    PML3D.cpp
    PML3D.h
)

add_library(OPS_Element_PML_f EXCLUDE_FROM_ALL)
target_sources(OPS_Element_PML_f
  PUBLIC
    pml_3d.f
    pml_2d.f
)

target_link_libraries(OPS_Element_PML PUBLIC OPS_Element_PML_f)
#target_link_libraries(OPS_Element PUBLIC OPS_Element_PML_f)
#target_include_directories(OPS_Element_PML PUBLIC ${CMAKE_CURRENT_LIST_DIR})

