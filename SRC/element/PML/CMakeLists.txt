add_library(OPS_Element_PML OBJECT)
add_library(OPS_Element_PML_f OBJECT)

target_sources(OPS_Element_PML
    PRIVATE
        PML2D.cpp
        PML3D.cpp
    PUBLIC
        PML2D.h
        PML3D.h
)

target_sources(OPS_Element_PML_f
    PRIVATE
        pml_3d.f
        pml_2d.f
)

target_link_libraries(OPS_Element_PML PUBLIC OPS_Element_PML_f)
target_include_directories(OPS_Element_PML PUBLIC ${CMAKE_CURRENT_LIST_DIR})

