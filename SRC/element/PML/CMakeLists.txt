add_library(OPS_Element_PML OBJECT)

target_sources(OPS_Element_PML
    PRIVATE
        PML2D.cpp
        pml_2d.f
        PML3D.cpp
        pml_3d.f
    PUBLIC
        PML2D.h
        PML3D.h
)
target_include_directories(OPS_Element_PML PUBLIC ${CMAKE_CURRENT_LIST_DIR})

