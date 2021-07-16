add_library(OPS_Element_beam3d OBJECT)

target_sources(OPS_Element_beam3d
    PRIVATE
        beam3d01.cpp
        beam3d02.cpp
    PUBLIC
        beam3d01.h
        beam3d02.h
)

target_include_directories(OPS_Element_beam3d PUBLIC ${CMAKE_CURRENT_LIST_DIR})

