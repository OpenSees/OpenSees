add_library(OPS_Element_mixedBeamColumn OBJECT)

target_sources(OPS_Element_mixedBeamColumn
    PRIVATE
        MixedBeamColumn2d.cpp
        MixedBeamColumn3d.cpp
    MixedBeamColumnAsym3d.cpp
    PUBLIC
        MixedBeamColumn2d.h
        MixedBeamColumn3d.h
    MixedBeamColumnAsym3d.h
)
target_include_directories(OPS_Element_mixedBeamColumn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

