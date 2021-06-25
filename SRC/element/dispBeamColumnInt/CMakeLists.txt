add_library(OPS_Element_dispBeamColumnInt OBJECT)

target_sources(OPS_Element_dispBeamColumnInt
    PRIVATE
        DispBeamColumn2dInt.cpp
        FiberSection2dInt.cpp
        LinearCrdTransf2dInt.cpp
        #TclDispBeamColumnIntCommand.cpp
    PUBLIC
        DispBeamColumn2dInt.h
        FiberSection2dInt.h
        LinearCrdTransf2dInt.h
)
target_include_directories(OPS_Element_dispBeamColumnInt PUBLIC $(CMAKE_CURRENT_LIST_DIR))

