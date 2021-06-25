add_library(OPS_Element_absorbentBoundaries OBJECT)

target_sources(OPS_Element_absorbentBoundaries
    PRIVATE
        LysmerTriangle.cpp
    PUBLIC
        LysmerTriangle.h
)
target_include_directories(OPS_Element_absorbentBoundaries PUBLIC $(CMAKE_CURRENT_LIST_DIR))

