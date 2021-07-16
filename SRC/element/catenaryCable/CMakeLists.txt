add_library(OPS_Element_catenaryCable OBJECT)

target_sources(OPS_Element_catenaryCable
    PRIVATE
        CatenaryCable.cpp
    PUBLIC
        CatenaryCable.h
)
target_include_directories(OPS_Element_catenaryCable PUBLIC ${CMAKE_CURRENT_LIST_DIR})

