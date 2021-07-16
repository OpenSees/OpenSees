add_library(OPS_Element_RockingBC OBJECT)

target_sources(OPS_Element_RockingBC
    PRIVATE
        RockingBC.cpp
    PUBLIC
        RockingBC.h
)
target_include_directories(OPS_Element_RockingBC PUBLIC ${CMAKE_CURRENT_LIST_DIR})

