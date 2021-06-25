add_library(OPS_String OBJECT)

target_sources(OPS_String
    PRIVATE
        G3string.cpp
    PUBLIC
        G3string.h
)
target_include_directories(OPS_String PUBLIC $(CMAKE_CURRENT_LIST_DIR))

