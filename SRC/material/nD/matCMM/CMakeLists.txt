target_sources(OPS_Material
    PRIVATE
    MaterialCMM.cpp
    PUBLIC
    MaterialCMM.h
)

target_include_directories(OPS_Material PUBLIC $(CMAKE_CURRENT_LIST_DIR))
