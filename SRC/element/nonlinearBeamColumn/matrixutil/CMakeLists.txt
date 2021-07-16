target_sources(G3_Matrix
    PRIVATE
    MatrixUtil.cpp
    PUBLIC
    MatrixUtil.h
)

target_include_directories(OPS_Element PUBLIC $(CMAKE_CURRENT_LIST_DIR))
target_include_directories(G3_Matrix PUBLIC $(CMAKE_CURRENT_LIST_DIR))


