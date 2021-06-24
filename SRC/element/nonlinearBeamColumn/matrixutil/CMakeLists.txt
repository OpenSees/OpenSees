target_sources(OPS_Element
    PRIVATE
	MatrixUtil.cpp
    PUBLIC
	MatrixUtil.h
)

target_include_directories(OPS_Element PUBLIC $(CMAKE_CURRENT_LIST_DIR))


