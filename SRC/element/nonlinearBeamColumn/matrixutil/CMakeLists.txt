target_sources(element
    PRIVATE
	MatrixUtil.cpp
    PUBLIC
	MatrixUtil.h
)

target_include_directories(element PUBLIC $(CMAKE_CURRENT_LIST_DIR))


