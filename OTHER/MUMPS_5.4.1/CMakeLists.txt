add_library(MUMPS)

add_subdirectory(PORD/lib)
add_subdirectory(src)

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
