
add_library(MUMPS_C)


add_subdirectory(src)
add_subdirectory(PORD/lib)

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

