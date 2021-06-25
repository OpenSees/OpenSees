add_library(OPS_Reliability OBJECT)

target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(./analysis/)
add_subdirectory(./domain/)
add_subdirectory(./FEsensitivity/)
add_subdirectory(./tcl/)
