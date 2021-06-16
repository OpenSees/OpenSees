target_sources(analysis
    PRIVATE
      FE_Element.cpp
    PUBLIC
      FE_Element.h
)

target_include_directories(analysis PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(penalty)
add_subdirectory(lagrange)
add_subdirectory(transformation)
