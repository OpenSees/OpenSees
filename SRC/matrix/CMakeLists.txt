add_library(matrix "")

target_sources(matrix
    PRIVATE
      Matrix.cpp
      Vector.cpp
      ID.cpp
    PUBLIC
      Matrix.h
      Vector.h
      ID.h
)


target_include_directories(matrix PUBLIC $(CMAKE_CURRENT_LIST_DIR))
