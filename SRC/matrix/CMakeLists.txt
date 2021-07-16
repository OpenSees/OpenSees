#add_library(matrix OBJECT)

target_sources(G3_Matrix
    PRIVATE
      Matrix.cpp
      Vector.cpp
      ID.cpp
    PUBLIC
      Matrix.h
      Vector.h
      ID.h
)


target_include_directories(G3_Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR})
