add_executable(matrix-vector-otherfact-test
  test_matrix_vector_otherfact.cpp
  ${PROJECT_SOURCE_DIR}/SRC/matrix/Matrix.cpp
  ${PROJECT_SOURCE_DIR}/SRC/matrix/Vector.cpp
  ${PROJECT_SOURCE_DIR}/SRC/matrix/ID.cpp
  ${PROJECT_SOURCE_DIR}/SRC/handler/OPS_Stream.cpp
  ${PROJECT_SOURCE_DIR}/SRC/handler/DummyStream.cpp
  ${PROJECT_SOURCE_DIR}/SRC/actor/actor/MovableObject.cpp
)

target_include_directories(matrix-vector-otherfact-test PRIVATE
  ${PROJECT_SOURCE_DIR}/SRC
  ${PROJECT_SOURCE_DIR}/SRC/matrix
  ${PROJECT_SOURCE_DIR}/SRC/handler
  ${PROJECT_SOURCE_DIR}/SRC/actor/actor
)

target_link_libraries(matrix-vector-otherfact-test PRIVATE
  ${LAPACK_LIBRARIES}
  ${CMAKE_DL_LIBS}
)

add_test(
  NAME matrix-vector-otherfact
  COMMAND matrix-vector-otherfact-test
)
