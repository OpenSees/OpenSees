cmake_minimum_required(VERSION 3.16)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(CFSTsteelRegression LANGUAGES CXX)
  enable_testing()
  if(NOT LAPACK_LIBRARIES)
    find_package(LAPACK REQUIRED)
  endif()
endif()
get_filename_component(OPS_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)

add_executable(cfststeel-test
  ../test_cfststeel.cpp
  ${OPS_ROOT}/SRC/material/uniaxial/CFSTsteel.cpp
  ${OPS_ROOT}/SRC/material/uniaxial/UniaxialMaterial.cpp
  ${OPS_ROOT}/SRC/material/Material.cpp
  ${OPS_ROOT}/SRC/tagged/TaggedObject.cpp
  ${OPS_ROOT}/SRC/tagged/storage/MapOfTaggedObjects.cpp
  ${OPS_ROOT}/SRC/tagged/storage/MapOfTaggedObjectsIter.cpp
  ${OPS_ROOT}/SRC/actor/actor/MovableObject.cpp
  ${OPS_ROOT}/SRC/actor/channel/Channel.cpp
  ${OPS_ROOT}/SRC/actor/objectBroker/FEM_ObjectBroker.cpp
  ${OPS_ROOT}/SRC/matrix/Vector.cpp
  ${OPS_ROOT}/SRC/matrix/Matrix.cpp
  ${OPS_ROOT}/SRC/matrix/ID.cpp
  ${OPS_ROOT}/SRC/handler/OPS_Stream.cpp
  ${OPS_ROOT}/SRC/handler/StandardStream.cpp
  ${OPS_ROOT}/SRC/recorder/response/Response.cpp
  ${OPS_ROOT}/SRC/recorder/response/MaterialResponse.cpp
  ${OPS_ROOT}/SRC/element/Information.cpp
)
target_compile_features(cfststeel-test PRIVATE cxx_std_17)
target_include_directories(cfststeel-test PRIVATE
  ${OPS_ROOT}/SRC
  ${OPS_ROOT}/SRC/api
  ${OPS_ROOT}/SRC/material
  ${OPS_ROOT}/SRC/material/uniaxial
  ${OPS_ROOT}/SRC/tagged
  ${OPS_ROOT}/SRC/tagged/storage
  ${OPS_ROOT}/SRC/actor/actor
  ${OPS_ROOT}/SRC/actor/channel
  ${OPS_ROOT}/SRC/actor/message
  ${OPS_ROOT}/SRC/actor/objectBroker
  ${OPS_ROOT}/SRC/domain/component
  ${OPS_ROOT}/SRC/matrix
  ${OPS_ROOT}/SRC/handler
  ${OPS_ROOT}/SRC/element
  ${OPS_ROOT}/SRC/recorder/response
)
if(MSVC)
  target_compile_options(cfststeel-test PRIVATE /utf-8)
  set_property(TARGET cfststeel-test PROPERTY MSVC_RUNTIME_LIBRARY MultiThreaded)
  target_compile_definitions(cfststeel-test PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()
target_link_libraries(cfststeel-test PRIVATE ${LAPACK_LIBRARIES} ${CMAKE_DL_LIBS})
add_test(NAME cfststeel-lifecycle COMMAND cfststeel-test)
