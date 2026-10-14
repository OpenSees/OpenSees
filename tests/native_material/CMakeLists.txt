cmake_minimum_required(VERSION 3.16)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(OpenSeesNativeMaterialTests LANGUAGES CXX)
  enable_testing()
endif()
get_filename_component(OPENSEES_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(HYSTERETICSM_SOURCE "${OPENSEES_SOURCE_ROOT}/SRC/material/uniaxial/HystereticSMMaterial.cpp"
    CACHE FILEPATH "Material source for native regression testing")

add_library(native-material-test-support STATIC
  ${OPENSEES_SOURCE_ROOT}/SRC/material/uniaxial/UniaxialMaterial.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/material/Material.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/matrix/Vector.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/matrix/Matrix.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/matrix/ID.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/handler/OPS_Stream.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/handler/DummyStream.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/actor/actor/MovableObject.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/tagged/TaggedObject.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/tagged/storage/MapOfTaggedObjects.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/tagged/storage/MapOfTaggedObjectsIter.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/recorder/response/Response.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/recorder/response/MaterialResponse.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/element/Information.cpp
  ${OPENSEES_SOURCE_ROOT}/SRC/domain/component/Parameter.cpp
  unused_entrypoints.cpp
)
foreach(directory IN ITEMS "" api material material/uniaxial matrix handler actor/actor
        actor/channel actor/address actor/objectBroker tagged tagged/storage recorder/response
        element domain/component)
  target_include_directories(native-material-test-support PUBLIC "${OPENSEES_SOURCE_ROOT}/SRC/${directory}")
endforeach()
target_compile_features(native-material-test-support PUBLIC cxx_std_17)

add_executable(hystereticsm-trial-state-test ../test_hystereticsm_trial_state.cpp ${HYSTERETICSM_SOURCE})
target_link_libraries(hystereticsm-trial-state-test PRIVATE native-material-test-support)
add_test(NAME hystereticsm-trial-state COMMAND hystereticsm-trial-state-test)

if(MSVC)
  foreach(target native-material-test-support hystereticsm-trial-state-test)
    if(TARGET ${target})
      set_property(TARGET ${target} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
      target_compile_options(${target} PRIVATE /fp:precise /utf-8)
    endif()
  endforeach()
endif()
