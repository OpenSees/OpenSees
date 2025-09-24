
target_compile_definitions(OpenSeesRT_Parallel PRIVATE 
  _PARALLEL_INTERPRETERS
  _PARALLEL_PROCESSING     # OpenSeesSP
)
target_link_libraries(OpenSeesRT_Parallel PRIVATE OPS_Parallel OPS_Actor)

target_sources(OpenSeesRT_Parallel PRIVATE
  parallel_commands.cpp
)

add_executable(OpenSeesSP    EXCLUDE_FROM_ALL )
add_executable(OpenSeesMP    EXCLUDE_FROM_ALL )

if(MKL_FOUND)
  message(STATUS "MKL was found.")
  set(SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
else()
  message(STATUS "MKL NOT found.")
  find_package(SCALAPACK)
  message(">>> ScaLAPACK: ${ScaLAPACK_LIBRARIES}")
endif()
  
target_link_libraries(OpenSeesRT_Parallel PRIVATE OPS_Parallel OPS_Actor
    ${MPI_CXX_LIBRARIES}
    /usr/lib/libscalapack.so
    ${MPI_Fortran_LIBRARIES}
    ${MPI_CXX_LINK_FLAGS}
)

target_sources(OpenSeesMP PRIVATE mpiMain.cpp)
target_link_libraries(OpenSeesMP PRIVATE OPS_Parallel OPS_Actor  OpenSeesRT ${MPI_LIBRARIES})

target_sources(OpenSeesSP PRIVATE mpiMain.cpp)
target_link_libraries(OpenSeesSP PRIVATE OPS_Parallel OPS_Actor  OpenSeesRT ${MPI_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    /usr/lib/libscalapack.so
    ${MPI_Fortran_LIBRARIES}       
    ${MPI_CXX_LINK_FLAGS}
)
