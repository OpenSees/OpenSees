#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

find_package(MPI)
if (NOT MPI_FOUND)
  return()
endif()


add_library(OpenSeesRT_Parallel OBJECT)

target_include_directories(OpenSeesRT_Parallel
  PRIVATE
  ${OPS_SRC_DIR}/domain/partitioner/
  ${OPS_SRC_DIR}/analysis/analysis/
)


target_link_libraries(OpenSeesRT_Parallel PRIVATE 
    OPS_Runtime 
    OPS_Parallel 
    OPS_Actor 
    MPI::MPI_CXX
)

target_sources(OpenSeesRT_Parallel PRIVATE
  machine.cpp
)

if (MKL_FOUND)
  message(STATUS "MKL was found.")
  set(SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
else()
  message(STATUS "MKL NOT found.")
  find_package(SCALAPACK)
  message(">>> ScaLAPACK: ${ScaLAPACK_LIBRARIES}")
endif()
  
target_link_libraries(OpenSeesRT_Parallel PRIVATE 
    OPS_Parallel 
    OPS_Actor
    MPI::MPI_CXX
    ${MPI_Fortran_LIBRARIES}
)

#
# OpenSeesMP
#
add_library(LibOpenSeesMP  SHARED  EXCLUDE_FROM_ALL )
set_property(TARGET LibOpenSeesMP PROPERTY POSITION_INDEPENDENT_CODE 1)

target_link_libraries(LibOpenSeesMP PRIVATE 
        OpenSeesRT_Parallel
	OPS_Parallel 
	OPS_Actor  
	OpenSeesRT 
	OPS_Runtime 
        MPI::MPI_CXX
)

target_sources(LibOpenSeesMP PRIVATE 
    communicate.cpp
    ${OPS_SRC_DIR}/parallel/OpenSeesMP.cpp
)

#
# OpenSeesSP
#
add_library(LibOpenSeesSP  SHARED  EXCLUDE_FROM_ALL )
set_property(TARGET LibOpenSeesSP PROPERTY POSITION_INDEPENDENT_CODE 1)

target_sources(LibOpenSeesSP PRIVATE 
    partition.cpp
    ${OPS_SRC_DIR}/parallel/OpenSeesSP.cpp
)

target_link_libraries(LibOpenSeesSP 
  PRIVATE 
    OpenSeesRT_Parallel
    OPS_Parallel 
    OPS_Partition
    OPS_Analysis
    OPS_Runtime 
    OPS_Actor  
    OpenSeesRT 
    MPI::MPI_CXX
    METIS
)
