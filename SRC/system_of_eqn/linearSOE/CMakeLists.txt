#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
  PRIVATE
    DomainSolver.cpp
    LinearSOE.cpp
    LinearSOESolver.cpp
  PUBLIC
    DomainSolver.h
    LinearSOE.h
    LinearSOESolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(bandGEN)
add_subdirectory(bandSPD)
add_subdirectory(diagonal)
add_subdirectory(fullGEN)
add_subdirectory(sparseGEN)
add_subdirectory(sparseSYM)
add_subdirectory(umfGEN)

add_subdirectory(profileSPD)
#add_subdirectory(cg)
if(DEFINED PETSC_FLAG)
  add_subdirectory(petsc)
endif()
message(STATUS "MUMPS_FLAG = ${MUMPS_FLAG}")
if(DEFINED MUMPS_FLAG)
endif()
#add_subdirectory(itpack)

