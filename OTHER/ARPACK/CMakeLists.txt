#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#------------------------------------------------------------------------------
#   arpack.caam.rice.edu
#==============================================================================
# Claudio Perez

project(ARPACK Fortran)

add_library(ARPACK)

# This checks whether the compiler has an "argument mismatch" flag;
# if it does, it needs to be set. This solution was inspired by 
# a CMake file in the OneLab repositories
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-fallow-argument-mismatch" ARGMISMATCH)
if (ARGMISMATCH)
    message("Setting -fallow-argument-mismatch")
    target_compile_options(ARPACK PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
endif()

target_sources(ARPACK PUBLIC
  dgetv0.f
  dlaqrb.f
  dstqrb.f
  dsortc.f
  dsortr.f
  dstatn.f
  dstats.f
  dnaitr.f
  dnapps.f
  dnaup2.f
  dnaupd.f
  dnconv.f
  dneigh.f
  dngets.f
  dsaitr.f
  dsapps.f
  dsaup2.f
  dsaupd.f
  dsconv.f
  dseigt.f
  dsgets.f
  
  dneupd.f
  dseupd.f
  dsesrt.f
  dvout.f
  dmout.f
  
  icnteq.f
  icopy.f
  iset.f
  iswap.f
  ivout.f
  second.f
)

target_link_libraries(ARPACK PRIVATE ${LAPACK_LIBRARIES})

# Required on Unix OS family to be able to be linked into shared libraries.
set_target_properties(${PROJECT_NAME}
                      PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Expose public includes to other subprojects through cache variable.
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

export(TARGETS ARPACK
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/ARPACK.cmake"
)


