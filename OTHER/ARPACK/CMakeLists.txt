############################################################################
#
#  Program:         ARPACK
#
#  Module:          Makefile
#
#  Purpose:         Sources Makefile
#
#  Creation date:   February 22, 1996
#
#  Modified:        2021
#
#  Send bug reports, comments or suggestions to arpack.caam.rice.edu
#
############################################################################

project(ARPACK Fortran)

add_library(ARPACK)

# Arg-mismatch check adapted from OneLab CMakeLists
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-fallow-argument-mismatch" ARGMISMATCH)
if (ARGMISMATCH)
    message("Setting -fallow-argument-mismatch")
    #add_compile_options("-fallow-argument-mismatch")
    target_compile_options(ARPACK PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
endif()

target_sources(ARPACK PUBLIC
    dgetv0.f dlaqrb.f dstqrb.f dsortc.f dsortr.f dstatn.f dstats.f
	dnaitr.f dnapps.f dnaup2.f dnaupd.f dnconv.f dneigh.f dngets.f
	dsaitr.f dsapps.f dsaup2.f dsaupd.f dsconv.f dseigt.f dsgets.f 
	dneupd.f dseupd.f dsesrt.f dvout.f  dmout.f 
	icnteq.f icopy.f iset.f iswap.f ivout.f second.f
)
target_link_libraries(ARPACK PUBLIC ${LAPACK_LIBRARIES})

# Required on Unix OS family to be able to be linked into shared libraries.
set_target_properties(${PROJECT_NAME}
                      PROPERTIES POSITION_INDEPENDENT_CODE ON)


# Expose public includes to other subprojects through cache variable.
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

export(TARGETS ARPACK
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/ARPACK.cmake"
       #EXPORT_LINK_INTERFACE_LIBRARIES #NAMESPACE MathFunctions::
)


