#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#------------------------------------------------------------------------------
#
# Major sections:
# - General Setup
# - External Libraries (Conan or User supplied paths)
# - Libraries
# - Executables
# - Option Application
#
#==============================================================================

cmake_minimum_required(VERSION 3.16)
project(OpenSeesFramework LANGUAGES C CXX Fortran)


set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(OPS_EXTERN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_BUNDLED_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_SRC_DIR "${PROJECT_SOURCE_DIR}/SRC/")

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE) # for the static libraries to be usable in OpenSeesPy
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake/cmake")
include(OpenSeesFunctions)

# include user config
include(${PROJECT_SOURCE_DIR}/Conf.cmake)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(LAPACK_VENDOR "OpenBLAS")
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(MKL_THREADING intel_thread)
endif()

include(FetchContent)

set(FETCHCONTENT_UPDATES_DISCONNECTED_MUMPSUPSTREAM true)
set(FETCHCONTENT_QUIET false)

# --- build MUMPS with options here
set(BUILD_SINGLE off)
set(BUILD_DOUBLE on)
set(BUILD_COMPLEX off)
set(BUILD_COMPLEX16 off)
set(MUMPS_parallel on)

find_package(MPI COMPONENTS C CXX Fortran REQUIRED)

# Ensure mumps uses its own FindScalapack module instead system config
set(OLD_CMAKE_FIND_PACKAGE_PREFER_CONFIG ${CMAKE_FIND_PACKAGE_PREFER_CONFIG})
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG FALSE)

FetchContent_Declare(mumps
  GIT_REPOSITORY https://github.com/scivision/mumps-superbuild.git
  GIT_TAG df62d7630760f3efd608f19c8e2002542897f7a9 # v5.8.2.1
)

FetchContent_MakeAvailable(mumps)

if(NOT DEFINED ENV{MKLROOT})
  # This means that the same FindLapack will be used for both OpenSees and MUMPS
  # but for some reason this doesn't work with MKL
  list(APPEND CMAKE_MODULE_PATH "${mumps_SOURCE_DIR}/cmake")
endif()

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ${OLD_CMAKE_FIND_PACKAGE_PREFER_CONFIG})

find_package(LAPACK REQUIRED)
find_package(Python COMPONENTS Interpreter Development)
find_package(HDF5 REQUIRED)
find_package(TCL REQUIRED)

include_directories(SYSTEM ${MPI_C_INCLUDE_DIRS})
include_directories(SYSTEM ${TCL_INCLUDE_DIRS})

#cmake -DCMAKE_RULE_MESSAGES:BOOL=OFF -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .

#
# external packages
#

add_library(OPS_External_packages INTERFACE)
add_library(OPS_OS_Specific_libs INTERFACE)

MESSAGE("${CMAKE_BINARY_DIR}")

# define user-selectable options for end target
set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS
    G3 OpenSees OpenSeesMP OpenSeesSP OpenSeesPy)

#
# OS Configuration:
#  defines: compile flags & user supplied paths
#

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  add_compile_definitions(_LINUX _UNIX _TCL85)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(OPS_Use_Graphics_Option None)
  add_compile_definitions(_LINUX _UNIX _TCL85)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    target_link_libraries(OPS_OS_Specific_libs INTERFACE wsock32 ws2_32)
    add_compile_definitions(_WIN32 _TCL85)
endif()

#==============================================================================
#                            General Setup
#
#==============================================================================
#----------------------------------------------------------------
# Compilers
#----------------------------------------------------------------

# Fortran
#--------------------------------------
enable_language(Fortran)

# C++
#--------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)
set(CMAKE_FORTRAN_STANDARD 08)

set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cqls <TARGET> <LINK_FLAGS> <OBJECTS>")
#add_compile_options(-g)

# Warnings
#opensees_add_cxx_flag(
#    GNU  -Wall
#    MSVC /W0
#)

# Floating-point
#opensees_add_cxx_flag(
#    GNU  -ffloat-store
#    MSVC /fp:precise
#)

#----------------------------------------------------------------
# Global Includes
#----------------------------------------------------------------
# include paths to main abstract classes

# NOTE BeamIntegration and MatrixUtil need to be removed from element/forceBEamColumn

include_directories(
  ${OPS_SRC_DIR}
  ${OPS_SRC_DIR}/matrix
  ${OPS_SRC_DIR}/handler
  ${OPS_SRC_DIR}/database
  ${OPS_SRC_DIR}/element
  ${OPS_SRC_DIR}/element/forceBeamColumn
  ${OPS_SRC_DIR}/element/nonlinearBeamColumn/matrixutil
  ${OPS_SRC_DIR}/coordTransformation
  ${OPS_SRC_DIR}/damping
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/tagged/storage
  ${OPS_SRC_DIR}/recorder
  ${OPS_SRC_DIR}/renderer
  ${OPS_SRC_DIR}/damage
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_SRC_DIR}/material
  ${OPS_SRC_DIR}/material/section
  ${OPS_SRC_DIR}/material/uniaxial
  ${OPS_SRC_DIR}/material/nD
  ${OPS_SRC_DIR}/graph/graph
  ${OPS_SRC_DIR}/graph/numberer
  ${OPS_SRC_DIR}/graph/partitioner
  ${OPS_SRC_DIR}/domain/component
  ${OPS_SRC_DIR}/domain/domain
  ${OPS_SRC_DIR}/domain/subdomain
  ${OPS_SRC_DIR}/domain/load
  ${OPS_SRC_DIR}/domain/loadBalancer
  ${OPS_SRC_DIR}/domain/pattern
  ${OPS_SRC_DIR}/domain/groundMotion
  ${OPS_SRC_DIR}/domain/node
  ${OPS_SRC_DIR}/domain/constraints
  ${OPS_SRC_DIR}/domain/region
  ${OPS_SRC_DIR}/analysis/algorithm
  ${OPS_SRC_DIR}/analysis/dof_grp
  ${OPS_SRC_DIR}/analysis/fe_ele
  ${OPS_SRC_DIR}/analysis/algorithm/equiSolnAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/eigenAlgo
  ${OPS_SRC_DIR}/analysis/algorithm/domainDecompAlgo
  ${OPS_SRC_DIR}/analysis/analysis
  ${OPS_SRC_DIR}/analysis/integrator
  ${OPS_SRC_DIR}/analysis/handler
  ${OPS_SRC_DIR}/analysis/numberer
  ${OPS_SRC_DIR}/analysis/model
  ${OPS_SRC_DIR}/convergenceTest
  ${OPS_SRC_DIR}/modelbuilder
  ${OPS_SRC_DIR}/system_of_eqn
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE
  ${OPS_SRC_DIR}/system_of_eqn/eigenSOE
  ${OPS_SRC_DIR}/actor/actor
  ${OPS_SRC_DIR}/actor/channel
  ${OPS_SRC_DIR}/actor/objectBroker
  ${OPS_SRC_DIR}/actor/message
  ${OPS_BUNDLED_DIR}/eigenAPI
)

# TODO: Temporary fix; include all directories with .h files through glob
file(GLOB_RECURSE ops_include_files CONFIGURE_DEPENDS
    ${OPS_SRC_DIR}/*.h
    ${OPS_BUNDLED_DIR}/*.h
)
set(OPS_INCLUDE_DIRS "")
foreach(filepath ${ops_include_files})
    get_filename_component(dir_path ${filepath} PATH)
    set(OPS_INCLUDE_DIRS ${OPS_INCLUDE_DIRS} ${dir_path})
endforeach()
list(REMOVE_DUPLICATES OPS_INCLUDE_DIRS)

include_directories(${OPS_INCLUDE_DIRS})
# end TODO

# ##############################################################
#
# NUMERICAL and OpenSees LIBRARIES
#
# ##############################################################

#
# NUMERICAL LIBRARIES
#

add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_5.1.1")
set (SUPERLU_LIBRARIES SUPERLU)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/UMFPACK")
set (UMFPACK_LIBRARIES UMFPACK)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/AMD")
set (AMD_LIBRARIES AMD)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/ARPACK")
set (ARPACK_LIBRARIES ARPACK)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/CSPARSE")
set (CSPARSE_LIBRARIES CSPARSE)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/tetgen1.4.3")
set (TETGEN_LIBRARIES tet)
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/Triangle")
set (TRIANGLE_LIBRARIES triangle)
# add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/eigenAPI")
# set (EIGENAPI_LIBRARIES EIGENAPI)


set (MUMPS_FLAG -D_MUMPS)

add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_DIST_4.3/SRC")
add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/METIS")


#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()



add_library(OPS_Numerics INTERFACE)

target_link_libraries(OPS_Numerics INTERFACE
       LAPACK::LAPACK
       ${ARPACK_LIBRARIES}
       ${CSPARSE_LIBRARIES}
       ${SUPERLU_LIBRARIES}
       ${UMFPACK_LIBRARIES}
       ${TETGEN_LIBRARIES}
       ${TRIANGLE_LIBRARIES}
       ${AMD_LIBRARIES}
       ${AMD_LIBRARIES}
       ${EIGENAPI_LIBRARIES}
)

set(TCL_LIBRARIES ${TCL_LIBRARY})


#
# OpenSees OBJECT LIBRARIES
#

add_library(OPS_Matrix             OBJECT)
add_library(OPS_Actor              OBJECT)
add_library(OPS_ObjectBroker       OBJECT)
add_library(OPS_Handler            OBJECT)
add_library(OPS_Recorder           OBJECT)
add_library(OPS_Reliability        OBJECT)
add_library(OPS_Tagged             OBJECT)
add_library(OPS_Utilities          OBJECT)
add_library(OPS_ModelBuilder       OBJECT)
add_library(OPS_Domain             OBJECT)
add_library(OPS_SysOfEqn           OBJECT)
if(Python_FOUND)
  add_library(OPS_SysOfEqnPython   OBJECT)
endif()
#add_library(OPS_SysOfEqnSP        OBJECT)
add_library(OPS_Analysis           OBJECT)
add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_Element            OBJECT)
add_library(OPS_ElementFortran     OBJECT)
#add_library(OPS_ElementSP         OBJECT)
add_library(OPS_Material           OBJECT)
add_library(OPS_MaterialFortran    OBJECT)
add_library(OPS_Damage             OBJECT)
add_library(OPS_Database           OBJECT)
add_library(OPS_INTERPRETER        OBJECT)



#
# Optional Extensions
#

add_library(OPS_Paraview           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Renderer           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_Default   OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Graphics_GL        OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_PFEM               OBJECT EXCLUDE_FROM_ALL)

#
# Tcl Interpreter Library & special include and link directives
#

add_library(OPS_InterpTcl          STATIC)

target_compile_definitions(OPS_InterpTcl PUBLIC _TCL85)
target_link_libraries(OPS_InterpTcl tcl::tcl)

#
# OpenSees G3 Library - a slimmed down OpenSees
#

add_library(G3)

target_link_libraries(G3
    coordTransformation
    damping
    OPS_Matrix
    OPS_Analysis
    OPS_ModelBuilder
    OPS_Domain
    OPS_ConvergenceTest
    OPS_Element
    OPS_Material
    OPS_Recorder
    OPS_Handler
    OPS_SysOfEqn
    OPS_Tagged
    OPS_Utilities
    graph
    OPS_Actor
    OPS_ObjectBroker
    OPS_Numerics
    OPS_INTERPRETER
)

#
# OpenSees Libray
#

add_library(OpenSeesLIB  EXCLUDE_FROM_ALL)

target_link_libraries(OpenSeesLIB
    ${OPS_Extension_List}
    OPS_Actor
    OPS_Analysis
    OPS_ConvergenceTest
    OPS_Damage
    OPS_Database
    OPS_Domain
    graph
    OPS_Element
    OPS_ElementFortran
    OPS_Handler
    OPS_INTERPRETER
    OPS_Material
    OPS_MaterialFortran
    OPS_Material_YieldSurface
    OPS_Matrix
    OPS_ModelBuilder
    OPS_ObjectBroker
    OPS_PFEM
    OPS_Recorder
    coordTransformation
    damping
    OPS_Renderer
    OPS_Section_Repres
    OPS_Section_YieldSurface
    OPS_SysOfEqn
    OPS_Thermal
    OPS_OS_Specific_libs
    OPS_Tagged
    OPS_Utilities
)


# ########################################################################
#
# EXECUTABLES
#
# ########################################################################

#
# OpenSees Tcl Interpreter
#

add_executable(OpenSees
  ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
  ${OPS_SRC_DIR}/tcl/commands.cpp
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
)

target_link_libraries(OpenSees
  OPS_InterpTcl
  coordTransformation
  OpenSeesLIB
  OPS_Reliability
  OPS_ReliabilityTcl
  OPS_Numerics
  OPS_Recorder
  ${CMAKE_DL_LIBS}
  HDF5::HDF5
)


#
# OpenSeesSP Tcl Parallel Interpreter
#


add_executable(OpenSeesSP
  ${OPS_SRC_DIR}/tcl/mpiMain.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
  ${OPS_SRC_DIR}/tcl/commands.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomain.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainEleIter.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainSubIter.cpp
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
  ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenColLinSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenRowLinSOE.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
)


if (DEFINED OPENMPI)
    target_compile_definitions(OpenSeesSP
    PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG} _OPENMPI)
else()
    target_compile_definitions(OpenSeesSP
    PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG})
endif()

target_link_libraries(OpenSeesSP
    OPS_InterpTcl
    OpenSeesLIB
    OPS_Reliability
    OPS_ReliabilityTcl
    OPS_Recorder
    METIS
    SUPERLU_DIST
    OPS_Numerics
    MUMPS::MUMPS
    ${CMAKE_DL_LIBS}
    HDF5::HDF5
    MPI::MPI_CXX
    MPI::MPI_Fortran
)


#
# OpenSeesMP Tcl MPI Interpreter
#

add_executable(OpenSeesMP
  ${OPS_SRC_DIR}/tcl/mpiParameterMain.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
  ${OPS_SRC_DIR}/tcl/commands.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomain.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainEleIter.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainSubIter.cpp
  ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
)

if (DEFINED OPENMPI)
    target_compile_definitions(OpenSeesMP
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG} _OPENMPI)
else()
    target_compile_definitions(OpenSeesMP
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG})
endif()

target_link_libraries(OpenSeesMP
    OPS_InterpTcl
    OpenSeesLIB
    OPS_Reliability
    OPS_ReliabilityTcl
    OPS_Recorder
    METIS
    SUPERLU_DIST
    OPS_Numerics
    MUMPS::MUMPS
    ${CMAKE_DL_LIBS}
    HDF5::HDF5
    MPI::MPI_CXX
    MPI::MPI_Fortran
)


#
# OpenSeesPy Module
#
if(Python_FOUND)
	add_library(OpenSeesPy SHARED
		${OPS_SRC_DIR}/interpreter/PythonModule.cpp
		${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
	)

	set_target_properties(OpenSeesPy PROPERTIES PREFIX "")

	target_link_libraries(OpenSeesPy
	   OpenSeesLIB
     OPS_SysOfEqnPython
	   OPS_Reliability
	   OPS_Recorder
	   OPS_Numerics
	   HDF5::HDF5
	   Python::Module
	)
endif()

#
# INSTALL
#

install(TARGETS OpenSees OpenSeesSP OpenSeesMP)
set (TCL_INIT_FILE "${TCL_INCLUDE_DIR}/../lib/tcl8.6/init.tcl")
install(FILES ${TCL_INIT_FILE} DESTINATION lib/tcl8.6)

if(WIN32)
  if (NOT DEFINED ENV{ONEAPI_ROOT})
    message(FATAL_ERROR "ONEAPI_ROOT environment variable not set.")
  endif()

  file(TO_CMAKE_PATH "$ENV{ONEAPI_ROOT}/bin" ONEAPI_BIN_DIR)

  list(APPEND SYSTEM_RUNTIME_LIBS_COMMON
      "${ONEAPI_BIN_DIR}/libmmd.dll"
      "${ONEAPI_BIN_DIR}/libiomp5md.dll"
      "${ONEAPI_BIN_DIR}/tbbmalloc.dll"
      "${ONEAPI_BIN_DIR}/mkl_intel_thread.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_sequential.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_core.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_avx2.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_avx512.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_def.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_mc3.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_vml_avx2.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_vml_avx512.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_vml_cmpt.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_vml_def.2.dll"
      "${ONEAPI_BIN_DIR}/mkl_vml_mc3.2.dll"
  )

  set(SYSTEM_RUNTIME_LIBS
      ${SYSTEM_RUNTIME_LIBS_COMMON}
      "${ONEAPI_BIN_DIR}/libifcoremd.dll"
  )

  set(SYSTEM_RUNTIME_LIBS_DEBUG
      ${SYSTEM_RUNTIME_LIBS_COMMON}
      "${ONEAPI_BIN_DIR}/libifcoremdd.dll"
  )
endif()

function(resolve_paths input_and_output_list)
  foreach(lib_file IN LISTS ${input_and_output_list})
    get_filename_component(real_file "${lib_file}" REALPATH)
    if(NOT EXISTS "${real_file}")
      message(FATAL_ERROR "Runtime library not found: ${real_file}")
    endif()
    list(APPEND resolved_list "${real_file}")
  endforeach()
  set(${input_and_output_list} ${resolved_list} PARENT_SCOPE)
endfunction()

resolve_paths(SYSTEM_RUNTIME_LIBS)
resolve_paths(SYSTEM_RUNTIME_LIBS_DEBUG)


if(WIN32)
  set(INSTALL_SYSTEM_RUNTIME_DESTINATION bin)
else()
  set(INSTALL_SYSTEM_RUNTIME_DESTINATION lib)
endif()

# Handle both multi-config generators (Visual Studio, Xcode) and single-config generators (Ninja, Unix Makefiles)
if(CMAKE_CONFIGURATION_TYPES)
  foreach(config IN CMAKE_CONFIGURATION_TYPES)
    if (config STREQUAL "Debug")
      install(PROGRAMS ${SYSTEM_RUNTIME_LIBS_DEBUG}
        DESTINATION ${INSTALL_SYSTEM_RUNTIME_DESTINATION}
        CONFIGURATIONS ${config}
      )
    else()
      install(PROGRAMS ${SYSTEM_RUNTIME_LIBS}
        DESTINATION ${INSTALL_SYSTEM_RUNTIME_DESTINATION}
        CONFIGURATIONS ${config}
      )
    endif()
  endforeach()
else()
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    install(PROGRAMS ${SYSTEM_RUNTIME_LIBS_DEBUG}
      DESTINATION ${INSTALL_SYSTEM_RUNTIME_DESTINATION}
    )
  else()
    install(PROGRAMS ${SYSTEM_RUNTIME_LIBS}
      DESTINATION ${INSTALL_SYSTEM_RUNTIME_DESTINATION}
    )
  endif()
endif()


if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  install(CODE [[
    set(EXES "$<TARGET_FILE:OpenSees>;$<TARGET_FILE:OpenSeesSP>;$<TARGET_FILE:OpenSeesMP>")
    # Collect dependencies of the installed targets
    file(GET_RUNTIME_DEPENDENCIES
        RESOLVED_DEPENDENCIES_VAR _deps
        UNRESOLVED_DEPENDENCIES_VAR _unresolved
        EXECUTABLES ${EXES}
        # Picked from https://github.com/AppImageCommunity/pkg2appimage/blob/19e30b276ffedf4d3b4b56bc6320f463625a74f8/excludelist
        PRE_EXCLUDE_REGEXES
            "ld-linux(-x86-64)?.so.2"
            "libdl.so.2"
            "libc.so.6"
            "libm.so.6"
            "libpthread.so.0"
            "libmpi" # Exclude MPI runtime libraries, which are expected to be provided by the user in their environment
    )

    if(_unresolved)
        message(WARNING "Unresolved runtime dependencies: ${_unresolved}")
    endif()

    set(_libdir "${CMAKE_INSTALL_PREFIX}/lib")
    file(MAKE_DIRECTORY "${_libdir}")

    set(_installed_libs)
    foreach(src IN LISTS _deps)
        get_filename_component(_name "${src}" NAME)
        file(INSTALL
            DESTINATION "${_libdir}"
            TYPE SHARED_LIBRARY
            FOLLOW_SYMLINK_CHAIN
            FILES "${src}"
        )
        list(APPEND _installed_libs "${_libdir}/${_name}")
    endforeach()

    find_program(PATCHELF_EXECUTABLE patchelf)
    if(NOT PATCHELF_EXECUTABLE)
        message(FATAL_ERROR "patchelf not found, cannot rewrite RPATH of third-party libs")
    endif()

    foreach(f IN LISTS _installed_libs)
        execute_process(
            COMMAND "${PATCHELF_EXECUTABLE}" --set-rpath "\$ORIGIN" "${f}"
            RESULT_VARIABLE _rc
        )
        if(NOT _rc EQUAL 0)
            message(FATAL_ERROR "patchelf failed on ${f}")
        endif()
    endforeach()
  ]])
endif()

# Remove extra files from mumps install
install(CODE [[
  file(REMOVE_RECURSE
       "${CMAKE_INSTALL_PREFIX}/cmake"
       "${CMAKE_INSTALL_PREFIX}/include"
       "${CMAKE_INSTALL_PREFIX}/share")

  set(libdir "${CMAKE_INSTALL_PREFIX}/lib")

  file(GLOB static_libs
       "${libdir}/*.a"
       "${libdir}/*.lib")

  if(static_libs)
    message(STATUS "Removing static libraries: ${static_libs}")
    file(REMOVE ${static_libs})
  endif()
]])


# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})




#==============================================================================
#                            Configure targets
#
#==============================================================================

#target_link_libraries(OPS_Numerics     INTERFACE ${OPS_Numlib_List})
#message("Numlibs: ${OPS_Numlib_List}")
#target_link_libraries(OPS_Element      PRIVATE ${OPS_Element_List})
#target_link_libraries(OPS_ObjectBroker PRIVATE ${OPS_Element_List})
#target_link_libraries(OPS_InterpTcl    PRIVATE ${OPS_Element_List})



#==============================================================================
#                            Apply Options
#
#==============================================================================
#----------------------------
# FMK
#----------------------------
if(FMK)
    add_compile_definitions(
        _HAVE_Damage2p
        _HAVE_PSUMAT
        _HAVE_PML
        _FILIP_LHNMYS
    )
endif()

#----------------------------
# Extensions
#----------------------------
message("OPS >>> Configuring OpenSees extensions")
foreach(extension IN LISTS OPS_SysOfEqn_List OPS_Element_List OPS_Extension_List)
    string(TOUPPER "${extension}" ext_flag)
    string(REGEX REPLACE "^OPS_" "OPSDEF_" ext_flag "${ext_flag}")
    add_compile_definitions(${ext_flag})
endforeach()
foreach(extension IN LISTS OPS_Exclude_List)
    string(TOUPPER "${extension}" ext_flag)
    string(REGEX REPLACE "^OPS_" "OPS_EXCLUDE_" ext_flag "${ext_flag}")
    message("    Adding macro definition '${ext_flag}'")
    add_compile_definitions(${ext_flag})
endforeach()


# Renderer
#----------------------------
if (${OPS_Use_Graphics_Option} STREQUAL "Base")
  target_link_libraries(${OPS_FINAL_TARGET} OPS_Graphics_Default OPS_Graphics)

elseif (${OPS_Use_Graphics_Option} STREQUAL "OpenGL")
  message("OPS >>> Including OpenGL graphics option")
  set_source_files_properties(
      "${OPS_SRC_DIR}/recorder/AlgorithmIncrements.cpp"
      "${OPS_SRC_DIR}/recorder/FilePlotter.cpp"
      "${OPS_SRC_DIR}/renderer/main.cpp"
      "${OPS_SRC_DIR}/renderer/OpenGlDevice.cpp"
      "${OPS_SRC_DIR}/renderer/OpenGlRenderer.cpp"
      "${OPS_SRC_DIR}/tcl/TclFeViewer.cpp"
      "${OPS_SRC_DIR}/tcl/TclVideoPlayer.cpp"

    PROPERTIES COMPILE_DEFINITIONS _AGL
  )
  target_link_libraries(${OPS_FINAL_TARGET} OPS_Graphics_GL OPS_Graphics)
else()
  add_compile_definitions(_NOGRAPHICS)
endif()


# Reliability
#----------------------------
#if ("OPS_Reliability" IN_LIST OPS_Extension_List)
    add_compile_definitions(_RELIABILITY)
    #target_link_libraries(${OPS_FINAL_TARGET} OPS_Reliability)
#endif()


#----------------------------
# HDF5
#----------------------------

add_compile_definitions(_H5DRM)
add_compile_definitions(_HDF5)
include_directories(${HDF5_INCLUDE_DIR})

#----------------------------
# Eigen3
#----------------------------
if(Eigen3_FOUND)
   include_directories(${Eigen3_INCLUDE_DIR})
   add_compile_definitions(_EIGEN3)
else()
   message(STATUS "OPS >>> Could not find Eigen3")
endif()

if (OPS_Use_Dev_Directories)
  add_subdirectory("${PROJECT_SOURCE_DIR}/DEVELOPER/")
endif()
