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
# - External Libreries (Conan or User supplied paths)
# - Libraries
# - Executables
# - Option Application
#
#==============================================================================

cmake_minimum_required(VERSION 3.23)
project(OpenSeesFramework LANGUAGES C CXX Fortran)


set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(OPS_EXTERN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_EXTERNALS_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_BUNDLED_DIR "${PROJECT_SOURCE_DIR}/OTHER/")
set(OPS_SRC_DIR "${PROJECT_SOURCE_DIR}/SRC/")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake/cmake")
#What is the purpose of "OpenSeesFunctions"... adds unnecessary steps to build
#include(OpenSeesFunctions)

# include user config
include(${PROJECT_SOURCE_DIR}/Conf.cmake)

#
# Compiler SPECIFIC
#

# On CMake 3.25 or older CXX_COMPILER_FRONTEND_VARIANT is not always set

if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "")
  set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "${CMAKE_CXX_COMPILER_ID}")
endif()
MESSAGE(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
MESSAGE(STATUS "CMAKE_CXX_SIMULATE_ID: ${CMAKE_CXX_SIMULATE_ID}")
MESSAGE(STATUS "CMAKE_CXX_COMPILER_FRONTEND_VARIANT: ${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}")

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(BUILD_SHARED_LIBS OFF)

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    if ("${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" MATCHES "MSVC")
      find_package(IntelSYCL REQUIRED)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
      set(BUILD_SHARED_LIBS OFF)
      set(CMAKE_POSITION_INDEPENDENT_CODE ON)
      set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")        
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MT /W0 /nologo")
    else()
      add_compile_options(-fPIC)
    endif()
    
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
    set(BUILD_SHARED_LIBS OFF)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")

else()
    MESSAGE("COMPILER: UKNOWN COMPILER ${CMAKE_CXX_COMPILER_ID}")
endif()
 
#
# external packages
#

add_library(OPS_External_packages INTERFACE)
add_library(OPS_OS_Specific_libs INTERFACE)

# 
# Conan (manages external dependencies for a typical build)
# 

MESSAGE(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")

   MESSAGE("USING CONAN")
   include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
   conan_basic_setup()
   set (HDF5_FOUND TRUE)
   set (HDF5_LIBRARIES ${CONAN_LIBS_HDF5} ${CONAN_LIBS_ZLIB})
   set (HDF5_VERSION "1.12.0")
   set(USING_CONAN TRUE)
   set(CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR} ${CMAKE_MODULE_PATH})
   set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR} ${CMAKE_PREFIX_PATH})

else()

  MESSAGE("NOT USING CONAN")
  include(OpenSeesFunctions)
  set(NOT_USING_CONAN TRUE)
  set (CONAN_LIBS )
  find_package(HDF5 REQUIRED)
  find_package(TCL REQUIRED)

endif()

# define user-selectable options for end target
set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS 
    G3 OpenSees OpenSeesMP OpenSeesSP OpenSeesPy)

#
# OS Configuration:
#  defines: compile flags & user supplied paths and includes if not using conan
#

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  add_compile_definitions(_LINUX _UNIX _TCL85)
  if (NOT_USING_CONAN)
    include(OpenSeesDependenciesUnix)
  endif()
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(OPS_Use_Graphics_Option None)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  add_compile_definitions(_LINUX _UNIX _TCL85)
  if (NOT_USING_CONAN)
    include(OpenSeesDependenciesUnix)
  endif()  
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  if(CMAKE_GENERATOR MATCHES "Ninja")
    set(CMAKE_Fortran_PREPROCESS OFF)
  endif()
  if (NOT_USING_CONAN)
    include(OpenSeesDependenciesWin)
  endif()      
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
#enable_language(Fortran)   # Redundant project(... Fortran ...) does this already

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
# include all directories with CMakeLists.txt files through glob
# Searching for *.h masks when CMakeLists.txt are not setup

file(GLOB_RECURSE ops_include_files CONFIGURE_DEPENDS 
    ${OPS_SRC_DIR}/CMakeLists.txt
    ${OPS_BUNDLED_DIR}/CMakeLists.txt
)
set(OPS_INCLUDE_DIRS "")
foreach(filepath ${ops_include_files})
    get_filename_component(dir_path ${filepath} PATH)
    set(OPS_INCLUDE_DIRS ${OPS_INCLUDE_DIRS} ${dir_path})
endforeach()
list(REMOVE_DUPLICATES OPS_INCLUDE_DIRS)

include_directories(${OPS_INCLUDE_DIRS})

# ##############################################################
#
# NUMERICAL and OpenSees LIBRARIES
#
# ##############################################################

#
# NUMERICAL LIBRARIES
#

###
# To allow the CMake build to work better with existing package or module build 
# workflows, check if the user has defined the location of any of the following
# external packages before adding the bundled "/OTHER/<package>" subdirectories 
# to the build.
# Even if a "FIND_<PACKAGE>" module does not currently exist in CMAKE_ROOT or 
# CMAKE_MODULE_PATH, add the command anyways (with "QUIET" designation) incase 
# one is added in future updates
###

#SUPERLU
MESSAGE(STATUS "Looking for SUPERLU...")
find_package(SUPERLU )
if(NOT SUPERLU_FOUND)
  MESSAGE(STATUS "SUPERLU Not Found... Use bundled version")
  set (SUPERLU_FOUND TRUE CACHE BOOL "providing own SUPERLU lib")
  unset(SUPERLU_LIBRARIES CACHE)
  unset(SUPERLU_INCLUDES CACHE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_5.1.1")
  set (SUPERLU_LIBRARIES SUPERLU)
endif()

#UMFPACK
MESSAGE(STATUS "Looking for User UMFPACK...${UMFPACK}")
find_package(UMFPACK QUIET)
if(NOT UMFPACK_FOUND)
  MESSAGE(STATUS "User UMFPACK Not Found... Use bundled version")
  set (UMFPACK_FOUND TRUE CACHE BOOL "providing own UMFPACK lib")
  unset(UMFPACK_DIR CACHE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/UMFPACK")
  set (UMFPACK_LIBRARIES UMFPACK)
endif()

#AMD
MESSAGE(STATUS "Looking for User AMD...${AMD}")
find_package(AMD )
if(NOT AMD_FOUND)
  MESSAGE(STATUS "User AMD Not Found... Use bundled version")
  set (AMD_FOUND TRUE CACHE BOOL "providing own AMD lib")
  unset(AMD_AMD_INCLUDE_DIR CACHE)
  unset(AMD_LIBRARY CACHE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/AMD")
  set (AMD_LIBRARIES AMD)
endif()

#ARPACK
MESSAGE(STATUS "Looking for User ARPACK...")
find_package(ARPACK )
if(NOT ARPACK_FOUND)
  MESSAGE(STATUS "User ARPACK Not Found... Use bundled version")
  set (ARPACK_FOUND TRUE CACHE BOOL "providing own ARPACK lib")
  unset(ARPACK_LIBRARIES CACHE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/ARPACK")
  set (ARPACK_LIBRARIES ARPACK)
endif()

#CSPARSE
MESSAGE(STATUS "Looking for User CSPARSE...")
find_package(CSPARSE QUIET)
if(NOT CSPARSE_FOUND)
  MESSAGE(STATUS "User CSPARSE Not Found... Use bundled version")
  set (CSPARSE_FOUND TRUE CACHE BOOL "providing own CSPARSE lib")
  unset(CSPARSE_DIR CACHE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/CSPARSE")
  set (CSPARSE_LIBRARIES CSPARSE)
endif()

#TETGEN
MESSAGE(STATUS "Looking for User TETGEN...")
find_package(TETGEN QUIET)
if(NOT TETGEN_FOUND)
  MESSAGE(STATUS "User TETGEN Not Found... Use bundled version")
  set (TETGEN_FOUND TRUE CACHE BOOL "providing own TETGEN lib")
  unset(TETGEN_DIR CACHE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/tetgen1.4.3")
  set (TETGEN_LIBRARIES tet)
endif()

#TRIANGLE
MESSAGE(STATUS "Looking for User TRIANGLE...")
find_package(TRIANGLE QUIET)
if(NOT TRIANGLE_FOUND)
  MESSAGE(STATUS "User TRIANGLE Not Found... Use bundled version")
  set (TRIANGLE_FOUND TRUE CACHE BOOL "providing own TRIANGLE lib")
  unset(TRIANGLE_DIR CACHE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/Triangle")
  set (TRIANGLE_LIBRARIES triangle)
endif()

#LAPACK
MESSAGE(STATUS "Looking for LAPACK...")
find_package(LAPACK)
if(NOT LAPACK_FOUND)  
  MESSAGE(STATUS "LAPACK Not Found...")
  set (LAPACK_FOUND TRUE CACHE BOOL "providing own LAPACK lib")
endif()

find_package (Python COMPONENTS Interpreter Development)


find_package(MPI)
find_package(MKL)


if(LAPACK_FOUND)
  message(STATUS "LAPACK was found.")
  message(STATUS "LAPACK_LINKER_FLAGS = ${LAPACK_LINKER_FLAGS}")
  message(STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}" )
else()
  MESSAGE(STATUS "LAPACK NOT FOUND ... Use bundled version")
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/BLAS") 
  add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/LAPACK")
  set(LAPACK_LIBRARIES LAPACK)
endif()

if(PYTHON_FOUND)
  message(STATUS "Python_FOUND:${Python_FOUND}")
  message(STATUS "Python_VERSION:${Python_VERSION}")
  message(STATUS "Python_Development_FOUND:${Python_Development_FOUND}")
  message(STATUS "Python_LIBRARIES:${Python_LIBRARIES}")
  message(STATUS "Python_INCLUDES:${Python_INCLUDE_DIRS}")  
else()
  message(STATUS "PYTHON NOT FOUND")
endif()


if(MPI_FOUND)
  message(STATUS "MPI was found .. path added ${MPI_C_INCLUDE_DIRS} ${MPI_FOUND}")
  include_directories(SYSTEM ${MPI_C_INCLUDE_DIRS})

  #MUMPS
  MESSAGE(STATUS "Looking for MUMPS...")
  set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${MUMPS_DIR}")
  find_package(MUMPS)
  if (NOT MUMPS_FOUND)
    # Add MUMPS as external project. Downloads source, builds it, and installs it in
    # top-level project build directory for linking against.
    message(STATUS "MUMPS NOT Found. Build it from source...")
    unset(MUMPS_DIR CACHE)
    include(ExternalProject)
      ExternalProject_Add(mumps
      GIT_REPOSITORY    https://github.com/OpenSees/mumps.git
      GIT_TAG           main
      GIT_SHALLOW       True
      PREFIX            "mumps-extern"
      CMAKE_CACHE_ARGS   "-Darith:STRING=d"
      BUILD_ALWAYS      True
      BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config <CONFIG>
      INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR> --prefix <INSTALL_DIR>
    )
    ExternalProject_Get_property(mumps BINARY_DIR)
    set (MUMPS_FLAG -D_MUMPS)  
    set(MUMPS_DIR ${BINARY_DIR})
    if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      set (MUMPS_LIBRARIES "${BINARY_DIR}/dmumps;${BINARY_DIR}/mumps_common;${BINARY_DIR}/pord")
    else()
      set (MUMPS_LIBRARIES "${BINARY_DIR}/libdmumps.a;${BINARY_DIR}/libmumps_common.a;${BINARY_DIR}/libpord.a")
    endif()
    message(STATUS "MUMPS_LIBRARIES: ${MUMPS_LIBRARIES}")
  else()
    set (MUMPS_FLAG -D_MUMPS)  
  endif()

  #SUPERLU_DIST
  MESSAGE(STATUS "Looking for SUPERLU_DIST...")
  find_package(SUPERLU_DIST QUIET)
  if(NOT SUPERLU_DIST_FOUND)
    MESSAGE(STATUS "SUPERLU_DIST Not Found... Use bundled version")
    set (SUPERLU_DIST_FOUND TRUE CACHE BOOL "providing own SUPERLU_DIST lib")
    unset(SUPERLU_DIST_DIR CACHE)
    add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/SuperLU_DIST_4.3")
  endif()
  
  #METIS
  MESSAGE(STATUS "Looking for METIS...")
  find_package(METIS QUIET)
  if(NOT METIS_FOUND)
    MESSAGE(STATUS "METIS Not Found... Use bundled version")
    set (METIS_FOUND TRUE CACHE BOOL "providing own METIS lib")
    unset(METIS_INCLUDES CACHE)
    unset(METIS_LIBRARIES CACHE)
    add_subdirectory("${PROJECT_SOURCE_DIR}/OTHER/METIS")
  endif()

else()
  message(STATUS "MPI was NOT found.")
endif()

if(MKL_FOUND)
  message(STATUS "MKL was found.")
  set (SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set (MKL_LPATH ${MKL_ROOT}/lib/intel64)
    set (SCALAPACK_LIBRARIES "${MKL_LPATH}/mkl_scalapack_ilp64.lib;${MKL_LPATH}/mkl_intel_ilp64.lib;${MKL_LPATH}/mkl_sequential.lib;${MKL_LPATH}/mkl_core.lib;${MKL_LPATH}/mkl_blacs_intelmpi_ilp64.lib")
  endif() 
  message(STATUS "SCALAPACK_LIBRARIES=${SCALAPACK_LIBRARIES}")     
else()
  message(STATUS "MKL NOT found .. user to provide -DSCALAPACK_LIBRARIES=")

  #SCALAPACK
  MESSAGE(STATUS "Looking for SCALAPACK...")
  find_package(SCALAPACK QUIET)
  if(SCALAPACK_FOUND)
    MESSAGE(STATUS "SCALAPACK Found")
    message(STATUS "SCALAPACK_LIBRARIES=${SCALAPACK_LIBRARIES}")
  endif()
  
endif()

add_library(OPS_Numerics INTERFACE)

target_link_libraries(OPS_Numerics INTERFACE
       ${ARPACK_LIBRARIES}
       ${CSPARSE_LIBRARIES} 
       ${SUPERLU_LIBRARIES}
       ${UMFPACK_LIBRARIES}
       ${TETGEN_LIBRARIES}
       ${TRIANGLE_LIBRARIES}                     
       ${AMD_LIBRARIES}   
       ${LAPACK_LIBRARIES}
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


include_directories(${TCL_INCLUDE_PATH})

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
target_include_directories(OPS_InterpTcl PUBLIC ${TCL_INCLUDE_PATH})
target_link_libraries(OPS_InterpTcl PRIVATE ${TCL_LIBRARIES})

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

# Adding EXCLUDE_FROM_ALL prevents these from being made by default.

#
# OpenSees Tcl Interpreter
#



add_executable(OpenSees EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/tcl/tclAppInit.cpp
  ${OPS_SRC_DIR}/tcl/tclMain.cpp
  ${OPS_SRC_DIR}/tcl/commands.cpp
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
)

if(NOT EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
  target_include_directories(OpenSees PUBLIC ${TCL_INCLUDE_PATH})

  FILE(COPY ${CONAN_LIB_DIRS_TCL}/tcl8.6
       DESTINATION ${PROJECT_BINARY_DIR}/lib/tcl8.6
       FILES_MATCHING PATTERN .tcl)
     
else ()
     FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6)
     add_custom_command(
      TARGET OpenSees POST_BUILD
  COMMENT("Copying init.tcl to ${CMAKE_CURRENT_BINARY_DIR}")
      COMMAND ${CMAKE_COMMAND} -E copy
        ${CONAN_LIB_DIRS_TCL}/tcl8.6/init.tcl
  ${CMAKE_CURRENT_BINARY_DIR}/lib/tcl8.6
     )
endif()


target_link_libraries(OpenSees
  OPS_InterpTcl 
  coordTransformation
  OpenSeesLIB
  OPS_Reliability
  OPS_ReliabilityTcl  
  OPS_Numerics
  OPS_Recorder
  ${CMAKE_DL_LIBS} 
  ${HDF5_LIBRARIES} 
  ${CONAN_LIBS}
)


#
# OpenSeesSP Tcl Parallel Interpreter
#


if(MPI_FOUND)

    add_executable(OpenSeesSP EXCLUDE_FROM_ALL 
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


    target_include_directories(OpenSeesSP PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})
    target_compile_options(OpenSeesSP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

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
       ${MUMPS_LIBRARIES}
       ${CMAKE_DL_LIBS} 
       ${HDF5_LIBRARIES} 
       ${CONAN_LIBS}
       ${MPI_CXX_LIBRARIES}
       ${SCALAPACK_LIBRARIES}
       ${MPI_Fortran_LIBRARIES}       
       ${MPI_CXX_LINK_FLAGS}
    )
    
endif()


#
# OpenSeesMP Tcl MPI Interpreter
#

if(MPI_FOUND)

    add_executable(OpenSeesMP EXCLUDE_FROM_ALL 
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

    target_include_directories(OpenSeesMP PRIVATE ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})

    target_compile_options(OpenSeesMP PRIVATE ${MPI_CXX_COMPILE_FLAGS})
    
    if (NOT DEFINED MUMPS_DIR)
        add_dependencies(OpenSeesMP mumps)
    endif()         

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
       ${MUMPS_LIBRARIES}
       ${CMAKE_DL_LIBS} 
       ${HDF5_LIBRARIES} 
       ${CONAN_LIBS}
       ${MPI_CXX_LIBRARIES}
       ${SCALAPACK_LIBRARIES}
       ${MPI_Fortran_LIBRARIES}       
       ${MPI_CXX_LINK_FLAGS}
    )
    
endif()


#
# OpenSeesPy Module
#

add_library(OpenSeesPy SHARED EXCLUDE_FROM_ALL
    ${OPS_SRC_DIR}/interpreter/PythonModule.cpp
    ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
)

set_target_properties(OpenSeesPy PROPERTIES PREFIX "")

target_include_directories(OpenSeesPy PUBLIC ${Python_INCLUDE_DIRS})

target_link_libraries(OpenSeesPy
   OpenSeesLIB 
   OPS_Reliability
   OPS_Recorder
   OPS_Numerics 
   ${HDF5_LIBRARIES} 
   ${CONAN_LIBS} 
   ${Python_LIBRARIES}
)


#
# INSTALL
#

install (TARGETS OpenSees DESTINATION bin)

# Set selected executable to be built by default
set_target_properties(${OPS_FINAL_TARGET} PROPERTIES EXCLUDE_FROM_ALL OFF)

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
add_compile_definitions(_RELIABILITY)


#----------------------------
# HDF5
#----------------------------

if(HDF5_FOUND)
   include_directories(${HDF5_INCLUDE_DIR})
   set(_hdf5_libs hdf5 hdf5_cpp)
   if (HDF5_VERSION VERSION_GREATER_EQUAL 1.12.0)
         add_compile_definitions(_H5DRM)
         add_compile_definitions(_HDF5)
         message(STATUS "OPS >>> Have HDF5 and VERSION >= 1.12.0")    
   endif()
else()
   message(STATUS "OPS >>> Could not find HDF5")
endif()


if (OPS_Use_Dev_Directories)
  add_subdirectory("${PROJECT_SOURCE_DIR}/DEVELOPER/")
endif()