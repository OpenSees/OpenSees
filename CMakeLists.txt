cmake_minimum_required(VERSION 3.18)
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
project(
                              OpenSees
)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/ETC/cmake)
set(OPS_EXTERN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/OTHER/)
set(OPS_EXTERNALS_DIR ${PROJECT_SOURCE_DIR}/OTHER/)
set(OPS_SRC_DIR ${PROJECT_SOURCE_DIR}/SRC/)
#include(OpenSeesFunctions)
include(${PROJECT_SOURCE_DIR}/Conf.cmake)
#------------------------------------------------------------------------------
#                            Basic Switches
#
#------------------------------------------------------------------------------
#option(FMK
#    "Special FMK Code"                                       OFF)
#
#option(OPS_THREADSAFE
#    "Only build thread safe components"                       ON)
#
## Component Libraries
##--------------------------------------
#option(OPS_MATERIAL_UNIAXIAL_PY 
#    "Include PY material library"                            OFF)
#
#option(OPS_MATERIAL_UNIAXIAL_SNAP 
#    "Include snap material library"                          OFF)
#
#option(OPS_OPTION_RELIABILITY   
#    "Include reliability"                                    OFF)
#
#option(OPS_OPTION_PFEM 
#    "Include PFEM library"                                   OFF)
#
#option(OPS_OPTION_HDF5
#    "HDF5 Dependent Code"                                    OFF)
#
#option(OPS_OPTION_THERMAL
#    "Include thermal components"                             OFF)
#
#option(OPS_OPTION_RENDERER
#     "Include renderer"                                      OFF)
#
##------------------------------------------------------------------------------
##                            Properties
##
##------------------------------------------------------------------------------
#
#define_property(TARGET
#    PROPERTY   OPS_INTERPRETER_GLOBAL #TODO
#    BRIEF_DOCS "Include functionality for using global interpreter"
#    FULL_DOCS  "..."
#)



#----------------------------------------------------------------
#                      External Libraries
#----------------------------------------------------------------
# Synopsis
# - opensees_load(<PACKAGE> [BUILD|FIND|SEARCH|PATHS] [<PATHS>])
#
# Options:
# - BUILD:  Build OpenSees provided library
# - FIND:   Use CMake to find library, fail if not found
# - SEARCH: Try finding library with CMake, build OpenSees
#           Version if not found.
# - PATHS:  Provide specific paths for library.
#
#----------------------------------------------------------------
#opensees_load(TCL      FIND)
#
#opensees_load(BLAS   SEARCH)
#
#opensees_load(LAPACK SEARCH)
#
#opensees_load(ARPACK SEARCH)
#
#opensees_load(METIS  SEARCH)


#----------------------------------------------------------------
# Compilers
#
#----------------------------------------------------------------

# Fortran
#--------------------------------------
enable_language(Fortran)

# C++
#--------------------------------------


#----------------------------------------------------------------
# Compilers
#
#----------------------------------------------------------------


#if(FMK)
#   add_compile_definitions(
#	_HAVE_Damage2p    
#        _HAVE_PSUMAT
#        _HAVE_PML
#        _FILIP_LHNMYS)    
#endif()
#
#if(OPS_OPTION_HDF5)
#   find_package(HDF5)
#    if(HDF5_FOUND)
#        include_directories(${HDF5_INCLUDE_DIR})
#        set(_hdf5_libs hdf5 hdf5_cpp)
#    add_compile_definitions(-D_H5DRM)
#    else()
#     message(STATUS ">>> Could not find HDF5")
#    endif()
#endif()
#
#if(APPLE)
# message(STATUS ">>> MacOS")
#endif()
#
#if(UNIX AND NOT APPLE)
#   message(STATUS ">>> LINUX")
#   add_compile_definitions(_LINUX _UNIX)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -ffloat-store")
#endif()
#
#if(WIN32)
# message(STATUS ">>> WIN32")
#endif()
#
# include paths to main abstract classes
#

# NOTE BeamIntegration and MatrixUtil need to be removed from element/forceBEamColumn

include_directories(src
#    src/base
#    src/matrix
#    src/handler
#    src/database
#    #src/element
    src/element/forceBeamColumn
    src/element/nonlinearBeamColumn/matrixutil
#    src/coordTransformation
#    src/tagged
#    src/tagged/storage
#    src/recorder
#    src/renderer
#    src/damage
#    src/recorder/response
#    src/material
#    src/material/section
#    src/material/uniaxial
#    #src/material/nD
#    src/graph/graph
#    src/graph/numberer
#    src/graph/partitioner
#    src/domain/component
#    src/domain/domain
#    src/domain/subdomain
#    src/domain/load
#    src/domain/loadBalancer
#    src/domain/pattern
#    src/domain/groundMotion
#    src/domain/node
#    src/domain/constraints
#    src/domain/region
#    src/analysis/algorithm
#    src/analysis/dof_grp
#    src/analysis/fe_ele
#    src/analysis/algorithm/equiSolnAlgo
#    src/analysis/algorithm/eigenAlgo
#    src/analysis/algorithm/domainDecompAlgo
#    src/analysis/analysis
#    src/analysis/integrator
#    src/analysis/handler
#    src/analysis/numberer
#    src/analysis/model
#    src/convergenceTest
#    src/modelbuilder
#    src/system_of_eqn
#    src/system_of_eqn/linearSOE
#    src/system_of_eqn/eigenSOE
#    src/actor/actor
#    src/actor/channel
#    src/actor/objectBroker
#    src/actor/message
)
#


#------------------------------------------------------------------------------
#                            Targets
#------------------------------------------------------------------------------
add_library(libg3)
add_library(libPFEM       EXCLUDE_FROM_ALL)
add_executable(OpenSees   EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclMain.cpp)
add_executable(OpenSeesSP EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclMain.cpp)
add_executable(OpenSeesMP EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/tcl/tclMain.cpp)
add_executable(openseespy EXCLUDE_FROM_ALL ${OPS_SRC_DIR}/interpreter/pythonMain.cpp)


# Set properties for targets
target_compile_definitions(OpenSeesSP 
    PUBLIC _OPS_PARALLEL_INTERPRETERS
)
target_compile_definitions(OpenSeesMP
    PUBLIC _OPS_PARALLEL_PROCESSING 
)

#include_directories(include)


#
# build
#
add_subdirectory(${OPS_SRC_DIR})


