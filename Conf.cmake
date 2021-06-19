#cmake_minimum_required(VERSION 3.18)
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
#project(
#                              OpenSees
#)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/etc/cmake)
#set(OPS_EXTERN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/)
#set(OPS_EXTERNALS_DIR ${PROJECT_SOURCE_DIR}/ext/)
#set(OPS_SRC_DIR ${PROJECT_SOURCE_DIR}/src/)
include(OpenSeesFunctions)
#------------------------------------------------------------------------------
#                            Basic Switches
#
#------------------------------------------------------------------------------
option(FMK
    "Special FMK Code"                                       OFF)

option(OPS_THREADSAFE
    "Only build thread safe components"                       ON)

# Component Libraries
#--------------------------------------
option(OPS_MATERIAL_UNIAXIAL_PY 
    "Include PY material library"                            OFF)

option(OPS_MATERIAL_UNIAXIAL_SNAP 
    "Include snap material library"                          OFF)

option(OPS_OPTION_RELIABILITY   
    "Include reliability"                                    OFF)

option(OPS_OPTION_PFEM 
    "Include PFEM library"                                   OFF)

option(OPS_OPTION_DRM
    "Include DRM library"                                    OFF)

option(OPS_OPTION_BROKER_ALL
    "Include broker for all classes"                         OFF)


option(OPS_OPTION_HDF5
    "HDF5 Dependent Code"                                    OFF)

option(OPS_OPTION_THERMAL
    "Include thermal components"                             OFF)

option(OPS_OPTION_RENDERER
     "Include renderer"                                      OFF)

#------------------------------------------------------------------------------
#                            Properties
#
#------------------------------------------------------------------------------

define_property(TARGET
    PROPERTY   OPS_INTERPRETER_GLOBAL #TODO
    BRIEF_DOCS "Include functionality for using global interpreter"
    FULL_DOCS  "..."
)



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
opensees_load(TCL      FIND)

opensees_load(BLAS   SEARCH)

opensees_load(LAPACK SEARCH)

opensees_load(ARPACK SEARCH)

opensees_load(METIS  SEARCH)


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


if(FMK)
   add_compile_definitions(
	_HAVE_Damage2p    
        _HAVE_PSUMAT
        _HAVE_PML
        _FILIP_LHNMYS)    
endif()

if(OPS_OPTION_HDF5)
   find_package(HDF5)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
        set(_hdf5_libs hdf5 hdf5_cpp)
    add_compile_definitions(-D_H5DRM)
    else()
     message(STATUS ">>> Could not find HDF5")
    endif()
endif()

if(APPLE)
 message(STATUS ">>> MacOS")
endif()

if(UNIX AND NOT APPLE)
   message(STATUS ">>> LINUX")
   add_compile_definitions(_LINUX _UNIX)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -fPIC -ffloat-store")
endif()

if(WIN32)
 message(STATUS ">>> WIN32")
endif()
#
# include paths to main abstract classes
#

# NOTE BeamIntegration and MatrixUtil need to be removed from element/forceBEamColumn

include_directories(${OPS_SRC_DIR}
#    ${OPS_SRC_DIR}/base
    ${OPS_SRC_DIR}/matrix
    ${OPS_SRC_DIR}/handler
    ${OPS_SRC_DIR}/database
    #${OPS_SRC_DIR}/element
    #${OPS_SRC_DIR}/element/forceBeamColumn
    #${OPS_SRC_DIR}/element/nonlinearBeamColumn/matrixutil
    ${OPS_SRC_DIR}/coordTransformation
    ${OPS_SRC_DIR}/tagged
    ${OPS_SRC_DIR}/tagged/storage
    ${OPS_SRC_DIR}/recorder
    ${OPS_SRC_DIR}/renderer
    ${OPS_SRC_DIR}/damage
    ${OPS_SRC_DIR}/recorder/response
    ${OPS_SRC_DIR}/material
    ${OPS_SRC_DIR}/material/section
    ${OPS_SRC_DIR}/material/uniaxial
    #${OPS_SRC_DIR}/material/nD
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
)



#------------------------------------------------------------------------------
#                            Targets
#------------------------------------------------------------------------------
#add_library(libg3)
#add_executable(OpenSees   EXCLUDE_FROM_ALL src/tcl/tclMain.cpp)
#add_executable(OpenSeesSP EXCLUDE_FROM_ALL src/tcl/tclMain.cpp)
#add_executable(OpenSeesMP EXCLUDE_FROM_ALL src/tcl/tclMain.cpp)
#add_executable(openseespy EXCLUDE_FROM_ALL src/interpreter/pythonMain.cpp)
#
#
## Set properties for targets
#target_compile_definitions(OpenSeesSP 
#    PUBLIC _OPS_PARALLEL_INTERPRETERS
#)
#target_compile_definitions(OpenSeesMP
#    PUBLIC _OPS_PARALLEL_PROCESSING
#)
#
include_directories(${PROJECT_SOURCE_DIR}/include)

#include_directories(ext/CSPARSE)
#include_directories(ext/AMD)
#include_directories(ext/UMFPACK)

#
# build
#
#add_subdirectory(${OPS_SRC_DIR})


