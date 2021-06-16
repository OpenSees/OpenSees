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
# libG3
#
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.13)
project(LG3)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/etc/cmake)
set(OPS_EXTERN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/)
include(OpenSeesFunctions)


option(FMK
    "Special FMK Code"                                       OFF)

option(OPS_THREADSAFE
    "Only build thread safe components"                       ON)

# Interpreter options
#----------------------------------------------------------
option(OPS_INTERPRETER_GLOBAL #TODO
    "Include functionality for using global interpreter"     OFF)

option(OPS_INTERPRETER_TCL
    "Use TCL interpreter"                                    OFF)

option(OPS_INTERPRETER_PYTHON
    "Use Python interpreter"                                 OFF)

# Parallel options
#----------------------------------------------------------
option(OPS_OPTION_PARALLEL   
    "Build OpenSeesSP"                                       OFF)

option(OPS_PARALLEL_INTERPRETER 
    "Build OpenSeesMP"                                       OFF)

# Misc. options
#----------------------------------------------------------
option(OPS_OPTION_RENDERER
    "Include renderer"                                       OFF)

option(OPS_OPTION_THERMAL
    "Include thermal components"                             OFF)

option(OPS_OPTION_HDF5
    "HDF5 Dependent Code"                                    OFF)

option(OPS_OPTION_SEQUENTIAL 
    "Build OpenSees"                                          ON)


option(OPS_OPTION_RELIABILITY   
    "Include reliability"                                    OFF)

option(OPS_OPTION_PFEM 
    "Include PFEM library"                                   OFF)

# Component Libraries
#--------------------------------------
option(OPS_MATERIAL_UNIAXIAL_PY 
    "Include PY material library"                            OFF)

option(OPS_MATERIAL_UNIAXIAL_SNAP 
    "Include snap material library"                          OFF)

# External Libraries
#--------------------------------------
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
opensees_load(TCL      FIND)

opensees_load(BLAS   SEARCH)

opensees_load(LAPACK SEARCH)

opensees_load(ARPACK SEARCH)

opensees_load(METIS  SEARCH)


#----------------------------------------------------------
# Compilers
#
#----------------------------------------------------------

# Fortran
#--------------------------------------
enable_language(Fortran)

# C++
#--------------------------------------


if(FMK)
   add_definitions(-D_HAVE_Damage2p    
           -D_HAVE_PSUMAT
        -D_HAVE_PML
        -D_FILIP_LHNMYS)    
endif()

if(OPS_OPTION_HDF5)
   find_package(HDF5)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
        set(_hdf5_libs hdf5 hdf5_cpp)
    add_definitions(-D_H5DRM)
    else()
     message(STATUS ">>> Could not find HDF5")
    endif()
endif()

if(OPS_OPTION_SEQUENTIAL)

elseif(OPS_OPTION_PARALLEL)
    add_definitions(-D_OPS_PARALLEL_PROCESSING)
else(OPS_PARALLEL_INTERPRETER)
    add_definitions(-D_OPS_PARALLEL_INTERPRETERS)
endif()



if(APPLE)
 message(STATUS ">>> MacOS")
endif()

if(UNIX AND NOT APPLE)
   message(STATUS ">>> LINUX")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -D_LINUX -D_UNIX -ffloat-store")

endif()

if(WIN32)
 message(STATUS ">>> WIN32")
endif()

#
# include paths to main abstract classes
#

# NOTE BeamIntegration and MatrixUtil need to be removed from element/forceBEamColumn

include_directories(src
    src/base
    src/matrix
    src/handler
    src/database
    #src/element
    #src/element/forceBeamColumn
    #src/element/nonlinearBeamColumn/matrixutil
    src/coordTransformation
    src/tagged
    src/tagged/storage
    src/recorder
    src/renderer
    src/damage
    src/recorder/response
    src/material
    src/material/section
    src/material/uniaxial
    #src/material/nD
    src/graph/graph
    src/graph/numberer
    src/graph/partitioner
    src/domain/component
    src/domain/domain
    src/domain/subdomain
    src/domain/load
    src/domain/loadBalancer
    src/domain/pattern
    src/domain/groundMotion
    src/domain/node
    src/domain/constraints
    src/domain/region
    src/analysis/algorithm
    src/analysis/dof_grp
    src/analysis/fe_ele
    src/analysis/algorithm/equiSolnAlgo
    src/analysis/algorithm/eigenAlgo
    src/analysis/algorithm/domainDecompAlgo
    src/analysis/analysis
    src/analysis/integrator
    src/analysis/handler
    src/analysis/numberer
    src/analysis/model
    src/convergenceTest
    src/modelbuilder
    src/system_of_eqn
    src/system_of_eqn/linearSOE
    src/system_of_eqn/eigenSOE
    src/actor/actor
    src/actor/channel
    src/actor/objectBroker
    src/actor/message)

include_directories(include)

include_directories(ext/CSPARSE)
include_directories(ext/AMD)
include_directories(ext/UMFPACK)

#
# build
#
add_subdirectory(src)


