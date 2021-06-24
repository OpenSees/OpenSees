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
#                           Select Executable
#
#------------------------------------------------------------------------------
set(OPS_FINAL_TARGET "OpenSeesTcl" 
    CACHE STRING "OpenSees final target"
)
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
opensees_load(TCL     
	LIBRARY /home/claudio/miniconda3/lib/libtcl8.6.so
	INCLUDE /home/claudio/miniconda3/include 
)
opensees_load(BLAS                                         SEARCH)

opensees_load(LAPACK                                       SEARCH)

opensees_load(ARPACK                                       SEARCH)

opensees_load(METIS                                        SEARCH)

#opensees_load(SUPERLU                                      SEARCH)

#opensees_load(HDF5                                         SEARCH)



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



include_directories(${PROJECT_SOURCE_DIR}/include)


#
# build
#
#add_subdirectory(${OPS_SRC_DIR})


