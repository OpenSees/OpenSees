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
project(
                              OpenSees
)
include(OpenSeesFunctions)
#==============================================================================
#                           Select Executable
#
#==============================================================================
set(OPS_FINAL_TARGET "OpenSeesTcl" 
    CACHE STRING "OpenSees final target"
)
#==============================================================================
#                            Basic Switches
#
#==============================================================================
option(FMK
    "Special FMK Code"                                       OFF)

option(OPS_THREADSAFE
    "Only build thread safe components"                       ON)

# Component Libraries
#--------------------------------------

option(OPS_Use_RELIABILITY   
    "Include reliability"                                    OFF)

option(OPS_Use_PFEM 
    "Include PFEM library"                                   OFF)

option(OPS_Use_DRM
    "DRM lib"                                                OFF)

option(OPS_Use_HDF5
    "HDF5 Dependent Code"                                    OFF)

option(OPS_Use_Thermal
    "Include thermal components"                             OFF)

option(OPS_Use_RENDERER
    "Include renderer"                                        ON)

option(OPS_MATERIAL_UNIAXIAL_PY 
    "Include PY material library"                            OFF)

option(OPS_MATERIAL_UNIAXIAL_SNAP 
    "Include snap material library"                          OFF)

#==============================================================================
#                            Properties
#
#==============================================================================

define_property(TARGET
    PROPERTY   OPS_INTERPRETER_GLOBAL #TODO
    BRIEF_DOCS "Include functionality for using global interpreter"
    FULL_DOCS  "..."
)

#==============================================================================
#                      External Libraries
#
#==============================================================================
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

opensees_load(BLAS                                         #SEARCH)
    LIBRARY /home/claudio/lib/libBlas.a
)
#opensees_load(CBLAS                                        SEARCH)

opensees_load(LAPACK                                       #SEARCH)
    LIBRARY /home/claudio/lib/libLapack.a
)
opensees_load(ARPACK                                       SEARCH)

opensees_load(METIS                                        SEARCH)

opensees_load(SUPERLU                                      #SEARCH)
    LIBRARY /home/claudio/lib/libSuperLU.a
    INCLUDE ${OPS_BUNDLED_DIR}/SuperLU_5.1.1/
)

opensees_load(HDF5                                         SEARCH)


#==============================================================================
#                           Select Element Libraries
#
# Each element in this list ows and associated macro definition
#==============================================================================
set(OPS_Element_List

    OPS_Element_absorbentBoundaries
    OPS_Element_adapter
    OPS_Element_beam3d
    OPS_Element_beamWithHinges
    OPS_Element_catenaryCable
    OPS_Element_componentElement
    OPS_Element_dispBeamColumnInt
    OPS_Element_elastomericBearing
    OPS_Element_feap
    OPS_Element_frictionBearing
    OPS_Element_generic
    OPS_Element_gradientInelasticBeamColumn
    OPS_Element_joint
    OPS_Element_LHMYS
    OPS_Element_mixedBeamColumn
    OPS_Element_mvlem
    OPS_Element_PFEMElement
    OPS_Element_PML
    OPS_Element_pyMacro
    OPS_Element_RockingBC
    OPS_Element_shell
    OPS_Element_surfaceLoad
    OPS_Element_updatedLagrangianBeamColumn
    OPS_Element_UP_ucsd
)

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


# Temporary fix
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/incl)



