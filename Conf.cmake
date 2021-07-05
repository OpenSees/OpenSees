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

option(OPS_Use_Reliability   
    "Include reliability"                                    OFF)

option(OPS_Use_Thermal
    "Include thermal components"                              ON)

option(OPS_Use_Graphics
    "Include graphics"                                       OFF)

option(OPS_Use_PFEM 
    "Include PFEM library"                                   OFF)

option(OPS_Use_ASDEA
    "Include ASDEA library"                                   ON)

option(OPS_Use_DRM
    "DRM lib"                                                 ON)

option(OPS_Use_HDF5
    "HDF5 Dependent Code"                                    OFF)


# TODO: Implement material options like elements
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
opensees_load(TCL                                          FIND
	#LIBRARY /home/claudio/miniconda3/lib/libtcl8.6.so
	#INCLUDE /home/claudio/miniconda3/include 
)

opensees_load(BLAS                                         SEARCH
	#LIBRARY /home/claudio/lib/libBlas.a
)

opensees_load(LAPACK                                       SEARCH
#LIBRARY /home/claudio/lib/libLapack.a
)

opensees_load(SUPERLU                                      SEARCH
#LIBRARY /home/claudio/lib/libSuperLU.a
#INCLUDE ${OPS_BUNDLED_DIR}/SuperLU_5.1.1/
)

opensees_load(ARPACK                                       SEARCH)

opensees_load(METIS                                        SEARCH)

opensees_load(HDF5                                           FIND)

opensees_load(MySQL                                          FIND)


#==============================================================================
#                           Select Element Libraries
#
# Each element in this list ows and associated macro definition
#==============================================================================
set(OPS_Element_List

    #OPS_Element_PFEMElement
    #OPS_Element_beamWithHinges
    #OPS_Element_feap
    OPS_Element_LHMYS
    OPS_Element_PML
    OPS_Element_RockingBC
    OPS_Element_UP_ucsd
    OPS_Element_absorbentBoundaries
    OPS_Element_adapter
    OPS_Element_beam3d
    #OPS_Element_beam2d
    OPS_Element_catenaryCable
    OPS_Element_componentElement
    OPS_Element_dispBeamColumnInt
    OPS_Element_forceBeamColumn
    OPS_Element_elastomericBearing
    OPS_Element_frictionBearing
    OPS_Element_generic
    OPS_Element_gradientInelasticBeamColumn
    OPS_Element_joint
    OPS_Element_mixedBeamColumn
    OPS_Element_mvlem
    OPS_Element_pyMacro
    OPS_Element_shell
    OPS_Element_surfaceLoad
    OPS_Element_truss
    OPS_Element_updatedLagrangianBeamColumn
)


