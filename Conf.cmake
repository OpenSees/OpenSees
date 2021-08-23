#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
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

# Optional Extensions
#--------------------------------------
option(OPS_Use_Reliability   
    "Include reliability"                                    OFF)


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

option(OPS_Use_Thermal
    "Include thermal components"                              ON)

#==============================================================================
#                           Select Auxiliary Libraries
#
# Each element in this list owns an associated macro definition that is the
# name of the lib converted to uppercase, and prepended with an underscore
# (e.g. using OPS_Element_truss defines the macro _OPS_ELEMENT_TRUSS)
#==============================================================================
set(OPS_Extension_List

    OPS_Reliability       # TODO: replace existing tests on '_RELIABILITY'

    OPS_NumLib_PETSC
    OPS_NumLib_METIS
    OPS_NumLib_UMFPACK

    #OPS_ExtLib_PFEM

    OPS_Graphics
    OPS_Renderer 
    OPS_Renderer_GLX      # TODO: replace existing tests on '_GLX'
    OPS_Renderer_X11
)

set(OPS_Element_List

    OPS_Element_truss
    #OPS_Element_beam2d
    OPS_Element_beam3d
    OPS_Element_dispBeamColumnInt
    OPS_Element_forceBeamColumn
    OPS_Element_mixedBeamColumn

    #OPS_Element_beamWithHinges
    OPS_Element_LHMYS
    #OPS_Element_Dmglib
    #OPS_Element_PML
    OPS_Element_RockingBC
    OPS_Element_UP_ucsd
    OPS_Element_absorbentBoundaries
    OPS_Element_adapter
    OPS_Element_catenaryCable
    OPS_Element_componentElement

    OPS_Element_elastomericBearing
    OPS_Element_frictionBearing

    OPS_Element_generic
    OPS_Element_gradientInelasticBeamColumn
    OPS_Element_joint
    OPS_Element_mvlem
    OPS_Element_pyMacro
    OPS_Element_shell
    OPS_Element_surfaceLoad
    OPS_Element_updatedLagrangianBeamColumn
    OPS_Element_masonry
    #OPS_Element_feap
    #OPS_Element_PFEMElement
)


