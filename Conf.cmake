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
set(OPS_FINAL_TARGET "OpenSeesTcl" CACHE STRING "OpenSees final target")
#==============================================================================
#                            Basic Switches
#
#==============================================================================

# Optional Extensions
#--------------------------------------
option(OPS_Use_Dev_Directories
  "Include files in DEVELOPER directory"                   OFF)

option(OPS_Use_DRM
  "DRM lib"                                                ON )

option(OPS_Use_HDF5
  "HDF5 Dependent Code"                                    OFF)

option(OPS_Use_PFEM
  "                   "                                    ON )

option(FMK
  "Special FMK Code"                                       OFF)

set(OPS_Use_Graphics_Option
  None
  # Base
  # OpenGL
)

#==============================================================================
#                           Select Auxiliary Components
#
# Each element in this list owns an associated macro definition that is the
# name of the lib converted to uppercase and prepended with "OPSDEF_"
# (e.g. using OPS_Element_truss defines the macro OPSDEF_ELEMENT_TRUSS)
#==============================================================================
set(OPS_Numlib_List
)

set(OPS_SysOfEqn_List
  OPS_SysOfEqn_UMF
  #OPS_SysOfEqn_ITPACK
)

set(OPS_Extension_List
  OPS_ASDEA
  OPS_Paraview
  #OPS_Reliability       # TODO: replace existing tests on '_RELIABILITY'


  #OPS_Graphics
  #OPS_Renderer 
  #OPS_Renderer_GLX      # TODO: replace existing tests on '_GLX'
  #OPS_Renderer_X11
)

set(OPS_Element_List
  OPS_Element_truss
  OPS_Element_InertiaTruss
  #OPS_Element_beam2d
  OPS_Element_beam3d
  OPS_Element_dispBeamColumnInt
  OPS_Element_forceBeamColumn
  OPS_Element_mixedBeamColumn
  #OPS_Element_beamWithHinges
  OPS_Element_LHMYS
  #OPS_Element_Dmglib
  OPS_Element_PML
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
  OPS_Element_PFEMElement
  OPS_Element_CEq
  OPS_Material_StressDensity
)

set(OPS_Exclude_List
  OPS_Element_feap
  OPS_Material_StressDensity
  OPS_Recorder_PVD
  OPS_Uniaxial_Fedeas
)

