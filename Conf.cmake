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
set(OPS_FINAL_TARGET "OpenSees" CACHE STRING "OpenSees final target")
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

