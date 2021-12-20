#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#add_library(OPS_Element_masonry)

target_sources(OPS_Element
  PUBLIC
    BeamGT.cpp
    MasonPan12.cpp
    MasonPan3D.cpp
  PRIVATE
    BeamGT.h
    MasonPan12.h
    MasonPan3D.h
)    

