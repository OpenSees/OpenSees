#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

#opensees_elements(OPS_Element_CEq
#  SOURCES

target_sources(OPS_Element
   PRIVATE
    ASDEmbeddedNodeElement.cpp
   PUBLIC    
    ASDEmbeddedNodeElement.h
)

