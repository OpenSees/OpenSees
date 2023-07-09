#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# opensees_elements(OPS_Element_beam3d

target_sources(OPS_Element
  PRIVATE
     beam3d01.cpp
     beam3d02.cpp
  PUBLIC     
     beam3d01.h
     beam3d02.h
)

#target_include_directories(OPS_Element_beam3d PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

