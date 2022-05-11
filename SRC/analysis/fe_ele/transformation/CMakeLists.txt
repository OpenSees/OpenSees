#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      TransformationFE.cpp
    PUBLIC
      TransformationFE.h
)

#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
