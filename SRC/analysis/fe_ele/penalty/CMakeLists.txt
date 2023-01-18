#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      PenaltySP_FE.cpp
      PenaltyMP_FE.cpp
    PUBLIC
      PenaltySP_FE.h
      PenaltyMP_FE.h
)

#target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
