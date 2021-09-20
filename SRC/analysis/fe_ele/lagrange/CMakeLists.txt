#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      LagrangeSP_FE.cpp
      LagrangeMP_FE.cpp
    PUBLIC
      LagrangeSP_FE.h
      LagrangeMP_FE.h
)

target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
