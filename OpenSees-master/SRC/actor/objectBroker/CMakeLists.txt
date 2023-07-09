#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_ObjectBroker
    PRIVATE
      FEM_ObjectBroker.cpp
    PUBLIC
      FEM_ObjectBroker.h
)

target_include_directories(OPS_ObjectBroker PUBLIC ${CMAKE_CURRENT_LIST_DIR})
