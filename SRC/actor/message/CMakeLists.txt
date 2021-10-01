#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Actor
    PRIVATE
      Message.cpp      
    PUBLIC
      Message.h
)

target_include_directories(OPS_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
