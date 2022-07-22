#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Actor
    PRIVATE
      Actor.cpp
      MovableObject.cpp
    PUBLIC
      Actor.h
      MovableObject.h
)

target_include_directories(OPS_Actor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
