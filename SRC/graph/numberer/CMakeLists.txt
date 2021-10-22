#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(graph
    PRIVATE
      RCM.cpp
      AMDNumberer.cpp
      SimpleNumberer.cpp
      GraphNumberer.cpp
      MyRCM.cpp
    PUBLIC
      RCM.h
      AMDNumberer.h
      SimpleNumberer.h
      GraphNumberer.h
      MyRCM.h
)


target_include_directories(graph PUBLIC ${CMAKE_CURRENT_LIST_DIR})
