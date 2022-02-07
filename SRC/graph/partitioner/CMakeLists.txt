#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(graph
    PRIVATE
      MetisWrapper.cpp
    PUBLIC
      MetisWrapper.h
      GraphPartitioner.h
)
