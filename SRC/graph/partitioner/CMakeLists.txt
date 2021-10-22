#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(graph
    PRIVATE
      Metis.cpp
    PUBLIC
      Metis.h
      GraphPartitioner.h
)
