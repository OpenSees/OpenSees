 #==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
  PRIVATE
    PartitionedDomainEleIter.cpp 
    PartitionedDomainSubIter.cpp
    PartitionedDomain.cpp
  PUBLIC
    PartitionedDomainEleIter.h
    PartitionedDomainSubIter.h
    PartitionedDomain.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
