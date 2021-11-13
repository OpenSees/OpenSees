#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
  PRIVATE
    DomainComponent.cpp
    Parameter.cpp
    ElementParameter.cpp
    NodeResponseParameter.cpp
    LoadFactorParameter.cpp
    MaterialStageParameter.cpp
    MatParameter.cpp
    ElementStateParameter.cpp
  PUBLIC
    DomainComponent.h
    Parameter.h
    ElementParameter.h
    NodeResponseParameter.h
    LoadFactorParameter.h
    MaterialStageParameter.h
    MatParameter.h
    ElementStateParameter.h
)

target_sources(OPS_Domain
  PRIVATE
    InitialStateParameter.cpp
  PUBLIC
    InitialStateParameter.h
)
        
#target_sources(OPS_Reliability
target_sources(OPS_Domain
  PRIVATE
    RVParameter.cpp
    DVParameter.cpp
  PUBLIC
    RVParameter.h
    DVParameter.h
)


target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

