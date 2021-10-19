#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
  PRIVATE
    Domain.cpp
    DomainModalProperties.cpp
  PUBLIC
    Domain.h
    DomainModalProperties.h
    ElementIter.h
    LoadCaseIter.h
    MP_ConstraintIter.h
    NodeIter.h
    ParameterIter.h
    Pressure_ConstraintIter.h
    SP_ConstraintIter.h
    SubdomainIter.h
)

add_subdirectory(single)
add_subdirectory(partitioned)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

