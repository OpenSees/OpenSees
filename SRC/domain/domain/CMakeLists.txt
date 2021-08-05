
target_sources(OPS_Domain
    PRIVATE
    Domain.cpp
    PUBLIC
    Domain.h
    ElementIter.h
    LoadCaseIter.h
    MP_ConstraintIter.h
    NodeIter.h
    ParameterIter.h
    Pressure_ConstraintIter.h
    SP_ConstraintIter.h
    SubdomainIter.h
)

target_sources(OPS_ASDEA 
    PRIVATE DomainModalProperties.cpp
    PUBLIC DomainModalProperties.h
)


target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(single)
add_subdirectory(partitioned)
