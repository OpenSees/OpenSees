target_sources(OPS_Reliability
    PRIVATE
        CorrelationCoefficient.cpp
        CorrelationCoefficientIter.cpp
        Cutset.cpp
        CutsetIter.cpp
        ParameterPositioner.cpp
        ParameterPositionerIter.cpp
        PrincipalPlane.cpp
        RandomVariable.cpp
        RandomVariableIter.cpp
        RandomVariablePositioner.cpp
        RandomVariablePositionerIter.cpp
        ReliabilityDomain.cpp
        ReliabilityDomainComponent.cpp
    PUBLIC
        CorrelationCoefficient.h
        CorrelationCoefficientIter.h
        Cutset.h
        CutsetIter.h
        ParameterPositioner.h
        ParameterPositionerIter.h
        PrincipalPlane.h
        RandomVariable.h
        RandomVariableIter.h
        RandomVariablePositioner.h
        RandomVariablePositionerIter.h
        ReliabilityDomain.h
        ReliabilityDomainComponent.h
)

target_include_directories(OPS_Reliability PUBLIC $(CMAKE_CURRENT_LIST_DIR))

