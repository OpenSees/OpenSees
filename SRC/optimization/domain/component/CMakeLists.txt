#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Optimization
    PRIVATE
        ConstraintFuntion.cpp
        DesignVariable.cpp
        DesignVariableIter.cpp
        DesignVariablePositioner.cpp
        DesignVariablePositionerIter.cpp
        ObjectiveFunction.cpp
    PUBLIC
        ConstraintFuntion.h
        DesignVariable.h
        DesignVariableIter.h
        DesignVariablePositioner.h
        DesignVariablePositionerIter.h
        ObjectiveFunction.h
)

target_include_directories(OPS_Optimization PUBLIC ${CMAKE_CURRENT_LIST_DIR})