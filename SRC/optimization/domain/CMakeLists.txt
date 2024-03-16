#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Optimization
    PRIVATE
        Optimization.cpp        
        OptimizationDomainComponent.cpp
        OptimizationDomain.cpp
    PUBLIC
        Optimization.h
        OptimizationDomainComponent.h
        OptimizationDomain.h        
)

target_include_directories(OPS_Optimization PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(component)
