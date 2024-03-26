#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Optimization
    PRIVATE
        TclOptimizationBuilder.cpp
    PUBLIC
        TclOptimizationBuilder.h
)

target_include_directories(OPS_Optimization PUBLIC ${CMAKE_CURRENT_LIST_DIR})