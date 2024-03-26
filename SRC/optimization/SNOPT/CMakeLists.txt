#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Optimization
    PRIVATE
        SNOPTAnalysis.cpp
        SNOPTClass.cpp
        toyOptFunction.cpp
    PUBLIC
        SNOPTAnalysis.h
        SNOPTClass.h
        toyOptFunction.h
)

target_include_directories(OPS_Optimization PUBLIC ${CMAKE_CURRENT_LIST_DIR})
