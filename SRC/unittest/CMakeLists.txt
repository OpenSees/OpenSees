#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Unittest OBJECT)

target_sources(OPS_Unittest
    PRIVATE
        unittest.cpp
    PUBLIC
        unittest.h
)
target_include_directories(OPS_Unittest PUBLIC ${CMAKE_CURRENT_LIST_DIR})

