#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Machine OBJECT)

target_sources(OPS_Machine
    PRIVATE
        etime.c
)
target_include_directories(OPS_Machine PUBLIC ${CMAKE_CURRENT_LIST_DIR})

