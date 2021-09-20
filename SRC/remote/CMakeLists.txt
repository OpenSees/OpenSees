#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Remote OBJECT)

target_sources(OPS_Remote
    PRIVATE
        remote.c
)
target_include_directories(OPS_Remote PUBLIC ${CMAKE_CURRENT_LIST_DIR})

