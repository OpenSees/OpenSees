#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Interp_Tcl
    PRIVATE
        Block2D.cpp
        Block3D.cpp
        myCommands.cpp
    PUBLIC
        Block2D.h
        Block3D.h
)
target_include_directories(OPS_Interp_Tcl PUBLIC ${CMAKE_CURRENT_LIST_DIR})

