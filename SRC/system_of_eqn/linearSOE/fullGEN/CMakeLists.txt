#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
    PRIVATE
        FullGenLinSOE.cpp
    FullGenLinSolver.cpp
        FullGenLinLapackSolver.cpp

    PUBLIC
        FullGenLinSOE.h
    FullGenLinSolver.h
        FullGenLinLapackSolver.h

)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})
