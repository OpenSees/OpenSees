#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_SysOfEqn
    PRIVATE
        UmfpackGenLinSOE.cpp
        UmfpackGenLinSolver.cpp

    PUBLIC
        UmfpackGenLinSOE.h
        UmfpackGenLinSolver.h

)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

