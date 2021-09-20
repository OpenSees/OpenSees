#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
    PRIVATE
    DiagonalSOE.cpp
    DiagonalSolver.cpp
    DiagonalDirectSolver.cpp
    DistributedDiagonalSOE.cpp
    DistributedDiagonalSolver.cpp
    PUBLIC 
        DiagonalSOE.h
    DiagonalSolver.h
    DiagonalDirectSolver.h
    DistributedDiagonalSOE.h
    DistributedDiagonalSolver.h
)



target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

