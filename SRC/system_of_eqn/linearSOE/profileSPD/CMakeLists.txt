#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
    PRIVATE
    ProfileSPDLinSOE.cpp
    ProfileSPDLinSolver.cpp
    ProfileSPDLinDirectSolver.cpp
    ProfileSPDLinSubstrSolver.cpp
    ProfileSPDLinDirectBlockSolver.cpp
    ProfileSPDLinDirectSkypackSolver.cpp
    #ProfileSPDLinSolverGather.cpp
    #ProfileSPDLinSOEGather.cpp
    SProfileSPDLinSOE.cpp
    SProfileSPDLinSolver.cpp
    PUBLIC
    ProfileSPDLinSOE.h
    ProfileSPDLinSolver.h
    ProfileSPDLinDirectSolver.h
    ProfileSPDLinSubstrSolver.h
    ProfileSPDLinDirectBlockSolver.h
    ProfileSPDLinDirectSkypackSolver.h
    #ProfileSPDLinSolverGather.h
    #ProfileSPDLinSOEGather.h
    SProfileSPDLinSOE.h
    SProfileSPDLinSolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

