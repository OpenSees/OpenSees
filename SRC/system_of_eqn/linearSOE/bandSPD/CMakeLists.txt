#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
    PRIVATE
    BandSPDLinSOE.cpp
    BandSPDLinSolver.cpp
    BandSPDLinLapackSolver.cpp
    DistributedBandSPDLinSOE.cpp
    PUBLIC
    BandSPDLinSOE.h
    BandSPDLinSolver.h
    BandSPDLinLapackSolver.h
    DistributedBandSPDLinSOE.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})
