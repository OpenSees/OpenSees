#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_SysOfEqn
    PRIVATE
        ArpackSOE.cpp
        ArpackSolver.cpp
        EigenSOE.cpp
        EigenSolver.cpp
        FullGenEigenSOE.cpp
        FullGenEigenSolver.cpp
        SymBandEigenSOE.cpp
        SymBandEigenSolver.cpp
    PUBLIC
        ArpackSOE.h
        ArpackSolver.h
        EigenSOE.h
        EigenSolver.h
        FullGenEigenSOE.h
        FullGenEigenSolver.h
        SymBandEigenSOE.h
        SymBandEigenSolver.h
)

target_sources(OPS_SysOfEqn
    PRIVATE
    BandArpackSOE.cpp
    BandArpackSolver.cpp
    SymArpackSOE.cpp
    SymArpackSolver.cpp
    PUBLIC
    BandArpackSOE.h
    BandArpackSolver.h
    SymArpackSOE.h
    SymArpackSolver.h
)

