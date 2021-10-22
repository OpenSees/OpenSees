#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
    PRIVATE
        ActorPetscSOE.cpp
        badPetscSOE.cpp
        main.cpp
        PetscSOE.cpp
        PetscSolver.cpp
        PetscSparseSeqSolver.cpp
        ShadowPetscSOE.cpp
    #TclPetsc.cpp

    PUBLIC
        ActorPetscSOE.h
        oldpetsc.h
        PetscSOE.h
        PetscSolver.h
        PetscSparseSeqSolver.h
        ShadowPetscSOE.h
)


