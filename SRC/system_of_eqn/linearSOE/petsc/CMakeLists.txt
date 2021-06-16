
target_sources(systemOfEqn
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


