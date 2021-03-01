add_library(systemOfEqn "")

target_sources(systemOfEqn
    PRIVATE
	ProfileSPDLinSOE.cpp
	ProfileSPDLinSolver.cpp
	ProfileSPDLinDirectSolver.cpp
	ProfileSPDLinSubstrSolver.cpp
	ProfileSPDLinDirectBlockSolver.cpp
	ProfileSPDLinDirectSkypackSolver.cpp
	ProfileSPDLinSolverGather.cpp
	ProfileSPDLinSOEGather.cpp
	DistributedProfileSPDLinSOE.cpp
	SProfileSPDLinSOE.cpp
	SProfileSPDLinSolver.cpp
    PUBLIC
	ProfileSPDLinSOE.h
	ProfileSPDLinSolver.h
	ProfileSPDLinDirectSolver.h
	ProfileSPDLinSubstrSolver.h
	ProfileSPDLinDirectBlockSolver.h
	ProfileSPDLinDirectSkypackSolver.h
	ProfileSPDLinSolverGather.h
	ProfileSPDLinSOEGather.h
	DistributedProfileSPDLinSOE.h
	SProfileSPDLinSOE.h
	SProfileSPDLinSolver.h
)

if(PARALLEL_INTERPRETERS)

elseif(PARALLEL)

endif

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))
