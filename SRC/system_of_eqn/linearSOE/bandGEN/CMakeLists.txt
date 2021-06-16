
target_sources(systemOfEqn
    PRIVATE
	BandGenLinSOE.cpp
	BandGenLinSolver.cpp
	BandGenLinLapackSolver.cpp
	DistributedBandGenLinSOE.cpp
	#BandGenLinSOE_Single.cpp
	#LinearSOESolver.cpp
    PUBLIC
	BandGenLinSOE.h
	BandGenLinSolver.h
	BandGenLinLapackSolver.h
	DistributedBandGenLinSOE.h
	#BandGenLinSOE_Single.h
	#LinearSOESolver.h
)

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))
