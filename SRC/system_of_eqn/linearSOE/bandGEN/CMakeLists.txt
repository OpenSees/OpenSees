
target_sources(libg3
    PRIVATE
	BandGenLinSOE.cpp
	BandGenLinSolver.cpp
	BandGenLinLapackSolver.cpp
	DistributedBandGenLinSOE.cpp
    PUBLIC
	BandGenLinSOE.h
	BandGenLinSolver.h
	BandGenLinLapackSolver.h
	DistributedBandGenLinSOE.h
)


target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))

