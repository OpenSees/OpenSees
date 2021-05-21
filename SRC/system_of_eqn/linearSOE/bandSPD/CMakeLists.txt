add_library(systemOfEqn "")

target_sources(systemOfEqn
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

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))
