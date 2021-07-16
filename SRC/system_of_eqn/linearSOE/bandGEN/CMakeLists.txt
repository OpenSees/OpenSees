
target_sources(G3_SysOfEqn
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


target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

