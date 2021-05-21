add_library(systemOfEqn "")

target_sources(systemOfEqn
    PRIVATE
	SymSparseLinSOE.cpp
	SymSparseLinSolver.cpp
	grcm.c
	nest.c
	nmat.c
	symbolic.c
	utility.c
	newordr.c
	nnsim.c
	genmmd.c
    PUBLIC
	SymSparseLinSOE.h
	SymSparseLinSolver.h
)

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))
