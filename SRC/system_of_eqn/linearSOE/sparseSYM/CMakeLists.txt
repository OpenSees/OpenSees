
target_sources(G3_SysOfEqn
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
	genmmd.f
    PUBLIC
	SymSparseLinSOE.h
	SymSparseLinSolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))
