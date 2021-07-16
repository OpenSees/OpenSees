
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
    PUBLIC
	    SymSparseLinSOE.h
	    SymSparseLinSolver.h
)

add_library(G3_SysOfEqn_f STATIC)
target_sources(G3_SysOfEqn_f PRIVATE genmmd.f)
target_link_libraries(G3_SysOfEqn PUBLIC G3_SysOfEqn_f)

#target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

