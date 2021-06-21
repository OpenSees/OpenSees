
target_sources(libg3
    PRIVATE
	DomainSolver.cpp
	LinearSOE.cpp
	LinearSOESolver.cpp
    PUBLIC
	DomainSolver.h
	LinearSOE.h
	LinearSOESolver.h
)

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))

add_subdirectory(bandGEN)
add_subdirectory(bandSPD)
add_subdirectory(diagonal)
add_subdirectory(fullGEN)
add_subdirectory(sparseGEN)
add_subdirectory(sparseSYM)
add_subdirectory(umfGEN)

#add_subdirectory(profileSPD)
#add_subdirectory(cg)
#add_subdirectory(petsc)
#add_subdirectory(mumps)
#add_subdirectory(itpack)
