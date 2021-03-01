add_library(systemOfEqn "")

target_sources(systemOfEqn
    PRIVATE
	DiagonalSOE.cpp
	DiagonalSolver.cpp
	DiagonalDirectSolver.cpp
	DistributedDiagonalSOE.cpp
	DistributedDiagonalSolver.cpp
    PUBLIC 
    	DiagonalSOE.h
	DiagonalSolver.h
	DiagonalDirectSolver.h
	DistributedDiagonalSOE.h
	DistributedDiagonalSolver.h
)



target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))
