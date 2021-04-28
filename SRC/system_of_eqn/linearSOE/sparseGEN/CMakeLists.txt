add_library(systemOfEqn "")

target_sources(systemOfEqn
    PRIVATE 
    	SparseGenColLinSOE.cpp
	SparseGenColLinSolver.cpp
	SparseGenRowLinSOE.cpp
	SparseGenRowLinSolver.cpp
	SuperLU.cpp
	PFEMSolver.cpp
	PFEMLinSOE.cpp
	PFEMCompressibleSolver.cpp
	PFEMCompressibleLinSOE.cpp
	PFEMSolver_LumpM.cpp
	PFEMSolver_Umfpack.cpp
	PFEMQuasiLinSOE.cpp
	PFEMQuasiSolver.cpp
	PFEMDiaLinSOE.cpp
	PFEMDiaSolver.cpp
	PFEMSolver_Laplace.cpp
	PFEMSolver_LumpM.cpp
    PUBLIC
    	SparseGenColLinSOE.h
	SparseGenColLinSolver.h
	SparseGenRowLinSOE.h
	SparseGenRowLinSolver.h
	SuperLU.h
	PFEMSolver.h
	PFEMLinSOE.h
	PFEMCompressibleSolver.h
	PFEMCompressibleLinSOE.h
	PFEMSolver_LumpM.h
	PFEMSolver_Umfpack.h
	PFEMQuasiLinSOE.h
	PFEMQuasiSolver.h
	PFEMDiaLinSOE.h
	PFEMDiaSolver.h
	PFEMSolver_Laplace.h
	PFEMSolver_LumpM.h
)


ifdef(PARALLEL)
target_sources(systemOfEqn
    PRIVATE 
	DistributedSuperLU.cpp
	DistributedSparseGenColLinSOE.cpp
	DistributedSparseGenRowLinSOE.cpp
	PFEMSolver_Mumps.cpp
	PFEMSolver_LumpM.cpp
	PFEMCompressibleSolver_Mumps.cpp
    PUBLIC
	DistributedSuperLU.h
	DistributedSparseGenColLinSOE.h
	DistributedSparseGenRowLinSOE.h
	PFEMSolver_Mumps.h
	PFEMCompressibleSolver_Mumps.h
)

target_include_directories(systemOfEqn PUBLIC $(CMAKE_CURRENT_LIST_DIR))
