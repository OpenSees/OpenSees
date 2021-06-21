
target_sources(libg3
    PRIVATE 
    	SparseGenColLinSOE.cpp
	SparseGenColLinSolver.cpp
	SparseGenRowLinSOE.cpp
	SparseGenRowLinSolver.cpp
	SuperLU.cpp
    PUBLIC
    	SparseGenColLinSOE.h
	SparseGenColLinSolver.h
	SparseGenRowLinSOE.h
	SparseGenRowLinSolver.h
	SuperLU.h
)

target_sources(libPFEM
    PRIVATE
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

target_sources(OpenSeesMP
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

