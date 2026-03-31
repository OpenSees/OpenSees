#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
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

target_sources(OPS_PFEM
    PRIVATE
        PFEMCompressibleLinSOE.cpp
        PFEMCompressibleSolver.cpp
        PFEMDiaLinSOE.cpp
        PFEMDiaSolver.cpp
        PFEMLinSOE.cpp
        PFEMQuasiLinSOE.cpp
        PFEMQuasiSolver.cpp
        PFEMSolver.cpp
        PFEMSolver_Laplace.cpp
        PFEMSolver_LumpM.cpp
        PFEMSolver_LumpM.cpp
        PFEMSolver_Umfpack.cpp
        PFEMSolver_Mumps.cpp	
    PUBLIC
        PFEMCompressibleLinSOE.h
        PFEMCompressibleSolver.h
        PFEMDiaLinSOE.h
        PFEMDiaSolver.h
        PFEMLinSOE.h
        PFEMQuasiLinSOE.h
        PFEMQuasiSolver.h
        PFEMSolver.h
        PFEMSolver_Laplace.h
        PFEMSolver_LumpM.h
        PFEMSolver_LumpM.h
        PFEMSolver_Umfpack.h
        PFEMSolver_Mumps.h	
)

if(MPI_FOUND)

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

target_sources(OpenSeesSP
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

endif()


target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

