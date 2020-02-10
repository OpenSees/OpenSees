!
!
!  Include file for Fortran use of the KSP package in PETSc
!
#if !defined (__PETSCKSPDEF_H)
#define __PETSCKSPDEF_H

#include "petsc/finclude/petscpc.h"

#define KSP type(tKSP)
#define KSPGuess type(tKSPGuess)

#define KSPType character*(80)
#define KSPGuessType character*(80)
#define KSPCGType PetscEnum
#define KSPFCDTruncationType PetscEnum
#define KSPConvergedReason PetscEnum
#define KSPNormType PetscEnum
#define KSPGMRESCGSRefinementType PetscEnum
#define MatSchurComplementAinvType PetscEnum
!
!  Various Krylov subspace methods
!
#define KSPRICHARDSON 'richardson'
#define KSPCHEBYSHEV 'chebyshev'
#define KSPCG 'cg'
#define KSPCGNE 'cgne'
#define KSPSTCG 'stcg'
#define KSPGLTR 'gltr'
#define KSPFCG 'fcg'
#define KSPGMRES 'gmres'
#define KSPFGMRES 'fgmres'
#define KSPLGMRES 'lgmres'
#define KSPDGMRES 'dgmres'
#define KSPPGMRES 'pgmres'
#define KSPTCQMR 'tcqmr'
#define KSPBCGS 'bcgs'
#define KSPIBCGS 'ibcgs'
#define KSPFBCGS  'fbcgs'
#define KSPFBCGSR 'fbcgsr'
#define KSPBCGSL 'bcgsl'
#define KSPCGS 'cgs'
#define KSPTFQMR 'tfqmr'
#define KSPCR 'cr'
#define KSPLSQR 'lsqr'
#define KSPPREONLY 'preonly'
#define KSPQCG 'qcg'
#define KSPBICG 'bicg'
#define KSPMINRES 'minres'
#define KSPSYMMLQ 'symmlq'
#define KSPLCD 'lcd'
#define KSPPYTHON 'python'
#define KSPGCR 'gcr'
#define KSPTSIRM 'tsirm'
#define KSPCGLS 'cgls'
#define KSPFETIDP 'fetidp'
!
!  Various Initial guesses for Krylov subspace methods
!
#define KSPGUESSFISCHER 'fischer'
#define KSPGUESSPOD 'pod'
#endif
