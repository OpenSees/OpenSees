!
!  This file is part of MUMPS 5.1.2, released
!  on Mon Oct  2 07:37:01 UTC 2017
!
!
!  Copyright 1991-2017 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
!  University of Bordeaux.
!
!  This version of MUMPS is provided to you free of charge. It is
!  released under the CeCILL-C license:
!  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
!
      INCLUDE 'zmumps_root.h'
      TYPE ZMUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters 
! for the interface to the user, plus internal
! information from the solver
!
! *****************
! INPUT PARAMETERS
! *****************
!    -----------------
!    MPI Communicator
!    -----------------
        INTEGER :: COMM
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite, 
!        SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
        INTEGER ::  SYM, PAR
        INTEGER ::  JOB 
!    --------------------
!    Order of Input matrix 
!    --------------------
        INTEGER ::  N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
        INTEGER    :: NZ  ! Standard integer input + bwd. compat.
        INTEGER(8) :: NNZ ! 64-bit integer input
        COMPLEX(kind=8), DIMENSION(:), POINTER :: A
        INTEGER, DIMENSION(:), POINTER :: IRN, JCN
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA, ROWSCA, pad0
!
!       ------------------------------------
!       Case of distributed assembled matrix
!       matrix on entry:
!       ------------------------------------
        INTEGER    :: NZ_loc  ! Standard integer input + bwd. compat.
        INTEGER    :: pad1
        INTEGER(8) :: NNZ_loc ! 64-bit integer input
        INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
        COMPLEX(kind=8), DIMENSION(:), POINTER :: A_loc, pad2
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
        INTEGER :: NELT, pad3
        INTEGER, DIMENSION(:), POINTER :: ELTPTR
        INTEGER, DIMENSION(:), POINTER :: ELTVAR
        COMPLEX(kind=8), DIMENSION(:), POINTER :: A_ELT, pad4
!
!    ---------------------------------------------
!    Symmetric permutation : 
!               PERM_IN if given by user (optional)
!    ---------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data 
! ******************
!    --------------------------------------------------------
!    RHS / SOL_loc
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
        COMPLEX(kind=8), DIMENSION(:), POINTER :: RHS, REDRHS
        COMPLEX(kind=8), DIMENSION(:), POINTER :: RHS_SPARSE
        COMPLEX(kind=8), DIMENSION(:), POINTER :: SOL_loc
        INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE
        INTEGER, DIMENSION(:), POINTER :: IRHS_PTR
        INTEGER, DIMENSION(:), POINTER :: ISOL_loc
        INTEGER ::  LRHS, NRHS, NZ_RHS, LSOL_loc, LREDRHS
        INTEGER ::  pad5
!    ----------------------------
!    Control parameters,
!    statistics and output data
!    ---------------------------
        INTEGER ::  ICNTL(40)
        INTEGER ::  INFO(40) 
        INTEGER :: INFOG(40)
        DOUBLE PRECISION ::  COST_SUBTREES
        DOUBLE PRECISION ::  CNTL(15)
        DOUBLE PRECISION ::  RINFO(40)
        DOUBLE PRECISION ::  RINFOG(40)
!    ---------------------------------------------------------
!    Permutations computed during analysis:
!       SYM_PERM: Symmetric permutation 
!       UNS_PERM: Column permutation (optional)
!    ---------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: SYM_PERM, UNS_PERM
! 
!    -----
!    Schur
!    -----
        INTEGER ::  NPROW, NPCOL, MBLOCK, NBLOCK
        INTEGER ::  SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
        INTEGER ::  SIZE_SCHUR
        COMPLEX(kind=8), DIMENSION(:), POINTER :: SCHUR
        COMPLEX(kind=8), DIMENSION(:), POINTER :: SCHUR_CINTERFACE
        INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
!    -------------------------------------
!    Case of distributed matrix on entry:
!    ZMUMPS potentially provides mapping
!    -------------------------------------
        INTEGER, DIMENSION(:), POINTER :: MAPPING
!    --------------
!    Version number
!    --------------
        CHARACTER(LEN=30) ::  VERSION_NUMBER
!    -----------
!    Out-of-core
!    -----------
        CHARACTER(LEN=255) :: OOC_TMPDIR
        CHARACTER(LEN=63) :: OOC_PREFIX
!    ------------------------------------------
!    To save the matrix in matrix market format
!    ------------------------------------------
        CHARACTER(LEN=255) ::  WRITE_PROBLEM
!    -----------
!    Save/Restore
!    -----------
        CHARACTER(LEN=255) :: SAVE_DIR
        CHARACTER(LEN=255)  :: SAVE_PREFIX
        CHARACTER(LEN=7)   ::  pad8  
!
!
! **********************
! INTERNAL Working data
! *********************
        INTEGER(8) :: KEEP8(150), MAX_SURF_MASTER
        INTEGER ::  INST_Number
!       For MPI
        INTEGER ::  COMM_NODES, MYID_NODES, COMM_LOAD
        INTEGER ::   MYID, NPROCS, NSLAVES
        INTEGER ::  ASS_IRECV
        INTEGER ::  LBUFR
        INTEGER ::  LBUFR_BYTES
        INTEGER, DIMENSION(:), POINTER ::  BUFR
!       IS is used for the factors + workspace for contrib. blocks
        INTEGER, DIMENSION(:), POINTER :: IS
!       IS1 (maxis1) contains working arrays computed 
!       and used only during analysis
        INTEGER, DIMENSION(:), POINTER :: IS1
!       For analysis/facto/solve phases
        INTEGER ::  MAXIS1, Deficiency
        INTEGER ::  KEEP(500)
!       The following data/arrays are computed during the analysis
!       phase and used during the factorization and solve phases.
        INTEGER ::  LNA
        INTEGER ::  NBSA
        INTEGER,POINTER,DIMENSION(:) :: STEP, NE_STEPS, ND_STEPS
!  Info for pruning tree 
        INTEGER,POINTER,DIMENSION(:) :: Step2node
!  ---------------------
        INTEGER,POINTER,DIMENSION(:) :: FRERE_STEPS, DAD_STEPS
        INTEGER,POINTER,DIMENSION(:) :: FILS, FRTPTR, FRTELT
        INTEGER(8),POINTER,DIMENSION(:) :: PTRAR
        INTEGER,POINTER,DIMENSION(:) :: NA, PROCNODE_STEPS
!       The two pointer arrays computed in facto and used by the solve
!          (except the factors) are PTLUST_S and PTRFAC. 
        INTEGER, DIMENSION(:), POINTER :: PTLUST_S
        INTEGER(8), DIMENSION(:), POINTER :: PTRFAC
!       main real working arrays for factorization/solve phases
        COMPLEX(kind=8), DIMENSION(:), POINTER :: S
!       Information on mapping
        INTEGER, DIMENSION(:), POINTER :: PROCNODE
!       Input matrix ready for numerical assembly 
!           -arrowhead format in case of assembled matrix
!           -element format otherwise
        INTEGER, DIMENSION(:), POINTER :: INTARR
        COMPLEX(kind=8), DIMENSION(:), POINTER :: DBLARR
!       Element entry: internal data
        INTEGER :: NELT_loc, LELTVAR
        INTEGER, DIMENSION(:), POINTER :: ELTPROC
!       Candidates and node partitionning
        INTEGER, DIMENSION(:,:), POINTER :: CANDIDATES
        INTEGER, DIMENSION(:),   POINTER :: ISTEP_TO_INIV2
        INTEGER, DIMENSION(:),   POINTER :: FUTURE_NIV2
        INTEGER, DIMENSION(:,:), POINTER :: TAB_POS_IN_PERE 
        LOGICAL, DIMENSION(:),   POINTER :: I_AM_CAND
!       For heterogeneous architecture
        INTEGER, DIMENSION(:), POINTER :: MEM_DIST
!       Compressed RHS
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP_ROW
        LOGICAL  :: POSINRHSCOMP_COL_ALLOC, pad11
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP_COL
        COMPLEX(kind=8), DIMENSION(:),   POINTER :: RHSCOMP
!       Info on the subtrees to be used during factorization
        DOUBLE PRECISION, DIMENSION(:), POINTER :: MEM_SUBTREE
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COST_TRAV
        INTEGER, DIMENSION(:),   POINTER :: MY_ROOT_SBTR
        INTEGER, DIMENSION(:),   POINTER :: MY_FIRST_LEAF
        INTEGER, DIMENSION(:),   POINTER :: MY_NB_LEAF
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST_SEQ
        INTEGER, DIMENSION(:),   POINTER :: SBTR_ID
        INTEGER, DIMENSION(:),   POINTER :: SCHED_DEP
        INTEGER, DIMENSION(:),   POINTER :: SCHED_GRP
        INTEGER, DIMENSION(:),   POINTER :: SCHED_SBTR
        INTEGER, DIMENSION(:),   POINTER :: CROIX_MANU
        COMPLEX(kind=8), DIMENSION(:), POINTER :: WK_USER
        INTEGER :: NBSA_LOCAL
        INTEGER :: LWK_USER
!    Internal control array
        DOUBLE PRECISION ::  DKEEP(230)
!    For simulating parallel out-of-core stack.
        DOUBLE PRECISION, DIMENSION(:),POINTER :: CB_SON_SIZE
!    Instance number used/managed by the C/F77 interface
        INTEGER ::  INSTANCE_NUMBER
!    OOC management data that must persist from factorization to solve.
        INTEGER ::  OOC_MAX_NB_NODES_FOR_ZONE
        INTEGER, DIMENSION(:,:),   POINTER :: OOC_INODE_SEQUENCE
        INTEGER(8),DIMENSION(:,:), POINTER :: OOC_SIZE_OF_BLOCK
        INTEGER(8), DIMENSION(:,:),   POINTER :: OOC_VADDR
        INTEGER,DIMENSION(:), POINTER :: OOC_TOTAL_NB_NODES
        INTEGER,DIMENSION(:), POINTER :: OOC_NB_FILES
        INTEGER :: OOC_NB_FILE_TYPE,pad12
        INTEGER,DIMENSION(:), POINTER :: OOC_FILE_NAME_LENGTH
        CHARACTER,DIMENSION(:,:), POINTER :: OOC_FILE_NAMES  
!    Indices of nul pivots
        INTEGER,DIMENSION(:), POINTER :: PIVNUL_LIST
!    Array needed to manage additionnal candidate processor 
        INTEGER, DIMENSION(:,:), POINTER :: SUP_PROC, pad14
!    Lists of nodes where processors work. Built/used in solve phase.
        INTEGER, DIMENSION(:), POINTER :: IPTR_WORKING, WORKING
!    Root structure(internal)
        TYPE (ZMUMPS_ROOT_STRUC) :: root
!    Low-rank
        INTEGER, POINTER, DIMENSION(:) :: LRGROUPS
        INTEGER :: NBGRP,pad13
!    Pointer encoding for FDM_F data
        CHARACTER(LEN=1), DIMENSION(:), POINTER :: FDM_F_ENCODING
!    Pointer array encoding BLR factors pointers
        CHARACTER(LEN=1), DIMENSION(:), POINTER :: BLRARRAY_ENCODING
!    Multicore
        INTEGER :: LPOOL_AFTER_L0_OMP, LPOOL_BEFORE_L0_OMP
        INTEGER :: L_PHYS_L0_OMP
        INTEGER :: L_VIRT_L0_OMP                                    
        INTEGER :: LL0_OMP_MAPPING,pad15
        INTEGER(8) :: THREAD_LA
! Pool before L0_OMP
        INTEGER, DIMENSION(:), POINTER :: IPOOL_BEFORE_L0_OMP
! Pool after L0_OMP
        INTEGER, DIMENSION(:), POINTER :: IPOOL_AFTER_L0_OMP
! Subtrees
        INTEGER, DIMENSION(:), POINTER :: PHYS_L0_OMP
! Amalgamated subtrees
        INTEGER, DIMENSION(:), POINTER :: VIRT_L0_OMP
! From heaviest to lowest subtree
        INTEGER, DIMENSION(:), POINTER :: PERM_L0_OMP
! To get leafs in global pool
        INTEGER, DIMENSION(:), POINTER :: PTR_LEAFS_L0_OMP
! Mapping of the subtrees
        INTEGER, DIMENSION(:), POINTER :: L0_OMP_MAPPING
! for RR on root
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SINGULAR_VALUES
        INTEGER ::  NB_SINGULAR_VALUES
! To know if OOC files are associated to a saved and so if they should be removed.
        LOGICAL :: ASSOCIATED_OOC_FILES
      END TYPE ZMUMPS_STRUC
