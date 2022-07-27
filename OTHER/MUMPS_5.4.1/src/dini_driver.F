C
C  This file is part of MUMPS 5.4.1, released
C  on Tue Aug  3 09:49:43 UTC 2021
C
C
C  Copyright 1991-2021 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      SUBROUTINE DMUMPS_INI_DRIVER( id )
      USE DMUMPS_STRUC_DEF
C
C  Purpose:
C  =======
C
C  Initialize an instance of the DMUMPS package.
C
      USE DMUMPS_BUF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_STRUC) id
      INTEGER MASTER, IERR,PAR_loc,SYM_loc
      PARAMETER( MASTER = 0 )
      INTEGER color
C     -----------------------------
C     Initialize MPI related data
C     -----------------------------
      CALL MPI_COMM_SIZE(id%COMM, id%NPROCS, IERR )
C     Now done in the main MUMPS driver:
C     CALL MPI_COMM_RANK(id%COMM, id%MYID, IERR )
C
      PAR_loc=id%PAR
      SYM_loc=id%SYM
C     Broadcasting PAR/SYM (KEEP(46)/KEEP(50)) in order to
C     have only one value available: the one from the master
      CALL MPI_BCAST(PAR_loc,1,MPI_INTEGER,MASTER,id%COMM,IERR)
      CALL MPI_BCAST(SYM_loc,1,MPI_INTEGER,MASTER,id%COMM,IERR)
C     Initialize a subcommunicator
C     for slave nodes
C
      IF ( PAR_loc .eq. 0 ) THEN
C       -------------------
C       Host is not working
C       -------------------
        IF ( id%MYID .eq. MASTER ) THEN
          color = MPI_UNDEFINED
        ELSE
          color = 0
        END IF
        CALL MPI_COMM_SPLIT( id%COMM, color, 0,
     &                       id%COMM_NODES, IERR )
        id%NSLAVES = id%NPROCS - 1
      ELSE
C       ----------------
C       Host is working
C       ----------------
        CALL MPI_COMM_DUP( id%COMM, id%COMM_NODES, IERR )
        id%NSLAVES = id%NPROCS
      END IF
C     ---------------------------
C     Use same slave communicator
C     for load information
C     ---------------------------
      IF (PAR_loc .ne. 0 .or. id%MYID .NE. MASTER) THEN
        CALL MPI_COMM_DUP( id%COMM_NODES, id%COMM_LOAD, IERR )
      ENDIF
C     ----------------------------------------------
C     Initialize default values for CNTL,ICNTL,KEEP,KEEP8
C     potentially depending on id%SYM and id%NSLAVES
C     ----------------------------------------------
      CALL DMUMPSID( id%NSLAVES, id%LWK_USER,
     &    id%CNTL(1), id%ICNTL(1),
     &    id%KEEP(1), id%KEEP8(1), id%INFO(1), id%INFOG(1),
     &    id%RINFO(1), id%RINFOG(1),
     &    SYM_loc, PAR_loc, id%DKEEP(1), id%MYID )
      id%WRITE_PROBLEM="NAME_NOT_INITIALIZED"
      CALL MUMPS_SET_VERSION( id%VERSION_NUMBER )
      id%OOC_TMPDIR="NAME_NOT_INITIALIZED"
      id%OOC_PREFIX="NAME_NOT_INITIALIZED"
      id%SAVE_DIR="NAME_NOT_INITIALIZED"
      id%SAVE_PREFIX="NAME_NOT_INITIALIZED"
C     Default value for NRHS is 1
      id%NRHS = 1
C     Leading dimension will be reset to id%N is DMUMPS_SOL_DRIVER
C     if id%NRHS remains equal to 1. Otherwise id%LRHS must be
C     set by user.
      id%LRHS = 0 ! Value will be checked in DMUMPS_CHECK_DENSE_RHS
                  ! Not accessed if id%NRHS=1
C     Similar behaviour for LREDRHS (value will
C     be checked in DMUMPS_CHECK_REDRHS)
      id%LREDRHS = 0
C
C     Module needs to know the size of an INTEGER
      CALL DMUMPS_BUF_INIT( id%KEEP( 34 ), id%KEEP(35) )
C
      id%INST_Number = -1
C
C     Define the options for Metis
C
      id%METIS_OPTIONS(:) = 0
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)      
#if defined(metis4) || defined(parmetis3)
C     Useful size is 8
C     set to default options
      id%METIS_OPTIONS(1) = 0
#else
C     Useful size is 40  
C     This sets the default values
      CALL METIS_SETDEFAULTOPTIONS(id%METIS_OPTIONS)
C     This number, 18, corresponds to METIS_OPTIONS_NUMBERING which
C     tells METIS to use fortran numbering and is found in metis.h
C     In Metis 5.0.3 and Parmetis 4.0.2, METIS_OPTIONS_NUMBERING 
C     was METIS_OPTIONS(17). MUMPS doesnot support those versions anymore.
C     To use them, just change METIS_OPTIONS(18) into METIS_OPTIONS(17)
C     like that: METIS_OPTIONS(17) = 1
      id%METIS_OPTIONS(18) = 1
#endif
#endif      
C
C  Nullify a few pointers and integers
C
      id%N = 0; id%NZ = 0; id%NNZ = 0_8
      NULLIFY(id%IRN)
      NULLIFY(id%JCN)
      NULLIFY(id%A)
      id%NZ_loc = 0; id%NNZ_loc = 0_8
      NULLIFY(id%IRN_loc)
      NULLIFY(id%JCN_loc)
      NULLIFY(id%A_loc)
      NULLIFY(id%MAPPING)
      NULLIFY(id%RHS)
      NULLIFY(id%REDRHS)
      id%NZ_RHS=0
      NULLIFY(id%RHS_SPARSE)
      NULLIFY(id%IRHS_SPARSE)
      NULLIFY(id%IRHS_PTR)
      NULLIFY(id%ISOL_loc)
      NULLIFY(id%IRHS_loc)
      id%LSOL_loc=0
      id%LRHS_loc=0
      id%Nloc_RHS=0
      NULLIFY(id%SOL_loc)
      NULLIFY(id%RHS_loc)
      NULLIFY(id%COLSCA)
      NULLIFY(id%ROWSCA)
      NULLIFY(id%PERM_IN)
      NULLIFY(id%IS)
      NULLIFY(id%STEP)
C     Info for analysis by block
      id%NBLK = 0
      NULLIFY(id%BLKPTR)
      NULLIFY(id%BLKVAR)
C     Info for pruning tree
      NULLIFY(id%Step2node)
      NULLIFY(id%DAD_STEPS)
      NULLIFY(id%NE_STEPS)
      NULLIFY(id%ND_STEPS)
      NULLIFY(id%FRERE_STEPS)
      NULLIFY(id%SYM_PERM)
      NULLIFY(id%UNS_PERM)
      NULLIFY(id%PIVNUL_LIST)
      NULLIFY(id%FILS)
      NULLIFY(id%PTRAR)
      NULLIFY(id%FRTPTR)
      NULLIFY(id%FRTELT)
      NULLIFY(id%NA)
      id%LNA=0
      NULLIFY(id%PROCNODE_STEPS)
      NULLIFY(id%S)
      NULLIFY(id%PTLUST_S)
      NULLIFY(id%PTRFAC)
      NULLIFY(id%INTARR) 
      NULLIFY(id%DBLARR)
      NULLIFY(id%DEPTH_FIRST)
      NULLIFY(id%DEPTH_FIRST_SEQ)
      NULLIFY(id%SBTR_ID)
      NULLIFY(id%SCHED_DEP)
      NULLIFY(id%SCHED_SBTR)
      NULLIFY(id%SCHED_GRP)
      NULLIFY(id%CROIX_MANU)
      NULLIFY(id%WK_USER)
      NULLIFY(id%MEM_SUBTREE)
      NULLIFY(id%MEM_SUBTREE)
      NULLIFY(id%MY_ROOT_SBTR)
      NULLIFY(id%MY_FIRST_LEAF)
      NULLIFY(id%MY_NB_LEAF)
      NULLIFY(id%COST_TRAV)
      NULLIFY(id%RHSCOMP)
      NULLIFY(id%POSINRHSCOMP_ROW)
      NULLIFY(id%POSINRHSCOMP_COL)
      id%POSINRHSCOMP_COL_ALLOC = .FALSE.
C
C     Out of Core management related data
C
      NULLIFY(id%OOC_INODE_SEQUENCE)
      NULLIFY(id%OOC_TOTAL_NB_NODES)
      NULLIFY(id%OOC_SIZE_OF_BLOCK)
      NULLIFY(id%OOC_FILE_NAME_LENGTH)
      NULLIFY(id%OOC_FILE_NAMES)
      NULLIFY(id%OOC_VADDR)
      NULLIFY(id%OOC_NB_FILES)
      NULLIFY(id%LRGROUPS)
      NULLIFY(id%FDM_F_ENCODING)
      NULLIFY(id%BLRARRAY_ENCODING)
      NULLIFY(id%MPITOOMP_PROCS_MAP)
C     Must be nullified because of routine
C     DMUMPS_SIZE_IN_STRUCT
      NULLIFY(id%CB_SON_SIZE)
C
C     Components of the root
C
      NULLIFY(id%root%RHS_CNTR_MASTER_ROOT)
      NULLIFY(id%root%RHS_ROOT)
      NULLIFY(id%root%RG2L_ROW)
      NULLIFY(id%root%RG2L_COL)
      NULLIFY(id%root%IPIV)
      NULLIFY(id%root%SCHUR_POINTER)
      NULLIFY(id%SCHUR_CINTERFACE)
C
C     Element-entry
C
      id%NELT=0
      NULLIFY(id%ELTPTR)
      NULLIFY(id%ELTVAR)
      NULLIFY(id%A_ELT)
      NULLIFY(id%ELTPROC)
C
C     Schur
C
      id%SIZE_SCHUR = 0
      NULLIFY( id%LISTVAR_SCHUR )
      NULLIFY( id%SCHUR )
C     -- Distributed Schur
      id%NPROW      = 0
      id%NPCOL      = 0
      id%MBLOCK     = 0
      id%NBLOCK     = 0
      id%SCHUR_MLOC = 0 ! Exit from analysis
      id%SCHUR_NLOC = 0 ! Exit from analysis
      id%SCHUR_LLD  = 0
C
C     Candidates and node partitionning
C
      NULLIFY(id%ISTEP_TO_INIV2)
      NULLIFY(id%I_AM_CAND)
      NULLIFY(id%FUTURE_NIV2)
      NULLIFY(id%TAB_POS_IN_PERE)
      NULLIFY(id%CANDIDATES)
      id%OOC_NB_FILE_TYPE=-123456
C
C     Initializations for L0_OMP mechanisms
C
      NULLIFY(id%IPOOL_B_L0_OMP)
      NULLIFY(id%IPOOL_A_L0_OMP)
      NULLIFY(id%PHYS_L0_OMP)
      NULLIFY(id%VIRT_L0_OMP)
      NULLIFY(id%VIRT_L0_OMP_MAPPING)
      NULLIFY(id%PERM_L0_OMP)
      NULLIFY(id%PTR_LEAFS_L0_OMP)
      NULLIFY(id%L0_OMP_MAPPING)
      NULLIFY(id%L0_OMP_FACTORS)
      NULLIFY(id%I4_L0_OMP)
      NULLIFY(id%I8_L0_OMP)
      id%LPOOL_B_L0_OMP = 0
      id%LPOOL_A_L0_OMP  = 0
      id%L_VIRT_L0_OMP       = 0
      id%L_PHYS_L0_OMP       = 0
      id%THREAD_LA           = 0
C
C     Mapping information used during solve.
C
      NULLIFY(id%IPTR_WORKING)
      NULLIFY(id%WORKING)
C
C     Initializations for Rank detection/null space
C
      NULLIFY(id%SINGULAR_VALUES)
      CALL DMUMPS_RR_INIT_POINTERS(id)
C     Architecture data
      NULLIFY(id%MEM_DIST)
C     Must be nullified because of routine
C     DMUMPS_SIZE_IN_STRUCT
      NULLIFY(id%SUP_PROC)
      id%Deficiency = 0
      id%root%LPIV = -1
      id%root%yes  = .FALSE.
      id%root%gridinit_done  = .FALSE.
C     NOT IN SAVE/RESTORE
      id%ASSOCIATED_OOC_FILES=.FALSE.
C
C     ----------------------------------------
C     Find MYID_NODES relatively to COMM_NODES
C     If  the calling processor is not inside
C     COMM_NODES, MYID_NODES will not be
C     significant / used anyway
C     ----------------------------------------
      IF ( id%KEEP( 46 ) .ne. 0  .OR.
     &     id%MYID .ne. MASTER ) THEN
          CALL MPI_COMM_RANK
     &         (id%COMM_NODES, id%MYID_NODES, IERR )
      ELSE
          id%MYID_NODES = -464646
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_INI_DRIVER
