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
C
      SUBROUTINE DMUMPS_ANA_DRIVER(id)
      USE DMUMPS_LOAD
      USE MUMPS_STATIC_MAPPING
      USE DMUMPS_STRUC_DEF
      USE MUMPS_MEMORY_MOD
      USE DMUMPS_PARALLEL_ANALYSIS
      USE DMUMPS_ANA_LR
      USE DMUMPS_LR_CORE
      USE DMUMPS_LR_STATS
      USE MUMPS_LR_COMMON
      USE DMUMPS_ANA_AUX_M
      USE MUMPS_ANA_BLK_M, ONLY: COMPACT_GRAPH_T, LMATRIX_T
      IMPLICIT NONE
C     
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
C
C     Purpose
C     =======
C
C     Performs analysis and (if required) Max-trans on the master, then
C     broadcasts information to the slaves. Also includes mapping.
C     
C     
C     Parameters
C     ==========
C     
      TYPE(DMUMPS_STRUC), TARGET :: id
C     
C     Local variables
C     ===============
C     
C     
C     Pointers inside integer array, various data
      INTEGER IKEEP, NE, NA
      INTEGER I, allocok
C     Other locals
      INTEGER NB_NIV2, IDEST
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER LOCAL_M, LOCAL_N
      INTEGER numroc
      EXTERNAL numroc
      INTEGER IRANK
      INTEGER MP, LP, MPG
      LOGICAL PROK, PROKG, LISTVAR_SCHUR_2BE_FREED, LPOK
      INTEGER SIZE_SCHUR_PASSED
      INTEGER SBUF_SEND_FR, SBUF_REC_FR
      INTEGER SBUF_SEND_LR, SBUF_REC_LR
      INTEGER TOTAL_MBYTES
      INTEGER(8) SBUF_RECOLD8, MIN_BUF_SIZE8
      INTEGER MIN_BUF_SIZE
      INTEGER(8) MAX_SIZE_FACTOR_TMP
      INTEGER LEAF, INODE, ISTEP, INN, LPTRAR
      INTEGER NBLEAF, NBROOT, MYROW_CHECK, INIV2
      DOUBLE PRECISION TIMEG
      INTEGER(8) ::  MAX_FRONT_SURFACE_LOCAL_L0, 
     &               MAX_SIZE_FACTOR_L0,
     &               ENTRIES_IN_FACTORS_UNDER_L0, 
     &               ENTRIES_IN_FACTORS_MASTERS_LO
      INTEGER             :: MAXFR_UNDER_L0
      DOUBLE PRECISION    :: COST_SUBTREES_UNDER_L0, OPSA_UNDER_L0
C     to store the size of the sequencial peak of stack
C     (or an estimation for not calling REORDER_TREE_N )
      DOUBLE PRECISION      :: PEAK
      INTEGER(8)::  SIZECB_UNDER_L0, SIZECB_UNDER_L0_IF_LRCB
      LOGICAL   :: ABOVE_L0
C     
C     INTEGER WORKSPACE 
C     
      INTEGER, ALLOCATABLE, DIMENSION(:):: IPOOL
      INTEGER                           :: LIPOOL
      INTEGER, ALLOCATABLE, TARGET, DIMENSION(:) :: PAR2_NODES
      INTEGER, DIMENSION(:), POINTER             :: PAR2_NODESPTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PROCNODE
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWtemp
      INTEGER, DIMENSION(:), ALLOCATABLE :: XNODEL, NODEL
      INTEGER, DIMENSION(:), POINTER :: SSARBR
C     Element matrix entry
      INTEGER, POINTER ::  NELT, LELTVAR
      INTEGER, DIMENSION(:), POINTER :: KEEP, INFO, INFOG
      INTEGER(8), DIMENSION(:), POINTER :: KEEP8
      INTEGER(8)                   :: ENTRIES_IN_FACTORS_LOC_MASTERS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFO
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFOG
      INTEGER, DIMENSION(:), POINTER :: ICNTL
      LOGICAL :: I_AM_SLAVE, PERLU_ON, COND
      INTEGER :: OOC_STRAT, BLR_STRAT
      INTEGER :: IDUMMY
      INTEGER, TARGET    :: IDUMMY_ARRAY(1)
      INTEGER, POINTER, DIMENSION(:) :: IRN_loc_PTR
      INTEGER, POINTER, DIMENSION(:) :: JCN_loc_PTR
      INTEGER, POINTER, DIMENSION(:) :: IRN_PTR
      INTEGER, POINTER, DIMENSION(:) :: JCN_PTR
      INTEGER, POINTER, DIMENSION(:) :: SIZEOFBLOCKS_PTR
      INTEGER, POINTER, DIMENSION(:) :: UNS_PERM_PTR
      LOGICAL :: BDUMMY
      INTEGER(8) :: K8_33relaxed, K8_34relaxed, K8_35relaxed,
     &              K8_50relaxed
      LOGICAL :: SUM_OF_PEAKS
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      INTEGER, EXTERNAL :: MUMPS_ENCODE_TPN_IPROC
      INTEGER :: PROCNODE_VALUE
      INTEGER K,J, IFS
      INTEGER SIZE_TEMP_MEM,SIZE_DEPTH_FIRST,SIZE_COST_TRAV
      LOGICAL IS_BUILD_LOAD_MEM_CALLED
      LOGICAL PRINT_MAXAVG
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: TEMP_MEM
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_ROOT
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_LEAF
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_SIZE
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH_FIRST
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH_FIRST_SEQ
      INTEGER, DIMENSION (:), ALLOCATABLE :: SBTR_ID
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: COST_TRAV_TMP
      INTEGER(8) :: TOTAL_BYTES, ITMP8
      INTEGER          :: SIZE_PAR2_NODESPTR
      INTEGER          :: LSIZEOFBLOCKS_PTR
      LOGICAL          :: READY_FOR_ANA_F
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAPCOL
      LOGICAL :: BLKPTR_ALLOCATED, BLKVAR_ALLOCATED
      INTEGER :: IB, BLKSIZE
      INTEGER :: IBcurrent, IPOS, IPOSB, II
C      Internal work arrays:
C      DOF2BLOCK(idof)=inode, idof in [1,N], inode in [1,NBLK]
C      SIZEBLOCK(1:NBLK) (for node valuation) 
      INTEGER, TARGET, DIMENSION(:), allocatable:: SIZEOFBLOCKS
      INTEGER, DIMENSION(:), allocatable:: DOF2BLOCK
      INTEGER    :: NBRECORDS
      INTEGER(8) :: NSEND8, NLOCAL8
C     LMAT_BLOCK: in case of centralized matrix, 
C                 to store on MASTER the cleaned Lmatrix
C                 used to compute GCOMP
C                 LMAT_BLOCK might also be saved to 
C                 be used during grouping
C     LUMAT     : in case of distributed matrix
C                 to store distributed the cleaned LU matrix
C                 LUMAT might also be saved to 
C                 be used for MPI based grouping
C     LUMAT_REMAP : in case of distributed matrix
C                 it is used to remap LUMAT
C 
C     GCOMP     : Graph "ready" to be called by orderings
C
      TYPE(LMATRIX_T)         :: LMAT_BLOCK, LUMAT, LUMAT_REMAP
      LOGICAL                 :: GCOMP_PROVIDED
      TYPE(COMPACT_GRAPH_T)   :: GCOMP
      TYPE(COMPACT_GRAPH_T)   :: GCOMP_DIST
      INTEGER, POINTER, DIMENSION(:) ::  
     &     NFSIZPTR,
     &     FILSPTR,
     &     FREREPTR, NE_STEPSPTR,
     &     IKEEP1, IKEEP2, IKEEP3,
     &     STEPPTR, LRGROUPSPTR
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET :: IKEEPALLOC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WORK2ALLOC
      ! Used because of multithreaded SIM_NP_
      INTEGER :: locMYID, locMYID_NODES
      LOGICAL, POINTER :: locI_AM_CAND(:)
      INTEGER(kind=8) :: NZ8, LIW8
C     NBLK : id%N or order of blocked matrix
      INTEGER :: NBLK
      INTEGER :: LIW_ELT
C
      INTERFACE
C     Explicit interface because of pointer arguments:
      SUBROUTINE DMUMPS_FREE_ID_DATA_MODULES(id_FDM_F_ENCODING,
     &  id_BLRARRAY_ENCODING, KEEP8)
      USE MUMPS_FRONT_DATA_MGT_M, only : MUMPS_FDM_STRUC_TO_MOD,
     &                                   MUMPS_FDM_END
      USE DMUMPS_LR_DATA_M, only : DMUMPS_BLR_STRUC_TO_MOD,
     &                             DMUMPS_BLR_END_MODULE
#     if defined(MUMPS_F2003)
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER, intent(inout) ::
     &                                            id_FDM_F_ENCODING
#     else
      CHARACTER, DIMENSION(:), POINTER :: id_BLRARRAY_ENCODING
      CHARACTER, DIMENSION(:), POINTER :: id_FDM_F_ENCODING
#     endif
      INTEGER(8), intent(inout) :: KEEP8(150)
      END SUBROUTINE DMUMPS_FREE_ID_DATA_MODULES
      END INTERFACE
C
C  Beginning of executable statements
C
      IS_BUILD_LOAD_MEM_CALLED=.FALSE.
      KEEP   => id%KEEP
      KEEP8  => id%KEEP8
      INFO   => id%INFO
      RINFO  => id%RINFO
      INFOG  => id%INFOG
      RINFOG => id%RINFOG
      ICNTL  => id%ICNTL
      NELT    => id%NELT
      LELTVAR => id%LELTVAR
      KEEP(264) = 0    ! reinitialise out-of-range status (0=yes)
      KEEP(265) = 0    ! reinitialise dupplicates (0=yes)
      PRINT_MAXAVG = .NOT.(id%NSLAVES.EQ.1 .AND. KEEP(46).EQ.1)
      NULLIFY ( NFSIZPTR,
     &     FILSPTR,
     &     FREREPTR, NE_STEPSPTR,
     &     IKEEP1, IKEEP2, IKEEP3, STEPPTR, LRGROUPSPTR, 
     &     SSARBR, SIZEOFBLOCKS_PTR, IRN_loc_PTR, JCN_loc_PTR,
     &     IRN_PTR, JCN_PTR,
     &     PAR2_NODESPTR )
      IF (associated(id%UNS_PERM)) DEALLOCATE(id%UNS_PERM)
      nullify(id%UNS_PERM)
      IDUMMY = 1
      BDUMMY = .FALSE.
C     Set default value that witl be reset in
C     case of blocked format matrices
      NBLK = id%N
      GCOMP_PROVIDED   = .FALSE.
      BLKPTR_ALLOCATED = .FALSE.
      BLKVAR_ALLOCATED = .FALSE.
C     -------------------------------------
C     Depending on the type of parallelism,
C     the master can now (soon) potentially
C     have the role of a slave
C     -------------------------------------
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &     ( id%MYID .eq. MASTER .AND.
     &     id%KEEP(46) .eq. 1 ) )
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
C     LP     : errors
C     MP     : INFO
      LPOK  = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK  = (( MP  .GT. 0 ).AND.(ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(ICNTL(4).GE.2))
      IF ( PROK ) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            WRITE(MP, '(A)') 'L U Solver for unsymmetric matrices'
         ELSE IF ( KEEP(50) .eq. 1 ) THEN
            WRITE(MP, '(A)') 
     & 'L D L^T Solver for symmetric positive definite matrices'
         ELSE
            WRITE(MP, '(A)') 
     &           'L D L^T Solver for general symmetric matrices'
         END IF
         IF ( KEEP(46) .eq. 1 ) THEN
            WRITE(MP, '(A)') 'Type of parallelism: Working host'
         ELSE
            WRITE(MP, '(A)') 'Type of parallelism: Host not working'
         END IF
      END IF
      IF ( PROKG .AND. (MP.NE.MPG)) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            WRITE(MPG, '(A)') 'L U Solver for unsymmetric matrices'
         ELSE IF ( KEEP(50) .eq. 1 ) THEN
            WRITE(MPG, '(A)') 
     & 'L D L^T Solver for symmetric positive definite matrices'
         ELSE
            WRITE(MPG, '(A)') 
     &           'L D L^T Solver for general symmetric matrices'
         END IF
         IF ( KEEP(46) .eq. 1 ) THEN
            WRITE(MPG, '(A)') 'Type of parallelism: Working host'
         ELSE
            WRITE(MPG, '(A)') 'Type of parallelism: Host not working'
         END IF
      END IF      
      IF (PROK) WRITE( MP, 110 )
      IF (PROKG .AND. (MPG.NE.MP)) WRITE( MPG, 110 )
C
C     BEGIN CASE OF ALLOCATED DATA FROM PREVIOUS CALLS
C     ----------------------------------------
C     Free some memory from factorization,
C     if allocated, at least large arrays.
C     This will also limit the amount of useless
C     data saved to disk in case of save-restore
C     ----------------------------------------
      IF (id%KEEP8(24).EQ.0_8) THEN
C       -- deallocate only when not provided/allocated by the user
         IF (associated(id%S)) THEN
            DEALLOCATE(id%S)
            id%KEEP8(23)=0_8
         ENDIF
      ENDIF
      NULLIFY(id%S)
      KEEP8(24) = 0_8  ! reinitialize last used size of WK_USER
      IF (associated(id%IS)) THEN
        DEALLOCATE(id%IS)
        NULLIFY(id%IS)
      ENDIF
C     also avoid keeping BLR factors allocated if analysis
C     called after a previous BLR factorization without
C     an intermediate JOB=-2 call.
      CALL DMUMPS_FREE_ID_DATA_MODULES(id%FDM_F_ENCODING,
     &            id%BLRARRAY_ENCODING, id%KEEP8(1))
      IF (associated(id%root%RG2L_ROW))THEN
        DEALLOCATE(id%root%RG2L_ROW)
        NULLIFY(id%root%RG2L_ROW)
      ENDIF
      IF (associated(id%root%RG2L_COL))THEN
        DEALLOCATE(id%root%RG2L_COL)
        NULLIFY(id%root%RG2L_COL)
      ENDIF
      IF (associated( id%PTLUST_S )) THEN
        DEALLOCATE(id%PTLUST_S)
        NULLIFY(id%PTLUST_S)
      ENDIF
      IF (associated(id%PTRFAC)) THEN
        DEALLOCATE(id%PTRFAC)
        NULLIFY(id%PTRFAC)
      END IF
      IF (associated(id%RHSCOMP)) THEN
        DEALLOCATE(id%RHSCOMP)
        NULLIFY(id%RHSCOMP)
        id%KEEP8(25)=0_8
      ENDIF
      IF (associated(id%POSINRHSCOMP_ROW)) THEN
        DEALLOCATE(id%POSINRHSCOMP_ROW)
        NULLIFY(id%POSINRHSCOMP_ROW)
      ENDIF
      IF (id%POSINRHSCOMP_COL_ALLOC) THEN
        DEALLOCATE(id%POSINRHSCOMP_COL)
        NULLIFY(id%POSINRHSCOMP_COL)
        id%POSINRHSCOMP_COL_ALLOC = .FALSE.
      ENDIF
C     --------------------------------------------
C     If analysis redone, suppress old,
C     meaningless, Step2node array.
C     This is necessary since we could otherwise
C     end up having a wrong Step2node during solve
C     --------------------------------------------
      IF (associated(id%Step2node)) THEN
        DEALLOCATE(id%Step2node)
        NULLIFY(id%Step2node)
      ENDIF
C     END CASE OF ALLOCATED DATA FROM PREVIOUS CALLS
C
C     Decode API (ICNTL parameters, mainly)
C     and check consistency of the KEEP array.
C     Note: DMUMPS_ANA_CHECK_KEEP also sets
C     some INFOG parameters
      CALL DMUMPS_ANA_CHECK_KEEP(id)
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
C     -------------------------------------------
C     Broadcast KEEP(60) since we need to broadcast
C     related information
C     ------------------------------------------
      CALL MPI_BCAST( KEEP(60), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C        broadcast also size of schur
      IF (id%KEEP(60) .NE. 0 ) THEN
       CALL MPI_BCAST( KEEP(116), 1, MPI_INTEGER, MASTER, 
     &                   id%COMM, IERR )
      ENDIF
      IF (id%KEEP(60) .EQ. 2 .or. id%KEEP(60). EQ. 3) THEN
         CALL MPI_BCAST( id%NPROW, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%NPCOL, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%MBLOCK, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%NBLOCK, 1,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
C     Note that DMUMPS_INIT_ROOT_ANA will
C     then use that information.
      ENDIF
C     ----------------------------------------------
C     Broadcast KEEP(54) now to know if the
C     structure of the graph is intially distributed
C     and should be assembled on the master
C     Broadcast KEEP(55) now to know if the
C     matrix is in assembled or elemental format
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(54), 2, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast KEEP(69) now to know if
C     we will need to communicate during analysis
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(69), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast Out of core strategy (used only on master so far)
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast analysis strategy (used only on master so far)
C     ----------------------------------------------
      CALL MPI_BCAST( KEEP(244), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ---------------------------
C     Fwd in facto
C     Broadcast KEEP(251,252,253) defined on master so far
      CALL MPI_BCAST( KEEP(251), 3, MPI_INTEGER,MASTER,id%COMM,IERR)
C
      CALL MPI_BCAST( id%KEEP(490), 5, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast N 
C     ----------------------------------------------
      CALL MPI_BCAST( id%N, 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C     ----------------------------------------------
C     Broadcast NZ for assembled entry
C     ----------------------------------------------
      IF ( KEEP(55) .EQ. 0) THEN
         IF ( KEEP(54) .eq. 3 ) THEN
C     Compute total number of non-zeros
          CALL MPI_ALLREDUCE( id%KEEP8(29), id%KEEP8(28), 1,
     &       MPI_INTEGER8, 
     &       MPI_SUM, id%COMM, IERR )
         ELSE
C     Broadcast NZ from the master node
            CALL MPI_BCAST( id%KEEP8(28), 1, MPI_INTEGER8, MASTER,
     &           id%COMM, IERR )
         END IF
      ELSE
C     Broadcast NA_ELT <=> KEEP8(30) for elemental entry
         CALL MPI_BCAST( id%KEEP8(30), 1, MPI_INTEGER8, MASTER,
     &        id%COMM, IERR )
      ENDIF
      IF( id%KEEP(54).EQ.3) THEN
C     test IRN_loc and JCN_loc allocated on working procs
       IF (I_AM_SLAVE .AND. id%KEEP8(29).GT.0 .AND.
     &     ( (.NOT. associated(id%IRN_loc)) .OR. 
     &       (.NOT. associated(id%JCN_loc)) )
     &   ) THEN
         id%INFO(1) = -22
         id%INFO(2) = 16
       ENDIF
      ENDIF
      IF ( associated(id%MEM_DIST) ) THEN
         DEALLOCATE( id%MEM_DIST )
      ENDIF
      allocate( id%MEM_DIST( 0:id%NSLAVES-1 ), STAT=IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO(1) = -7
         INFO(2) = id%NSLAVES
         IF ( LPOK ) THEN
            WRITE(LP, 150) 'MEM_DIST'
         END IF
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
      id%MEM_DIST(0:id%NSLAVES-1) = 0
      CALL MUMPS_INIT_ARCH_PARAMETERS(
     &     id%COMM,id%COMM_NODES,KEEP(69),KEEP(46),
     &     id%NSLAVES,id%MEM_DIST,INFO)
C     ========================
C     Write problem to a file,
C     if requested by the user
C     ========================
      CALL DMUMPS_DUMP_PROBLEM(id)
C     =================
C     ANALYSIS BY BLOCK 
C     =================
      IF ( id%MYID .EQ. MASTER ) THEN
       IF (KEEP(13).NE.0) THEN
C      Analysis by block with block data provided by user
C
C       Check if block structure is centralized or distributed
        IF (.NOT.associated(id%BLKVAR)) THEN
C        BLKVAR is identity and implicitly centralized
         KEEP(14) = 0
        ELSE
         IF (size(id%BLKVAR).EQ.id%N) THEN
C         Centralized block stucture
          KEEP(14) = 0
         ELSE
C         Distributed block stucture
          KEEP(14) = 1
          IF ( LPOK ) THEN
             WRITE(LP,'(A,A,I8)') 
     &       " ERROR with centralized matrix. Size of id%BLKVAR ", 
     &       "should be equal to id%N instead of ", 
     &       size(id%BLKVAR)
          ENDIF
          id%INFO(1) = -57
          id%INFO(2) = 3
         ENDIF
        ENDIF
        IF (KEEP(13).GE.1) THEN
C        BLKPTR provided by user
C        check input data
         IF ( .NOT.associated(id%BLKPTR)) THEN
          IF ( LPOK ) THEN
             WRITE(LP,'(A,I8)') 
     &       " id%BLKPTR should be provided by user on host "
          ENDIF
          id%INFO(1) = -57
          id%INFO(2) = 2
         ENDIF 
         IF ( (id%NBLK.LE.0).OR.(id%NBLK.GT.id%N)
     &       .OR. (id%NBLK+1.NE.size(id%BLKPTR))
     &      ) THEN
          IF ( LPOK ) THEN
             WRITE(LP,'(A,I8)') 
     &       " ERROR incorrect value of id%NBLK:", id%NBLK
          ENDIF
          id%INFO(1) = -57
          id%INFO(2) = 1
         ENDIF 
         NBLK=id%NBLK
         IF (id%BLKPTR(id%NBLK+1)-1.NE.id%N) THEN
           IF ( LPOK ) THEN
             WRITE(LP,'(A,A,I8)') 
     &       " ERROR id%BLKPTR(id%NBLK+1)-1 ", 
     &       "should be equal to id%N instead of ", 
     &       id%BLKPTR(id%NBLK+1)-1
           ENDIF
           id%INFO(1) = -57
           id%INFO(2) = 2
         ENDIF
         IF (id%BLKPTR(1).NE.1) THEN
           IF ( LPOK ) THEN
             WRITE(LP,'(A,A,I8)') 
     &       " ERROR id%BLKPTR(1)", 
     &       "should be equal to 1 instead of ", 
     &       id%BLKPTR(1)
           ENDIF
           id%INFO(1) = -57
           id%INFO(2) = 2
         ENDIF
        ELSE IF (KEEP(13).LT.0) THEN
C        regular blocks in BLKVAR of size -KEEP(13)
C        mod(id%N,-KEEP(13)) has already been checked
         NBLK = id%N/(-KEEP(13))
        ENDIF
C      end of KEEP(13).NE.0
       ENDIF 
C      end of  id%MYID .EQ. MASTER 
      ENDIF 
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 500
C
C     Broadcast KEEP(13-14), NBLK
      CALL MPI_BCAST( KEEP(13), 2, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( NBLK, 1, MPI_INTEGER, MASTER, id%COMM, IERR )
C
C     ===========================
      IF (KEEP(13).NE.0) THEN
C     { BEGIN preparation ANA_BLK
C     ===========================
       IF ( ( (KEEP(54).NE.3).AND.(id%MYID.EQ.MASTER) ) 
     &    .OR.  (KEEP(54).EQ.3) ) THEN
C         ----------------------------------------
C         Allocate SIZEOFBLOCKS, DOF2BLOCK
C         ----------------------------------------
          IF (allocated(SIZEOFBLOCKS)) DEALLOCATE(SIZEOFBLOCKS)
          IF (allocated(DOF2BLOCK)) DEALLOCATE(DOF2BLOCK)
          allocate(SIZEOFBLOCKS(NBLK), DOF2BLOCK(id%N), 
     &                  STAT=allocok)
C
          IF (allocok.NE.0) THEN
           id%INFO( 1 ) = -7
           id%INFO( 2 ) = id%N+NBLK
           IF ( LPOK ) WRITE(LP, 150) ' SIZEOFBLOCKS, DOF2BLOCK'
          ENDIF
C
          IF (id%MYID.EQ.MASTER.AND.allocok.EQ.0) THEN
C          BLKPTR and BLKVAR needed for DMUMPS_EXPAND_TREE
C          allocate then if not associated
           IF (.NOT.associated(id%BLKPTR)) THEN
            BLKPTR_ALLOCATED = .TRUE.
            allocate(id%BLKPTR(NBLK+1), STAT=allocok)
            IF (allocok.NE.0) THEN
             BLKPTR_ALLOCATED = .TRUE.
             id%INFO( 1 ) = -7
             id%INFO( 2 ) = NBLK+1
             IF ( LPOK ) WRITE(LP, 150) ' id%BLKPTR '
            ENDIF
           ENDIF
           IF (.NOT.associated(id%BLKVAR).AND.allocok.EQ.0) THEN
            allocate(id%BLKVAR(id%N), STAT=allocok)
            BLKVAR_ALLOCATED = .TRUE.
            IF (allocok.NE.0) THEN
             BLKVAR_ALLOCATED = .FALSE.
             id%INFO( 1 ) = -7
             id%INFO( 2 ) = id%N
             IF ( LPOK ) WRITE(LP, 150) ' id%BLKVAR '
            ENDIF
           ENDIF
          ENDIF
       ENDIF
       CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
       IF (INFO(1).LT.0) GOTO 500
       IF ( id%MYID .EQ. MASTER ) THEN
C      -----------------------------------------
C      Compute SIZEOFBLOCKS, DOF2BLOCK on MASTER
C      based on id%BLKPTR and id%BLKVAR 
C      and compute id%BLKPTR and id%BLKVAR if not 
C      provided by user
C      -----------------------------------------
           IF (BLKVAR_ALLOCATED) THEN
C           implicitly id%BLKVAR(I)=I
            DO I=1, id%N
             id%BLKVAR(I)=I
            ENDDO
           ENDIF
           IF (BLKPTR_ALLOCATED) THEN
            IB=0
            BLKSIZE=-KEEP(13)
            DO I=1, id%N, BLKSIZE
              IB=IB+1
              id%BLKPTR(IB) = I
            ENDDO
            id%BLKPTR(NBLK+1) = id%N+1
           ENDIF
C
           CALL MUMPS_AB_COMPUTE_SIZEOFBLOCK (
     &          NBLK, id%N, id%BLKPTR(1), id%BLKVAR(1),
     &          SIZEOFBLOCKS, DOF2BLOCK)
       ENDIF
C      =======================
       IF (KEEP(54).NE.3) THEN
C      =======================
C      ---------------------
C      Matrix structure available on host
C      ---------------------
        KEEP(14) = 0 
        IF (id%MYID.EQ.MASTER) THEN
C         Store input matrix (IRN/JCN) as a cleaned blocked Lmatrix 
C         of nodes (indices \in [1,NBLK])
          IF (id%KEEP8(28) .EQ. 0_8) THEN
            IRN_PTR => IDUMMY_ARRAY
            JCN_PTR => IDUMMY_ARRAY
          ELSE
            IRN_PTR => id%IRN
            JCN_PTR => id%JCN
          ENDIF
          CALL MUMPS_AB_COORD_TO_LMAT ( id%MYID,
     &          NBLK, id%N, id%KEEP8(28), IRN_PTR(1), JCN_PTR(1),
     &          DOF2BLOCK, 
     &          INFO(1), INFO(2), LP, LPOK, 
     &          LMAT_BLOCK )
        ENDIF
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
C
        IF (id%MYID.EQ.MASTER) THEN
C         From LMAT_BLOCK build GCOMP format wich requires
C         symmetrizing the Lmatrix
          CALL MUMPS_AB_LMAT_TO_CLEAN_G ( id%MYID, .TRUE.,
     &         .TRUE., ! not relevant because unfold is true
     &         LMAT_BLOCK, GCOMP, 
     &         INFO(1), ICNTL(1))
          GCOMP_PROVIDED = .TRUE.
          IF (KEEP(494).EQ.0) THEN
            CALL MUMPS_AB_FREE_LMAT(LMAT_BLOCK)
          ENDIF
        ENDIF
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
        IF ( INFO(1) .LT. 0 ) GOTO 500
C      ====
       ELSE
C      ====
C      -------------------------------
C      Matrix structure is distributed
C      and since KEEP(13).NE.0 then
C      ordering is centralized since
C      -------------------------------
C 
         IF (.NOT. I_AM_SLAVE .OR. ! non-working master
     &       id%KEEP8(29) .EQ. 0_8) THEN ! NNZ_loc or NZ_loc
C          Master non-working
           IRN_loc_PTR => IDUMMY_ARRAY
           JCN_loc_PTR => IDUMMY_ARRAY
           id%KEEP8(29) = 0_8
         ELSE
           IRN_loc_PTR => id%IRN_loc
           JCN_loc_PTR => id%JCN_loc
         ENDIF
C
C        Given distributed matrix  IRN_loc_PTR, JCN_loc_PTR
C        build distributed cleaned graph GCOMP and
C        save  distributed LUMAT in case of grouping
C
         IF (id%NPROCS.EQ.1) THEN
C          Centralized cleaned graph is ready
C          call directly with GCOMP
           READY_FOR_ANA_F = .TRUE.
           CALL MUMPS_AB_DCOORD_TO_DCOMPG ( 
     &          id%MYID, id%NPROCS, id%COMM, 
     &          NBLK, id%N,
     &          id%KEEP8(29),  ! => NNZ_loc or NZ_loc
     &          IRN_loc_PTR(1), JCN_loc_PTR(1), 
     &          DOF2BLOCK(1), 
     &          id%ICNTL(1), id%INFO(1), id%KEEP(1),
     &          LUMAT, GCOMP, READY_FOR_ANA_F)
           GCOMP_PROVIDED = .TRUE.
         ELSE
           READY_FOR_ANA_F = .FALSE.
           CALL MUMPS_AB_DCOORD_TO_DCOMPG ( 
     &          id%MYID, id%NPROCS, id%COMM, 
     &          NBLK, id%N,
     &          id%KEEP8(29),  ! => NNZ_loc or NZ_loc
     &          IRN_loc_PTR(1), JCN_loc_PTR(1), 
     &          DOF2BLOCK(1), 
     &          id%ICNTL(1), id%INFO(1), id%KEEP(1),
     &          LUMAT, GCOMP_DIST, READY_FOR_ANA_F)
         ENDIF
C
C
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
          IF ( INFO(1) .LT. 0 ) GOTO 500
C      =====
       ENDIF
C      =====
       IF (allocated(DOF2BLOCK)) THEN
C        DOF2BLOCK reused on master if pivot order given by user
         IF ( (id%MYID.EQ.MASTER).AND. (KEEP(256) .NE. 1)) THEN
           DEALLOCATE(DOF2BLOCK)
         ENDIF
       ENDIF
C     ========================
      ENDIF
C     } END preparation ANA_BLK
C     =========================
C     ====================================================
C     TEST FOR SEQUENTIAL OR PARALLEL ANALYSIS (KEEP(244))
C     ====================================================
      IF ( (KEEP(244).EQ.1) .AND. (KEEP(54) .eq. 3) ) THEN
C     -----------------------------------------------
C     Sequential analysis: 
C     Collect on the host -- if matrix is distributed
C     at analysis -- all integer information needed
C     to perform ordering
C     -----------------------------------------------
          IF (KEEP(13).NE.0) THEN
           IF (id%NPROCS.NE.1) THEN
             CALL MUMPS_AB_GATHER_GRAPH(
     &       id%ICNTL(1), KEEP(1), id%COMM, id%MYID, id%NPROCS,
     &       id%INFO(1), 
     &       GCOMP_DIST, GCOMP)
             GCOMP_PROVIDED = .TRUE.
C
             CALL MUMPS_AB_FREE_GCOMP(GCOMP_DIST)
           ENDIF
          ELSE
            CALL DMUMPS_GATHER_MATRIX(id)
            CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &            id%COMM, id%MYID )
          ENDIF
          IF ( INFO(1) .LT. 0 ) GOTO 500
      ENDIF
 1234 CONTINUE
      IF (KEEP(244) .EQ. 1) THEN
C     Sequential analysis : Schur
        IF ( id%MYID .eq. MASTER ) THEN
C           Prepare arguments for call to DMUMPS_ANA_F and
C           DMUMPS_ANA_F_ELT in case id%SCHUR was not allocated
C           by user. The objective is to avoid passing a null
C           pointer. 
C FIXME  Block fomat for Schur
            IF ( .NOT. associated( id%LISTVAR_SCHUR ) ) THEN
               SIZE_SCHUR_PASSED = 1
               LISTVAR_SCHUR_2BE_FREED=.TRUE.
               allocate( id%LISTVAR_SCHUR( 1 ), STAT=allocok )
               IF ( allocok .GT. 0 ) THEN
                  WRITE(*,*)
     &                 'PB allocating an array of size 1 for Schur!! '
                  INFO(1)=-7
                  INFO(2)=1
               END IF
            ELSE
               SIZE_SCHUR_PASSED=id%SIZE_SCHUR
               LISTVAR_SCHUR_2BE_FREED = .FALSE.
            END IF
        ENDIF
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &            id%COMM, id%MYID )
        IF ( INFO(1) < 0 ) GOTO 500
      ENDIF
C
      IF ((id%MYID.EQ.MASTER).AND.(KEEP(244) .EQ. 1) 
     &     .AND. (id%N.EQ.NBLK) 
     &   ) THEN
C     Sequential analysis : maximum transversal on master
         IF  ((KEEP(50).NE.1).AND.
     &        .NOT.((KEEP(23).EQ.7).AND.KEEP(50).EQ.0) 
     &       ) THEN
C             (KEEP(23).EQ.7).AND.KEEP(50).EQ.0) :
C             For unsymmetric matrix, if automatic setting is requested
C             default setting of Maximum Transversal is decided during
C             DMUMPS_ANA_F and is based on matrix unsymmetry. 
C             Thus in this case we skip DMUMPS_ANA_O
              IF ( ( KEEP(23) .NE. 0 )  .OR.
C                Automatic choice for scaling does not force Maxtrans
C                Only when scaling is explicitly asked during analysis
C                (KEEP(52)=-2) DMUMPS_ANA_O is called
     &           KEEP(52) .EQ. -2 ) THEN
C
C                Maximum Trans. algorithm called on original matrix.
C                We compute a permutation of the original matrix to
C                have a zero free diagonal
C                KEEP(23)=7 means that automatic choice
C                      of max trans value will be done during analysis
C                Permutation is held in UNS_PERM(1, ...,N).  
C                Maximum transversal is not available for element 
C                entry format
C              UNS_PERM that might be set to
C              to permutation computed during Max transversal
               ALLOCATE(id%UNS_PERM(id%N),IKEEPALLOC(3*id%N),
     &                  WORK2ALLOC(id%N), stat=IERR)
               IF (IERR.GT.0) THEN
                INFO(1)=-7
                INFO(2)=5*id%N
               ELSE
                CALL DMUMPS_ANA_O(id%N, id%KEEP8(28), KEEP(23),
     &              id%UNS_PERM, IKEEPALLOC, id%IRN, id%JCN, id%A,
     &              id%ROWSCA, id%COLSCA,
     &              WORK2ALLOC, id%KEEP, id%ICNTL, id%INFO, id%INFOG)
                IF (allocated(WORK2ALLOC)) DEALLOCATE(WORK2ALLOC)
                IF (KEEP(23).EQ.0) THEN
C                 Maximum tranversal did not produce a permutation
                  IF (associated( id%UNS_PERM )) 
     &                    DEALLOCATE(id%UNS_PERM)
                  NULLIFY(id%UNS_PERM)
                ENDIF
C               Check if IKEEPALLOC needed for ANA_F
                IF (KEEP(23).EQ.0.AND.(KEEP(95).EQ.1)) THEN
                  IF (allocated(IKEEPALLOC)) DEALLOCATE(IKEEPALLOC)
                ENDIF
               ENDIF
               IF (INFO(1) .LT. 0) THEN
C              Fatal error
C              Permutation was not computed; reset keep(23)
                  KEEP(23) = 0
               ELSE
               ENDIF
              ELSE
                KEEP(23)    = 0
C               Switch off
C               compressed/contrained ordering
                id%KEEP(95) = 1
              END IF
         ENDIF
C     END OF MAX-TRANS ON THE MASTER
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) GOTO 500
C
      IF ( KEEP(244) .EQ. 1) THEN
C        Sequential analysis: allocate data for ordering on MASTER
         IF (id%MYID.EQ.MASTER) THEN
C          allocate IKEEPALLOC and TREE related pointers
C          IKEEPALLOC might have been allocated in DMUMPS_ANA_O
C          and IKEEPALLOC(1:N) might hold information to 
C          be given to ANA_F.
           IF (allocated(IKEEPALLOC)) THEN
             ALLOCATE( FILSPTR(NBLK), FREREPTR(NBLK),
     &               NFSIZPTR(NBLK), stat=IERR)
             IF (IERR.GT.0) THEN
              INFO(1)=-7
              INFO(2)=3*NBLK
             ENDIF
           ELSE
             ALLOCATE(IKEEPALLOC(NBLK+2*id%N), 
     &               FILSPTR(NBLK), FREREPTR(NBLK),
     &               NFSIZPTR(NBLK), stat=IERR)
            IF (IERR.GT.0) THEN
             INFO(1)=-7
             INFO(2)=4*NBLK+2*id%N
            ENDIF
           ENDIF
         ENDIF
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
         IF ( INFO(1) < 0 ) GOTO 500
      ENDIF
C
      IF (KEEP(244) .EQ. 1) THEN
C     Sequential analysis
        IF ( id%MYID .eq. MASTER ) THEN
C       BEGINNING OF ANALYSIS ON THE MASTER
C       ------------------------------------------------------
C       For element entry (KEEP(55).ne.0), we do not know NZ, 
C       and so the whole allocation of IW cannot be done at this 
C       point and more workspace is declared/allocated/used
C       inside DMUMPS_ANA_F_ELT.
C       ------------------------------------------------------
C
            IF (KEEP(55) .EQ. 0) THEN
C              ----------------
C              Assembled format
C              ----------------
               NZ8=id%KEEP8(28)
C              Compute LIW8:
C              For local orderings a contiguous space IW 
C              of size LIW8 must be provided. 
C              IW must hold the graph (with double adjacency 
C              list) and and extra space of size the number of 
C              nodes in the graph: 
C                ==>   LIW8 = 2_8 * NZ8 +  int(NBLK,8) + 1_8
C              In case of analysis by block and 
C              However, when GCOMP is provided directly then 
C              IW is not allocated 
C                ==>   LIW8 = 0
C              In this case 
C                size(LCOMP%ADJ)>= 2_8*NZ8+int(NBLK,8)+1_8 
C                should hold
               IF (KEEP(13).NE.0) THEN
C                Compact graph is provided on entry to DMUMPS_ANA_F
                 NZ8=0_8 ! GCOMP is provided on entry
               ENDIF
               IF (NZ8.EQ.0_8) THEN 
                 LIW8 = 0_8
               ELSE
                 LIW8 = 2_8 * NZ8 +  int(NBLK,8) + 1_8
               ENDIF
C        
            ELSE
C              ----------------
C              Elemental format
C              ----------------
C              Only available for AMD, METIS, and given ordering
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
               COND = (KEEP(60) .NE. 0) .OR. (KEEP(256) .EQ. 5)
#else
               COND = (KEEP(60) .NE. 0)
#endif
               IF( COND ) THEN
C
C
C                 we suppress supervariable detection when Schur
C                 is active or when METIS is applied
C                 Workspaces for FLAG(N), and either LEN(N) or some pointers(N+1)
                  LIW_ELT = id%N + id%N + 1
               ELSE
C                 Spaces FLAG(N), LEN(N), N+3, SVAR(0:N),
                  LIW_ELT =  id%N + id%N + id%N + 3 + id%N + 1
               ENDIF
C     
            ENDIF
C           We must ensure that an array of order
C           3*N is available for DMUMPS_ANA_LNEW
            IF (KEEP(55) .EQ. 0) THEN
              IF (LIW8.LT.3_8*int(NBLK,8)) LIW8 = 3_8*int(NBLK,8)
            ELSE
              IF (LIW_ELT.LT.3*id%N) LIW_ELT = 3*id%N
            ENDIF
C
            IF ( KEEP(256) .EQ. 1 ) THEN
C            It has been checked that id%PERM_IN is associated but
C            values of pivot order will be checked later and
C            should be checked here too
C            PERM_IN( I ) = position of I in the pivot order
             IKEEP2 => IKEEPALLOC(NBLK+1:NBLK+id%N)
C            Build inverse permutation and check PERM_IN 
             DO I = 1, id%N
                IKEEP2(I) = 0
             ENDDO
             DO I = 1, id%N
                IF ( id%PERM_IN(I) .LT.1 .OR.
     &               id%PERM_IN(I) .GT. id%N ) THEN
C                 PERM_IN entry is out-of-range
                  INFO(1) = -4
                  INFO(2) = I
                  GOTO 10
                ELSE IF ( IKEEP2(id%PERM_IN(I)) .NE. 0 ) THEN
C                 Duplicate entry in PERM_IN was found
                  INFO(1) = -4
                  INFO(2) = I
                  GOTO 10
                ELSE
C                 Store entry in inverse permutation
                  IKEEP2(id%PERM_IN( I )) = I
                ENDIF
             ENDDO
             IF ((KEEP(55) .EQ. 0).AND.(KEEP(13).NE.0)
     &           .AND.(KEEP(13).NE.-1)
     &         ) THEN
C             Build blocked permutation:
C              IKEEPALLOC(IB)= IBPos where IB, IBPos \in [1:NBLK]
C             IKEEP2 holds inverse permutation
              IPOSB = 0
              IPOS  = 1
              DO WHILE (IPOS.LE.id%N)
                IPOSB     = IPOSB+1
                I         = IKEEP2(IPOS)
                IBcurrent = DOF2BLOCK(I)
                BLKSIZE   = SIZEOFBLOCKS(IBcurrent)
                IKEEPALLOC(IBcurrent) = IPOSB
                IF (BLKSIZE.GT.1) THEN
                 DO II = 1, BLKSIZE-1
                  IPOS    = IPOS+1
                  I       = IKEEP2(IPOS)
                  IB      = DOF2BLOCK(I)
                  IF (IB.NE.IBcurrent) THEN
                   INFO(1)= -4
                   INFO(2)= I
                   GOTO 10
                  ENDIF
                 ENDDO
                ENDIF
                IPOS = IPOS+1
              ENDDO
C             IF PERM_IN is correct then
C             on exit last position should be NBLK
              IF (IPOSB.NE.NBLK) THEN
                   INFO(1)= -4
C                  N+1 to indicate "global" error
                   INFO(2)= id%N+1
                   GOTO 10
              ENDIF
             ELSE
               DO I = 1, id%N
                  IKEEPALLOC( I ) = id%PERM_IN( I )
               END DO
             ENDIF
             IF (allocated(DOF2BLOCK)) DEALLOCATE(DOF2BLOCK)
            END IF
            INFOG(1) = 0
            INFOG(2) = 0
C           Initialize structural symmetry value to not yet computed.
            INFOG(8) = -1
            IF (KEEP(55) .EQ. 0) THEN
               IKEEP1 => IKEEPALLOC(1:NBLK)
               IKEEP2 => IKEEPALLOC(NBLK+1:NBLK+id%N)
               IKEEP3 => IKEEPALLOC(NBLK+id%N+1:NBLK+2*id%N)
C              id%UNS_PERM corresponds to argument PIV
C              in DMUMPS_ANA_F, it should be an assumed-shape
C              array rather than a possibly null pointer:
               IF (associated(id%UNS_PERM)) THEN
                 UNS_PERM_PTR => id%UNS_PERM
               ELSE
                 UNS_PERM_PTR => IDUMMY_ARRAY
               ENDIF
               IF (KEEP(13).EQ.0) THEN
                CALL DMUMPS_ANA_F(id%N, NZ8,
     &              id%IRN, id%JCN,
     &              LIW8, IKEEP1, IKEEP2, IKEEP3,
     &              KEEP(256), NFSIZPTR,
     &              FILSPTR, FREREPTR,
     &              id%LISTVAR_SCHUR, SIZE_SCHUR_PASSED,
     &              id%ICNTL, id%INFOG, id%KEEP,id%KEEP8,id%NSLAVES, 
     &              UNS_PERM_PTR,
     &              id%CNTL(4), id%COLSCA,  id%ROWSCA
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)         
     &              , id%METIS_OPTIONS(1)
#endif               
     &          )               
               ELSE
                IRN_loc_PTR => IDUMMY_ARRAY
                JCN_loc_PTR => IDUMMY_ARRAY
                CALL DMUMPS_ANA_F(NBLK, NZ8,
     &              IRN_loc_PTR, JCN_loc_PTR,
     &              LIW8, IKEEP1, IKEEP2, IKEEP3,
     &              KEEP(256), NFSIZPTR,
     &              FILSPTR, FREREPTR,
     &              id%LISTVAR_SCHUR, SIZE_SCHUR_PASSED,
     &              id%ICNTL, id%INFOG, id%KEEP,id%KEEP8,id%NSLAVES, 
     &              UNS_PERM_PTR,
     &              id%CNTL(4), id%COLSCA,  id%ROWSCA
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)         
     &              , id%METIS_OPTIONS(1)
#endif               
     &              , id%N, SIZEOFBLOCKS, GCOMP_PROVIDED, GCOMP
     &          )              
                IF (GCOMP_PROVIDED) CALL MUMPS_AB_FREE_GCOMP(GCOMP)
C
               ENDIF
               INFOG(7)     = KEEP(256)
C              UNS_PERM_PTR was only used locally
C              for the call to DMUMPS_ANA_F
               NULLIFY(UNS_PERM_PTR)
            ELSE
               allocate( XNODEL ( id%N+1 ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = id%N + 1
                  IF ( LPOK ) THEN
                     WRITE(LP, 150) 'XNODEL'
                  END IF
                  GOTO 10
               ENDIF
               IF (LELTVAR.ne.id%ELTPTR(NELT+1)-1)  THEN
C     -- internal error
                  INFO(1) = -2002
                  INFO(2) = id%ELTPTR(NELT+1)-1
                  GOTO 10
               ENDIF
               allocate( NODEL ( LELTVAR ), stat = IERR )
               IF ( IERR .GT. 0 ) THEN
                  INFO( 1 ) = -7
                  INFO( 2 ) = LELTVAR
                  IF ( LPOK ) THEN
                     WRITE(LP, 150) 'NODEL'
                  END IF
                  GOTO 10
               ENDIF
               CALL DMUMPS_ANA_F_ELT(id%N, NELT,
     &              id%ELTPTR(1), id%ELTVAR(1), LIW_ELT,
     &              IKEEPALLOC(1),
     &              KEEP(256), NFSIZPTR(1), FILSPTR(1),
     &              FREREPTR(1), id%LISTVAR_SCHUR(1),
     &              SIZE_SCHUR_PASSED,
     &              ICNTL(1), INFOG(1), KEEP(1),KEEP8(1),
     &              id%NSLAVES,
     &              XNODEL(1), NODEL(1)
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)         
     &              , id%METIS_OPTIONS(1)
#endif               
     &          )
               INFOG(7)=KEEP(256)
C     
C              XNODEL and NODEL as output to DMUMPS_ANA_F_ELT 
C              be used in DMUMPS_FRTELT and thus 
C              cannot be deallocated at this point
C     
            ENDIF
            IF ( LISTVAR_SCHUR_2BE_FREED ) THEN
C              We do not want to have LISTVAR_SCHUR
C              allocated of size 1 if Schur is off. 
               DEALLOCATE( id%LISTVAR_SCHUR )
               NULLIFY   ( id%LISTVAR_SCHUR )
               LISTVAR_SCHUR_2BE_FREED = .TRUE.
            ENDIF
C           ------------------------------
C           Significant error codes should
C           always be in INFO(1/2)
C           ------------------------------
            INFO(1)=INFOG(1)
            INFO(2)=INFOG(2)
C           save statistics in KEEP array.
            KEEP(28) = INFOG(6)
            IKEEP = 1
            NA      = IKEEP +     id%N
            NE      = IKEEP + 2 * id%N
C       -- if (id%myid.eq.master)
        ENDIF
C       -- if sequential analysis
      ENDIF
C
  10  CONTINUE
      IF (KEEP(244).EQ.1) THEN
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
         IF ( INFO(1) < 0 ) GOTO 500
      ENDIF
      IF ((KEEP(244).EQ.1).AND.(KEEP(55).EQ.0)) THEN
C     Sequential analysis on assembled matrix
C     check if  max transversal should be called
        CALL MPI_BCAST(KEEP(23),1,MPI_INTEGER,MASTER,id%COMM,IERR)
        IF ( (KEEP(23).LE.-1).AND.(KEEP(23).GE.-6) ) THEN
C          -- Perform max transversal
           KEEP(23) = -KEEP(23)
           IF (id%MYID.EQ.MASTER) THEN
             IF (.NOT. associated(id%A)) KEEP(23) = 1
             IF (associated(id%UNS_PERM)) DEALLOCATE(id%UNS_PERM)
             NULLIFY(id%UNS_PERM)
             IF (allocated(IKEEPALLOC)) DEALLOCATE(IKEEPALLOC)
             IF (associated(FILSPTR) ) THEN
               DEALLOCATE(FILSPTR)
               NULLIFY(FILSPTR)
             ENDIF
             IF (associated(FREREPTR) ) THEN
               DEALLOCATE(FREREPTR)
               NULLIFY(FREREPTR)
             ENDIF
             IF (associated(NFSIZPTR) ) THEN
               DEALLOCATE(NFSIZPTR)
               NULLIFY(NFSIZPTR)
             ENDIF
           ENDIF
        GOTO 1234
        ENDIF
      ENDIF
      IF (id%MYID.EQ.MASTER) THEN
        IF ((KEEP(244).EQ.1).AND. (KEEP(55).EQ.0))  THEN
C              Sequential ordering on assembled matrix
               IF ((KEEP(54).EQ.3).AND.KEEP(494).EQ.0) THEN
                IF (associated(id%IRN)) THEN
                  DEALLOCATE(id%IRN)
                  NULLIFY(id%IRN)
                ENDIF
                IF (associated(id%JCN)) THEN
                  DEALLOCATE(id%JCN)
                  NULLIFY(id%JCN)
                ENDIF
               ENDIF
        ENDIF
      ENDIF
      IF (KEEP(244).NE.1) THEN
C     Parallel analysis
         IKEEP   = 1
         NA      = IKEEP +     id%N
         NE      = IKEEP + 2 * id%N
         IF (id%MYID .EQ. MASTER) THEN
            ALLOCATE( IKEEPALLOC(3*id%N), WORK2ALLOC(4*id%N), 
     &                FILSPTR(id%N), FREREPTR(id%N), NFSIZPTR(id%N),
     &      stat=IERR)
         ELSE
C           Because our purpose is to minimize the peak memory consumption,
C           we can afford to allocate on processes other than host
            ALLOCATE(IKEEPALLOC(3*id%N),WORK2ALLOC(4*id%N), stat=IERR )
         ENDIF
         IF (IERR.GT.0) THEN
           INFO(1) = -7
           IF (id%MYID .EQ. MASTER) THEN
             INFO( 2 ) = 10*id%N 
           ELSE
             INFO( 2 ) = 7*id%N 
          ENDIF
         ENDIF
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
         IF ( INFO(1) < 0 ) GOTO 500
         CALL DMUMPS_ANA_F_PAR(id,
     &        IKEEPALLOC,
     &        WORK2ALLOC,
     &        NFSIZPTR,
     &        FILSPTR,
     &        FREREPTR)
         DEALLOCATE(WORK2ALLOC) 
         IF(id%MYID .NE. MASTER) THEN
           DEALLOCATE(IKEEPALLOC)
         ENDIF
         KEEP(28) = INFOG(6)
      END IF
C     Allocated PROCNODE on MASTER
      IF (id%MYID.EQ.MASTER) THEN
       allocok = 0
       allocate(PROCNODE(NBLK), STAT=allocok)
       IF (allocok .ne. 0) THEN
            INFO(1) = -7
            INFO(2) = NBLK
       ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) GOTO 500
      IF(id%MYID .EQ. MASTER) THEN
C        Save ICNTL(14) value into KEEP(12)
         CALL MUMPS_GET_PERLU(KEEP(12),ICNTL(14),
     &        KEEP(50),KEEP(54),ICNTL(6),KEEP(52))
         CALL DMUMPS_ANA_R(NBLK, FILSPTR(1), FREREPTR(1),
     &        IKEEPALLOC(NE), IKEEPALLOC(NA))
C      **********************************************************
C      Continue with CALL to MAPPING routine
C        *********************
C        BEGIN SEQUENTIAL CODE
C        No mapping computed
C        *********************
C
C        In sequential, if no special root
C        reset KEEP(20) and KEEP(38) to 0
C
         IF (id%NSLAVES .EQ. 1
     &      ) THEN
            id%NBSA = 0
            IF ( (id%KEEP(60).EQ.0).
     &           AND.(id%KEEP(53).EQ.0))  THEN 
C     If Schur is on (keep(60).ne.0)
C     or if RR is on (keep (53) > 0 
C     then we keep root numbers
C              root node number in seq  
               id%KEEP(20)=0    
C              root node number in paral  
               id%KEEP(38)=0   
            ENDIF
C     No type 2 nodes:
            id%KEEP(56)=0
C     All mapped on MPI process 0, and of type TPN=0
C     (treated as if they were all root of subtree)
            PROCNODE_VALUE = MUMPS_ENCODE_TPN_IPROC(0, 0, KEEP(199))
            DO I = 1, NBLK
              PROCNODE(I) = PROCNODE_VALUE
            END DO
C     It may also happen that KEEP(38) has already been set,
C     in the case of a distributed Schur complement (KEEP(60)=2 or 3).
C     In that case, PROCNODE should be set accordingly and KEEP(38) is
C     not modified.
            IF (id%KEEP(60) .EQ. 2 .OR. id%KEEP(60).EQ.3) THEN
               PROCNODE_VALUE = MUMPS_ENCODE_TPN_IPROC(3, 0, KEEP(199))
               CALL DMUMPS_SET_PROCNODE(id%KEEP(38), PROCNODE(1),
     &              PROCNODE_VALUE, FILSPTR(1), NBLK)
            ENDIF
C        *******************
C        END SEQUENTIAL CODE
C        *******************
         ELSE
C        *****************************
C        BEGIN MAPPING WITH CANDIDATES
C        (NSLAVES > 1)
C        *****************************
C     
C     
C      peak is set by default to 1 largest front + One largest CB
       PEAK = dble(id%INFOG(5))*dble(id%INFOG(5)) + ! front matrix
     &        dble(id%KEEP(2))*dble(id%KEEP(2))     ! cb bloc
C     IKEEP(1:N,1) can be used as a work space since it is set
C     to its final state by the SORT_PERM subroutine below.
            SSARBR => IKEEPALLOC(IKEEP:IKEEP+NBLK-1)
C     ======================================================
C     Map nodes and assign candidates for dynamic scheduling
C     ======================================================
      IF ((KEEP(13).NE.0).AND.(NBLK.NE.id%N)) THEN
       SIZEOFBLOCKS_PTR => SIZEOFBLOCKS(1:NBLK)
       LSIZEOFBLOCKS_PTR = NBLK
      ELSE
       SIZEOFBLOCKS_PTR => IDUMMY_ARRAY
       LSIZEOFBLOCKS_PTR = 1
       IDUMMY_ARRAY(1) = -1
      ENDIF
            CALL DMUMPS_DIST_AVOID_COPIES(
     &           NBLK,id%NSLAVES,ICNTL(1),
     &           INFOG(1),
     &           IKEEPALLOC(NE),
     &           NFSIZPTR(1),
     &           FREREPTR(1),
     &           FILSPTR(1),
     &           KEEP(1),KEEP8(1),PROCNODE(1),
     &           SSARBR(1),id%NBSA,PEAK,IERR
     &           , SIZEOFBLOCKS_PTR(1), LSIZEOFBLOCKS_PTR 
     &           )
            NULLIFY(SSARBR)
            if(IERR.eq.-999) then 
               write(6,*) ' Internal error during static mapping '
               INFO(1) = IERR
               GOTO 11
            ENDIF
            IF(IERR.NE.0) THEN 
               INFO(1) = -135
               INFO(2) = IERR
               GOTO 11
            ENDIF
            CALL DMUMPS_ANA_R(NBLK, FILSPTR(1),
     &           FREREPTR(1), IKEEPALLOC(NE),
     &           IKEEPALLOC(NA))
         ENDIF
 11      CONTINUE
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) GOTO 500
C     The following part is done in parallel
      CALL MPI_BCAST( id%NELT, 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      IF (KEEP(55) .EQ. 0) THEN
C     Assembled matrix format. Fill up the id%PTRAR array
C     Broadcast id%SYM_PERM needed to fill up id%PTRAR
C     postpone to after computation  of id%SYM_PERM 
C     computed after id%DAD_STEPS
         if (associated(id%FRTPTR)) DEALLOCATE(id%FRTPTR)
         if (associated(id%FRTELT)) DEALLOCATE(id%FRTELT)
         allocate( id%FRTPTR(1), id%FRTELT(1) ,STAT=allocok)
         IF (allocok .GT. 0) THEN
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'FRTPTR,FRTELT'
            END IF
            INFO(1)= -7
            INFO(2)= 2
         END IF
      ELSE
C     Element Entry: 
C     -------------------------------
C     COMPUTE THE LIST OF ELEMENTS THAT WILL BE ASSEMBLED
C     AT EACH NODE OF THE ELIMINATION TREE. ALSO COMPUTE
C     FOR EACH ELEMENT THE TREE NODE TO WHICH IT IS ASSIGNED.
C     
C     FRTPTR is an INTEGER array of length N+1 which need not be set by
C     the user. On output, FRTPTR(I) points in FRTELT to first element 
C     in the list of elements assigned to node I in the elimination tree.
C     
C     FRTELT is an INTEGER array of length NELT which need not be set by
C     the user. On output, positions FRTELT(FRTPTR(I)) to
C     FRTELT(FRTPTR(I+1)-1) contain the list of elements assigned to 
C     node I in the elimination tree.
C     
         LPTRAR = id%NELT+id%NELT+2
         CALL MUMPS_I8REALLOC(id%PTRAR, LPTRAR, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%PTRAR (Analysis)', ERRCODE=-7)
         CALL MUMPS_REALLOC(id%FRTPTR, id%N+1, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%FRTPTR (Analysis)', ERRCODE=-7)
         CALL MUMPS_REALLOC(id%FRTELT, id%NELT, id%INFO, LP,
     &        FORCE=.TRUE., STRING='id%FRTELT (Analysis)', ERRCODE=-7)
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
         IF ( INFO(1) < 0 ) GOTO 500
         IF(id%MYID .EQ. MASTER) THEN
C     In the elemental format case, PTRAR&friends are still
C     computed sequentially and then broadcasted
            CALL DMUMPS_FRTELT(
     &           id%N, NELT, id%ELTPTR(NELT+1)-1, FREREPTR(1),
     &           FILSPTR(1),
     &           IKEEPALLOC(NA), IKEEPALLOC(NE), XNODEL, 
     &           NODEL, id%FRTPTR(1), id%FRTELT(1), id%ELTPROC(1))
            DO I=1, id%NELT+1
C              PTRAR declared 64-bit
               id%PTRAR(id%NELT+I+1)=int(id%ELTPTR(I),8)
            ENDDO
            DEALLOCATE(XNODEL)
            DEALLOCATE(NODEL)
         END IF
         CALL MPI_BCAST( id%PTRAR(id%NELT+2), id%NELT+1, MPI_INTEGER8,
     &        MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%FRTPTR(1), id%N+1, MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
         CALL MPI_BCAST( id%FRTELT(1), id%NELT,  MPI_INTEGER,
     &        MASTER, id%COMM, IERR )
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) GOTO 500
C     We switch again to sequential computations on the master node
      IF(id%MYID .EQ. MASTER) THEN
         IF ( INFO( 1 ) .LT. 0 ) GOTO 12
         IF ( KEEP(55) .ne. 0 ) THEN
C        ---------------------------------------
C        Build ELTPROC: correspondance between elements and slave ranks
C        in COMM_NODES with special values -1 (all procs) and -2 and -3
C        (no procs). This is used later to distribute the elements on
C        the processes at the beginning of the factorisation phase
C        ---------------------------------------
            CALL DMUMPS_ELTPROC(NBLK, NELT, id%ELTPROC(1),id%NSLAVES,
     &           PROCNODE(1), id%KEEP(1))
         END IF
         NB_NIV2 = KEEP(56)
         IF ( NB_NIV2.GT.0 ) THEN
C     
            allocate(PAR2_NODES(NB_NIV2),
     &           STAT=allocok)
            IF (allocok .GT.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'PAR2_NODES'
               END IF
               GOTO 12
            END IF
         ENDIF
         IF ((NB_NIV2.GT.0) .AND. (KEEP(24).EQ.0)) THEN
            INIV2 = 0
            DO 777 INODE = 1, NBLK
               IF ( ( FREREPTR(INODE) .NE. NBLK ) .AND.
     &              ( MUMPS_TYPENODE(PROCNODE(INODE),id%KEEP(199))
     &              .eq. 2) ) THEN
                  INIV2 = INIV2 + 1
                  PAR2_NODES(INIV2) = INODE
               END IF
 777        CONTINUE
            IF ( INIV2 .NE. NB_NIV2 ) THEN
               WRITE(*,*) "Internal Error 2 in DMUMPS_ANA_DRIVER",
     &              INIV2, NB_NIV2
               CALL MUMPS_ABORT()
            ENDIF
         ENDIF
         IF ( (KEEP(24) .NE. 0) .AND. (NB_NIV2.GT.0) ) THEN
C           allocate array to store cadidates stategy
C           for each level two nodes
            IF ( associated(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
            allocate( id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     &           stat=allocok)
            if (allocok .gt.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2*(id%NSLAVES+1)
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'CANDIDATES'
               END IF
               GOTO 12
            END IF
            CALL MUMPS_RETURN_CANDIDATES
     &           (PAR2_NODES,id%CANDIDATES,
     &            IERR)
            IF(IERR.NE.0)  THEN
               INFO(1) = -2002
               GOTO 12
            ENDIF
C     deallocation of variables of module mumps_static_mapping
            CALL MUMPS_END_ARCH_CV()
            IF(IERR.NE.0)  THEN
               INFO(1) = -2002
               GOTO 12
            ENDIF
         ELSE
            IF (associated(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
            allocate(id%CANDIDATES(1,1), stat=allocok)
            IF (allocok .NE. 0) THEN
               INFO(1)= -7
               INFO(2)= 1
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'CANDIDATES'
               END IF
               GOTO 12
            ENDIF
         ENDIF
C*******************************************************************
C     ---------------
 12      CONTINUE
C     ---------------
*     
*     ===============================
*     End of analysis phase on master
*     ===============================
*     
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) GOTO 500
C     
C     We now allocate and compute arrays in NSTEPS
C     on the master, as this makes more sense.
C     
C     Broadcast KEEP8(101) to be used in MUMPS_ANA_L0_OMP
      CALL MPI_BCAST( id%KEEP8(101), 1, MPI_INTEGER8, MASTER,
     &     id%COMM, IERR )
C
C     ==============================
C     PREPARE DATA FOR FACTORIZATION
C     ==============================
C     ------------------
      CALL MPI_BCAST( id%KEEP(1), 110, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     We also need to broadcast KEEP8(21) 
      CALL MPI_BCAST( id%KEEP8(21), 1, MPI_INTEGER8, MASTER,
     &     id%COMM, IERR )
C     --------------------------------------------------
C     Broadcast KEEP(205) which is outside the first 110
C     KEEP entries but is needed for factorization.
C     --------------------------------------------------
      CALL MPI_BCAST( id%KEEP(205), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     --------------
C     Broadcast NBSA 
      CALL MPI_BCAST( id%NBSA, 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     -----------------
C     Global MAXFRT (computed in DMUMPS_ANA_M)
C     is needed on all the procs during DMUMPS_ANA_DISTM
C     to evaluate workspace for solve. 
C     We could also recompute it in DMUMPS_ANA_DISTM
      IF (id%MYID==MASTER) KEEP(127)=INFOG(5)
      CALL MPI_BCAST( id%KEEP(127), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     -----------------
C     Global max panel size KEEP(226)
      CALL MPI_BCAST( id%KEEP(226), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     -----------------
      CALL MPI_BCAST( id%KEEP(464), 2, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(471), 2, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(475), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(482), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(487), 2, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     Number of leaves not belonging to L0 KEEP(262)
C              and KEEP(263) : inner or outer sends for blocked facto
      CALL MPI_BCAST( id%KEEP(262), 2, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
C     ----------------------------------------
C     Allocate new workspace on all processors
C     ----------------------------------------
      IF (id%MYID.EQ.MASTER) THEN
C        id%STEP is of size NBLK because it 
C                is computed on compressed graph and then extended
C                and broadcasted on all procs
        CALL MUMPS_REALLOC(id%STEP, NBLK, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%STEP (Analysis)', ERRCODE=-7)
      ELSE
C        id%STEP is of size id%N because it 
C                is received in extended form
        CALL MUMPS_REALLOC(id%STEP, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%STEP (Analysis)', ERRCODE=-7)
      ENDIF
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%PROCNODE_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE.,
     &     STRING='id%PROCNODE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%NE_STEPS, id%KEEP(28), id%INFO, LP, 
     &     FORCE=.TRUE., 
     &     STRING='id%NE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%ND_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE., 
     &     STRING='id%ND_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%FRERE_STEPS, id%KEEP(28), id%INFO, LP,
     &     FORCE=.TRUE., 
     &     STRING='id%FRERE_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
      CALL MUMPS_REALLOC(id%DAD_STEPS, id%KEEP(28), id%INFO, LP, 
     &     FORCE=.TRUE., 
     &     STRING='id%DAD_STEPS (Analysis)', ERRCODE=-7)
      IF(INFO(1).LT.0) GOTO 94
C     id%FILS is allocated before expand tree
      IF (KEEP(55) .EQ. 0) THEN
        LPTRAR = id%N+id%N
        CALL MUMPS_I8REALLOC(id%PTRAR, LPTRAR, id%INFO, LP,
     &       FORCE=.TRUE., STRING='id%PTRAR (Analysis)', ERRCODE=-7)
        IF(INFO(1).LT.0) GOTO 94
      ENDIF
      IF (id%MYID.EQ.MASTER) THEN
       CALL MUMPS_REALLOC(id%LRGROUPS, NBLK, id%INFO, LP, 
     &     FORCE=.TRUE.
     &     ,STRING='id%LRGROUPS (Analysis)', ERRCODE=-7)
      ELSE
       CALL MUMPS_REALLOC(id%LRGROUPS, id%N, id%INFO, LP, 
     &     FORCE=.TRUE.
     &     ,STRING='id%LRGROUPS (Analysis)', ERRCODE=-7)
      ENDIF
      IF(INFO(1).LT.0) GOTO 94
C     Copy data for factorization and/or solve.
C     ================================
C     COMPUTE ON THE MASTER, BROADCAST
C     TO OTHER PROCESSES
C     ================================
      IF ( id%MYID .NE. MASTER .OR. id%KEEP(23) .EQ. 0 ) THEN
       IF ( associated( id%UNS_PERM ) ) THEN
         DEALLOCATE(id%UNS_PERM)
       ENDIF
      ENDIF
 94   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%MYID .EQ. MASTER ) THEN
C     NA -> compressed NA containing only list
C     of leaves of the elimination tree and list of roots 
C     (the two useful informations for factorization/solve).
         IF (NBLK.eq.1) THEN
            NBROOT = 1
            NBLEAF = 1
         ELSE IF (IKEEPALLOC(NA+NBLK-1) .LT.0) THEN
            NBLEAF= NBLK
            NBROOT= NBLK
         ELSE IF (IKEEPALLOC(NA+NBLK-2) .LT.0) THEN
            NBLEAF = NBLK-1
            NBROOT = IKEEPALLOC(NA+NBLK-1)
         ELSE
            NBLEAF = IKEEPALLOC(NA+NBLK-2)
            NBROOT = IKEEPALLOC(NA+NBLK-1)
         ENDIF
         id%LNA = 2+NBLEAF+NBROOT
      ENDIF
      CALL MPI_BCAST( id%LNA, 1, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MUMPS_REALLOC(id%NA, id%LNA, id%INFO, LP, FORCE=.TRUE., 
     &     STRING='id%NA (Analysis)', ERRCODE=-7)
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF (id%MYID .EQ.MASTER ) THEN
C{    The structure of NA is the following:
C        NA(1) is the number of leaves.
C        NA(2) is the number of roots.
C        NA(3:2+NA(1)) are the leaves.
C        NA(3+NA(1):2+NA(1)+NA(2)) are the roots.
         id%NA(1) = NBLEAF
         id%NA(2) = NBROOT
C     
C        Initialize NA with the leaves and roots
         LEAF = 3
         IF ( NBLK == 1 ) THEN
            id%NA(LEAF) = 1
            LEAF = LEAF + 1
         ELSE IF (IKEEPALLOC(NA+NBLK-1) < 0) THEN
            id%NA(LEAF) = - IKEEPALLOC(NA+NBLK-1)-1
            LEAF = LEAF + 1
            DO I = 1, NBLEAF - 1
               id%NA(LEAF) = IKEEPALLOC(NA+I-1)
               LEAF = LEAF + 1
            ENDDO
         ELSE IF (IKEEPALLOC(NA+NBLK-2) < 0 ) THEN
            INODE = - IKEEPALLOC(NA+NBLK-2) - 1
            id%NA(LEAF) = INODE
            LEAF =LEAF + 1
            IF ( NBLEAF > 1 ) THEN
               DO I = 1, NBLEAF - 1
                  id%NA(LEAF) = IKEEPALLOC(NA+I-1)
                  LEAF = LEAF + 1
               ENDDO
            ENDIF
         ELSE
            DO I = 1, NBLEAF
               id%NA(LEAF) = IKEEPALLOC(NA+I-1)
               LEAF = LEAF + 1
            ENDDO
         END IF
C
C        Build array STEP(1:id%N) to hold step numbers in
C        range 1..id%KEEP(28), allowing compression of
C        other arrays from id%N to id%KEEP(28)
C        (the number of nodes/steps in the assembly tree)
         ISTEP = 0
         DO I = 1, NBLK
            IF ( FREREPTR(I) .ne. NBLK + 1 ) THEN
C        New node in the tree.
c        (Set step( inode_n ) = inode_nsteps for principal
C        variables and -inode_nsteps for internal variables
C        of the node)
               ISTEP = ISTEP + 1
               id%STEP(I)=ISTEP
               INN = FILSPTR(I)
               DO WHILE ( INN .GT. 0 )
                  id%STEP(INN) = - ISTEP
                  INN = FILSPTR(INN)
               END DO
               IF (FREREPTR(I) .eq. 0) THEN
C              Keep root nodes list in NA
                  id%NA(LEAF) = I
                  LEAF = LEAF + 1
               ENDIF
            ENDIF
         END DO
         IF ( LEAF - 1 .NE. 2+NBROOT + NBLEAF ) THEN
            WRITE(*,*) 'Internal error 2 in DMUMPS_ANA_DRIVER'
            CALL MUMPS_ABORT()
         ENDIF
         IF ( ISTEP .NE. id%KEEP(28) ) THEN
            write(*,*) 'Internal error 3 in DMUMPS_ANA_DRIVER', 
     &        ISTEP, id%KEEP(28)
            CALL MUMPS_ABORT()
         ENDIF
C        ============
C        SET PROCNODE, FRERE, NE
C        ============
C        copies to NSTEP array should be ok
         DO I = 1, NBLK
            IF (FREREPTR(I) .NE. NBLK+1) THEN
               id%PROCNODE_STEPS(id%STEP(I)) = PROCNODE( I )
               id%FRERE_STEPS(id%STEP(I))    = FREREPTR(I)
               id%NE_STEPS(id%STEP(I))    = IKEEPALLOC(NE+I-1)
               id%ND_STEPS(id%STEP(I))    = NFSIZPTR(I)
            ENDIF
         ENDDO
C        ===============================
C        Algorithm to compute array DAD_STEPS:
C        ----
C        For each node set dad for all of its sons
C        plus, for root nodes set dad to zero.
C     
C        ===============================
         DO I = 1, NBLK
C     -- skip non principal nodes
            IF ( id%STEP(I) .LE. 0) CYCLE
C     -- (I) is a principal node
            IF (FREREPTR(I) .eq. 0) THEN
C     -- I is a root node and has no father
               id%DAD_STEPS(id%STEP(I)) = 0
            ENDIF
C     -- Find first son node (IFS)
            IFS = FILSPTR(I)
            DO WHILE ( IFS .GT. 0 )
               IFS= FILSPTR(IFS)
            END DO
C     -- IFS > 0 if I is not a leave node
C     -- Go through list of brothers of IFS if any
            IFS = -IFS
            DO WHILE (IFS.GT.0) 
C     -- I is not a leave node and has a son node IFS
               id%DAD_STEPS(id%STEP(IFS)) = I
               IFS   = FREREPTR(IFS)
            ENDDO
         END DO
C
C     
C        Following arrays (PROCNODE and IKEEPALLOC) not used anymore 
C        during analysis
         IF (allocated(PROCNODE)) DEALLOCATE(PROCNODE)
         IF (allocated(IKEEPALLOC)) DEALLOCATE(IKEEPALLOC)
         IF (associated(FREREPTR)) DEALLOCATE(FREREPTR)
         NULLIFY(FREREPTR)
         IF (associated(NFSIZPTR)) DEALLOCATE(NFSIZPTR)
         NULLIFY(NFSIZPTR)
      ENDIF 
      IF (KEEP(494).NE.0) THEN
C{
         IF (id%MYID.EQ.MASTER) THEN
            IF (PROKG) THEN
                  CALL MUMPS_SECDEB(TIMEG)
            END IF
         ENDIF
C     =======================================================
C     Compute a grouping of variables for LR approximations.
C     Grouping may be performed on a distributed matrix
C     =======================================================
C     
C     I/ Prepare data before call to grouping
        IF ((KEEP(54).EQ.3).AND.(KEEP(13).NE.0)) THEN
C         Matrix is distributed on entry and compression computed
          IF (KEEP(487).NE.1) CALL MUMPS_ABORT()
          ALLOCATE(MAPCOL(id%KEEP(28)), stat=allocok)
          IF (allocok .ne.0) then
           INFO(1)= -7
           INFO(2)= id%KEEP(28)
          ENDIF
C         Broadcast errors
          CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
          IF ( INFO(1).LT.0 ) GOTO 500
C
          CALL MUMPS_INIALIZE_REDIST_LUMAT ( 
     &     id%INFO, id%ICNTL, id%KEEP, id%COMM, id%MYID, NBLK,
     &     LUMAT, id%PROCNODE_STEPS(1), id%KEEP(28), MAPCOL, 
     &     LUMAT_REMAP, NBRECORDS, id%STEP(1))
C         INFO(1) has been broadcasted already in routine
          IF ( id%INFO(1).LT.0 ) GOTO 500
C
C         -- Redistribute LUMAT into LU_REMAP relying on procnode
          CALL MUMPS_AB_DIST_LMAT_TO_LUMAT ( 
     &      .FALSE., ! do not UNFOLD
     &      .TRUE.,  ! MAPCOL in NSTEPS=> STEP array needed 
     &      id%INFO, id%ICNTL, id%COMM, id%MYID, NBLK, id%NPROCS,
     &      LUMAT, MAPCOL, id%KEEP(28), id%STEP(1), NBLK,
     &      LUMAT_REMAP, NBRECORDS, NSEND8, NLOCAL8
     &     )
          CALL MUMPS_AB_FREE_LMAT(LUMAT)
C         Distribute SIZEOFBLOCKS that was defined only on master
          CALL MPI_BCAST( SIZEOFBLOCKS, NBLK, MPI_INTEGER, MASTER, 
     &     id%COMM, IERR )
C
        ELSE IF ((KEEP(54).NE.3).AND.(KEEP(13).NE.0)
     &         .AND. (KEEP(487).EQ.1) ) THEN
C         Centralized matrix and LMAT_BLOCK available
C         ---> build LUMAT_REMAP on MASTER
          IF (id%MYID.EQ.MASTER) THEN
           CALL MUMPS_AB_LMAT_TO_LUMAT ( 
     &              LMAT_BLOCK, LUMAT_REMAP,
     &              INFO(1), ICNTL(1))
C            --- LMAT_BLOCK not needed anymore
           CALL MUMPS_AB_FREE_LMAT(LMAT_BLOCK)
          ENDIF
C         Broadcast errors
          CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
          IF ( INFO(1).LT.0 ) GOTO 500
C
        ELSE IF ((KEEP(54).EQ.3).AND.(KEEP(13).EQ.0)
     &           .AND. KEEP(487).EQ.1) THEN
C        Matrix is distributed on entry and compression not requested
C        (this will be the case when ICNTL(15).EQ.0 and 
C         // analysis, or Schur, etc...)
C        note that with distributed matrix and centralized ordering
C        compression is forced to limit memory peak)
C        Free centralized matrix before grouping to 
C        limit memory peak
         IF (associated(id%IRN)) THEN
                  DEALLOCATE(id%IRN)
                  NULLIFY(id%IRN)
         ENDIF
         IF (associated(id%JCN)) THEN
                  DEALLOCATE(id%JCN)
                  NULLIFY(id%JCN)
         ENDIF
         IF (.NOT. I_AM_SLAVE .OR. ! non-working master
     &       id%KEEP8(29) .EQ. 0_8) THEN ! NNZ_loc or NZ_loc
C          Master non-working
           IRN_loc_PTR => IDUMMY_ARRAY
           JCN_loc_PTR => IDUMMY_ARRAY
         ELSE
           IRN_loc_PTR => id%IRN_loc
           JCN_loc_PTR => id%JCN_loc
         ENDIF
         ALLOCATE(MAPCOL(id%KEEP(28)), stat=allocok)
         IF (allocok .ne.0) then
           INFO(1)= -7
           INFO(2)= id%KEEP(28)
         ENDIF
C        Broadcast errors
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
         IF ( INFO(1).LT.0 ) GOTO 500
C
C        Build MAPCOL and LUMAT_REMAP mapped according 
C        to MAPCOL (outputs available on all MPI procs).
         CALL MUMPS_AB_DCOORD_TO_DTREE_LUMAT ( 
     &          id%MYID, id%NPROCS, id%COMM, 
     &          NBLK, id%N,
     &          id%KEEP8(29),  ! => NNZ_loc or NZ_loc
     &          IRN_loc_PTR(1), JCN_loc_PTR(1), 
     &          id%PROCNODE_STEPS(1), id%KEEP(28), id%STEP(1),
     &          id%ICNTL(1), id%INFO(1), id%KEEP(1),
     &          MAPCOL, LUMAT_REMAP )
         IF (INFO(1).GE.0) THEN
C         SIZEOFBLOCKS needed on all procs during MPI grouping
          ALLOCATE(SIZEOFBLOCKS(NBLK), stat=allocok)
          IF (allocok .ne.0) then
           INFO(1)= -7
           INFO(2)= NBLK
          ENDIF
          DO I=1, NBLK
            SIZEOFBLOCKS(I) = 1
          ENDDO
         ENDIF
C        Broadcast errors
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
         IF ( INFO(1).LT.0 ) GOTO 500
        ELSE IF ((KEEP(54).EQ.3) .AND. (KEEP(244).EQ.2)
     &       .AND.  (KEEP(487).NE.1) 
     &     ) THEN
C        Grouping preparation on slaves:
C        If the input matrix is distributed and the parallel analysis is
C        chosen, the graph used to be centralized in order to compute the
C        clustering. 
C
            CALL DMUMPS_GATHER_MATRIX(id)
        ENDIF
C     ============
C     ============
C     II/ GROUPING
C     ============
        IF ((KEEP(54).EQ.3).AND.(KEEP(487).EQ.1)) THEN
C         Matrix is distributed on entry and halo of size 1
C         Distributed memory based grouping is used
          IF (id%MYID.NE.MASTER) THEN
            ALLOCATE(FILSPTR(NBLK), stat=IERR)
            IF (IERR.GT.0) THEN
             INFO(1)=-7
             INFO(2)=NBLK
            ENDIF
          ENDIF
C         Broadcast errors
          CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
          IF ( INFO(1).LT.0 ) GOTO 500
C         Distribute SIZEOFBLOCKS that was defined only on master
C          CALL MPI_BCAST( SIZEOFBLOCKS, NBLK, MPI_INTEGER, MASTER, 
C     &     id%COMM, IERR )
          CALL DMUMPS_AB_LR_MPI_GROUPING(NBLK, 
     &           MAPCOL, id%KEEP(28), 
     &           id%KEEP(28), LUMAT_REMAP, FILSPTR, 
     &           id%FRERE_STEPS,
     &           id%DAD_STEPS, id%STEP, id%NA,
     &           id%LNA, id%LRGROUPS, SIZEOFBLOCKS(1), id%KEEP(50),
     &           id%ICNTL(1), id%KEEP(487), id%KEEP(488),
     &           id%KEEP(490), id%KEEP(38), id%KEEP(20), id%KEEP(60), 
     &           id%INFO(1), id%INFO(2),
     &           id%KEEP(264), id%KEEP(265), id%KEEP(482), id%KEEP(472),
     &           id%KEEP(127), id%KEEP(469), id%KEEP(10), 
     &           id%KEEP(54), LPOK, LP, id%COMM, id%MYID, id%NPROCS)
          IF (allocated(MAPCOL)) DEALLOCATE(MAPCOL)
          IF (id%MYID.NE.MASTER) THEN
            DEALLOCATE(FILSPTR)
            NULLIFY(FILSPTR)
          ENDIF
C
        ELSE IF (id%MYID.EQ.MASTER) THEN
            IF ((KEEP(54).NE.3).AND.(KEEP(13).NE.0)
     &         .AND. (KEEP(487).EQ.1) ) THEN
C            Centralized matrix and LMAT_BLOCK available
C            --- build LUMAT
C            -- LR grouping exploiting LUMAT
C            -- centralized => MAPCOL not needed
C            FIXME 5.4: call to DMUMPS_AB_LR_GROUPING "ready" to be 
C                   replaced by call to DMUMPS_AB_LR_MPI_GROUPING 
C
             IDUMMY_ARRAY(1) = -1
             CALL DMUMPS_AB_LR_GROUPING(NBLK, 
     &           IDUMMY_ARRAY, 1,
     &           id%KEEP(28), LUMAT_REMAP, FILSPTR, 
     &           id%FRERE_STEPS,
     &           id%DAD_STEPS, id%STEP, id%NA,
     &           id%LNA, id%LRGROUPS, SIZEOFBLOCKS(1), id%KEEP(50),
     &           id%ICNTL(1), id%KEEP(487), id%KEEP(488),
     &           id%KEEP(490), id%KEEP(38), id%KEEP(20), id%KEEP(60), 
     &           id%INFO(1), id%INFO(2),
     &           id%KEEP(264), id%KEEP(265), id%KEEP(482), id%KEEP(472),
     &           id%KEEP(127), id%KEEP(469), id%KEEP(10), 
     &           id%KEEP(54),
     &           LPOK, LP, id%MYID, id%COMM)
            ELSE
C            grouping based on centralized matrix
             IF (KEEP(469).EQ.0) THEN
               CALL DMUMPS_LR_GROUPING(id%N, id%KEEP8(28), id%KEEP(28),
     &           id%IRN,
     &           id%JCN, FILSPTR, id%FRERE_STEPS,
     &           id%DAD_STEPS, id%NE_STEPS, id%STEP, id%NA,
     &           id%LNA, id%LRGROUPS,
     &           id%KEEP(50),
     &           id%ICNTL(1), id%KEEP(487), id%KEEP(488), 
     &           id%KEEP(490), id%KEEP(38), id%KEEP(20), id%KEEP(60),
     &           id%INFO(1), id%INFO(2),
     &           id%KEEP(264), id%KEEP(265), id%KEEP(482), id%KEEP(472),
     &           id%KEEP(127), id%KEEP(10), 
     &           id%KEEP(54),
     &           LPOK, LP)
             ELSE
               CALL DMUMPS_LR_GROUPING_NEW(id%N, id%KEEP8(28), 
     &           id%KEEP(28), id%IRN,
     &           id%JCN, FILSPTR, id%FRERE_STEPS,
     &           id%DAD_STEPS, id%STEP, id%NA,
     &           id%LNA, id%LRGROUPS, id%KEEP(50),
     &           id%ICNTL(1), id%KEEP(487), id%KEEP(488),
     &           id%KEEP(490), id%KEEP(38), id%KEEP(20), id%KEEP(60), 
     &           id%INFO(1), id%INFO(2),
     &           id%KEEP(264), id%KEEP(265), id%KEEP(482), id%KEEP(472),
     &           id%KEEP(127), id%KEEP(469), id%KEEP(10), 
     &           id%KEEP(54),
     &           LPOK, LP)
             ENDIF
            ENDIF
        ENDIF
C     ============
C     III/ CLEANUP 
C     ============
C       Free LUMAT_REMAP is allocated
        CALL MUMPS_AB_FREE_LMAT(LUMAT_REMAP)
        IF (allocated(MAPCOL)) DEALLOCATE(MAPCOL)
        IF (allocated(SIZEOFBLOCKS)) DEALLOCATE(SIZEOFBLOCKS)
        IF ( (KEEP(54).EQ.3) .AND. (KEEP(244).EQ.2).AND.
     &       (KEEP(487).NE.1) ) THEN
C           Cleanup the irn and jcn arrays filled up by the
C           cmumps_gather_matrix above. It might have been done
C           during grouping
                IF (associated(id%IRN)) THEN
                  DEALLOCATE(id%IRN)
                  NULLIFY(id%IRN)
                ENDIF
                IF (associated(id%JCN)) THEN
                  DEALLOCATE(id%JCN)
                  NULLIFY(id%JCN)
                ENDIF
        END IF
        IF (PROKG) THEN
          CALL MUMPS_SECFIN(TIMEG)
          WRITE(MPG,145) TIMEG
        END IF
C}    Grouping: KEEP(494) .NE. 0
      ENDIF
      IF (id%MYID.NE. MASTER) THEN
        CALL MUMPS_REALLOC(id%FILS, id%N, id%INFO, LP, FORCE=.TRUE.,
     &     STRING='id%FILS (Analysis)', ERRCODE=-7)
        IF(INFO(1).LT.0) GOTO 97
      ENDIF
C
      IF ((id%MYID.EQ.MASTER) .AND.(KEEP(13).NE.0)) THEN
C{      ===========
C       Expand tree  
C       ===========
C       Current tree is relative to the analysis by block.
C       Expand the tree on the master if compression is effective
C       (in all cases, grouping done or not)
          IF (NBLK.LT.id%N.OR.(.NOT.BLKVAR_ALLOCATED)) THEN
C         even if NBLK.EQ.N BLKVAR provided by user might hold
C         a permutation of the variables and this expand_tree_steps 
C         should also be called
C           Expand FILSPTR, id%STEP into id%FILS, STEPPTR
C           and update arrays of size NSTEPS
            ALLOCATE(STEPPTR(id%N), LRGROUPSPTR(id%N), stat=IERR)
            IF (IERR.GT.0) THEN
             INFO(1)=-7
             INFO(2)=id%N
             GOTO 97
            ENDIF
            IF (NB_NIV2.EQ.0) THEN
               IDUMMY_ARRAY(1)   = -9999
               PAR2_NODESPTR     => IDUMMY_ARRAY(1:1)
               SIZE_PAR2_NODESPTR=1
            ELSE
               PAR2_NODESPTR => PAR2_NODES(1:NB_NIV2)
               SIZE_PAR2_NODESPTR=NB_NIV2
            ENDIF
            CALL MUMPS_REALLOC(id%FILS, id%N, id%INFO, LP, 
     &           FORCE=.TRUE.,
     &           STRING='id%FILS (Analysis)', ERRCODE=-7)
            IF(INFO(1).LT.0) GOTO 97
            CALL DMUMPS_EXPAND_TREE_STEPS (id%ICNTL,
     &        id%N, NBLK, id%BLKPTR(1), id%BLKVAR(1), 
     &        FILSPTR(1), id%FILS(1), id%KEEP(28),
     &        id%STEP(1), STEPPTR(1),
     &        PAR2_NODESPTR(1), SIZE_PAR2_NODESPTR,
     &        id%DAD_STEPS(1), id%FRERE_STEPS(1), 
     &        id%NA(1), id%LNA, id%LRGROUPS(1), LRGROUPSPTR(1), 
     &        id%KEEP(20), id%KEEP(38)
     &           ) 
            NULLIFY(PAR2_NODESPTR)
            DEALLOCATE(id%STEP)
            id%STEP=>STEPPTR
            NULLIFY(STEPPTR)
            DEALLOCATE(id%LRGROUPS)
            id%LRGROUPS=>LRGROUPSPTR
            NULLIFY(LRGROUPSPTR)
            DEALLOCATE(FILSPTR)
            NULLIFY(FILSPTR)
          ELSE
           if (associated(id%FILS)) DEALLOCATE(id%FILS)
           id%FILS=>FILSPTR
           NULLIFY(FILSPTR)
          ENDIF
C}
      ENDIF
      IF ((id%N.EQ.NBLK).AND.associated(FILSPTR)) THEN
C          id%FILS has not been initialized
           if (associated(id%FILS)) DEALLOCATE(id%FILS)
           id%FILS=>FILSPTR
           NULLIFY(FILSPTR)
      ENDIF
 97   CONTINUE
       CALL MUMPS_REALLOC(id%SYM_PERM, id%N, id%INFO, LP, 
     &     FORCE=.TRUE.,
     &     STRING='id%SYM_PERM (Analysis)', ERRCODE=-7)
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1), id%COMM, id%MYID )
          IF ( INFO(1) < 0 ) GOTO 500
      IF (id%MYID.EQ.MASTER) THEN
C     =================================================================
C     Reorder the tree using a variant of Liu's algorithm. Note that
C     REORDER_TREE MUST always be called since it sorts NA (the list of
C     leaves) in a valid order in the sense of a depth-first traversal.
C     =================================================================
               CALL DMUMPS_REORDER_TREE(id%N, id%FRERE_STEPS(1),
     &              id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &              id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1), 
     &              id%KEEP(28), .TRUE., id%KEEP(28), id%KEEP(70),
     &              id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(215),
     &              id%KEEP(234), id%KEEP(55), id%KEEP(199),
     &              id%PROCNODE_STEPS(1),id%NSLAVES,PEAK,id%KEEP(90)
     &              )
            IF(id%KEEP(261).EQ.1)THEN
               CALL MUMPS_SORT_STEP(id%N, id%FRERE_STEPS(1),
     &              id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &              id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1), 
     &              id%KEEP(28), .TRUE., id%KEEP(28), id%INFO(1),
     &              id%ICNTL(1),id%PROCNODE_STEPS(1),id%NSLAVES
     &              )
            ENDIF
C     Compute and export some global information on the tree needed by
C     dynamic schedulers during the factorization. The type of
C     information depends on the selected strategy.
         IF ((id%KEEP(76).GE.4).OR.(id%KEEP(76).GE.6).OR.
     &              (id%KEEP(47).EQ.4).OR.((id%KEEP(81).GT.0)
     &              .AND.(id%KEEP(47).GE.2)))THEN
            IS_BUILD_LOAD_MEM_CALLED=.TRUE.
            IF ((id%KEEP(47) .EQ. 4).OR.
     &           (( id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
               IF(id%NSLAVES.GT.1) THEN
C                 NBSA is the total number of subtrees  and
C                 is an upperbound of the local number of
C                 subtrees
                  SIZE_TEMP_MEM = id%NBSA
               ELSE
C                 Only one processor, NA(2) is the number of leaves
                  SIZE_TEMP_MEM = id%NA(2)
               ENDIF
            ELSE
               SIZE_TEMP_MEM = 1
            ENDIF
            IF((id%KEEP(76).EQ.4).OR.(id%KEEP(76).EQ.6))THEN
               SIZE_DEPTH_FIRST=id%KEEP(28)
            ELSE
               SIZE_DEPTH_FIRST=1
            ENDIF
            allocate(TEMP_MEM(SIZE_TEMP_MEM,id%NSLAVES),STAT=allocok) 
            IF (allocok .NE.0) THEN
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_MEM'
               END IF
               GOTO 80 !! FIXME propagate error
            END IF
            allocate(TEMP_LEAF(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80 !! FIXME propagate error
            end if
            allocate(TEMP_SIZE(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_SIZE'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(TEMP_ROOT(SIZE_TEMP_MEM,id%NSLAVES),
     &           stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'TEMP_ROOT'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_TEMP_MEM*id%NSLAVES
               GOTO 80
            end if
            allocate(DEPTH_FIRST(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'DEPTH_FIRST'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            ALLOCATE(DEPTH_FIRST_SEQ(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'DEPTH_FIRST_SEQ'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            ALLOCATE(SBTR_ID(SIZE_DEPTH_FIRST),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'SBTR_ID'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_DEPTH_FIRST
               GOTO 80
            end if
            IF(id%KEEP(76).EQ.5)THEN
C     We reuse the same variable as before
               SIZE_COST_TRAV=id%KEEP(28)
            ELSE
               SIZE_COST_TRAV=1
            ENDIF
            allocate(COST_TRAV_TMP(SIZE_COST_TRAV),stat=allocok) 
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'COST_TRAV_TMP'
               END IF
               INFO(1)= -7
               INFO(2)= SIZE_COST_TRAV
               GOTO 80
            END IF
            IF(id%KEEP(76).EQ.5)THEN
               IF(id%KEEP(70).EQ.0)THEN
                  id%KEEP(70)=5
               ENDIF
               IF(id%KEEP(70).EQ.1)THEN
                  id%KEEP(70)=6
               ENDIF
            ENDIF
            IF(id%KEEP(76).EQ.4)THEN
               IF(id%KEEP(70).EQ.0)THEN
                  id%KEEP(70)=3
               ENDIF
               IF(id%KEEP(70).EQ.1)THEN
                  id%KEEP(70)=4
               ENDIF
            ENDIF
            CALL DMUMPS_BUILD_LOAD_MEM_INFO(id%N, id%FRERE_STEPS(1),
     &           id%STEP(1),id%FILS(1), id%NA(1), id%LNA,
     &           id%NE_STEPS(1), id%ND_STEPS(1), id%DAD_STEPS(1),
     &           id%KEEP(28), .TRUE., id%KEEP(28), id%KEEP(70),
     &           id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(47),
     &           id%KEEP(81),id%KEEP(76),id%KEEP(215),
     &           id%KEEP(234), id%KEEP(55), id%KEEP(199),
     &           id%PROCNODE_STEPS(1),TEMP_MEM,id%NSLAVES,
     &           SIZE_TEMP_MEM, PEAK,id%KEEP(90),SIZE_DEPTH_FIRST,
     &           SIZE_COST_TRAV,DEPTH_FIRST(1),DEPTH_FIRST_SEQ(1),
     &           COST_TRAV_TMP(1),
     &           TEMP_LEAF,TEMP_SIZE,TEMP_ROOT,SBTR_ID(1)
     &              )
         END IF
      ENDIF
      IF (id%MYID.EQ.MASTER) THEN
         CALL DMUMPS_SORT_PERM(id%N, id%NA(1), id%LNA,
     &        id%NE_STEPS(1), id%SYM_PERM(1),
     &        id%FILS(1), id%DAD_STEPS(1),
     &        id%STEP(1), id%KEEP(28), 
     &        id%KEEP(60), id%KEEP(20), id%KEEP(38),
     &        id%INFO(1) )
      ENDIF
C     Root principal variable
C     for scalapack (KEEP(38)) might have been updated
C     since root variables might have been permuted
C     and/or expanded (MUMPS_EXPAND_TREE) in case of compressed graph
C     It should thus be redistributed to all procs
      IF(((KEEP(494).NE.0).OR.KEEP(13).NE.0)
     &      .AND.(id%KEEP(38).GT.0)) 
     &             THEN  ! grouping at analysis (1 => LR
       CALL MPI_BCAST( id%KEEP(38), 1, MPI_INTEGER, MASTER,
     &     id%COMM, IERR )
      ENDIF
 80   CONTINUE
C     Broadcast errors
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
C     ---------------------------------------------------
C     Broadcast information computed on the master to
C     the slaves.
C     The matrix itself with numerical values and
C     integer data for the arrowhead/element description 
C     will be received at the beginning of FACTO.
C     ---------------------------------------------------
      CALL MPI_BCAST( id%FILS(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%NA(1), id%LNA, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%STEP(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%PROCNODE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%DAD_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%FRERE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR)
      CALL MPI_BCAST( id%NE_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%ND_STEPS(1), id%KEEP(28), MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%SYM_PERM(1), id%N, MPI_INTEGER,
     &     MASTER, id%COMM, IERR )
      IF(KEEP(494).NE.0) THEN
            CALL MPI_BCAST( id%LRGROUPS(1), id%N, MPI_INTEGER,
     &           MASTER, id%COMM, IERR )
      END IF
      IF (KEEP(55) .EQ. 0) THEN
C     Assembled matrix format. Fill up the id%PTRAR array
C     Broadcast id%SYM_PERM needed to fill up id%PTRAR
C     At the end of ANA_N_DIST, id%PTRAR is already on every processor
C     because it is computed in a distributed way.
C     No need to broadcast it again
         CALL DMUMPS_ANA_N_DIST(id, id%PTRAR)
         IF(id%MYID .EQ. MASTER) THEN
C           -----------------------------------
C           For distributed structure on entry,
C           we can now deallocate the complete
C           structure IRN / JCN.
C           -----------------------------------
            IF ( (KEEP(244) .EQ. 1) .AND. (KEEP(54) .EQ. 3) ) THEN
C               IRN and JCN might have already been deallocated
                IF (associated(id%IRN)) THEN
                  DEALLOCATE(id%IRN)
                  NULLIFY(id%IRN)
                ENDIF
                IF (associated(id%JCN)) THEN
                  DEALLOCATE(id%JCN)
                  NULLIFY(id%JCN)
                ENDIF
            END IF
         END IF
      ENDIF
C     
C     Store size of the stack memory for each
C     of the sequential subtree.
      IF((id%KEEP(76).EQ.4).OR.(id%KEEP(76).EQ.6))THEN
         IF(associated(id%DEPTH_FIRST)) THEN
            DEALLOCATE(id%DEPTH_FIRST)
         ENDIF
         allocate(id%DEPTH_FIRST(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         IF(associated(id%DEPTH_FIRST_SEQ)) THEN
            DEALLOCATE(id%DEPTH_FIRST_SEQ)
         ENDIF
         ALLOCATE(id%DEPTH_FIRST_SEQ(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(associated(id%SBTR_ID)) THEN
            DEALLOCATE(id%SBTR_ID)
         ENDIF
         ALLOCATE(id%SBTR_ID(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%DEPTH_FIRST(1:id%KEEP(28))=DEPTH_FIRST(1:id%KEEP(28))
            id%DEPTH_FIRST_SEQ(1:id%KEEP(28))=
     &           DEPTH_FIRST_SEQ(1:id%KEEP(28))
            id%SBTR_ID(1:KEEP(28))=SBTR_ID(1:KEEP(28))
         ENDIF
         CALL MPI_BCAST( id%DEPTH_FIRST(1), id%KEEP(28), MPI_INTEGER,
     &           MASTER, id%COMM, IERR )         
         CALL MPI_BCAST( id%DEPTH_FIRST_SEQ(1), id%KEEP(28),
     &           MPI_INTEGER,MASTER, id%COMM, IERR )  
         CALL MPI_BCAST( id%SBTR_ID(1), id%KEEP(28),
     &           MPI_INTEGER,MASTER, id%COMM, IERR )  
      ELSE
         IF(associated(id%DEPTH_FIRST)) THEN
            DEALLOCATE(id%DEPTH_FIRST)
         ENDIF
         allocate(id%DEPTH_FIRST(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         IF(associated(id%DEPTH_FIRST_SEQ)) THEN
            DEALLOCATE(id%DEPTH_FIRST_SEQ)
         ENDIF
         ALLOCATE(id%DEPTH_FIRST_SEQ(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         IF(associated(id%SBTR_ID)) THEN
            DEALLOCATE(id%SBTR_ID)
         ENDIF
         ALLOCATE(id%SBTR_ID(1),stat=allocok)
         IF (allocok .ne.0) then
            INFO(1)= -7
            INFO(2)= 1
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%DEPTH_FIRST_SEQ'
            END IF
            GOTO 87
         END IF
         id%SBTR_ID(1)=0
         id%DEPTH_FIRST(1)=0
         id%DEPTH_FIRST_SEQ(1)=0
      ENDIF
      IF(id%KEEP(76).EQ.5)THEN
         IF(associated(id%COST_TRAV)) THEN
            DEALLOCATE(id%COST_TRAV)
         ENDIF
         allocate(id%COST_TRAV(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%COST_TRAV'
            END IF
            INFO(1)= -7
            INFO(2)= id%KEEP(28)
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%COST_TRAV(1:id%KEEP(28))=
     &      dble(COST_TRAV_TMP(1:id%KEEP(28)))
         ENDIF
         CALL MPI_BCAST( id%COST_TRAV(1), id%KEEP(28),
     &        MPI_DOUBLE_PRECISION,MASTER, id%COMM, IERR )         
      ELSE
         IF(associated(id%COST_TRAV)) THEN
            DEALLOCATE(id%COST_TRAV)
         ENDIF
         allocate(id%COST_TRAV(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%COST_TRAV(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         id%COST_TRAV(1)=0.0d0
      ENDIF
      IF (id%KEEP(47) .EQ. 4 .OR.
     &     ((id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
         IF(id%MYID .EQ. MASTER)THEN
            DO K=1,id%NSLAVES
               DO J=1,SIZE_TEMP_MEM
                  IF(TEMP_MEM(J,K) < 0.0D0) GOTO 666 
               ENDDO
 666           CONTINUE
               J=J-1
               IF (id%KEEP(46) == 1) THEN
                  IDEST = K - 1
               ELSE
                  IDEST = K
               ENDIF
               IF (IDEST .NE. MASTER) THEN
                  CALL MPI_SEND(J,1,MPI_INTEGER,IDEST,0,
     &                 id%COMM,IERR)
                  CALL MPI_SEND(TEMP_MEM(1,K),J,MPI_DOUBLE_PRECISION,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_LEAF(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_SIZE(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)
                  CALL MPI_SEND(TEMP_ROOT(1,K),J,MPI_INTEGER,
     &                 IDEST, 0, id%COMM,IERR)             
               ELSE
                  IF(associated(id%MEM_SUBTREE)) THEN
                     DEALLOCATE(id%MEM_SUBTREE)
                  ENDIF
                  allocate(id%MEM_SUBTREE(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MEM_SUBTREE'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%NBSA_LOCAL = J
                  id%MEM_SUBTREE(1:J)=TEMP_MEM(1:J,1)
                  IF(associated(id%MY_ROOT_SBTR)) THEN
                     DEALLOCATE(id%MY_ROOT_SBTR)
                  ENDIF
                  allocate(id%MY_ROOT_SBTR(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MY_ROOT_SBTR'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_ROOT_SBTR(1:J)=TEMP_ROOT(1:J,1)
                  IF(associated(id%MY_FIRST_LEAF)) THEN
                     DEALLOCATE(id%MY_FIRST_LEAF)
                  ENDIF
                  allocate(id%MY_FIRST_LEAF(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MY_FIRST_LEAF'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_FIRST_LEAF(1:J)=TEMP_LEAF(1:J,1)
                  IF(associated(id%MY_NB_LEAF)) THEN
                     DEALLOCATE(id%MY_NB_LEAF)
                  ENDIF
                  allocate(id%MY_NB_LEAF(J),stat=allocok)
                  IF (allocok .ne.0) then
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'id%MY_NB_LEAF'
                     END IF
                     INFO(1)= -7
                     INFO(2)= J
                     GOTO 87
                  END IF
                  id%MY_NB_LEAF(1:J)=TEMP_SIZE(1:J,1)
               ENDIF
            ENDDO
         ELSE
            CALL MPI_RECV(id%NBSA_LOCAL,1,MPI_INTEGER,
     &           MASTER,0,id%COMM,STATUS, IERR)
            IF(associated(id%MEM_SUBTREE)) THEN
               DEALLOCATE(id%MEM_SUBTREE)
            ENDIF
            allocate(id%MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'id%MEM_SUBTREE'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_ROOT_SBTR)) THEN
               DEALLOCATE(id%MY_ROOT_SBTR)
            ENDIF
            allocate(id%MY_ROOT_SBTR(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'id%MY_ROOT_SBTR'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_FIRST_LEAF)) THEN
               DEALLOCATE(id%MY_FIRST_LEAF)
            ENDIF
            allocate(id%MY_FIRST_LEAF(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'MY_FIRST_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            IF(associated(id%MY_NB_LEAF)) THEN
               DEALLOCATE(id%MY_NB_LEAF)
            ENDIF
            allocate(id%MY_NB_LEAF(id%NBSA_LOCAL),stat=allocok)
            IF (allocok .ne.0) then
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'MY_NB_LEAF'
               END IF
               INFO(1)= -7
               INFO(2)= id%NBSA_LOCAL
               GOTO 87
            END IF
            CALL MPI_RECV(id%MEM_SUBTREE(1),id%NBSA_LOCAL,
     &           MPI_DOUBLE_PRECISION,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_FIRST_LEAF(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_NB_LEAF(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
            CALL MPI_RECV(id%MY_ROOT_SBTR(1),id%NBSA_LOCAL,
     &           MPI_INTEGER,MASTER,0,
     &           id%COMM,STATUS,IERR)
         ENDIF
      ELSE
         id%NBSA_LOCAL = -999999
         IF(associated(id%MEM_SUBTREE)) THEN
            DEALLOCATE(id%MEM_SUBTREE)
         ENDIF
         allocate(id%MEM_SUBTREE(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MEM_SUBTREE(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_ROOT_SBTR)) THEN
            DEALLOCATE(id%MY_ROOT_SBTR)
         ENDIF
         allocate(id%MY_ROOT_SBTR(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MY_ROOT_SBTR(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_FIRST_LEAF)) THEN
            DEALLOCATE(id%MY_FIRST_LEAF)
         ENDIF
         allocate(id%MY_FIRST_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MY_FIRST_LEAF(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
         IF(associated(id%MY_NB_LEAF)) THEN
            DEALLOCATE(id%MY_NB_LEAF)
         ENDIF
         allocate(id%MY_NB_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MY_NB_LEAF(1)'
            END IF
            INFO(1)= -7
            INFO(2)= 1
            GOTO 87
         END IF
      ENDIF
      IF(id%MYID.EQ.MASTER)THEN
         IF(IS_BUILD_LOAD_MEM_CALLED)THEN 
            DEALLOCATE(TEMP_MEM)
            DEALLOCATE(TEMP_SIZE)
            DEALLOCATE(TEMP_ROOT)
            DEALLOCATE(TEMP_LEAF)
            DEALLOCATE(COST_TRAV_TMP)
            DEALLOCATE(DEPTH_FIRST)
            DEALLOCATE(DEPTH_FIRST_SEQ)
            DEALLOCATE(SBTR_ID)
         ENDIF
      ENDIF
 87   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
C     
      NB_NIV2 = KEEP(56)        ! KEEP(1:110) was broadcast earlier
C     NB_NIV2 is now available on all processors.
      IF (  NB_NIV2.GT.0  ) THEN
C        Allocate arrays on slaves
         if (id%MYID.ne.MASTER) then
            IF (associated(id%CANDIDATES)) THEN
               DEALLOCATE(id%CANDIDATES)
            ENDIF
            allocate(PAR2_NODES(NB_NIV2),
     &           id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     &           STAT=allocok)
            IF (allocok .ne.0) then
               INFO(1)= -7
               INFO(2)= NB_NIV2*(id%NSLAVES+1)
               IF ( LPOK ) THEN
                  WRITE(LP, 150) 'PAR2_NODES/id%CANDIDATES'
               END IF
            end if
         end if
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF ( INFO(1).LT.0 ) GOTO 500
         CALL MPI_BCAST(PAR2_NODES(1),NB_NIV2,
     &        MPI_INTEGER, MASTER, id%COMM, IERR )
         IF (KEEP(24) .NE.0 ) THEN
            CALL MPI_BCAST(id%CANDIDATES(1,1),
     &           (NB_NIV2*(id%NSLAVES+1)),
     &           MPI_INTEGER, MASTER, id%COMM, IERR )
         ENDIF
      ENDIF
      IF ( associated(id%ISTEP_TO_INIV2)) THEN
         DEALLOCATE(id%ISTEP_TO_INIV2)
         NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
      IF ( associated(id%I_AM_CAND)) THEN
         DEALLOCATE(id%I_AM_CAND)
         NULLIFY(id%I_AM_CAND)
      ENDIF
      IF (NB_NIV2.EQ.0) THEN 
C     allocate dummy arrays
C     ISTEP_TO_INIV2 will never be used 
C     Add a parameter SIZE_ISTEP_TO_INIV2 and make
C     it always available in a keep(71)
         id%KEEP(71) = 1
      ELSE
         id%KEEP(71) = id%KEEP(28)
      ENDIF
      allocate(id%ISTEP_TO_INIV2(id%KEEP(71)),
     &     id%I_AM_CAND(max(NB_NIV2,1)),
     &     stat=allocok)
      IF (allocok .gt.0) THEN
         IF ( LPOK ) THEN
            WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
            WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
         END IF
         INFO(1)= -7
         IF (NB_NIV2.EQ.0) THEN
            INFO(2)= 2
         ELSE
            INFO(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
         END IF
         GOTO 321
      ENDIF
      IF ( NB_NIV2 .GT.0 ) THEN
C   If BLR grouping was performed then PAR2_NODES(INIV2) 
C   might then point to a non principal variable
C   for which STEP might be negative
C
         id%ISTEP_TO_INIV2 = -9999
         DO INIV2 = 1, NB_NIV2
            INN = PAR2_NODES(INIV2)
            id%ISTEP_TO_INIV2(abs(id%STEP(INN))) = INIV2
         END DO 
         CALL DMUMPS_BUILD_I_AM_CAND( id%NSLAVES, KEEP(79),
     &        NB_NIV2, id%MYID_NODES,
     &        id%CANDIDATES(1,1), id%I_AM_CAND(1) )
      ENDIF
      IF ( I_AM_SLAVE ) THEN
         IF (associated(id%FUTURE_NIV2)) THEN
            DEALLOCATE(id%FUTURE_NIV2)
            NULLIFY(id%FUTURE_NIV2)
         ENDIF
         allocate(id%FUTURE_NIV2(id%NSLAVES), stat=allocok)
         IF (allocok .gt.0) THEN
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'FUTURE_NIV2'
            END IF
            INFO(1)= -7
            INFO(2)= id%NSLAVES
            GOTO 321
         ENDIF
         id%FUTURE_NIV2=0
         DO INIV2 = 1, NB_NIV2
            IDEST = MUMPS_PROCNODE(
     &           id%PROCNODE_STEPS(abs(id%STEP(PAR2_NODES(INIV2)))),
     &           id%KEEP(199))
            id%FUTURE_NIV2(IDEST+1)=id%FUTURE_NIV2(IDEST+1)+1
         ENDDO
C     Allocate id%TAB_POS_IN_PERE, 
C     TAB_POS_IN_PERE is an array of size (id%NSLAVES+2,NB_NIV2)
C     where NB_NIV2 is the number of type 2 nodes in the tree.
         IF ( associated(id%TAB_POS_IN_PERE)) THEN
            DEALLOCATE(id%TAB_POS_IN_PERE)
            NULLIFY(id%TAB_POS_IN_PERE)
         ENDIF
         allocate(id%TAB_POS_IN_PERE(id%NSLAVES+2,max(NB_NIV2,1)),
     &        stat=allocok)
         IF (allocok .gt.0) THEN
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
               WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
            END IF
            INFO(1)= -7
            IF (NB_NIV2.EQ.0) THEN
               INFO(2)= 2
            ELSE
               INFO(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
            END IF
            GOTO 321
         ENDIF
      END IF
C     deallocate PAR2_NODES  that was computed
C     on master and broadcasted on all slaves
      IF (NB_NIV2.GT.0) DEALLOCATE (PAR2_NODES)
 321  CONTINUE
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
C
      IF ( KEEP(38) .NE. 0 ) THEN
C     -------------------------
C     Initialize root structure
C     -------------------------
         CALL DMUMPS_INIT_ROOT_ANA( id%MYID,
     &        id%NSLAVES, id%N, id%root,
     &        id%COMM_NODES, KEEP( 38 ), id%FILS(1),
     &        id%KEEP(50), id%KEEP(46),
     &        id%KEEP(51)
     &        , id%KEEP(60), id%NPROW, id%NPCOL, id%MBLOCK, id%NBLOCK
     &        )
      ELSE
         id%root%yes = .FALSE.
      END IF
      IF ( KEEP(38) .NE. 0 .and. I_AM_SLAVE ) THEN
C     -----------------------------------------------
C     Check if at least one processor belongs to the
C     root. In the case where all of them have MYROW
C     equal to -1, this could be a problem due to the
C     BLACS. (mpxlf90_r and IBM BLACS).
C     -----------------------------------------------
         CALL MPI_ALLREDUCE(id%root%MYROW, MYROW_CHECK, 1,
     &        MPI_INTEGER, MPI_MAX, id%COMM_NODES, IERR)
         IF ( MYROW_CHECK .eq. -1) THEN
            INFO(1) = -25
            INFO(2) = 0
         END IF
         IF ( id%root%MYROW .LT. -1 .OR.
     &        id%root%MYCOL .LT. -1 ) THEN
            INFO(1) = -25
            INFO(2) = 0
         END IF
         IF ( LPOK .AND. INFO(1) == -25 ) THEN
            WRITE(LP, '(A)')
     &           'Problem with your version of the BLACS.'
            WRITE(LP, '(A)') 'Try using a BLACS version from netlib.'
         ENDIF
      END IF
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( I_AM_SLAVE ) THEN
C{
C     
C     
         IF (KEEP(55) .EQ. 0) THEN
            CALL DMUMPS_ANA_DIST_ARROWHEADS( id%MYID,
     &           id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     &           id%STEP(1), id%PTRAR(1),
     &           id%PTRAR(id%N +1),
     &           id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &           KEEP(1),KEEP8(1), ICNTL(1), id )
         ELSE
            CALL DMUMPS_ANA_DIST_ELEMENTS( id%MYID,
     &           id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     &           id%STEP(1),
     &           id%PTRAR(1),
     &           id%PTRAR(id%NELT+2 ),
     &           id%NELT, 
     &           id%FRTPTR(1), id%FRTELT(1),
     &           KEEP(1), KEEP8(1), ICNTL(1), id%KEEP(50) )
         ENDIF
C}
      ENDIF
C     -----------------------------------------
C     Perform some local analysis on the slaves
C     to estimate the size of the working space
C     for factorization
C     -----------------------------------------
      IF ( I_AM_SLAVE ) THEN
C{
         locI_AM_CAND => id%I_AM_CAND
         locMYID_NODES = id%MYID_NODES
         locMYID       = id%MYID
C           ===================================================
C           Precompute estimates of local_m,local_n
C           (number of rows/columns mapped on each processor)
C           in case of parallel root node.
C           and allocate CANDIDATES
C           ===================================================
C     
            IF ( id%root%yes ) THEN
              LOCAL_M = numroc( id%ND_STEPS(id%STEP(KEEP(38))),
     &              id%root%MBLOCK, id%root%MYROW, 0,
     &              id%root%NPROW )
              LOCAL_M = max(1, LOCAL_M)
              LOCAL_N = numroc( id%ND_STEPS(id%STEP(KEEP(38))),
     &              id%root%NBLOCK, id%root%MYCOL, 0,
     &              id%root%NPCOL )
            ELSE
              LOCAL_M = 0
              LOCAL_N = 0
            END IF
            IF  ( KEEP(60) .EQ. 2 .OR. KEEP(60) .EQ. 3 ) THEN
C           Return minimum nb rows/cols to user
               id%SCHUR_MLOC=LOCAL_M
               id%SCHUR_NLOC=LOCAL_N
C           Also store them in root structure for convenience
               id%root%SCHUR_MLOC=LOCAL_M
               id%root%SCHUR_NLOC=LOCAL_N
            ENDIF
            IF ( .NOT. associated(id%CANDIDATES)) THEN
               ALLOCATE(id%CANDIDATES(id%NSLAVES+1,1), stat=allocok)
               IF (allocok .gt.0) THEN
                     IF ( LPOK ) THEN
                        WRITE(LP, 150) 'CANDIDATES'
                     END IF
                     INFO(1)= -7
                     INFO(2)= id%NSLAVES+1
               ENDIF
            ENDIF
C}
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &             id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
C     -- Allocate and initialise IPOOL with leaves 
C     -- on which stats are performed
      IF ( I_AM_SLAVE ) THEN
C{
           LIPOOL = id%NA(1)
C         LIPOOL is number of leaf nodes and can be 0 
C         (for ex AboveL0 with nbthreads is 1)
          ALLOCATE( IPOOL(max(LIPOOL,1)), 
     &                   stat=allocok)
          IF (allocok .gt.0) THEN
                   IF ( LPOK ) THEN
                      WRITE(LP, 150) 'Allocation IPOOL'
                   END IF
                   INFO(1)= -7
                   INFO(2)=  LIPOOL
          ENDIF
C}
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
C
      IF ( I_AM_SLAVE ) THEN
C{
C            Initialize IPOOL with leaves of complete tree
             DO I=1, LIPOOL
              IPOOL(I) = id%NA(3+I-1)
             ENDDO
             ABOVE_L0                   =.FALSE.
             SIZECB_UNDER_L0            = 0_8
             SIZECB_UNDER_L0_IF_LRCB    = 0_8
             MAX_FRONT_SURFACE_LOCAL_L0 = 0_8
             MAX_SIZE_FACTOR_L0         = 0_8
             ENTRIES_IN_FACTORS_UNDER_L0= 0_8
             ENTRIES_IN_FACTORS_MASTERS_LO = 0_8
             MAXFR_UNDER_L0             = 0
             COST_SUBTREES_UNDER_L0     = 0.0D0
             OPSA_UNDER_L0              = 0.0D0
C
             NE_STEPSPTR => id%NE_STEPS
           KEEP(139) = MAXFR_UNDER_L0
           CALL DMUMPS_ANA_DISTM( locMYID_NODES, id%N, id%STEP(1),
     & id%FRERE_STEPS(1), id%FILS(1), IPOOL, LIPOOL, NE_STEPSPTR(1),
     & id%DAD_STEPS(1), id%ND_STEPS(1), id%PROCNODE_STEPS(1),
     & id%NSLAVES, ABOVE_L0,SIZECB_UNDER_L0,SIZECB_UNDER_L0_IF_LRCB,
     & MAXFR_UNDER_L0, MAX_FRONT_SURFACE_LOCAL_L0, MAX_SIZE_FACTOR_L0, 
     & ENTRIES_IN_FACTORS_UNDER_L0, ENTRIES_IN_FACTORS_MASTERS_LO,
     & COST_SUBTREES_UNDER_L0, OPSA_UNDER_L0, KEEP8(53), KEEP8(54),
     & KEEP8(11), KEEP(26), KEEP(15), KEEP8(12),  KEEP8(14),  
     & KEEP8(32), KEEP8(33), KEEP8(34), KEEP8(35), KEEP8(50), 
     & KEEP8(36), KEEP8(47), KEEP8(37), KEEP8(38), KEEP8(39),
     & KEEP8(40), KEEP8(41), KEEP8(42), KEEP8(43), KEEP8(44), KEEP8(45),
     & KEEP8(46), KEEP8(51), KEEP8(52), KEEP(224),KEEP(225),KEEP(27),
     & RINFO(1),id%CNTL(1), KEEP(1), KEEP8(1), LOCAL_M, LOCAL_N,
     & SBUF_RECOLD8, SBUF_SEND_FR, SBUF_REC_FR, SBUF_SEND_LR,
     & SBUF_REC_LR, id%COST_SUBTREES, KEEP(28), locI_AM_CAND(1),
     & max(KEEP(56),1), id%ISTEP_TO_INIV2(1), id%CANDIDATES(1,1), 
     & INFO(1), INFO(2), KEEP8(15),MAX_SIZE_FACTOR_TMP, 
     & KEEP8(9), ENTRIES_IN_FACTORS_LOC_MASTERS,
     & id%root%yes, id%root%NPROW, id%root%NPCOL
     &           )
           IF (ALLOCATED(IPOOL)) DEALLOCATE(IPOOL)
           NULLIFY(NE_STEPSPTR)
C            SUM_NIRNEC under L0 OMP
             KEEP(137)=0
C            SUM_NIRNEC_OOC under L0 OMP
             KEEP(138)=0
C           DKEEP(15) is used for dynamic load balancing only
C           it corresponds to the number of local operations
C           (in Millions)
            id%DKEEP(15)    = RINFO(1)/1000000.0
            IF(ASSOCIATED(locI_AM_CAND)) NULLIFY(locI_AM_CAND)
            id%MAX_SURF_MASTER = KEEP8(15)
C
            KEEP8(19)=MAX_SIZE_FACTOR_TMP
            KEEP( 29 ) = KEEP(15) + 3* max(KEEP(12),10)
     &           * ( KEEP(15) / 100 + 1)
C     Relaxed value of size of IS is not needed internally;
C     we save it directly in INFO(19)
            INFO( 19 ) = KEEP(225) + 3* max(KEEP(12),10)
     &           * ( KEEP(225) / 100 + 1)
C     =================================
C     Size of S (relaxed with ICNTL(14)
C     ===========================
C     size of S relaxed (FR, IC)
C     ===========================
            KEEP8(13)  = KEEP8(12) + int(KEEP(12),8) *
     &           ( KEEP8(12) / 100_8 + 1_8 )
C     size of S relaxed (FR or LR LU, OOC)
            KEEP8(17)  = KEEP8(14) + int(KEEP(12),8) *
     &           ( KEEP8(14) /100_8 +1_8)
C      size of S relaxed (LR LU, IC)
            K8_33relaxed  = KEEP8(33) + int(KEEP(12),8) *
     &           ( KEEP8(33) /100_8 +1_8)
C     size of S relaxed (LR LU+CB, OOC)
            K8_34relaxed  = KEEP8(34) + int(KEEP(12),8) *
     &           ( KEEP8(34) /100_8 +1_8)
C     size of S relaxed (LR LU+CB, OOC)
            K8_35relaxed  = KEEP8(35) + int(KEEP(12),8) *
     &           ( KEEP8(35) /100_8 +1_8)
C     size of S relaxed (LR CB, IC)
            K8_50relaxed  = KEEP8(50) + int(KEEP(12),8) *
     &           ( KEEP8(50) /100_8 +1_8)
C     KEEP8( 22 ) is the OLD maximum size of receive buffer 
C     that includes CB related communications.
C     KEEP( 43 ) : min size for send buffer
C     KEEP( 44 ) : min size for receive buffer
C     KEEP(43-44) kept for allocating buffers during
C                 factorization phase
         CALL MUMPS_ALLREDUCEI8 ( SBUF_RECOLD8, KEEP8(22), MPI_MAX,
     &                            id%COMM_NODES )
C     We do a max with KEEP(27)=maxfront because for small
C     buffers, we need at least one row of cb to be sent/
C     received.
         SBUF_SEND_FR = max(SBUF_SEND_FR,KEEP(27))
         SBUF_SEND_LR = max(SBUF_SEND_LR,KEEP(27))
         SBUF_REC_FR  = max(SBUF_REC_FR ,KEEP(27))
         SBUF_REC_LR  = max(SBUF_REC_LR ,KEEP(27))
         CALL MPI_ALLREDUCE (SBUF_REC_FR, KEEP(44), 1, 
     &        MPI_INTEGER, MPI_MAX,
     &        id%COMM_NODES, IERR)
         CALL MPI_ALLREDUCE (SBUF_REC_LR, KEEP(380), 1, 
     &        MPI_INTEGER, MPI_MAX,
     &        id%COMM_NODES, IERR)
         IF (KEEP(48)==5) THEN
            KEEP(43)=KEEP(44)
            KEEP(379)=KEEP(380)
         ELSE
            KEEP(43)=SBUF_SEND_FR
            KEEP(379)=SBUF_SEND_LR
         ENDIF
C     
         MIN_BUF_SIZE8 = KEEP8(22) / int(KEEP(238),8)
         MIN_BUF_SIZE8 = min( MIN_BUF_SIZE8, int(huge (KEEP(43)),8))
         MIN_BUF_SIZE  = int( MIN_BUF_SIZE8 )
C
         KEEP(44)  = max(KEEP(44), MIN_BUF_SIZE)
         KEEP(380) = max(KEEP(380), MIN_BUF_SIZE)
         KEEP(43)  = max(KEEP(43), MIN_BUF_SIZE)
         KEEP(379) = max(KEEP(379), MIN_BUF_SIZE)
            IF ( PROK ) THEN
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated INTEGER space for factors         :',
     &              KEEP(26)
               WRITE(MP,'(A,I16) ') 
     &              ' INFO(3), est. real space to store factors   :',
     &              KEEP8(11)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated number of entries in factors      :',
     &              KEEP8(9)
               WRITE(MP,'(A,I16) ') 
     &              ' Current value of space relaxation parameter :',
     &              KEEP(12)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated size of IS (In Core factorization):',
     &              KEEP(29)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated size of S  (In Core factorization):',
     &              KEEP8(13)
               WRITE(MP,'(A,I16) ') 
     &              ' Estimated size of S  (OOC factorization)    :',
     &              KEEP8(17)
            END IF
C}
      ELSE
C     ---------------------
C     Master is not working
C     ---------------------
         ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
         KEEP8(13) = 0_8
         KEEP(29) = 0
         KEEP8(17)= 0_8
         INFO(19) = 0
         KEEP8(11) = 0_8
         KEEP8(81) = 0_8
         KEEP8(82) = 0_8
         KEEP(26) = 0
         KEEP(27) = 0
         RINFO(1) = 0.0D0
         K8_33relaxed = 0_8
         K8_34relaxed = 0_8
         K8_35relaxed = 0_8
         K8_50relaxed = 0_8
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
C     --------------------------------------
C     KEEP8( 26 ) : Real arrowhead size
C     KEEP8( 27 ) : Integer arrowhead size
C     INFO(3)/KEEP8( 11 ) : Estimated real space needed for factors
C     INFO(4)/KEEP( 26 )  : Estimated integer space needed for factors
C     INFO(5)/KEEP( 27 )  : Estimated max front size
C     KEEP8(109)          : Estimated number of entries in factor
C                         (based on ENTRIES_IN_FACTORS_LOC_MASTERS computed 
C                          during DMUMPS_ANA_DISTM, where we assume 
C                          that each master of a node computes
C                          the complete factor size.
C     --------------------------------------
C     note that summing ENTRIES_IN_FACTORS_LOC_MASTERS or 
C     ENTRIES_IN_FACTORS_LOC_MASTERS should lead to the same result
      CALL MUMPS_ALLREDUCEI8( ENTRIES_IN_FACTORS_LOC_MASTERS, 
     &     KEEP8(109), MPI_SUM, id%COMM)
      CALL MUMPS_ALLREDUCEI8( KEEP8(19), KEEP8(119),
     &     MPI_MAX, id%COMM)
      CALL MPI_ALLREDUCE( KEEP(27), KEEP(127), 1,
     &     MPI_INTEGER, MPI_MAX,
     &     id%COMM, IERR)
      CALL MPI_ALLREDUCE( KEEP(26), KEEP(126), 1,
     &     MPI_INTEGER, MPI_SUM,
     &     id%COMM, IERR)
C     NRLADU related: KEEP8(11) holds factors above and under L0
      CALL MUMPS_REDUCEI8( KEEP8(11),
     &   KEEP8(111), MPI_SUM,
     &   MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4( KEEP8(111), INFOG(3) )
C     NRLADU_if_LR_LU related: KEEP8(32) holds factors above 
C                                        and under L0
C     convert it in Megabytes
      RINFO(5) = dble(KEEP8(32)
     &               *int(KEEP(35),8))/1D6
      CALL MUMPS_REDUCEI8( KEEP8(32),
     &                   ITMP8, MPI_SUM,
     &                   MASTER, id%COMM )
C     in Megabytes
      IF (id%MYID.EQ.MASTER) THEN
       RINFOG(15) = dble(ITMP8*int(KEEP(35),8))/1D6
      ENDIF
C     --------------
C     Flops estimate
C     --------------
      CALL MPI_ALLREDUCE( RINFO(1), RINFOG(1), 1,
     &     MPI_DOUBLE_PRECISION, MPI_SUM,
     &     id%COMM, IERR)
C
      CALL MUMPS_SETI8TOI4( KEEP8(11), INFO(3) )
      INFO ( 4 ) = KEEP(  26 )
      INFO ( 5 ) = KEEP(  27 )
      INFO ( 7 ) = KEEP(  29 )
      CALL MUMPS_SETI8TOI4( KEEP8(13), INFO(8) )
      CALL MUMPS_SETI8TOI4( KEEP8(17), INFO(20) )
      CALL MUMPS_SETI8TOI4( KEEP8(9), INFO(24) )
C
      CALL MUMPS_SETI8TOI4( K8_33relaxed, INFO(29) )
      CALL MUMPS_SETI8TOI4( K8_34relaxed, INFO(32) )
      CALL MUMPS_SETI8TOI4( K8_35relaxed, INFO(33) )
      CALL MUMPS_SETI8TOI4( K8_50relaxed, INFO(36) )
      INFOG( 4 ) = KEEP( 126 )
      INFOG( 5 ) = KEEP( 127 )
      CALL MUMPS_SETI8TOI4( KEEP8(109), INFOG(20) )
      CALL DMUMPS_DIAG_ANA(id%MYID, id%COMM, KEEP(1), KEEP8(1),
     &     INFO(1), INFOG(1), RINFO(1), RINFOG(1), ICNTL(1))
C     --------------------------
C     COMPUTE MEMORY ESTIMATIONS
      IF (PROK) WRITE( MP, 112 )
      IF (PROKG .AND. (MPG.NE.MP)) WRITE( MPG, 112 )
C     --------------------------
C     =========================
C     IN-CORE MEMORY STATISTICS
C     =========================
C
         OOC_STRAT = KEEP(201)
         BLR_STRAT = 0  ! no BLR compression
         IF (KEEP(201) .NE. -1) OOC_STRAT=0 ! We want in-core statistics
         PERLU_ON = .FALSE.     ! switch off PERLU to compute KEEP8(2)
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. ! UNDER_L0_OMP
     &        )
         KEEP8(2) = TOTAL_BYTES    
C
C
         PERLU_ON  = .TRUE.
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. ! UNDER_L0_OMP 
     &        )
         IF ( PROK ) THEN
            WRITE(MP,'(A,I12) ')
     & ' Estimated space in MBytes for IC factorization   (INFO(15)):',
     &           TOTAL_MBYTES
         END IF
         id%INFO(15) = TOTAL_MBYTES
C     
C     Centralize memory statistics on the host
C     
C     INFOG(16) = after analysis, est. mem size in Mbytes for facto,
C     for the processor using largest memory
C     INFOG(17) = after analysis, est. mem size in Mbytes for facto,
C     sum over all processors
C     INFOG(18/19) = idem at facto.
C     
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &     id%INFO(15), id%INFOG(16), IRANK )
      IF ( PROKG ) THEN
       IF (PRINT_MAXAVG) THEN
         WRITE( MPG,'(A,I12) ')
     & '    Maximum estim. space in Mbytes, IC facto.    (INFOG(16)):',
     &        id%INFOG(16)
       ENDIF
       WRITE(MPG,'(A,I12) ')
     & '    Total space in MBytes, IC factorization      (INFOG(17)):'
     &        ,id%INFOG(17)
      END IF
C        =========================================
C        NOW COMPUTE OUT-OF-CORE MEMORY STATISTICS
C        (except when OOC_STRAT is equal to -1 in
C        which case IC and OOC statistics are
C        identical)
C        =========================================
         OOC_STRAT = KEEP(201)
         BLR_STRAT = 0  ! no BLR compression
#if defined(OLD_OOC_NOPANEL)
         IF (OOC_STRAT .NE. -1) OOC_STRAT=2
#else
         IF (OOC_STRAT .NE. -1) OOC_STRAT=1
#endif
         PERLU_ON = .FALSE.     ! PERLU NOT taken into account
C        Used to compute KEEP8(3) (minimum number of bytes for OOC)
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES, 
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. ! UNDER_L0_OMP 
     &         )
         KEEP8(3) = TOTAL_BYTES
C
         PERLU_ON  = .TRUE.     ! PERLU taken into account
         CALL DMUMPS_MAX_MEM( KEEP(1), KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30),
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     &        OOC_STRAT, BLR_STRAT, PERLU_ON, TOTAL_BYTES,
     &        IDUMMY, BDUMMY, .FALSE., 
     &        .FALSE. ! UNDER_L0_OMP 
     &         )
         id%INFO(17) = TOTAL_MBYTES
C
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &     id%INFO(17), id%INFOG(26), IRANK )
      IF ( PROKG  ) THEN
       IF (PRINT_MAXAVG) THEN
         WRITE( MPG,'(A,I12) ')
     & '    Maximum estim. space in Mbytes, OOC facto.   (INFOG(26)):',
     &        id%INFOG(26)
       ENDIF
       WRITE(MPG,'(A,I12) ')
     & '    Total space in MBytes,  OOC factorization    (INFOG(27)):'
     &        ,id%INFOG(27)
      END IF
      IF (KEEP(494).NE.0) THEN
C        =========================================
C        NOW COMPUTE BLR statistics 
C        =========================================
         SUM_OF_PEAKS = .TRUE.
         CALL  DMUMPS_MEM_ESTIM_BLR_ALL( SUM_OF_PEAKS, 
     &        KEEP(1), KEEP8(1),
     &        id%MYID, id%COMM, 
     &        id%N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &        id%KEEP8(30), id%NSLAVES,
     &        id%INFO, id%INFOG, PROK, MP, PROKG, MPG 
     &        )
C
      END IF
C     -------------------------
C     Define a specific mapping
C     for the user
C     -------------------------
      IF ( id%MYID. eq. MASTER .AND. KEEP(54) .eq. 1 ) THEN
         IF (associated( id%MAPPING)) THEN
            DEALLOCATE( id%MAPPING)
         ENDIF
         allocate( id%MAPPING(id%KEEP8(28)), stat=allocok)
         IF ( allocok .GT. 0 ) THEN
            INFO(1) = -7
            CALL MUMPS_SETI8TOI4(id%KEEP8(28), INFO(2))
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'id%MAPPING'
            END IF
            GOTO 92
         END IF
         allocate(IWtemp( id%N ), stat=allocok)
         IF ( allocok .GT. 0 ) THEN
            INFO(1)=-7
            INFO(2)=id%N
            IF ( LPOK ) THEN
               WRITE(LP, 150) 'IWtemp(N)'
            END IF
            GOTO 92
         END IF
         IF ( id%KEEP8(28) .EQ. 0_8 ) THEN
           IRN_PTR => IDUMMY_ARRAY
           JCN_PTR => IDUMMY_ARRAY
         ELSE
           IRN_PTR => id%IRN
           JCN_PTR => id%JCN
         ENDIF
         CALL DMUMPS_BUILD_MAPPING(
     &        id%N, id%MAPPING(1), id%KEEP8(28),
     &        IRN_PTR(1),JCN_PTR(1), id%PROCNODE_STEPS(1),
     &        id%STEP(1),
     &        id%NSLAVES, id%SYM_PERM(1),
     &        id%FILS(1), IWtemp, id%KEEP(1),id%KEEP8(1),
     &        id%root%MBLOCK, id%root%NBLOCK,
     &        id%root%NPROW, id%root%NPCOL )
         DEALLOCATE( IWtemp )
 92      CONTINUE
      END IF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 500
C 
  500 CONTINUE
C     Deallocate allocated working space
      IF (allocated(PROCNODE)) DEALLOCATE(PROCNODE)
      IF (allocated(XNODEL)) DEALLOCATE(XNODEL)
      IF (allocated(NODEL)) DEALLOCATE(NODEL)
      IF (allocated(IPOOL)) DEALLOCATE(IPOOL)
      IF (allocated(SIZEOFBLOCKS)) DEALLOCATE(SIZEOFBLOCKS)
      IF (allocated(DOF2BLOCK)) DEALLOCATE(DOF2BLOCK)
      CALL MUMPS_AB_FREE_LMAT(LMAT_BLOCK)
      CALL MUMPS_AB_FREE_LMAT(LUMAT)
      CALL MUMPS_AB_FREE_LMAT(LUMAT_REMAP)
      CALL MUMPS_AB_FREE_GCOMP(GCOMP)
      CALL MUMPS_AB_FREE_GCOMP(GCOMP_DIST)
C     Standard deallocations (error or not)
      IF (associated(NFSIZPTR)) DEALLOCATE(NFSIZPTR)
      IF (associated(FREREPTR)) DEALLOCATE(FREREPTR)
      IF (associated(FILSPTR)) DEALLOCATE(FILSPTR)
      IF (associated(id%BLKPTR).AND.BLKPTR_ALLOCATED) THEN
        DEALLOCATE(id%BLKPTR)
        nullify(id%BLKPTR)
      ENDIF
      IF (associated(id%BLKVAR).AND.BLKVAR_ALLOCATED) THEN
        DEALLOCATE(id%BLKVAR)
        nullify(id%BLKVAR)
      ENDIF
      KEEP8(26)=max(1_8,KEEP8(26))
      KEEP8(27)=max(1_8,KEEP8(27))
      RETURN
 110  FORMAT(/' ****** ANALYSIS STEP ********'/)
 112  FORMAT(/' MEMORY ESTIMATIONS ... '/
     &   ' Estimations with standard Full-Rank (FR) factorization:')
 145  FORMAT(' ELAPSED TIME SPENT IN BLR CLUSTERING    =',F12.4)
 150  FORMAT(
     & /' ** FAILURE DURING DMUMPS_ANA_DRIVER, DYNAMIC ALLOCATION OF',
     &     A30)
      END SUBROUTINE DMUMPS_ANA_DRIVER
      SUBROUTINE DMUMPS_ANA_CHECK_KEEP(id)
C     This subroutine decodes the control parameters,
C     stores them in the KEEP array, and performs a
C     consistency check on the KEEP array.
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC)  :: id
C     internal variables
      INTEGER   :: LP, MP, MPG, I
      INTEGER   :: MASTER
      LOGICAL   :: PROK, PROKG, LPOK
      PARAMETER( MASTER = 0 )
      LP  = id%ICNTL( 1 )
      MP  = id%ICNTL( 2 )
      MPG = id%ICNTL( 3 )
C     LP     : errors
C     MP     : INFO
      LPOK  = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK  = (( MP  .GT. 0 ).AND.(id%ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(id%ICNTL(4).GE.2))
C     Re-intialize few KEEPs entries corresponding
C     to stat that are incremented such
C     the number of split nodes:
      id%KEEP(61)=0
      IF (id%MYID.eq.MASTER) THEN
        id%KEEP(256) = id%ICNTL(7) ! copy ordering option
        id%KEEP(252) = id%ICNTL(32)
        IF (id%KEEP(252) < 0 .OR. id%KEEP(252) > 1 ) THEN
          id%KEEP(252) = 0
        ENDIF
C       Which factors to store
        id%KEEP(251) = id%ICNTL(31)
        IF (id%KEEP(251) < 0 .OR. id%KEEP(251) > 2 ) THEN
          id%KEEP(251)=0
        ENDIF
C       For unsymmetric matrices, if forward solve 
C       performed during facto,
C       no reason to store L factors at all. Reset
C       KEEP(251) accordingly... except if the user
C       tells that no solve is needed.
        IF (id%KEEP(50) .EQ. 0 .AND. id%KEEP(252).EQ.1) THEN
          IF (id%KEEP(251) .NE. 1) id%KEEP(251) = 2
        ENDIF
C       Symmetric case, even if no backward needed,
C       store all factors
        IF (id%KEEP(50) .NE.0 .AND. id%KEEP(251) .EQ. 2) THEN
          id%KEEP(251) = 0
        ENDIF
C       Case of solve not needed:
        IF (id%KEEP(251) .EQ. 1) THEN
          id%KEEP(201) = -1
C         In that case, id%ICNTL(22) will
C         be ignored in future phases
        ELSE
C         Reset id%KEEP(201) -- typically for the case
C         of a previous analysis with KEEP(201)=-1
          id%KEEP(201) = 0
        ENDIF
        IF (id%KEEP(252).EQ.1) THEN
          id%KEEP(253) = id%NRHS
          IF (id%KEEP(253) .LE. 0) THEN
            id%INFO(1)=-42
            id%INFO(2)=id%NRHS
            RETURN
          ENDIF
        ELSE
          id%KEEP(253) = 0
        ENDIF
      ENDIF
      IF ( (id%KEEP(24).NE.0) .AND.
     &     id%NSLAVES.eq.1 ) THEN
         id%KEEP(24) = 0
      END IF
      IF ( (id%KEEP(24).EQ.0) .AND.
     &     id%NSLAVES.GT.1 ) THEN
         id%KEEP(24) = 8
      ENDIF
      IF ( (id%KEEP(24).NE.0)  .AND. (id%KEEP(24).NE.1)  .AND.
     &     (id%KEEP(24).NE.8)  .AND. (id%KEEP(24).NE.10) .AND.
     &     (id%KEEP(24).NE.12) .AND. (id%KEEP(24).NE.14) .AND.
     &     (id%KEEP(24).NE.16) .AND. (id%KEEP(24).NE.18)) THEN
         id%KEEP(24) = 8
      END IF
C****************************************************
C     
C     The master is doing most of the work
C     
C     NOTE:  Treatment of the errors on the master=
C     Go to the next SPMD part of the code in which
C     the first statement must be a call to PROPINFO
C     
C****************************************************
C     =========================================
C     Check (raise error or modify) some input
C     parameters or KEEP values on the master.
C     =========================================
      id%KEEP8(21) = int(id%KEEP(85),8)
      IF ( id%MYID .EQ. MASTER ) THEN
C     -- OOC/Incore strategy 
         IF (id%KEEP(201).NE.-1) THEN
           id%KEEP(201)=id%ICNTL(22)
           IF (id%KEEP(201) .GT. 0) THEN
#if defined(OLD_OOC_NOPANEL)
             id%KEEP(201)=2
#else
             id%KEEP(201)=1
#endif
           ENDIF
         ENDIF
C        ----------------------------
C        Save id%ICNTL(18) (distributed
C        matrix on entry) in id%KEEP(54)
C        ----------------------------
         id%KEEP(54) = id%ICNTL(18)
         IF ( id%KEEP(54) .LT. 0 .or. id%KEEP(54).GT.3 ) THEN
            IF ( PROKG ) THEN
               WRITE(MPG, *) ' Out-of-range value for id%ICNTL(18).'
               WRITE(MPG, *) ' Used 0 ie matrix not distributed'
            END IF
            id%KEEP(54) = 0
         END IF
         IF ( id%KEEP(54) .EQ. 1 ) THEN
            IF ( PROKG ) THEN
               WRITE(MPG, *) ' Option id%ICNTL(18)=1 is obsolete.'
               WRITE(MPG, *) ' We recommend not to use it.'
               WRITE(MPG, *) ' It will disappear in a future release'
            END IF
         END IF
C     -----------------------------------------
C     Save id%ICNTL(5) (matrix format) in id%KEEP(55)
C     -----------------------------------------
         id%KEEP(55) = id%ICNTL(5)
         IF ( id%KEEP(55) .LT. 0 .OR. id%KEEP(55) .GT. 1 ) THEN
            IF ( PROKG ) THEN
               WRITE(MPG, *) ' Out-of-range value for id%ICNTL(5).'
               WRITE(MPG, *) ' Used 0 ie matrix is assembled'
            END IF
            id%KEEP(55) = 0
         END IF
         id%KEEP(60) = id%ICNTL(19)
         IF ( id%KEEP( 60 ) .LE. 0 ) id%KEEP( 60 ) = 0
         IF ( id%KEEP( 60 ) .GT. 3 ) id%KEEP( 60 ) = 0
         IF (id%KEEP(60) .NE. 0 .AND. id%SIZE_SCHUR == 0 ) THEN
            IF (PROKG) THEN
              WRITE(MPG,'(A)')
     &        ' ** Schur option ignored because SIZE_SCHUR=0'
            ENDIF
            id%KEEP(60)=0
         END IF
C        ---------------------------------------
C        Save SIZE_SCHUR in a KEEP, for possible
C        check at factorization and solve phases
C        ---------------------------------------
         IF ( id%KEEP(60) .NE.0 ) THEN
            id%KEEP(116) = id%SIZE_SCHUR
            IF (id%SIZE_SCHUR .LT. 0 .OR. id%SIZE_SCHUR .GE. id%N) THEN
              id%INFO(1)=-49
              id%INFO(2)=id%SIZE_SCHUR
              RETURN
            ENDIF
C           List of Schur variables provided by user.
            IF ( .NOT. associated( id%LISTVAR_SCHUR ) ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 8
               RETURN
            ELSE IF (size(id%LISTVAR_SCHUR)<id%SIZE_SCHUR) THEN
               id%INFO(1) = -22
               id%INFO(2) = 8
               RETURN
            END IF
         ENDIF
         IF (id%KEEP(60) .EQ. 3 .AND. id%KEEP(50).NE.0) THEN
            IF (id%MBLOCK > 0 .AND. id%NBLOCK > 0 .AND.
     &           id%NPROW > 0 .AND. id%NPCOL > 0 ) THEN
               IF (id%NPROW *id%NPCOL .LE. id%NSLAVES) THEN
C     We will eventually have to "symmetrize the
C     Schur complement. For that NBLOCK and MBLOCK
C     must be equal.
                  IF (id%MBLOCK .NE. id%NBLOCK ) THEN
                     id%INFO(1)=-31
                     id%INFO(2)=id%MBLOCK - id%NBLOCK
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C     Check the ordering strategy and compatibility with
C     other control parameters
         id%KEEP(244) = id%ICNTL(28)
         id%KEEP(245) = id%ICNTL(29)
#if ! defined(parmetis) && ! defined(parmetis3)        
         IF ((id%KEEP(244) .EQ. 2) .AND. (id%KEEP(245) .EQ. 2)) THEN
            id%INFO(1)  = -38
            IF ( LPOK ) THEN
               WRITE(LP,'("ParMETIS not available.")')
            END IF
            RETURN
         END IF
#endif
#if ! defined(ptscotch)         
         IF ((id%KEEP(244) .EQ. 2) .AND. (id%KEEP(245) .EQ. 1)) THEN
            id%INFO(1)  = -38
            IF ( LPOK ) THEN
               WRITE(LP,'("PT-SCOTCH not available.")')
            END IF
            RETURN
         END IF
#endif
C        Analysis strategy is set to automatic in case of out-of-range values.
         IF((id%KEEP(244) .GT. 2) .OR.
     &        (id%KEEP(244) .LT. 0)) id%KEEP(244)=0
         IF(id%KEEP(244) .EQ. 0) THEN ! Automatic
C           One could check for availability of parallel ordering
C           tools, or for possible options incompatible with //
C           analysis to decide (e.g. avoid returning an error if
C           // analysis not compatible with some option but user
C           lets MUMPS decide to choose sequential or paralllel
C           analysis)
C           Current strategy for automatic is sequential analysis
            id%KEEP(244) = 1
         ELSE IF (id%KEEP(244) .EQ. 2) THEN
            IF(id%KEEP(55) .NE. 0) THEN
               id%INFO(1)  = -39
               IF (LPOK) THEN
               WRITE(LP,
     &              '("Incompatible values for ICNTL(5), ICNTL(28)")')
               WRITE(LP,
     &              '("Parallel analysis is not possible if the")')
               WRITE(LP,
     &              '("matrix is not assembled")')
               ENDIF
               RETURN
            ELSE IF(id%KEEP(60) .NE. 0) THEN
               id%INFO(1)  = -39
               IF (LPOK) THEN
               WRITE(LP,
     &              '("Incompatible values for ICNTL(19), ICNTL(28)")')
               WRITE(LP,
     &              '("Parallel analysis is not possible if SCHUR")')
               WRITE(LP,
     &              '("complement must be returned")')
               ENDIF
               RETURN
            END IF
C     In the case where there are too few processes to do
C     the parallel analysis we simply revert to sequential version
            IF(id%NSLAVES .LT. 2) THEN
               id%KEEP(244) = 1
               IF(PROKG) WRITE(MPG,
     &              '("Too few processes.
     & Reverting to sequential analysis")',advance='no')
               IF(id%KEEP(245) .EQ. 1) THEN
C                 Scotch necessarily available because pt-scotch
C                 is, otherwise an error would have occurred
                  IF(PROKG) WRITE(MPG, '(" with SCOTCH.")')
                  id%KEEP(256) = 3
               ELSE IF(id%KEEP(245) .EQ. 2) THEN
C                 Metis necessarily available because parmetis
C                 is, otherwise an error would have occurred
                  IF(PROKG) WRITE(MPG, '(" with Metis.")')
                  id%KEEP(256) = 5
               ELSE
                  IF(PROKG) WRITE(MPG, '(".")')
                  id%KEEP(256) = 7
               END IF
            END IF
C     In the case where there the input matrix is too small to do
C     the parallel analysis we simply revert to sequential version
            IF(id%N .LE. 50) THEN
               id%KEEP(244) = 1
               IF(PROKG) WRITE(MPG,
     &            '("Input matrix is too small for the parallel
     & analysis. Reverting to sequential analysis")',advance='no')
               IF(id%KEEP(245) .EQ. 1) THEN
                  IF(PROKG) WRITE(MPG, '(" with SCOTCH.")')
                  id%KEEP(256) = 3
               ELSE IF(id%KEEP(245) .EQ. 2) THEN
                  IF(PROKG) WRITE(MPG, '(" with Metis.")')
                  id%KEEP(256) = 5
               ELSE
                  IF(PROKG) WRITE(MPG, '(".")')
                  id%KEEP(256) = 7
               END IF
            END IF
         END IF
         id%INFOG(32) = id%KEEP(244)
         IF ( (id%KEEP(244) .EQ. 1) .AND.
     &        (id%KEEP(256) .EQ. 1) ) THEN
C     ordering given, PERM_IN must be of size N
            IF ( .NOT. associated( id%PERM_IN ) ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 3
               RETURN
            ELSE IF ( size( id%PERM_IN ) < id%N ) THEN
               id%INFO(1) = -22
               id%INFO(2) = 3
               RETURN
            END IF
         ENDIF
C     Check KEEP(9-10) for level 2
         IF (id%KEEP(9) .LE. 1 ) id%KEEP(9) = 500
         IF ( id%KEEP8(21) .GT. 0_8 ) THEN 
            IF ((id%KEEP8(21).LE.1_8) .OR.
     &          (id%KEEP8(21).GT.int(id%KEEP(9),8)))
     &         id%KEEP8(21) = int(min(id%KEEP(9),100),8)
         ENDIF
C     
         IF (id%KEEP(48). EQ. 1 ) id%KEEP(48) = -12345
C     
         IF ( (id%KEEP(48).LT.0) .OR. (id%KEEP(48).GT.5) ) THEN
            id%KEEP(48)=5
         ENDIF
C     Schur 
C     Given ordering must be compatible with Schur variables.
         IF ( (id%KEEP(60) .NE. 0) .AND. (id%KEEP(256) .EQ. 1) ) THEN
            DO I = 1, id%SIZE_SCHUR
               IF (id%PERM_IN(id%LISTVAR_SCHUR(I))
     &              .EQ. id%N-id%SIZE_SCHUR+I)
     &              CYCLE
C              -------------------------------
C              Problem with PERM_IN: -22/3
C              Above constrained explained in
C              doc of PERM_IN in user guide.
C              -------------------------------
               id%INFO(1) = -4
               id%INFO(2) = id%LISTVAR_SCHUR(I)
               RETURN
               IF (PROKG) THEN
                  WRITE(MPG,'(A)')
     & ' ** Ignoring user-ordering, because incompatible with Schur.'
                  WRITE(MPG,'(A)') ' ** id%ICNTL(7) treated as 0.'
               END IF
               EXIT
            ENDDO
         END IF
C     
C     Note that schur is not compatible with
C     
C     1/Max-trans DONE
C     2/Null space
C     3/Ordering given DONE
C     4/Scaling
C     5/Iterative Refinement
C     6/Error analysis
C     7/Parallel Analysis
C     
C     Graph modification prior to ordering (id%ICNTL(12) option)
C     id%KEEP (95) will hold the eventually modified value of id%ICNTL(12)
C     
         id%KEEP(95) = id%ICNTL(12)
C        reset to usual ordering (KEEP(95)=1)
C          - when matrix is not general symmetric 
C          - for out-of-range values
         IF (id%KEEP(50).NE.2) id%KEEP(95) = 1
         IF ((id%KEEP(95).GT.3).OR.(id%KEEP(95).LT.0)) id%KEEP(95) = 1
C     MAX-TRANS
C     
C     id%KEEP (23) will hold the eventually modified value of id%ICNTL(6)
C     (maximum transversal if >= 1)
C     
         id%KEEP(23) = id%ICNTL(6)
C     
C     
C     --------------------------------------------
C     Avoid max-trans unsymmetric permutation in case of
C     matrix is symmetric with SYM=1 or 
C     ordering is given,
C     or matrix is in element form, or Schur is asked
C     or initial matrix is distributed
C     --------------------------------------------
         IF (id%KEEP(23).LT.0.OR.id%KEEP(23).GT.7) id%KEEP(23) = 0
C        still forbid max trans for SYM=1 case
         IF ( id%KEEP(50) .EQ. 1 ) THEN
            IF (id%KEEP(23) .NE. 0) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A)')
     & ' ** Max-trans not needed with SYM=1 factorization'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(95) .GT. 1) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A)')
     & ' ** ICNTL(12) ignored: not needed with SYM=1 factorization'
               END IF
            ENDIF
            id%KEEP(95) = 1
         END IF
C     
         IF  (id%KEEP(60) .GT. 0) THEN
            IF (id%KEEP(23) .NE. 0) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A)')
     &                 ' ** Max-trans not allowed because of Schur'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(52).EQ.-2) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A)')
     & ' ** Scaling during analysis not allowed because of Schur'
               ENDIF
               id%KEEP(52) = 0
            ENDIF
C     also forbid compressed/constrained ordering...
            IF (id%KEEP(95) .GT. 1) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A)')
     & ' ** ICNTL(12) option not allowed because of Schur'
               END IF
            ENDIF
            id%KEEP(95) = 1
         END IF
         IF ( (id%KEEP(23) .NE. 0) .AND. (id%KEEP(256).EQ.1)) THEN
            id%KEEP(23) = 0
            IF (PROKG) THEN
               WRITE(MPG,'(A,A)')
     &          ' ** Maximum transversal (ICNTL(6)) not allowed ',
     &         'because ordering is given'
            END IF
         END IF
         IF ( id%KEEP(256) .EQ. 1 ) THEN
            IF (id%KEEP(95) > 1 .AND. PROKG) THEN
               WRITE(MPG,'(A)')
     &      ' ** ICNTL(12) option incompatible with given ordering'
            END IF
            id%KEEP(95) = 1
         END IF
         IF (id%KEEP(54) .NE. 0) THEN
            IF( id%KEEP(23) .NE. 0 ) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A,A)')
     &            ' ** Maximum transversal (ICNTL(6)) not allowed ',
     &            'because matrix is distributed'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (id%KEEP(52).EQ.-2) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A,A)')
     &            ' ** Scaling (ICNTL(8)) during analysis not ', 
     &            'allowed because matrix is distributed)'
               ENDIF
            ENDIF
            id%KEEP(52) = 0
            IF (id%KEEP(95) .GT. 1 .AND. MPG.GT.0) THEN
               WRITE(MPG,'(A,A)')
     &         ' ** ICNTL(12) option not allowed because matrix is ', 
     &         'distributed'
            ENDIF
            id%KEEP(95) = 1
         END IF
         IF ( id%KEEP(55) .NE. 0 ) THEN
            IF( id%KEEP(23) .NE. 0 ) THEN
               IF (PROKG) THEN
                  WRITE(MPG,'(A,A)')
     &            ' ** Maximum transversal (ICNTL(6)) not allowed ',
     &            'for matrices in elemental format'
               END IF
               id%KEEP(23) = 0
            ENDIF
            IF (PROKG .AND. id%KEEP(52).EQ.-2) THEN
               WRITE(MPG,'(A)')
     &            ' ** Scaling (ICNTL(8)) not allowed ',
     &            'for matrices in elemental format'
            ENDIF
            id%KEEP(52) = 0
            id%KEEP(95) = 1
         ENDIF
C     In the case where parallel analysis is done, column permutation
C     is not allowed
         IF(id%KEEP(244) .EQ. 2) THEN
            IF(id%KEEP(23) .EQ. 7) THEN
C     Automatic hoice: set it to 0
               id%KEEP(23) = 0
            ELSE IF (id%KEEP(23) .GT. 0) THEN
               id%INFO(1)  = -39
               id%KEEP(23) = 0
               IF (LPOK) THEN
               WRITE(LP,
     &              '("Incompatible values for ICNTL(6), ICNTL(28)")')
               WRITE(LP,
     &              '("Maximum transversal not allowed
     &                 in parallel analysis")')
               ENDIF
               RETURN
            END IF
         END IF
C     --------------------------------------------
C     Avoid distributed entry for element matrix.
C     --------------------------------------------
         IF ( id%KEEP(54) .NE. 0 .AND. id%KEEP(55) .NE. 0 ) THEN
            id%KEEP(54) = 0
            IF (PROKG) THEN
               WRITE(MPG,'(A)')
     & ' ** Distributed entry not available for element matrix'
            END IF
         ENDIF
C     ----------------------------------
C     Choice of symbolic analysis option
C     ----------------------------------
         IF (id%ICNTL(58).NE.1 .and. id%ICNTL(58).NE.2
     &        .and. id%ICNTL(58).NE.3 ) THEN
            id%KEEP(106)=1
C     Automatic choice leads to new symbolic
C     factorization except(see below) if KEEP(256)==1.
         ELSE
            id%KEEP(106)=id%ICNTL(58)
            IF (id%KEEP(106).EQ.3) THEN
C            option not available 
             id%KEEP(106)=1
            ENDIF
         ENDIF
C     modify input parameters to avoid incompatible
C     input data between ordering, scaling and maxtrans
C     note that if id%ICNTL(12)/id%KEEP(95) = 0 then
C     the automatic choice will be done in ANA_O
         IF(id%KEEP(50) .EQ. 2) THEN
C     LDLT case
            IF( .NOT. associated(id%A) ) THEN
C     constraint ordering can be computed only if values are
C     given to analysis
               IF(id%KEEP(95) .EQ. 3) THEN
                  id%KEEP(95) = 2
               ENDIF
            ENDIF
            IF(id%KEEP(95) .EQ. 3 .AND. id%KEEP(256) .NE. 2) THEN
C     if constraint and ordering is not AMF then use compress
               IF (PROK) WRITE(MP,*)
     &              'WARNING: DMUMPS_ANA_O constrained ordering not ', 
     &              'available with selected ordering'
               id%KEEP(95) = 2
            ENDIF 
            IF(id%KEEP(95) .EQ. 3) THEN
C     if constraint ordering required then we need to compute scaling
C     and max trans
C     NOTE that if we enter this condition then
C     id%A is associated because of the test above:
C     (IF( .NOT. associated(id%A) ) THEN)
               id%KEEP(23) = 5
               id%KEEP(52) = -2
            ELSE IF(id%KEEP(95) .EQ. 2 .AND. 
     &              (id%KEEP(23) .EQ. 0 .OR. id%KEEP(23) .EQ. 7) ) THEN
C     compressed ordering requires max trans but not necessary scaling
               IF( associated(id%A) ) THEN
                  id%KEEP(23) = 5
               ELSE
C     we can do compressed ordering without
C     information on the numerical values:
C     a maximum transversal already provides
C     information on the location of off-diagonal
C     nonzeros which can be candidates for 2x2
C     pivots
                  id%KEEP(23) = 1
               ENDIF
            ELSE IF(id%KEEP(95) .EQ. 1) THEN
               id%KEEP(23) = 0
            ELSE IF(id%KEEP(95) .EQ. 0 .AND. id%KEEP(23) .EQ. 0) THEN
C     if max trans desactivated then the automatic choice for type of ord
C     is set to 1, which means that we will use usual ordering
C     (no constraints or compression)
               id%KEEP(95) = 1
            ENDIF
         ELSE
            id%KEEP(95) = 1
         ENDIF
C     --------------------------------
C     Save ICNTL(56) (QR) in KEEP(53)
C     Will be broadcasted to all other
C     nodes in routine DMUMPS_BDCAST
C     --------------------------------
         id%KEEP(53)=0
         IF(id%KEEP(86).EQ.1)THEN
C     Force the exchange of both the memory and flops information during 
C     the factorization
            IF(id%KEEP(47).LT.2) id%KEEP(47)=2
         ENDIF
         IF(id%KEEP(48).EQ.5)THEN
            IF(id%KEEP(50).EQ.0)THEN
               id%KEEP(87)=50
               id%KEEP(88)=50
            ELSE
               id%KEEP(87)=70
               id%KEEP(88)=70
            ENDIF
         ENDIF
         IF((id%NSLAVES.EQ.1).AND.(id%KEEP(76).GT.3))THEN
            id%KEEP(76)=2
         ENDIF
         IF(id%KEEP(81).GT.0)THEN
            IF(id%KEEP(47).LT.2) id%KEEP(47)=2
         ENDIF
C
C        -- Save Block Low Rank input parameter 
         id%KEEP(494) = id%ICNTL(35)
         IF (id%KEEP(494).EQ.1) THEN
C        -- Automatic BLR option setting
           id%KEEP(494)= 2
         ENDIF
         IF ( id%KEEP(494).EQ.4) id%KEEP(494)=0
         IF ((id%KEEP(494).LT.0).OR.(id%KEEP(494).GT.4)) THEN
C          Out of range values treated as 0
           id%KEEP(494) = 0
         ENDIF
         IF(id%KEEP(494).NE.0) THEN
C        test BLR incompatibilities
C
          id%KEEP(464) = id%ICNTL(38)
          IF (id%KEEP(464).LT.0.OR.(id%KEEP(464).GT.1000)) THEN
C          Out of range values treated as 0
           id%KEEP(464) = 0
          ENDIF
C         LR is incompatible with elemental matrices, forbid it at analysis
          IF (id%KEEP(55).NE.0) THEN
            IF (LPOK) WRITE(LP,*)
     &           " *** BLR feature currently incompatible "
     &           ,"with elemental matrices"
C           BLR for elt entry might be developed in the future
            id%INFO(1)=-800
            id%INFO(2)=5
            RETURN
          ENDIF
C
C         LR incompatible with forward in facto 
          IF (id%KEEP(252).NE.0) THEN
            IF (LPOK) WRITE(LP,*)
     &        " *** BLR feature currently incompatible"
     &       ," with forward during factorization"
             id%INFO(1) = -43
             id%INFO(2) = 35
             RETURN
          ENDIF
C
         ENDIF
C
         IF(id%KEEP(494).NE.0) THEN
C         id%KEEP(469)=0,1,2,3,4
          IF ((id%KEEP(469).GT.4).OR.(id%KEEP(469).LT.0)) THEN
               id%KEEP(469)=0
          ENDIF
C         Not implemented yet               
          IF (id%KEEP(469).EQ.4) id%KEEP(469)=0
C         id%KEEP(471)=-1,0,1
          IF ((id%KEEP(471).LT.-1).AND.(id%KEEP(471).GT.1)) THEN
               id%KEEP(471)=-1
          ENDIF
C         id%KEEP(472)=0 or 1
          IF ((id%KEEP(472).NE.0).AND.(id%KEEP(472).NE.1)) THEN
               id%KEEP(472)=1
          ENDIF
C         id%KEEP(475)=0,1,2,3
          IF ((id%KEEP(475).GT.3).OR.(id%KEEP(475).LT.0)) THEN
               id%KEEP(475)=0
          ENDIF
C         id%KEEP(482)=0,1,2,3
          IF ((id%KEEP(482).GT.3).OR.(id%KEEP(482).LT.0)) THEN
               id%KEEP(482)=0
          ENDIF
          IF((id%KEEP(487).LT.0)) THEN
             id%KEEP(487)= 2 ! default value
          ENDIF
C         id%KEEP(488)>0
          IF((id%KEEP(488).LE.0)) THEN
              id%KEEP(488)=  8*id%KEEP(6) 
          ENDIF
C         id%KEEP(490)>0
          IF((id%KEEP(490).LE.0)) THEN
             id%KEEP(490) = 128
          ENDIF
C         KEEP(491)>0
          IF((id%KEEP(491).LE.0)) THEN
            id%KEEP(491) = 1000
          ENDIF
         ENDIF
C 
      id%KEEP(13) = 0
C     Analysis by Blocks
      id%KEEP(13) = id%ICNTL(15)
      IF (id%KEEP(13).GT.1) THEN
CV0   out-of range values
        id%KEEP(13) = 0
      ENDIF
      IF (id%KEEP(13).LT.0) THEN
        IF (mod(id%N,-id%KEEP(13)) .NE.0) THEN
          IF ( LPOK ) THEN
             WRITE(LP,'(A,I8)') 
     &       " ICNTL(15)=", id%ICNTL(15), 
     &       " is incompatible with N=", id%N
          ENDIF
          id%INFO(1) = -57
          id%INFO(2) = 1
        ENDIF
        IF (associated(id%BLKPTR)) THEN
          IF ( LPOK ) THEN
             WRITE(LP,'(A,I8)') 
     &       " ICNTL(15)=", id%ICNTL(15), 
     &       " is incompatible with BLKPTR provided by user"
          ENDIF
          id%INFO(1) = -57
          id%INFO(2) = 4
        ENDIF
      ENDIF
      IF (  (id%KEEP(13).EQ.0)            .AND.
     &      (.NOT. associated(id%BLKPTR)) .AND.
     &      (.NOT. associated(id%BLKVAR))
     &   )
     &  THEN
       IF ((id%KEEP(54).EQ.3).AND.(id%KEEP(244).NE.2)) THEN
          id%KEEP(13)=-1
       ENDIF
      ENDIF
      IF ( (id%KEEP(13).EQ.0           ) .AND.
     &      (.NOT. associated(id%BLKPTR)) .AND.
     &      (.NOT. associated(id%BLKVAR)) .AND.
     &      (id%KEEP(244).NE.2)
     &   )
     &  THEN
C       unsymmetic assembled matrices with or without BLR,
C       also in case of centralized matrix (if
C       matrix is distributed, then KEEP(13) has
C       been set to -1 in the block above)
        IF (id%KEEP(50).EQ.0.AND. id%KEEP(55).EQ.0) THEN
C         Respect decision taken for Maxtrans
C         since it will be switch off because 
C         if one activates the analysis by block
          IF ( (id%KEEP(23).LT.0) .OR. (id%KEEP(23).GT.7)
     &       ) THEN
             id%KEEP(13)=-1
          ENDIF
        ENDIF
      ENDIF
      IF ( (id%KEEP(13).EQ.0)  .AND.
     &     (id%KEEP(55).NE.0) 
     &   ) THEN
         IF (PROKG) WRITE(MPG,'(A,A)')
     &           " ** Analysis by block is incompatible ",
     &           "with elemental matrices"
C           switch off analysis by block
            id%KEEP(13)= 0
      ENDIF
      IF ( (id%KEEP(13).NE.0) .AND.
     &     (id%KEEP(106).NE.1) 
     &   ) THEN
         IF (PROKG) WRITE(MPG,'(A,A,I4)')
     &           " ** Analysis by block compatible ",
     &           "ONLY with SYMQAMD based symbolic factorization ", 
     &           id%KEEP(106)
C           switch off analysis by block
            id%KEEP(13)= 0
      ENDIF
      IF ( (id%KEEP(13).NE.0) .AND.
     &     (id%KEEP(244).EQ.2) 
     &   ) THEN
         IF (PROKG) WRITE(MPG,'(A,A)')
     &           " ** Analysis by block is incompatible ",
     &           "with parallel ordering "
C           switch off analysis by block
            id%KEEP(13)= 0
      ENDIF
      IF ( (id%KEEP(13).NE.0) .AND.
     &     (id%KEEP(60).NE.0) 
     &   ) THEN
         IF (PROKG) WRITE(MPG,'(A,A)')
     &           " ** Analysis by block is incompatible ",
     &           "with Schur "
C           switch off analysis by block
            id%KEEP(13)= 0
      ENDIF
      IF (id%KEEP(13).NE.0) THEN
C      Maximum transversal not compatible with analysis by block
       IF (id%KEEP(23).NE.0) THEN
C        in case of automatic choice (id%KEEP(27).EQ.7)
C        do not print message
         IF (PROKG.AND.id%KEEP(23).NE.7) WRITE(MPG,'(A,A)')
     &      " ** Maximum transversal (ICNTL(6)) ",
     &      "not compatible with analysis by block"
C        switch off max transversal
         id%KEEP(23)= 0
       ENDIF
C      - compression for LDLT
       IF (id%KEEP(95).NE.1) THEN
C        in case of automatic choice (id%KEEP(95).EQ.0)
C        do not print message
         IF (PROKG.AND.id%KEEP(95).NE.0) WRITE(MPG,'(A,A)')
     &      " ** ICNTL(12) not compatible with ",
     &      " analysis by block"
C        switch off 2x2 preprocessing for symmetric matrices
         id%KEEP(95) = 1 
       ENDIF
      ENDIF
C
C     end id%MYID.EQ.MASTER 
      END IF 
      RETURN
      END SUBROUTINE DMUMPS_ANA_CHECK_KEEP
      SUBROUTINE DMUMPS_GATHER_MATRIX(id)
C     This subroutine gathers a distributed matrix
C     on the host node 
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      TYPE(DMUMPS_STRUC)  :: id
C     local variables
      INTEGER, ALLOCATABLE :: REQPTR(:,:)
      INTEGER(8), ALLOCATABLE :: MATPTR(:)
      INTEGER(8), ALLOCATABLE :: MATPTR_cp(:)
      INTEGER(8)           :: IBEG8, IEND8
      INTEGER              :: INDX
      INTEGER              :: LP, MP, MPG, I, K
      INTEGER(8)           :: I8
      LOGICAL              :: PROK, PROKG
C     
C     messages are split into blocks of size BLOCKSIZE
C               (smaller than IOVFLO (=2^31-1))
C     on all processors
      INTEGER(4)           :: IOVFLO
      INTEGER              :: BLOCKSIZE 
      INTEGER              :: MAX_NBBLOCK_loc, NBBLOCK_loc
      INTEGER              :: SIZE_SENT, NRECV
      LOGICAL              :: OMP_FLAG, I_AM_SLAVE
      INTEGER(8)           :: NZ_loc8
C     for validation only:
      INTEGER              :: NB_BLOCKS, NB_BLOCK_SENT
      LP  = id%ICNTL( 1 )
      MP  = id%ICNTL( 2 )
      MPG = id%ICNTL( 3 )
C     LP     : errors
C     MP     : INFO
      PROK  = (( MP  .GT. 0 ).AND.(id%ICNTL(4).GE.2))
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG = (PROKG.AND.(id%ICNTL(4).GE.2))
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &     ( id%MYID .eq. MASTER .AND.
     &     id%KEEP(46) .eq. 1 ) )
C     iovflo = huge(INTEGER, kind=4)
      IOVFLO = huge(IOVFLO)  
C     we do not want too large messages
      BLOCKSIZE = int(max(100000_8,int(IOVFLO,8)/200_8))
      IF ( id%KEEP(46) .EQ. 0 .AND. id%MYID .EQ. MASTER ) THEN
C     host-node mode: master has no entries.
         id%KEEP8(29) = 0_8
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
C     -----------------------------------
C     Allocate small arrays for pointers
C     into arrays IRN/JCN
C     -----------------------------------
         ALLOCATE( MATPTR( id%NPROCS ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) =  id%NPROCS
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) ' array MATPTR'
            END IF
            GOTO 13
         END IF
         ALLOCATE( MATPTR_cp( id%NPROCS ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) =  id%NPROCS
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) ' array MATPTR'
            END IF
            GOTO 13
         END IF
C     -----------------------------------
C     Allocate a small array for requests
C     -----------------------------------
         ALLOCATE( REQPTR( id%NPROCS-1, 2 ), STAT = IERR )         
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            id%INFO(2) = 2 * (id%NPROCS-1)
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'array REQPTR'
            END IF
            GOTO 13
         END IF
C     --------------------
C     Allocate now IRN/JCN
C     --------------------
         ALLOCATE( id%IRN( id%KEEP8(28) ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            CALL MUMPS_SETI8TOI4(id%KEEP8(28),id%INFO(2))
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'array IRN'
            END IF
            GOTO 13
         END IF
         ALLOCATE( id%JCN( id%KEEP8(28) ), STAT = IERR )
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -7
            CALL MUMPS_SETI8TOI4(id%KEEP8(28),id%INFO(2))
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'array JCN'
            END IF
            GOTO 13
         END IF
      END IF
 13   CONTINUE
C     Propagate errors
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF ( id%INFO(1) < 0 ) RETURN
C     -------------------------------------
C     Get numbers of non-zeros for everyone
C     and count total and maximum 
C     nb of blocks of size BLOCKSIZE
C     that slaves will sent
C     -------------------------------------
      IF ( id%MYID .EQ. MASTER ) THEN
C        each block will correspond to 2 messages (IRN_LOC,JCN_LOC)
         NB_BLOCK_SENT = 0
         MAX_NBBLOCK_loc  = 0
         DO I = 1, id%NPROCS - 1
            CALL MPI_RECV( MATPTR( I+1 ), 1, 
     &           MPI_INTEGER8, I,
     &           COLLECT_NZ, id%COMM, STATUS, IERR )
            NBBLOCK_loc = ceiling(dble(MATPTR(I+1))/dble(BLOCKSIZE))
            MAX_NBBLOCK_loc = max(MAX_NBBLOCK_loc, NBBLOCK_loc)
            NB_BLOCK_SENT = NB_BLOCK_SENT + NBBLOCK_loc
         END DO
         IF ( id%KEEP(46) .eq. 0 ) THEN
            MATPTR( 1 ) = 1_8
         ELSE
            NZ_loc8=id%KEEP8(29)
            MATPTR( 1 ) = NZ_loc8 + 1_8
         END IF
C     --------------
C     Build pointers
C     --------------
         DO I = 2, id%NPROCS
            MATPTR( I ) = MATPTR( I ) + MATPTR( I-1 )
         END DO
      ELSE
         NZ_loc8=id%KEEP8(29)
         CALL MPI_SEND( NZ_loc8, 1, MPI_INTEGER8, MASTER,
     &        COLLECT_NZ, id%COMM, IERR )
      END IF
      IF ( id%MYID .eq. MASTER  ) THEN
C     -----------------------------------------------
C     Bottleneck is here master; use synchronous send
C     for slaves, but asynchronous receives on master
C     Then while master receives indices do the local
C     copies for better overlap. 
C     (If master has other things to do, he could try
C     to do them here.)
C     ------------------------------------
C       copy pointers to position in IRN/JCN 
        MATPTR_cp = MATPTR
        IF ( id%KEEP8(29) .NE. 0_8 ) THEN
            OMP_FLAG = ( id%KEEP8(29).GE.50000_8 )
!$OMP PARALLEL DO PRIVATE(I8)
!$OMP&         IF(OMP_FLAG)
            DO I8=1,id%KEEP8(29)
               id%IRN(I8) = id%IRN_loc(I8)
               id%JCN(I8) = id%JCN_loc(I8)
            ENDDO
!$OMP END PARALLEL DO
        ENDIF
C
C     Compute position for each block to be received
C     and store it.
        NB_BLOCKS = 0
C       at least one slave will send MAX_NBBLOCK_loc 
C       couple of messages (IRN_loc/JCN_loc)
        DO K = 1, MAX_NBBLOCK_loc
C        Post irecv for all messages from proc I 
C        that have been sent
         NRECV = 0
         DO I = 1, id%NPROCS - 1
C           Check if message was sent
            IBEG8     = MATPTR_cp( I )
            IF ( IBEG8 .LT. MATPTR(I+1))  THEN
C             Count number of request in NRECV
              NRECV = NRECV + 2
              IEND8 = min(IBEG8+int(BLOCKSIZE,8)-1_8, 
     &                    MATPTR(I+1)-1_8)
C             update pointer for receiving messages
C             from proc I in MATPTR_cp:
              MATPTR_cp( I ) = IEND8 + 1_8
              SIZE_SENT   = int(IEND8 -  IBEG8 + 1_8)
              NB_BLOCKS   = NB_BLOCKS + 1
C
              CALL MPI_IRECV( id%IRN(IBEG8), SIZE_SENT, MPI_INTEGER,
     &           I, COLLECT_IRN, id%COMM, REQPTR(I,1), IERR )
C
              CALL MPI_IRECV( id%JCN(IBEG8), SIZE_SENT, MPI_INTEGER,
     &           I, COLLECT_JCN, id%COMM, REQPTR(I,2), IERR )
            ELSE
             REQPTR( I,1 ) = MPI_REQUEST_NULL
             REQPTR( I,2 ) = MPI_REQUEST_NULL
            ENDIF
         END DO
C        Wait set of messages corresponding to current block
C        ( we dont exploit the fact that 
C            messages are not overtaking 
C            (if sent by one source to the same destination)  )
C
C        Loop on only non MPI_REQUEST_NULL requests
         DO I = 1, NRECV
             CALL MPI_WAITANY
     &           ( 2 * (id%NPROCS-1), REQPTR( 1, 1 ), INDX, 
     &           STATUS, IERR )
         ENDDO
C
C       process next block
      END DO
        DEALLOCATE( REQPTR )
        DEALLOCATE( MATPTR )
        DEALLOCATE( MATPTR_cp )
C     end of reception by master
      ELSE
C     -----------------------------
C     Send only if size is not zero
C     -----------------------------
         IF ( id%KEEP8(29) .NE. 0_8 ) THEN
           NZ_loc8=id%KEEP8(29) 
C          send by blocks of size BLOCKSIZE
           DO I8=1_8, NZ_loc8, int(BLOCKSIZE,8)
            SIZE_SENT = BLOCKSIZE
            IF (NZ_loc8-I8+1_8.LT.int(BLOCKSIZE,8)) THEN
              SIZE_SENT = int(NZ_loc8-I8+1_8)
            ENDIF
            CALL MPI_SEND( id%IRN_loc(I8), SIZE_SENT,
     &           MPI_INTEGER, MASTER,
     &           COLLECT_IRN, id%COMM, IERR )
            CALL MPI_SEND( id%JCN_loc(I8), SIZE_SENT,
     &           MPI_INTEGER, MASTER,
     &           COLLECT_JCN, id%COMM, IERR )
           END DO
         END IF
      END IF
      RETURN
 150  FORMAT(
     &/' ** FAILURE DURING DMUMPS_GATHER_MATRIX, DYNAMIC ALLOCATION OF',
     &     A30)
      END SUBROUTINE DMUMPS_GATHER_MATRIX
      SUBROUTINE DMUMPS_DUMP_PROBLEM(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     If id%WRITE_PROBLEM has been set by the user,
C     possibly on all processors in case of distributed
C     matrix, open a file and dumps the matrix and/or
C     the right hand side. In case the last characters
C     of id.WRITE_PROBLEM are "bin" (uppercase letters
C     are also accepted), then the matrix is written
C     in binary stream format (a C routine is called to
C     avoid depending on the access='stream' mode that
C     is only available since Fortran 2003). In that case,
C     a small header file is also written.
C     Otherwise, this subroutine calls
C     DMUMPS_DUMP_MATRIX (to write the matrix in
C     matrix-market format) and DMUMPS_DUMP_RHS.
C     The routine should be called on all MPI processes.
C
C     Examples:
C     1/ WRITE_PROBLEM='mymatrix.txt', centralized matrix
C       mymatrix.txt contains the matrix in matrix-market format
C     2/ WRITE_PROBLEM='mymatrix.txt', distributed matrix
C       mymatrix.txt<i> contains the portion of the matrix
C       on process <i>, in matrix-market format
C     3/ WRITE_PROBLEM='mymatrix.bin', centralized matrix
C       mymatrix.bin contains the matrix in binary format
C       mymatrix.header contains a short description in text format,
C                       with the first line identical to the one of
C                       a matrix-market format
C     4/ WRITE_PROBLEM='mymatrix.bin', distributed matrix
C       mymatrix.bin<i> contains the portion of the matrix
C                       on process <i>, in binary format
C
C       mymatrix.header contains a short description in text format,
C          with the first line identical to matrix-market format
C
C     If a centralized, dense, RHS is available, it is also written,
C     either in matrix-market or binary format (if WRITE_PROBLEM
C     has a .bin extension). In that case the filename for the RHS
C     is WRITE_PROBLEM//".rhs". If written in binary form, information
C     on the RHS is also provided in the header file.
C
      INCLUDE 'mpif.h'
C
C     Arguments
C     =========
C
      TYPE(DMUMPS_STRUC)  :: id
C
C     Local variables
C     ===============
C
      INTEGER              :: MASTER, IERR, I
      INTEGER              :: IUNIT
      LOGICAL              :: IS_ELEMENTAL
      LOGICAL              :: IS_DISTRIBUTED
      LOGICAL              :: NAME_INITIALIZED
      INTEGER              :: DO_WRITE, DO_WRITE_CHECK
      CHARACTER(LEN=20)    :: IDSTR
      LOGICAL              :: I_AM_SLAVE, I_AM_MASTER
      INTEGER              :: L
      LOGICAL              :: BINARY_FORMAT, DUMP_RHS,
     &                        DUMP_BLKPTR, DUMP_BLKVAR
      INTEGER              :: IS_A_PROVIDED, IS_A_PROVIDED_GLOB
      DOUBLE PRECISION, TARGET                :: A_DUMMY(1)
      INTEGER, TARGET                :: IRN_DUMMY(1), JCN_DUMMY(1)
      INTEGER, POINTER, DIMENSION(:) :: IRN_PASSED, JCN_PASSED
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: A_PASSED
      PARAMETER( MASTER = 0 )
      IUNIT = 69
      I_AM_SLAVE = ( id%MYID .NE. MASTER  .OR.
     &     ( id%MYID .EQ. MASTER .AND.
     &     id%KEEP(46) .EQ. 1 ) )
      I_AM_MASTER = (id%MYID.EQ.MASTER)
      NAME_INITIALIZED = id%WRITE_PROBLEM(1:20)
     &     .NE. "NAME_NOT_INITIALIZED"
      BINARY_FORMAT = .FALSE.
      L=len_trim(id%WRITE_PROBLEM)
      IF (L.GT.4) THEN
        IF ( id%WRITE_PROBLEM(L-3:L-3) .EQ. '.' .AND.
     &       ( id%WRITE_PROBLEM(L-2:L-2) .EQ. 'b' .OR.
     &         id%WRITE_PROBLEM(L-2:L-2) .EQ. 'B' ) .AND.
     &       ( id%WRITE_PROBLEM(L-1:L-1) .EQ. 'i' .OR.
     &         id%WRITE_PROBLEM(L-1:L-1) .EQ. 'I' ) .AND.
     &       ( id%WRITE_PROBLEM(L:L) .EQ. 'n' .OR.
     &         id%WRITE_PROBLEM(L:L) .EQ. 'N' ) ) THEN
          BINARY_FORMAT = .TRUE.
        ENDIF
      ENDIF
C     Check if RHS should also be dumped
      DUMP_RHS = id%MYID.EQ.MASTER .AND.
     &     associated(id%RHS) .AND. NAME_INITIALIZED
      DUMP_RHS = DUMP_RHS .AND. id%NRHS .GE. 1
      DUMP_RHS = DUMP_RHS .AND. id%N .GE. 1
      DUMP_RHS = DUMP_RHS .AND. id%ICNTL(20) .EQ. 0
C     Check if BLKPTR and/or BLKVAR should also be dumped
      DUMP_BLKPTR = .FALSE.
      DUMP_BLKVAR = .FALSE.
C     Remark: if id%KEEP(54) = 1 or 2, the structure
C     is centralized at analysis. Since DMUMPS_DUMP_PROBLEM
C     is called at analysis phase, we define IS_DISTRIBUTED
C     as below, which implies that the structure of the problem
C     is distributed in IRN_loc/JCN_loc at analysis.
C     equal to 
      IS_DISTRIBUTED = (id%KEEP(54) .EQ. 3)
      IS_ELEMENTAL   = (id%KEEP(55) .NE. 0)
      IF (id%MYID.EQ.MASTER .AND. .NOT. IS_DISTRIBUTED) THEN
C        ====================
C        Matrix is assembled
C        and centralized
C        ====================
        IF (NAME_INITIALIZED) THEN
          IF ( BINARY_FORMAT ) THEN
            IF (id%KEEP8(28) .EQ. 0_8) THEN
C             Special case of empty matrix
              A_PASSED   => A_DUMMY
              IRN_PASSED => IRN_DUMMY
              JCN_PASSED => JCN_DUMMY
              IS_A_PROVIDED = 1
            ELSE IF (associated(id%A)) THEN
              A_PASSED=>id%A
              IRN_PASSED => id%IRN
              JCN_PASSED => id%JCN
              IS_A_PROVIDED = 1
            ELSE
              A_PASSED => A_DUMMY
              IRN_PASSED => id%IRN
              JCN_PASSED => id%JCN
              IS_A_PROVIDED = 0
            ENDIF
            OPEN( IUNIT, FILE=id%WRITE_PROBLEM(1:L-4)//'.header' )
            CALL DMUMPS_DUMP_HEADER( IUNIT, id%N,
     &      IS_A_PROVIDED, id%KEEP(50), IS_DISTRIBUTED,
     &      id%NSLAVES, id%KEEP8(28), DUMP_RHS, id%NRHS,
     &      DUMP_BLKPTR, DUMP_BLKVAR, id%NBLK, id%ICNTL(15) )
            CLOSE( IUNIT )
            CALL MUMPS_DUMPMATBINARY_C( id%N, id%KEEP8(28),
     &      id%KEEP(35),
     &      IRN_PASSED(1), JCN_PASSED(1), A_PASSED(1),
     &      IS_A_PROVIDED,
     &      trim(id%WRITE_PROBLEM)//char(0) )
          ELSE
            OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM))
            CALL DMUMPS_DUMP_MATRIX( id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &         IS_DISTRIBUTED,  ! = .FALSE., centralized
     &         IS_ELEMENTAL,    ! Elemental or not
     &         .FALSE.)
            CLOSE(IUNIT)
          ENDIF
        ENDIF
      ELSE IF (id%KEEP(54).EQ.3) THEN
C        =====================
C        Matrix is distributed
C        =====================
         IF ( .NOT.NAME_INITIALIZED
     &        .OR. .NOT. I_AM_SLAVE )THEN
            DO_WRITE = 0
         ELSE
            DO_WRITE = 1
         ENDIF
         CALL MPI_ALLREDUCE(DO_WRITE, DO_WRITE_CHECK, 1,
     &        MPI_INTEGER, MPI_SUM, id%COMM, IERR)
C        -----------------------------------------
C        If yes, each processor writes its share
C        of the matrix in a file in matrix market
C        format (otherwise nothing written). We
C        append the process id to the filename.
C        Safer in case all filenames are the
C        same if all processors share the same
C        file system.
C        -----------------------------------------
         IF (DO_WRITE_CHECK.EQ.id%NSLAVES .AND. I_AM_SLAVE) THEN
            WRITE(IDSTR,'(I9)') id%MYID_NODES
            IF (BINARY_FORMAT) THEN
              IF (id%KEEP8(29) .EQ. 0_8) THEN
C               Special case of empty matrix
                A_PASSED   => A_DUMMY
                IRN_PASSED => IRN_DUMMY
                JCN_PASSED => JCN_DUMMY
C               (consider that A is provided when NNZ_loc=0)
                IS_A_PROVIDED = 1
              ELSE IF (associated(id%A_loc)) THEN
                A_PASSED=>id%A_loc
                IRN_PASSED => id%IRN_loc
                JCN_PASSED => id%JCN_loc
                IS_A_PROVIDED = 1
              ELSE
                A_PASSED => A_DUMMY
                IRN_PASSED => id%IRN_loc
                JCN_PASSED => id%JCN_loc
                IS_A_PROVIDED = 0
              ENDIF
              CALL MPI_ALLREDUCE( IS_A_PROVIDED,
     &             IS_A_PROVIDED_GLOB, 1,
     &             MPI_INTEGER, MPI_PROD, id%COMM_NODES, IERR )
C             IS_A_PROVIDED_GLOB = 1 => dump numerical values
C             IS_A_PROVIDED_GLOB = 0 => some processes did not provide
C                                   numerical values, dump only pattern,
C                                   and indicate this in the header
              IF ( id%MYID_NODES.EQ.0) THEN
C               Print header on first MPI worker (only one global header
C               file in case of distributed matrix), replacing the .bin
C               extension by a .header extension
                OPEN( IUNIT, FILE=id%WRITE_PROBLEM(1:L-4)//'.header' )
                CALL DMUMPS_DUMP_HEADER( IUNIT, id%N,
     &          IS_A_PROVIDED_GLOB, id%KEEP(50), IS_DISTRIBUTED,
     &          id%NSLAVES, id%KEEP8(28), DUMP_RHS, id%NRHS,
     &          DUMP_BLKPTR, DUMP_BLKVAR, id%NBLK, id%ICNTL(15) )
                CLOSE( IUNIT )
              ENDIF
              CALL MUMPS_DUMPMATBINARY_C( id%N, id%KEEP8(29),
     &        id%KEEP(35),
     &        IRN_PASSED(1), JCN_PASSED(1), A_PASSED(1),
     &        IS_A_PROVIDED_GLOB,
     &        trim(id%WRITE_PROBLEM)//trim(adjustl(IDSTR))//char(0) )
            ELSE
              OPEN(IUNIT,
     &             FILE=trim(id%WRITE_PROBLEM)//trim(adjustl(IDSTR)))
              CALL DMUMPS_DUMP_MATRIX(id,
     &           IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &           IS_DISTRIBUTED,           ! =.TRUE., distributed
     &           IS_ELEMENTAL,             ! Elemental or not
     &           .FALSE.)
              CLOSE(IUNIT)
            ENDIF
         ENDIF
C     ELSE ...
C     Nothing written in other cases.
      ENDIF
C     ===============
C     Right-hand side
C     ===============
      IF ( DUMP_RHS ) THEN
        IF (BINARY_FORMAT) THEN
C         dump RHS in binary format
          CALL MUMPS_DUMPRHSBINARY_C( id%N, id%NRHS, id%LRHS, id%RHS(1),
     &    id%KEEP(35),
     &    trim(id%WRITE_PROBLEM)//'.rhs'//char(0) )
        ELSE
C         dump RHS in matrix-market format
          OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM) //".rhs")
          CALL DMUMPS_DUMP_RHS(IUNIT, id)
          CLOSE(IUNIT)
        ENDIF
      ENDIF
      IF ( DUMP_BLKPTR ) THEN
          IF (BINARY_FORMAT) THEN
!           suppress trailing '.bin' and use '.blkptr'
            OPEN( IUNIT, FILE=id%WRITE_PROBLEM(1:L-4)//'.blkptr' )
          ELSE
!           just append '.blkptr'
            OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM)//".blkptr")
          ENDIF
          WRITE(IUNIT,'(I9)') id%NBLK
          DO I=1,id%NBLK+1
            WRITE(IUNIT,'(I9)') id%BLKPTR(I)
          ENDDO
          CLOSE(IUNIT)
      ENDIF
      IF ( DUMP_BLKVAR ) THEN
          IF (BINARY_FORMAT) THEN
!           suppress trailing '.bin' and use '.blkvar'
            OPEN( IUNIT, FILE=id%WRITE_PROBLEM(1:L-4)//'.blkvar' )
          ELSE
!           just append '.blkvar'
            OPEN(IUNIT,FILE=trim(id%WRITE_PROBLEM)//".blkvar")
          ENDIF
          DO I=1,id%N
            WRITE(IUNIT,'(I9)') id%BLKVAR(I)
          ENDDO
          CLOSE(IUNIT)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_DUMP_PROBLEM
      SUBROUTINE DMUMPS_DUMP_HEADER( IUNIT, N, IS_A_PROVIDED_GLOB,
     &  SYM, IS_DISTRIBUTED, NSLAVES, NNZTOT, DUMP_RHS, NRHS,
     &  DUMP_BLKPTR, DUMP_BLKVAR, NBLK, ICNTL15 )
C
C  Purpose:
C  =======
C
C  Write a small header file, similar to matrix-market headers,
C  to accompany a matrix written in binary format.
C
      INTEGER, INTENT(IN) :: IUNIT, N, IS_A_PROVIDED_GLOB , SYM, NSLAVES
      INTEGER(8), INTENT(IN) :: NNZTOT
      LOGICAL, INTENT(IN) :: IS_DISTRIBUTED, DUMP_RHS
      INTEGER, INTENT(IN) :: NRHS
      LOGICAL, INTENT(IN) :: DUMP_BLKPTR, DUMP_BLKVAR
      INTEGER, INTENT(IN) :: NBLK
      INTEGER, INTENT(IN) :: ICNTL15
C
C  Local declarations:
C  ==================
C
      CHARACTER (LEN=10)   :: SYMM
      CHARACTER (LEN=8)    :: ARITH
C     1/ write a line identical to first line of matrix-market header
      IF ( IS_A_PROVIDED_GLOB .EQ. 1 ) THEN
          ARITH='real'
      ELSE
        ARITH='pattern'
      ENDIF
      IF (SYM .eq. 0) THEN
        SYMM="general"
      ELSE
        SYMM="symmetric"
      END IF
      WRITE(IUNIT,'(A,A,A,A)') '%%MatrixMarket matrix coordinate ',
     &           trim(ARITH)," ",trim(SYMM)
C     2/ indicate if matrix is distributed or centralized,
C     then describe binary file content and format
      IF ( IS_DISTRIBUTED ) THEN
        WRITE(IUNIT,FMT='(A,I5,A)')
     &  '% Matrix is distributed (MPI ranks=',NSLAVES,')'
      ELSE
        WRITE(IUNIT,FMT='(A)')
     &  '% Matrix is centralized'
      ENDIF
      WRITE(IUNIT,FMT='(A)')
     &    '% Unformatted stream IO (no record boundaries):'
      IF (ARITH(1:7).EQ.'pattern') THEN
        IF (IS_DISTRIBUTED) THEN
          WRITE(IUNIT,'(A)')
     &    '%    N,NNZ_loc,IRN_loc(1:NNZ_loc),JCN_loc(1:NNZ_loc)'
        ELSE
          WRITE(IUNIT,'(A)')
     &    '%    N,NNZ,IRN(1:NNZ),JCN(1:NNZ)'
        ENDIF
        WRITE(IUNIT,'(A)') '%    (numerical values not provided)'
      ELSE
        IF (IS_DISTRIBUTED) THEN
          WRITE(IUNIT,'(A)')
     &    '%    N,NNZ_loc,IRN_loc(1:NNZ_loc),JCN_loc(1:NNZ_loc),'//
     &    'A_loc(1:NNZ_loc)'
        ELSE
          WRITE(IUNIT,'(A)') '%    N/NNZ/IRN(1:NNZ),JCN(1:NNZ),A(1:NNZ)'
        ENDIF
        WRITE(IUNIT,'(A)') '%    Double precision storage'
      ENDIF
      IF ( IS_DISTRIBUTED ) THEN
        WRITE(IUNIT,'(A,/,A)')
     &  '%    N,IRN_loc(i),JCN_loc(i): 32 bits',
     &  '%    NNZ_loc: 64 bits'
      ELSE
        WRITE(IUNIT,'(A,/,A)')
     &  '%    N,IRN(i),JCN(i): 32 bits',
     &  '%    NNZ: 64 bits'
      ENDIF
      WRITE(IUNIT,FMT='(A,I12)') '% Matrix order: N=',N
      WRITE(IUNIT,FMT='(A,I12)') '% Matrix nonzeros: NNZ=',NNZTOT
      IF (DUMP_RHS) THEN
        WRITE(IUNIT,FMT='(A)') '%'
        WRITE(IUNIT,FMT='(A,/,A,I10,A,I5)')
     &  '% A RHS was also written to disk by columns in binary form.',
     &  '%    Size: N rows x NRHS columns with N=',N,'  NRHS=',NRHS
        WRITE(IUNIT,FMT='(A,I12,A)')
     &  '%    Total:',int(N,8)*int(NRHS,8),' scalar values.'
        WRITE(IUNIT,'(A)') '%    Double precision storage'
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_DUMP_HEADER
      SUBROUTINE DMUMPS_DUMP_MATRIX
     & (id, IUNIT, I_AM_SLAVE, I_AM_MASTER,
     &  IS_DISTRIBUTED, IS_ELEMENTAL, PATTERN_ONLY )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose:
C  =======
C     This subroutine dumps a routine in matrix-market format
C     if the matrix is assembled, and in "MUMPS" format (see
C     example in the MUMPS users'guide, if the matrix is
C     centralized and elemental).
C     The routine can be called on all processors. In case of
C     distributed assembled matrix, each processor writes its
C     share as a matrix market file on IUNIT (IUNIT may have
C     different values on different processors).
C
C
C
C  Arguments (input parameters)
C  ============================
C
C     IUNIT: should be set to the Fortran unit where
C            data should be written.
C     I_AM_SLAVE: .TRUE. except on a non working master
C     IS_DISTRIBUTED: .TRUE. if matrix is distributed,
C                     i.e., if IRN_loc/JCN_loc are provided.
C     IS_ELEMENTAL  : .TRUE. if matrix is elemental
C     id            : main MUMPS structure
C
      LOGICAL, intent(in) :: I_AM_SLAVE,
     &                       I_AM_MASTER,
     &                       IS_DISTRIBUTED,
     &                       IS_ELEMENTAL,
     &                       PATTERN_ONLY
      INTEGER, intent(in) :: IUNIT
      TYPE(DMUMPS_STRUC), intent(in)  :: id
C
C  Local variables:
C  ===============
C
      CHARACTER (LEN=10)   :: SYMM
      CHARACTER (LEN=8)    :: ARITH
      INTEGER(8)           :: I8, NNZ_i
C
C  Executable statements:
C  =====================
      IF (I_AM_MASTER .AND. .NOT. IS_DISTRIBUTED .AND.
     &     .NOT. IS_ELEMENTAL) THEN
C        ==================
C        CENTRALIZED MATRIX
C        ==================
         IF (id%KEEP8(28) .EQ. 0_8) THEN
           CALL MUMPS_GET_NNZ_INTERNAL(id%NNZ, id%NZ, NNZ_i)
         ELSE
           NNZ_i=id%KEEP8(28)
         ENDIF
         IF ((associated(id%A)).AND.(.NOT.PATTERN_ONLY)) THEN
C     Write header line:
               ARITH='real'
         ELSE
            ARITH='pattern '
         ENDIF
         IF (id%KEEP(50) .eq. 0) THEN
            SYMM="general"
         ELSE
            SYMM="symmetric"
         END IF
         WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix coordinate ',
     &           trim(ARITH)," ",trim(SYMM)
         WRITE(IUNIT,*) id%N, id%N, NNZ_i
         IF ((associated(id%A)).AND.(.NOT.PATTERN_ONLY)) THEN
            DO I8=1_8,NNZ_i
               IF (id%KEEP(50).NE.0 .AND. id%IRN(I8).LT.id%JCN(I8)) THEN
C              permute upper diag entry
                     WRITE(IUNIT,*) id%JCN(I8), id%IRN(I8), id%A(I8)
               ELSE
                     WRITE(IUNIT,*) id%IRN(I8), id%JCN(I8), id%A(I8)
               ENDIF
            ENDDO
         ELSE
C           pattern only
            DO I8=1_8,id%KEEP8(28)
               IF (id%KEEP(50).NE.0 .AND. id%IRN(I8).LT.id%JCN(I8)) THEN
C                 permute upper diag entry
                  WRITE(IUNIT,*) id%JCN(I8), id%IRN(I8)
               ELSE
                     WRITE(IUNIT,*) id%IRN(I8), id%JCN(I8)
               ENDIF
            ENDDO
         ENDIF
      ELSE IF ( IS_DISTRIBUTED .AND. I_AM_SLAVE ) THEN
C        ==================
C        DISTRIBUTED MATRIX
C        ==================
         IF (id%KEEP8(29) .EQ. 0_8) THEN
           CALL MUMPS_GET_NNZ_INTERNAL(id%NNZ_loc, id%NZ_loc, NNZ_i)
         ELSE
           NNZ_i=id%KEEP8(29)
         ENDIF
         IF ((associated(id%A_loc)).AND.(.NOT.PATTERN_ONLY)) THEN
               ARITH='real'
         ELSE
               ARITH='pattern '
         ENDIF
         IF (id%KEEP(50) .eq. 0) THEN
            SYMM="general"
         ELSE
            SYMM="symmetric"
         END IF
         WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix coordinate ',
     &           trim(ARITH)," ",trim(SYMM)
         WRITE(IUNIT,*) id%N, id%N, NNZ_i
         IF ((associated(id%A_loc)).AND.(.NOT.PATTERN_ONLY)) THEN
            DO I8=1_8,NNZ_i
               IF (id%KEEP(50).NE.0 .AND.
     &             id%IRN_loc(I8).LT.id%JCN_loc(I8)) THEN
                     WRITE(IUNIT,*) id%JCN_loc(I8), id%IRN_loc(I8),
     &                    id%A_loc(I8)
               ELSE
                     WRITE(IUNIT,*) id%IRN_loc(I8), id%JCN_loc(I8),
     &                    id%A_loc(I8)
               ENDIF
            ENDDO
         ELSE
            DO I8=1_8,NNZ_i
               IF (id%KEEP(50).NE.0 .AND. 
     &            id%IRN_loc(I8).LT.id%JCN_loc(I8)) THEN
C                 permute upper diag entry
                  WRITE(IUNIT,*) id%JCN_loc(I8), id%IRN_loc(I8)
               ELSE
                  WRITE(IUNIT,*) id%IRN_loc(I8), id%JCN_loc(I8)
               ENDIF
            ENDDO
         ENDIF
      ELSE IF (IS_ELEMENTAL .AND. I_AM_MASTER) THEN
C        ==================
C        ELEMENTAL MATRIX
C        ==================         
         WRITE(IUNIT,*) id%N," :: N"
         WRITE(IUNIT,*) id%NELT," :: NELT"
         WRITE(IUNIT,*) size(id%ELTVAR)," :: NELTVAR"
         WRITE(IUNIT,*) size(id%A_ELT)," :: NELTVL"
         WRITE(IUNIT,*) id%ELTPTR(:)," ::ELTPTR"
         WRITE(IUNIT,*) id%ELTVAR(:)," ::ELTVAR"
         IF(.NOT.PATTERN_ONLY) THEN
            WRITE(IUNIT,*) id%A_ELT(:)         
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_DUMP_MATRIX
      SUBROUTINE DMUMPS_DUMP_RHS(IUNIT, id)
C
C  Purpose:
C  =======
C     Dumps a dense, centralized,
C     right-hand side in matrix market format on unit
C     IUNIT. Should be called on the host only.
C
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
C  Arguments
C  =========
      TYPE(DMUMPS_STRUC), intent(in)  :: id
      INTEGER, intent(in)             :: IUNIT
C
C  Local variables
C  ===============
C
      CHARACTER (LEN=8)    :: ARITH
      INTEGER              :: I, J
      INTEGER(8)           :: LD_RHS8, K8
C
C  Executable statements
C  =====================
C
      IF (associated(id%RHS)) THEN
               ARITH='real'
        WRITE(IUNIT,FMT=*)'%%MatrixMarket matrix array ',
     &           trim(ARITH),
     &           ' general'
        WRITE(IUNIT,*) id%N, id%NRHS
        IF ( id%NRHS .EQ. 1 ) THEN
           LD_RHS8 = int(id%N,8)
        ELSE
           LD_RHS8 = int(id%LRHS,8)
        ENDIF
        DO J = 1, id%NRHS
           DO I = 1, id%N
              K8=int(J-1,8)*LD_RHS8+int(I,8)
                 WRITE(IUNIT,*) id%RHS(K8)
        ENDDO
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_DUMP_RHS
      SUBROUTINE DMUMPS_BUILD_I_AM_CAND( NSLAVES, K79, 
     &     NB_NIV2, MYID_NODES,
     &     CANDIDATES, I_AM_CAND )
      IMPLICIT NONE
C
C     Purpose:
C     =======
C     Given  a list of candidate processors per node,
C     returns an array of booleans telling whether the
C     processor is candidate or not for a given node.
C
C     K79 holds splitting strategy (KEEP(79)). If K79>1 then
C     TPYE4,5,6 nodes might have been introduced and 
C     in this case "hidden" slaves should be taken 
C     into account to enable dynamic redistribution 
C     of the hidden slaves while climbing the chain of 
C     split nodes. The master of the first node in the 
C     chain requires a special treatment and is thus here
C     not considered as a slave. 
C     
      INTEGER, intent(in) :: NSLAVES, NB_NIV2, MYID_NODES, K79
      INTEGER, intent(in) :: CANDIDATES( NSLAVES+1, NB_NIV2 )
      LOGICAL, intent(out):: I_AM_CAND( NB_NIV2 )
      INTEGER I, INIV2, NCAND
      IF (K79.GT.0) THEN
C      Because of potential restarting the number of
C      candidates that will be used to distribute 
C      arrowheads have to include all possible candidates.
       DO INIV2=1, NB_NIV2
         I_AM_CAND(INIV2)=.FALSE.
         NCAND = CANDIDATES(NSLAVES+1,INIV2)
C        check if some hidden slaves are there
C        Note that if hidden candidates exists (type 5 or 6 nodes) then
C        in position CANDIDATES (NCAND+1,INIV2) must be the master 
C        of the first node in the chain (type 4) that we skip here because
C        a special treatment (it has to be "considered as a master" for all 
C        nodes in the list) is needed.
         DO I=1, NSLAVES
            IF (CANDIDATES(I,INIV2).LT.0) EXIT ! end of extra slaves
            IF (I.EQ.NCAND+1) CYCLE 
!     skip master of associated TYPE 4 node 
            IF (CANDIDATES(I,INIV2).EQ.MYID_NODES) THEN
               I_AM_CAND(INIV2)=.TRUE.
               EXIT
            ENDIF
         ENDDO
       END DO
      ELSE
       DO INIV2=1, NB_NIV2
         I_AM_CAND(INIV2)=.FALSE.
         NCAND = CANDIDATES(NSLAVES+1,INIV2)
         DO I=1, NCAND
            IF (CANDIDATES(I,INIV2).EQ.MYID_NODES) THEN
               I_AM_CAND(INIV2)=.TRUE.
               EXIT
            ENDIF
         ENDDO
       END DO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_BUILD_I_AM_CAND
