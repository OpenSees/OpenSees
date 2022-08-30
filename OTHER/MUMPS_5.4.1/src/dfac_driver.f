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
      SUBROUTINE DMUMPS_FAC_DRIVER( id)
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      USE DMUMPS_STRUC_DEF
      USE DMUMPS_LR_STATS
      USE DMUMPS_LR_DATA_M, only: DMUMPS_BLR_INIT_MODULE, 
     &                            DMUMPS_BLR_END_MODULE
     &                          , DMUMPS_BLR_STRUC_TO_MOD
     &                          , DMUMPS_BLR_MOD_TO_STRUC
      USE MUMPS_FRONT_DATA_MGT_M
#if ! defined(NO_FDM_DESCBAND)
      USE MUMPS_FAC_DESCBAND_DATA_M
#endif
#if ! defined(NO_FDM_MAPROW)
      USE MUMPS_FAC_MAPROW_DATA_M
#endif
!$    USE OMP_LIB
C     Derived datatype to pass pointers with implicit interfaces
      USE DMUMPS_FAC_S_IS_POINTERS_M, ONLY : S_IS_POINTERS_T
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C  Performs scaling, sorting in arrowhead, then
C  distributes the matrix, and perform
C  factorization.
C
C
      INTERFACE
      SUBROUTINE DMUMPS_ANORMINF(id, ANORMINF, LSCAL)
      USE DMUMPS_STRUC_DEF
      TYPE (DMUMPS_STRUC), TARGET :: id
      DOUBLE PRECISION, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      END SUBROUTINE DMUMPS_ANORMINF
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
C  Parameters
C  ==========
C
      TYPE(DMUMPS_STRUC), TARGET :: id
C
C  MPI
C  ===
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
C
C  Local variables
C  ===============
C
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: NSEND8, NSEND_TOT8
      INTEGER(8) :: NLOCAL8, NLOCAL_TOT8
      INTEGER :: LDPTRAR, NELT_arg, NBRECORDS
      INTEGER :: ITMP
      INTEGER :: KEEP464COPY, KEEP465COPY
      INTEGER(8) :: KEEP826_SAVE
      INTEGER(8) :: K67, K68, K70, K74, K75
      INTEGER(8) ITMP8
      INTEGER  MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER MP, LP, MPG, allocok
      LOGICAL PROK, PROKG, LSCAL, LPOK, COMPUTE_ANORMINF
C     Reception buffer
      INTEGER    :: DMUMPS_LBUFR, DMUMPS_LBUFR_BYTES
      INTEGER(8) :: DMUMPS_LBUFR_BYTES8 ! for intermediate computation
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFR
C     Size of send buffers (in bytes)
      INTEGER    :: DMUMPS_LBUF, DMUMPS_LBUF_INT
      INTEGER(8) :: DMUMPS_LBUF8 ! for intermediate computation
C
      INTEGER PTRIST, PTRWB, MAXELT_SIZE,
     &     ITLOC, IPOOL, K28, LPOOL
      INTEGER IRANK, ID_ROOT
      INTEGER KKKK
      INTEGER(8) :: NZ_locMAX8
      INTEGER(8) MEMORY_MD_ARG
      INTEGER(8) MAXS_BASE8, MAXS_BASE_RELAXED8
      DOUBLE PRECISION CNTL4, AVG_FLOPS
      INTEGER MIN_PERLU, MAXIS_ESTIM
C
      TYPE (S_IS_POINTERS_T) :: S_IS_POINTERS
      INTEGER   MAXIS
      INTEGER(8) :: MAXS
C     For S argument to arrowhead routines:
      INTEGER(8) :: MAXS_ARG
      DOUBLE PRECISION, TARGET :: S_DUMMY_ARG(1)
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: S_PTR_ARG
      INTEGER NPIV_CRITICAL_PATH
      DOUBLE PRECISION TIME, TIMEET
      DOUBLE PRECISION ZERO, ONE, MONE
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0, MONE = -1.0D0)
      DOUBLE PRECISION CZERO
      PARAMETER( CZERO = 0.0D0 )
      INTEGER PERLU, TOTAL_MBYTES, K231, K232, K233, BLR_STRAT
      INTEGER, PARAMETER :: IDUMMY = -9999
      LOGICAL, PARAMETER :: BDUMMY =.FALSE.
      INTEGER COLOUR, COMM_FOR_SCALING ! For Simultaneous scaling
      INTEGER LIWK, LWK_REAL
      INTEGER(8) :: LWK
C     I_AM_SLAVE: used to determine if proc has the role of a slave
C     WK_USER_PROVIDED is set to true when WK_USER is provided by user
      LOGICAL I_AM_SLAVE, PERLU_ON, WK_USER_PROVIDED, EARLYT3ROOTINS
      LOGICAL PRINT_MAXAVG
      DOUBLE PRECISION :: ANORMINF, SEUIL, SEUIL_LDLT_NIV2, Thresh_Seuil
      DOUBLE PRECISION :: CNTL1, CNTL3, CNTL5, CNTL6, EPS
      INTEGER N, LPN_LIST,POSBUF
      INTEGER, DIMENSION (:), ALLOCATABLE :: ITMP2
      INTEGER I,K
      INTEGER(8) :: ITEMP8
      INTEGER    :: PARPIV_T1
      INTEGER FRONTWISE
C temporary variables for collecting stats from all processors
      DOUBLE PRECISION :: TMP_MRY_LU_FR
      DOUBLE PRECISION :: TMP_MRY_LU_LRGAIN
      DOUBLE PRECISION :: TMP_MRY_CB_FR
      DOUBLE PRECISION :: TMP_MRY_CB_LRGAIN
      DOUBLE PRECISION :: TMP_FLOP_LRGAIN
      DOUBLE PRECISION :: TMP_FLOP_TRSM
      DOUBLE PRECISION :: TMP_FLOP_PANEL
      DOUBLE PRECISION :: TMP_FLOP_FRFRONTS
      DOUBLE PRECISION :: TMP_FLOP_TRSM_FR
      DOUBLE PRECISION :: TMP_FLOP_TRSM_LR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_FR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_LR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_LRLR3
      DOUBLE PRECISION :: TMP_FLOP_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_DECOMPRESS
      DOUBLE PRECISION :: TMP_FLOP_MIDBLK_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_FRSWAP_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_ACCUM_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_CB_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_CB_DECOMPRESS
      DOUBLE PRECISION :: TMP_FLOP_FACTO_FR
      DOUBLE PRECISION :: TMP_FLOP_SOLFWD_FR
      DOUBLE PRECISION :: TMP_FLOP_SOLFWD_LR
      INTEGER :: TMP_CNT_NODES
      DOUBLE PRECISION :: TMP_TIME_UPDATE
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR1
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR2
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR3
      DOUBLE PRECISION :: TMP_TIME_UPDATE_FRLR
      DOUBLE PRECISION :: TMP_TIME_UPDATE_FRFR
      DOUBLE PRECISION :: TMP_TIME_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_MIDBLK_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_FRSWAP_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_CB_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_PANEL
      DOUBLE PRECISION :: TMP_TIME_FAC_I
      DOUBLE PRECISION :: TMP_TIME_FAC_MQ
      DOUBLE PRECISION :: TMP_TIME_FAC_SQ
      DOUBLE PRECISION :: TMP_TIME_LRTRSM
      DOUBLE PRECISION :: TMP_TIME_FRTRSM
      DOUBLE PRECISION :: TMP_TIME_FRFRONTS
      DOUBLE PRECISION :: TMP_TIME_LR_MODULE
      DOUBLE PRECISION :: TMP_TIME_DIAGCOPY
      DOUBLE PRECISION :: TMP_TIME_DECOMP
      DOUBLE PRECISION :: TMP_TIME_DECOMP_UCFS
      DOUBLE PRECISION :: TMP_TIME_DECOMP_ASM1
      DOUBLE PRECISION :: TMP_TIME_DECOMP_LOCASM2
      DOUBLE PRECISION :: TMP_TIME_DECOMP_MAPLIG1
      DOUBLE PRECISION :: TMP_TIME_DECOMP_ASMS2S
      DOUBLE PRECISION :: TMP_TIME_DECOMP_ASMS2M
C
C  Workspace.
C
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WK_REAL
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IWK8
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURS
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCS
      INTEGER BUREGISTRE(12)
      INTEGER BUINTSZ, BURESZ, BUJOB
      INTEGER BUMAXMN, M, SCMYID, SCNPROCS
      DOUBLE PRECISION    SCONEERR, SCINFERR
C
C  Parameters arising from the structure
C  =====================================
C
      INTEGER, POINTER ::  JOB
*     Control parameters: see description in DMUMPSID
      DOUBLE PRECISION,DIMENSION(:),POINTER::RINFO, RINFOG
      DOUBLE PRECISION,DIMENSION(:),POINTER::    CNTL
      INTEGER,DIMENSION(:),POINTER:: INFOG, KEEP
      INTEGER, DIMENSION(:), POINTER :: MYIRN_loc, MYJCN_loc
      DOUBLE PRECISION, DIMENSION(:), POINTER :: MYA_loc
      INTEGER, TARGET :: DUMMYIRN_loc(1), DUMMYJCN_loc(1)
      DOUBLE PRECISION, TARGET :: DUMMYA_loc(1)
      INTEGER,DIMENSION(:),POINTER::ICNTL
      EXTERNAL MUMPS_GET_POOL_LENGTH
      INTEGER MUMPS_GET_POOL_LENGTH
      INTEGER(8) :: TOTAL_BYTES
      INTEGER(8) :: I8TMP, LWK_USER_SUM8
C
C  External references
C  ===================
      INTEGER numroc
      EXTERNAL numroc
      INTEGER:: NWORKING
      LOGICAL:: MEM_EFF_ALLOCATED
C  Fwd in facto:
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_MUMPS
      LOGICAL :: RHS_MUMPS_ALLOCATED
      INTEGER :: NB_ACTIVE_FRONTS_ESTIM
      INTEGER :: NB_FRONTS_F_ESTIM
C
C 
      JOB=>id%JOB
      RINFO=>id%RINFO
      RINFOG=>id%RINFOG
      CNTL=>id%CNTL
      INFOG=>id%INFOG
      KEEP=>id%KEEP
      ICNTL=>id%ICNTL
      IF (id%KEEP8(29) .NE. 0) THEN
        MYIRN_loc=>id%IRN_loc
        MYJCN_loc=>id%JCN_loc
        MYA_loc=>id%A_loc
      ELSE
        MYIRN_loc=>DUMMYIRN_loc
        MYJCN_loc=>DUMMYJCN_loc
        MYA_loc=>DUMMYA_loc
      ENDIF
      N = id%N
      EPS = epsilon ( ZERO )
C     TIMINGS: reset to 0
      id%DKEEP(92)=0.0D0
      id%DKEEP(93)=0.0D0
      id%DKEEP(94)=0.0D0
      id%DKEEP(97)=0.0D0
      id%DKEEP(98)=0.0D0
      id%DKEEP(56)=0.0D0
C     Count of MPI messages: reset to 0
      id%KEEP(266)=0
      id%KEEP(267)=0
C     MIN/MAX pivots reset to 0
      id%DKEEP(19)=huge(0.0D0)
      id%DKEEP(20)=huge(0.0D0)
      id%DKEEP(21)=0.0D0
C     Number of symmetric swaps
      id%KEEP8(80)=0_8
C     Largest increase of internal panel size
      id%KEEP(425) =0
C
      PRINT_MAXAVG = .NOT.(id%NSLAVES.EQ.1 .AND. KEEP(46).EQ.1)
C
C     BEGIN CASE OF ALLOCATED DATA FROM PREVIOUS CALLS
C     Data from factorization is now freed asap
C     id%S, id%IS
      IF (id%KEEP8(24).EQ.0_8) THEN
C       -- deallocate only when not provided/allocated by the user
         IF (associated(id%S)) THEN
            DEALLOCATE(id%S)
            id%KEEP8(23)=0_8
            NULLIFY(id%S)
         ENDIF
      ENDIF
      IF (associated(id%IS)) THEN
        DEALLOCATE(id%IS)
        NULLIFY(id%IS)
      ENDIF
C     Free BLR factors, if any
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
C
C     END CASE OF ALLOCATED DATA FROM PREVIOUS CALLS
C
C     Related to forward in facto functionality (referred to as "Fwd in facto")
      NULLIFY(RHS_MUMPS)
      RHS_MUMPS_ALLOCATED = .FALSE.
C     -----------------------------------------------------------------------
C     Set WK_USER_PROVIDED to true when workspace WK_USER is provided by user
C     We can accept WK_USER to be provided on only one proc and 
C     different values of WK_USER per processor
C     
      IF (id%KEEP8(24).GT.0_8) THEN 
C                We nullify S so that later when we test
C                if (associated(S) we can free space and reallocate it).
           NULLIFY(id%S)
      ENDIF
C
C     --  KEEP8(24) can now then be reset safely
      WK_USER_PROVIDED = (id%LWK_USER.NE.0)
      IF (WK_USER_PROVIDED) THEN
          IF (id%LWK_USER.GT.0) THEN
            id%KEEP8(24) = int(id%LWK_USER,8)
          ELSE
            id%KEEP8(24) = -int(id%LWK_USER,8)* 1000000_8 
          ENDIF
      ELSE
          id%KEEP8(24) = 0_8
      ENDIF
C     Compute sum of LWK_USER provided by user 
      LWK_USER_SUM8 = 0_8
      CALL MPI_REDUCE ( id%KEEP8(24), LWK_USER_SUM8, 1, MPI_INTEGER8,
     &                  MPI_SUM, MASTER, id%COMM, IERR )
C
C     KEEP8(26) might be modified
C       (element entry format) 
C       but need be restore for 
C       future factorisation
C       with different scaling option
C 
      KEEP826_SAVE = id%KEEP8(26)
C     In case of loop on factorization with
C     different scaling options, initialize
C     DKEEP(4:5) to 0.
      id%DKEEP(4)=-1.0D0
      id%DKEEP(5)=-1.0D0
C  Mapping information used during solve. In case of several facto+solve
C  it has to be recomputed. In case of several solves with the same
C  facto, it is not recomputed.
      IF (associated(id%IPTR_WORKING)) THEN
        DEALLOCATE(id%IPTR_WORKING)
        NULLIFY(id%IPTR_WORKING)
      END IF
      IF (associated(id%WORKING)) THEN 
        DEALLOCATE(id%WORKING)
        NULLIFY(id%WORKING)
      END IF
C
C  Units for printing
C  MP: diagnostics
C  LP: errors
C
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LPOK    = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.2))
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG   = (PROKG.AND.(id%ICNTL(4).GE.2))
      IF ( PROK ) WRITE( MP, 130 )
      IF ( PROKG ) WRITE( MPG, 130 )
C     -------------------------------------
C     Depending on the type of parallelism,
C     the master can now (soon) potentially
C     have the role of a slave
C     -------------------------------------
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
C
C  Prepare work for out-of-core
C
      IF (id%MYID .EQ. MASTER .AND. KEEP(201) .NE. -1) THEN
C       Note that if KEEP(201)=-1, then we have decided
C       at analysis phase that factors will not be stored
C       (neither in memory nor on disk). In that case,
C       ICNTL(22) is ignored.
C       -- ICNTL(22) must be set before facto phase 
C          (=1 OOC on; =0 OOC off)
C          and cannot be changed for subsequent solve phases.
        KEEP(201)=id%ICNTL(22)
        IF (KEEP(201) .NE. 0) THEN
#         if defined(OLD_OOC_NOPANEL)
            KEEP(201)=2
#         else
            KEEP(201)=1
#         endif
        ENDIF
      ENDIF
C     ----------------------
C     Broadcast KEEP options
C     defined for facto:
C     ----------------------
      CALL MPI_BCAST( KEEP(12), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(19), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(21), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      PERLU = KEEP(12)
      IF (id%MYID.EQ.MASTER) THEN
C       KEEP(50)  case
C       ==============
C
C       KEEP(50)  = 0 : matrix is unsymmetric
C       KEEP(50) /= 0 : matrix is symmetric
C       KEEP(50) = 1 : Ask L L^T on the root. Matrix is PSD.
C       KEEP(50) = 2 : Ask for L U on the root
C       KEEP(50) = 3 ... L D L^T ??
C     
        CNTL1 = id%CNTL(1)
C       ---------------------------------------
C       For symmetric (non general) matrices 
C       set (directly) CNTL1 = 0.0
C       ---------------------------------------
        KEEP(17)=0
        IF ( KEEP(50) .eq. 1 ) THEN
          IF (CNTL1 .ne. ZERO ) THEN
            IF ( PROKG ) THEN
              WRITE(MPG,'(A)')
     & '** Warning : SPD solver called, resetting CNTL(1) to 0.0D0'
            END IF
          END IF
          CNTL1 = ZERO
        END IF
C       CNTL1 threshold value must be between 
C       0.0 and 1.0 (for SYM=0) and 0.5 (for SYM=1,2)
        IF (CNTL1.GT.ONE)   CNTL1=ONE
        IF (CNTL1.LT.ZERO)  CNTL1=ZERO
        IF (KEEP(50).NE.0.AND.CNTL1.GT.0.5D0) THEN
          CNTL1 = 0.5D0
        ENDIF
        PARPIV_T1 = id%KEEP(268)
        IF (PARPIV_T1.EQ.77) THEN
         PARPIV_T1 = 0 
        ENDIF
        IF (PARPIV_T1.EQ.-3) THEN
          PARPIV_T1 = 0
        ENDIF
        IF ((PARPIV_T1.LT.-3).OR.(PARPIV_T1.GT.1)) THEN
C        out of range values
         PARPIV_T1 =0
        ENDIF
C       note that KEEP(50).EQ.1 => CNTL1=0.0
        IF (CNTL1.EQ.0.0.OR.(KEEP(50).eq.1)) PARPIV_T1 = 0 
C
        IF (PARPIV_T1.EQ.-2) THEN
         IF (KEEP(19).NE.0) THEN
C         switch off PARPIV_T1 if RR activated
C         but do NOT switch off PARPIV_1 with null pivot detection
          PARPIV_T1 = 0
         ENDIF
        ENDIF
        id%KEEP(269) = PARPIV_T1
      ENDIF
      CALL MPI_BCAST(CNTL1, 1, MPI_DOUBLE_PRECISION,
     &             MASTER, id%COMM, IERR)
        CALL MPI_BCAST( KEEP(269), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (id%MYID.EQ.MASTER) THEN
C         -----------------------------------------------------
C         Decoding of ICNTL(35) for factorization: same as
C         at analysis except that we store a copy of ICNTL(35)
C         in KEEP(486) instead of KEEP(494) and need to check
C         compatibility of KEEP(486) and KEEP(494): If LR was
C         not activated during analysis, it cannot be activated 
C         at factorization.
C         ------------------------------------------------------
          id%KEEP(486) = id%ICNTL(35)
          IF (id%KEEP(486).EQ.1) THEN
C         -- Automatic BLR option setting
           id%KEEP(486)= 2
          ENDIF
          IF ( id%KEEP(486).EQ.4) id%KEEP(486)=0
          IF ((id%KEEP(486).LT.0).OR.(id%KEEP(486).GT.4)) THEN
C           Out of range values treated as 0
            id%KEEP(486) = 0
          ENDIF
          IF ((KEEP(486).NE.0).AND.(KEEP(494).EQ.0)) THEN
C           To activate BLR during factorization,
C           ICNTL(35) must have been set at analysis.
            IF (LPOK)  THEN
              WRITE(LP,'(A)') 
     &      " *** Error with BLR setting "
              WRITE(LP,'(A)') " *** BLR was not activated during ",
     &      " analysis but is requested during factorization."
            ENDIF
            id%INFO(1)=-54
            id%INFO(2)=0
            GOTO 105
          ENDIF
          KEEP464COPY  = id%ICNTL(38)
          IF (KEEP464COPY.LT.0.OR.KEEP464COPY.GT.1000) THEN
C          Out of range values treated as 0
           KEEP464COPY = 0
          ENDIF
          IF (id%KEEP(461).LT.1) THEN
            id%KEEP(461) = 10
          ENDIF
          KEEP465COPY=0
          IF (id%ICNTL(36).EQ.1.OR.id%ICNTL(36).EQ.3) THEN
            IF (CNTL1.EQ.ZERO .OR. KEEP(468).LE.1) THEN
              KEEP(475) = 3
            ELSE IF ( (KEEP(269).GT.0).OR. (KEEP(269).EQ.-2)) THEN
              KEEP(475) = 2
            ELSE IF (KEEP(468).EQ.2) THEN
              KEEP(475) = 2
            ELSE
              KEEP(475) = 1
            ENDIF
          ELSE
            KEEP(475) = 0
          ENDIF
          KEEP(481)=0
          IF (id%ICNTL(36).LT.0 .OR. id%ICNTL(36).GE.2) THEN
C           Only options 1 and 2 are allowed
            KEEP(475) = 0
          ENDIF
C         K489 is set according to ICNTL(37)
          IF (id%ICNTL(37).EQ.0.OR.id%ICNTL(37).EQ.1) THEN
            KEEP(489) = id%ICNTL(37)
          ELSE
C           Other values treated as zero
            KEEP(489) = 0
          ENDIF
          IF (KEEP(79).GE.1) THEN
C          CompressCB incompatible with type4,5,6 nodes
           KEEP(489)=0
          ENDIF
            KEEP(489)=0
C         id%KEEP(476) \in [1,100] 
          IF ((id%KEEP(476).GT.100).OR.(id%KEEP(476).LT.1)) THEN
              id%KEEP(476)=  50
          ENDIF
C         id%KEEP(477) \in [1,100] 
          IF ((id%KEEP(477).GT.100).OR.(id%KEEP(477).LT.1)) THEN
              id%KEEP(477)=  100
          ENDIF
C         id%KEEP(483) \in [1,100] 
          IF ((id%KEEP(483).GT.100).OR.(id%KEEP(483).LT.1)) THEN
              id%KEEP(483)=  50
          ENDIF
C         id%KEEP(484) \in [1,100] 
          IF ((id%KEEP(484).GT.100).OR.(id%KEEP(484).LT.1)) THEN
              id%KEEP(484)=  50
          ENDIF
C         id%KEEP(480)=0,2,3,4,5,6
          IF ((id%KEEP(480).GT.6).OR.(id%KEEP(480).LT.0)
     &                           .OR.(id%KEEP(480).EQ.1)) THEN
               id%KEEP(480)=0
          ENDIF
C         id%KEEP(473)=0 or 1
          IF ((id%KEEP(473).NE.0).AND.(id%KEEP(473).NE.1)) THEN
               id%KEEP(473)=0
          ENDIF
C         id%KEEP(474)=0,1,2,3
          IF ((id%KEEP(474).GT.3).OR.(id%KEEP(474).LT.0)) THEN
               id%KEEP(474)=0
          ENDIF
C         id%KEEP(479)>0
          IF (id%KEEP(479).LE.0) THEN
               id%KEEP(479)=1
          ENDIF
         IF (id%KEEP(474).NE.0.AND.id%KEEP(480).EQ.0) THEN
              id%KEEP(474) = 0
          ENDIF
          IF (id%KEEP(478).NE.0.AND.id%KEEP(480).LT.4) THEN
              id%KEEP(478) = 0
          ENDIF
          IF (id%KEEP(480).GE.5 .OR.
     &           (id%KEEP(480).NE.0.AND.id%KEEP(474).EQ.3)) THEN
            IF (id%KEEP(475).LT.2) THEN
C             Reset to 3 if 5 or to 4 if 6
              id%KEEP(480) = id%KEEP(480) - 2
              write(*,*) ' Resetting KEEP(480) to ', id%KEEP(480)
            ENDIF
          ENDIF
 105    CONTINUE
        ENDIF  ! id%MYID .EQ. MASTER
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C
        IF (id%INFO(1).LT.0) GOTO 530
        CALL MPI_BCAST( KEEP(473), 14, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (KEEP(486).NE.0) THEN
          CALL MPI_BCAST( KEEP(489), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP464COPY, 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP465COPY, 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
         ENDIF
        IF (id%MYID.EQ.MASTER) THEN
          IF (KEEP(217).GT.2.OR.KEEP(217).LT.0) THEN
            KEEP(217)=0
          ENDIF
          KEEP(214)=KEEP(217)
          IF (KEEP(214).EQ.0) THEN
            IF (KEEP(201).NE.0) THEN ! OOC or no factors
              KEEP(214)=1
            ELSE
              KEEP(214)=2
            ENDIF
            IF (KEEP(486).EQ.2) THEN
              KEEP(214)=1
            ENDIF
          ENDIF
        ENDIF
        CALL MPI_BCAST( KEEP(214), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (KEEP(201).NE.0) THEN
C         -- Low Level I/O strategy
          CALL MPI_BCAST( KEEP(99), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(205), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(211), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        ENDIF
C     Fwd in facto: explicitly forbid
C     sparse RHS and A-1 computation
      IF (id%KEEP(252).EQ.1 .AND. id%MYID.EQ.MASTER) THEN
        IF (id%ICNTL(20).EQ.1) THEN ! out-of-range => 0
C         NB: in doc ICNTL(20) only accessed during solve
C         In practice, will have failed earlier if RHS not allocated.
C         Still it looks safer to keep this test.
          id%INFO(1)=-43
          id%INFO(2)=20
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: Sparse RHS is incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(30).NE.0) THEN ! out-of-range => 1
          id%INFO(1)=-43
          id%INFO(2)=30
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: A-1 functionality incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(9) .NE. 1) THEN
          id%INFO(1)=-43
          id%INFO(2)=9
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: Transpose system (ICNTL(9).NE.0) not ',
     &       ' compatible with forward performed during',
     &       ' factorization (ICNTL(32)=1)'
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C
      IF (id%INFO(1).LT.0) GOTO 530
C
C     The memory allowed is given by ICNTL(23) in Mbytes
C     0 means that nothing is provided.
C     Save memory available, ICNTL(23) in KEEP8(4)
C
      IF ( id%MYID.EQ.MASTER ) THEN
        ITMP = ICNTL(23)
      END IF
      CALL MPI_BCAST( ITMP, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C    
C     Ignore ICNTL(23) when WK_USER is provided
c     by resetting ITMP to zero on each proc where WK_USER is provided
      IF (WK_USER_PROVIDED) ITMP = 0
      ITMP8 = int(ITMP, 8)
      id%KEEP8(4) = ITMP8 * 1000000_8   ! convert to nb of bytes
      IF ( PROKG ) THEN
          NWORKING = id%NSLAVES
          WRITE( MPG, 172 ) NWORKING, id%ICNTL(22), KEEP(486),
     &    KEEP(12), 
     &    id%KEEP8(111), KEEP(126), KEEP(127), KEEP(28),
     &    id%KEEP8(4)/1000000_8, LWK_USER_SUM8, CNTL1
          IF (KEEP(252).GT.0) 
     &    WRITE(MPG,173) KEEP(253)
          IF (KEEP(269).NE.0) 
     &    WRITE(MPG,174) KEEP(269)
      ENDIF
      IF (KEEP(201).LE.0) THEN
C       In-core version or no factors
        KEEP(IXSZ)=XSIZE_IC
      ELSE IF (KEEP(201).EQ.2) THEN
C       OOC version, no panels
        KEEP(IXSZ)=XSIZE_OOC_NOPANEL
      ELSE IF (KEEP(201).EQ.1) THEN
C     Panel versions:
        IF (KEEP(50).EQ.0) THEN
          KEEP(IXSZ)=XSIZE_OOC_UNSYM
        ELSE
          KEEP(IXSZ)=XSIZE_OOC_SYM
        ENDIF
      ENDIF
      IF ( KEEP(486) .NE. 0 ) THEN !LR is activated
C       Stats initialization for LR
        CALL INIT_STATS_GLOBAL(id)
       END IF
C
*     **********************************
*     Begin intializations regarding the
*     computation of the determinant
*     **********************************
      IF (id%MYID.EQ.MASTER) KEEP(258)=ICNTL(33)
      CALL MPI_BCAST(KEEP(258), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      IF (KEEP(258) .NE. 0) THEN
        KEEP(259) = 0      ! Initial exponent of the local determinant
        KEEP(260) = 1      ! Number of permutations
        id%DKEEP(6)  = 1.0D0  ! real part of the local determinant
      ENDIF
*     ********************************
*     End intializations regarding the
*     computation of the determinant
*     ********************************
C
*     **********************
*     Begin of Scaling phase
*     **********************
C
C     SCALING MANAGEMENT
C     * Options 1, 3, 4 centralized only
C  
C     * Options 7, 8  : also works for distributed matrix
C
C     At this point, we have the scaling arrays allocated
C     on the master. They have been allocated on the master
C     inside the main MUMPS driver.
C
      CALL MPI_BCAST(KEEP(52), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      LSCAL = ((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 8))
      IF (LSCAL) THEN
C
        IF ( id%MYID.EQ.MASTER ) THEN
          CALL MUMPS_SECDEB(TIMEET)
        ENDIF
C       -----------------------
C       Retrieve parameters for
C       simultaneous scaling
C       -----------------------
        IF (KEEP(52) .EQ. 7) THEN 
C       -- Cheap setting of SIMSCALING (it is the default in 4.8.4)
           K231= KEEP(231)
           K232= KEEP(232)
           K233= KEEP(233)
        ELSEIF (KEEP(52) .EQ. 8) THEN
C       -- More expensive setting of SIMSCALING (it was the default in 4.8.1,2,3)
           K231= KEEP(239)
           K232= KEEP(240)
           K233= KEEP(241)
        ENDIF
        CALL MPI_BCAST(id%DKEEP(3),1,MPI_DOUBLE_PRECISION,MASTER,
     &       id%COMM,IERR)
C
        IF ( ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) .AND. 
     &       KEEP(54).NE.0 ) THEN
C         ------------------------------
C         Scaling for distributed matrix
C         We need to allocate scaling
C         arrays on all processors, not
C         only the master.
C         ------------------------------
           IF ( id%MYID .NE. MASTER ) THEN
              IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
              IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
         ENDIF
         M = N
         BUMAXMN=M
         IF(N > BUMAXMN) BUMAXMN = N
         LIWK = 4*BUMAXMN
         ALLOCATE (IWK(LIWK),BURP(M),BUCP(N),
     &            BURS(2* (id%NPROCS)),BUCS(2* (id%NPROCS)),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            id%INFO(2)=LIWK+M+N+4* (id%NPROCS)
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
C        -- estimation of memory and construction of partvecs
         BUJOB = 1
C        -- LWK not used
         LWK_REAL   = 1
         ALLOCATE(WK_REAL(LWK_REAL),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            id%INFO(2)=LWK_REAL
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
         CALL DMUMPS_SIMSCALEABS(
     &        MYIRN_loc(1), MYJCN_loc(1), MYA_loc(1),
     &        id%KEEP8(29),
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         IF(LIWK < BUINTSZ) THEN
            DEALLOCATE(IWK)
            LIWK = BUINTSZ
            ALLOCATE(IWK(LIWK), stat=allocok)
            IF (allocok > 0) THEN
               id%INFO(1)=-13
               id%INFO(2)=LIWK
            ENDIF
         ENDIF
         LWK_REAL = BURESZ
         DEALLOCATE(WK_REAL)
         ALLOCATE (WK_REAL(LWK_REAL), stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            id%INFO(2)=LWK_REAL
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
C        -- estimation of memory and construction of partvecs
         BUJOB = 2
         CALL DMUMPS_SIMSCALEABS(
     &        MYIRN_loc(1), MYJCN_loc(1), MYA_loc(1),
     &        id%KEEP8(29),
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         id%DKEEP(4) = SCONEERR
         id%DKEEP(5) = SCINFERR
CXXXX 
         DEALLOCATE(IWK, WK_REAL,BURP,BUCP,BURS, BUCS)
        ELSE IF ( KEEP(54) .EQ. 0 ) THEN
C         ------------------
C         Centralized matrix
C         ------------------
          IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8))  THEN
C             -------------------------------
C             Create a communicator of size 1
C             -------------------------------
              IF (id%MYID.EQ.MASTER) THEN
                COLOUR = 0
              ELSE
                COLOUR = MPI_UNDEFINED
              ENDIF
              CALL MPI_COMM_SPLIT( id%COMM, COLOUR, 0,
     &             COMM_FOR_SCALING, IERR )
              IF (id%MYID.EQ.MASTER) THEN
                 M = N
                 BUMAXMN=N
CXXXX 
                 IF(N > BUMAXMN) BUMAXMN = N
                 LIWK = 1
                 ALLOCATE (IWK(LIWK),BURP(1),BUCP(1),
     &                BURS(1),BUCS(1),
     &                stat=allocok)
                 IF (allocok > 0) THEN
                    id%INFO(1)=-13
                    id%INFO(2)=LIWK+1+1+1+1
                 ENDIF
                 LWK_REAL = M + N  
                 ALLOCATE (WK_REAL(LWK_REAL), stat=allocok)
                 IF (allocok > 0) THEN
                    id%INFO(1)=-13
                    id%INFO(2)=1
                 ENDIF
                 IF (id%INFO(1) .LT. 0) GOTO 400
                 CALL MPI_COMM_RANK(COMM_FOR_SCALING, SCMYID, IERR)
                 CALL MPI_COMM_SIZE(COMM_FOR_SCALING, SCNPROCS, IERR)
                 BUJOB = 1
                 CALL DMUMPS_SIMSCALEABS(
     &                id%IRN(1), id%JCN(1), id%A(1),
     &                id%KEEP8(28),
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 IF(LWK_REAL < BURESZ) THEN
                    ! internal error since LWK_REAL=BURESZ=M+N
                    id%INFO(1) = -136
                    GOTO 400
                 ENDIF
                 BUJOB = 2
                 CALL DMUMPS_SIMSCALEABS(id%IRN(1),
     &                id%JCN(1), id%A(1),
     &                id%KEEP8(28),
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 id%DKEEP(4) = SCONEERR
                 id%DKEEP(5) = SCINFERR
CXXXX 
                 DEALLOCATE(WK_REAL)                 
                 DEALLOCATE (IWK,BURP,BUCP,
     &                BURS,BUCS)
              ENDIF
C             Centralized matrix: make DKEEP(4:5) available to all processors
              CALL MPI_BCAST( id%DKEEP(4),2,MPI_DOUBLE_PRECISION,
     &                        MASTER, id%COMM, IERR )
  400         CONTINUE
              IF (id%MYID.EQ.MASTER) THEN
C               Communicator should only be
C               freed on the master process
                CALL MPI_COMM_FREE(COMM_FOR_SCALING, IERR)
              ENDIF
              CALL MUMPS_PROPINFO(ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID)
              IF (id%INFO(1).LT.0) GOTO 517
          ELSE IF (id%MYID.EQ.MASTER) THEN
C           ----------------------------------
C           Centralized scaling, options 1 to 6
C           ----------------------------------
            IF (KEEP(52).GT.0 .AND. KEEP(52).LE.6) THEN
C             ---------------------
C             Allocate temporary
C             workspace for scaling
C             ---------------------
              IF ( KEEP(52) .eq. 5 .or. 
     &          KEEP(52) .eq. 6 ) THEN
C               We have an explicit copy of the original
C               matrix in complex format which should probably
C               be avoided (but do we want to keep all
C               those old scaling options ?)
                LWK = id%KEEP8(28)
              ELSE
                LWK = 1_8
              END IF
              LWK_REAL = 5 * N
              ALLOCATE( WK_REAL( LWK_REAL ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                id%INFO(2) = LWK_REAL
                GOTO 137
              END IF
              ALLOCATE( WK( LWK ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
                GOTO 137
              END IF
              CALL DMUMPS_FAC_A(N, id%KEEP8(28), KEEP(52), id%A(1),
     &             id%IRN(1), id%JCN(1),
     &             id%COLSCA(1), id%ROWSCA(1),
     &             WK, LWK, WK_REAL, LWK_REAL, ICNTL(1), id%INFO(1) )
              DEALLOCATE( WK_REAL )
              DEALLOCATE( WK )
            ENDIF
          ENDIF
        ENDIF ! Scaling distributed matrices or centralized
        IF (id%MYID.EQ.MASTER) THEN
            CALL MUMPS_SECFIN(TIMEET)
            id%DKEEP(92)=TIMEET
C         Print inf-norm after last KEEP(233) iterations of
C         scaling option KEEP(52)=7 or 8 (SimScale)
C
          IF (PROKG.AND.(KEEP(52).EQ.7.OR.KEEP(52).EQ.8) 
     &             .AND. (K233+K231+K232).GT.0) THEN
           IF (K232.GT.0) WRITE(MPG, 166) id%DKEEP(4)
          ENDIF
        ENDIF
      ENDIF ! LSCAL
C
C       scaling might also be provided by the user
        LSCAL = (LSCAL .OR. (KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        IF (LSCAL .AND. KEEP(258).NE.0 .AND. id%MYID .EQ. MASTER) THEN
          DO I = 1, id%N
            CALL DMUMPS_UPDATEDETER_SCALING(id%ROWSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
          ENDDO
          IF (KEEP(50) .EQ. 0) THEN ! unsymmetric
            DO I = 1, id%N
              CALL DMUMPS_UPDATEDETER_SCALING(id%COLSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
            ENDDO
          ELSE
C           -----------------------------------------
C           In this case COLSCA = ROWSCA
C           Since determinant was initialized to 1,
C           compute square of the current determinant
C           rather than going through COLSCA.
C           -----------------------------------------
            CALL DMUMPS_DETER_SQUARE(id%DKEEP(6), KEEP(259))
          ENDIF
C         Now we should have taken the
C         inverse of the scaling vectors
          CALL DMUMPS_DETER_SCALING_INVERSE(id%DKEEP(6), KEEP(259))
        ENDIF
C
C       ********************
C       End of Scaling phase
C       At this point: either (matrix is distributed and KEEP(52)=7 or 8)
C       in which case scaling arrays are allocated on all processors,
C       or scaling arrays are only on the host processor.
C       In case of distributed matrix input, we will free the scaling
C       arrays on procs with MYID .NE. 0 after the all-to-all distribution
C       of the original matrix.
C       ********************
C
 137  CONTINUE
C     Fwd in facto: in case of repeated factorizations
C     with different Schur options we prefer to free
C     systematically this array now than waiting for
C     the root node. We rely on the fact that it is
C     allocated or not during the solve phase so if
C     it was allocated in a 1st call to facto and not
C     in a second, we don't want the solve to think
C     it was allocated in the second call.
      IF (associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE (id%root%RHS_CNTR_MASTER_ROOT)
        NULLIFY (id%root%RHS_CNTR_MASTER_ROOT)
      ENDIF
C     Fwd in facto: check that id%NRHS has not changed
      IF ( id%MYID.EQ.MASTER.AND. KEEP(252).EQ.1 .AND.
     &      id%NRHS .NE. id%KEEP(253) ) THEN
C         Error: NRHS should not have
C         changed since the analysis
          id%INFO(1)=-42
          id%INFO(2)=id%KEEP(253)
      ENDIF
C     Fwd in facto: allocate and broadcast RHS_MUMPS
C     to make it available on all processors.
      IF (id%KEEP(252) .EQ. 1) THEN
          IF ( id%MYID.NE.MASTER ) THEN
            id%KEEP(254) = N              ! Leading dimension
            id%KEEP(255) = N*id%KEEP(253) ! Tot size
            ALLOCATE(RHS_MUMPS(id%KEEP(255)),stat=IERR)
            IF (IERR > 0) THEN
               id%INFO(1)=-13
               id%INFO(2)=id%KEEP(255)
               IF (LPOK)
     &         WRITE(LP,*) 'ERROR while allocating RHS on a slave'
               NULLIFY(RHS_MUMPS)
            ENDIF
            RHS_MUMPS_ALLOCATED = .TRUE.
          ELSE 
C           Case of non working master
            id%KEEP(254)=id%LRHS              ! Leading dimension
            id%KEEP(255)=id%LRHS*(id%KEEP(253)-1)+id%N ! Tot size
            RHS_MUMPS=>id%RHS
            RHS_MUMPS_ALLOCATED = .FALSE.
            IF (LSCAL) THEN
C             Scale before broadcast: apply row
C             scaling (remark that we assume no
C             transpose).
              DO K=1, id%KEEP(253)
                DO I=1, N
                  RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          = RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          * id%ROWSCA(I)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
      ELSE
          id%KEEP(255)=1
          ALLOCATE(RHS_MUMPS(1),stat=IERR)
          IF (IERR > 0) THEN
             id%INFO(1)=-13
             id%INFO(2)=1
             IF (LPOK)
     &            WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
             NULLIFY(RHS_MUMPS)
          ENDIF
          RHS_MUMPS_ALLOCATED = .TRUE.
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 517
      IF (KEEP(252) .EQ. 1) THEN
C
C         Broadcast the columns of the right-hand side
C         one by one. Leading dimension is keep(254)=N
C         on procs with MYID > 0 but may be larger on
C         the master processor.
          DO I= 1, id%KEEP(253)
            CALL MPI_BCAST(RHS_MUMPS((I-1)*id%KEEP(254)+1), N,
     &           MPI_DOUBLE_PRECISION, MASTER,id%COMM,IERR)
          END DO
      ENDIF
C     Keep a copy of ICNTL(24) and make it
C     available on all working processors.
      KEEP(110)=id%ICNTL(24)
      CALL MPI_BCAST(KEEP(110), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
C     KEEP(110) defaults to 0 for out of range values
      IF (KEEP(110).NE.1) KEEP(110)=0
      IF (KEEP(219).NE.0) THEN
       CALL DMUMPS_BUF_MAX_ARRAY_MINSIZE(max(KEEP(108),1),IERR)
       IF (IERR .NE. 0) THEN
C      ------------------------
C      Error allocating DMUMPS_BUF
C      ------------------------
          id%INFO(1) = -13
          id%INFO(2) = max(KEEP(108),1)
       END IF
      ENDIF
C     -----------------------------------------------
C     Depending on the option used for 
C       -detecting null pivots (ICNTL(24)/KEEP(110))
C         CNTL(3) is used to set DKEEP(1)
C               ( A row is considered as null if ||row|| < DKEEP(1) )
C         CNTL(5) is then used to define if a large 
C                 value is set on the diagonal or if a 1 is set
C                 and other values in the row are reset to zeros.
C         SEUIL* corresponds to the minimum required 
C                absolute value of pivot.
C         SEUIL_LDLT_NIV2 is used only in the 
C                case of SYM=2 within a niv2 node for which 
C                we have only a partial view of the fully summed rows.
      IF (id%MYID .EQ. MASTER) CNTL3 = id%CNTL(3)
      CALL MPI_BCAST(CNTL3, 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL5 = id%CNTL(5)
      CALL MPI_BCAST(CNTL5, 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL6 = id%CNTL(6)
      CALL MPI_BCAST(CNTL6, 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) id%DKEEP(8) = id%CNTL(7)
      CALL MPI_BCAST(id%DKEEP(8), 1, MPI_DOUBLE_PRECISION,
     &               MASTER, id%COMM, IERR)
      id%DKEEP(11) = id%DKEEP(8)/id%KEEP(461)
      id%DKEEP(12) = id%DKEEP(8)/id%KEEP(462)
      IF (KEEP(486).EQ.0) id%DKEEP(8) = ZERO
      COMPUTE_ANORMINF = .FALSE.
      IF ( (KEEP(486) .NE. 0).AND. (id%DKEEP(8).LT.ZERO)) THEN
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (KEEP(19).NE.0) THEN
C       Rank revealing factorisation
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (KEEP(110).NE.0) THEN
C       Null pivot detection       
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
C     -------------------------------------------------------
C        We compute ANORMINF, when needed, based on
C        the infinite norm of Rowsca *A*Colsca
C        and make it available on all working processes.
      IF (COMPUTE_ANORMINF) THEN
         CALL DMUMPS_ANORMINF(  id , ANORMINF, LSCAL )
      ELSE
         ANORMINF = ZERO
      ENDIF
C
C     Set BLR threshold
      IF (id%DKEEP(8).LT.ZERO) THEN
        id%DKEEP(8) = abs(id%DKEEP(8))*ANORMINF
      ENDIF
       IF ((KEEP(19).NE.0).OR.(KEEP(110).NE.0)) THEN 
         IF (PROKG) THEN 
           WRITE(MPG,'(A,1PD16.4)')
     &    ' Effective value of CNTL(3)                 =',CNTL3      
         ENDIF
       ENDIF
       IF (KEEP(19).EQ.0) THEN 
C        -- RR is off
         SEUIL = ZERO
         id%DKEEP(9) = ZERO
       ELSE
C        -- RR is on
C
C      CNTL(3) is the threshold used in the following to compute  
C      DKEEP(9) the threshold under which the sing val. are considered
C      as null and from which we start to look for a gap between two
C      sing val.
         IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(9) = abs(CNTL(3))
         ELSE IF  (CNTL3 .GT. ZERO) THEN
           id%DKEEP(9) = CNTL3*ANORMINF
         ELSE  !  (CNTL(3) .EQ. ZERO) THEN
         ENDIF
         IF (PROKG) THEN 
           WRITE(MPG, '(A,I10)')
     &    'ICNTL(56) rank revealing effective value   =',KEEP(19)   
           WRITE(MPG,'(A,1PD10.3)')
     &    ' ...Threshold for singularities on the root =',id%DKEEP(9)
         ENDIF
C        RR postponing considers that pivot rows with norm smaller 
C        than SEUIL should be postponed. 
C        SEUIL should be bigger than DKEEP(9), this means that 
C        DKEEP(13) should be bigger than 1. 
         Thresh_Seuil = id%DKEEP(13)
         IF (id%DKEEP(13).LT.1)  Thresh_Seuil = 10
         SEUIL = id%DKEEP(9)*Thresh_Seuil
         IF (PROKG) WRITE(MPG,'(A,1PD10.3)')
     &   ' ...Threshold for postponing                =',SEUIL
       ENDIF !end KEEP(19)
       SEUIL_LDLT_NIV2 = SEUIL
C     -------------------------------
C     -- Null pivot row detection
C     -------------------------------
       IF (KEEP(110).EQ.0) THEN 
C        -- Null pivot is off
C        Initialize DKEEP(1) to a negative value
C        in order to avoid detection of null pivots
C        (test max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL
C        in DMUMPS_FAC_I, where PIVNUL=DKEEP(1))
         id%DKEEP(1) = -1.0D0
         id%DKEEP(2) = ZERO
       ELSE
C        -- Null pivot is on
        IF (KEEP(19).NE.0) THEN
C        -- RR is on
C      RR postponing considers that pivot rows of norm smaller that SEUIL 
C      should be postponed, but pivot rows smaller than DKEEP(1) are 
C      directly added to null space and thus considered as null pivot rows. 
         IF ((id%DKEEP(10).LE.0).OR.(id%DKEEP(10).GT.1)) THEN
C           DKEEP(10) is out of range, set to the default value 10-1
           id%DKEEP(1) = id%DKEEP(9)*1D-1
         ELSE
           id%DKEEP(1) = id%DKEEP(9)*id%DKEEP(10) 
         ENDIF    
        ELSE
C        -- RR is off
C        -- only Null pivot detection
C         We keep strategy currently used in MUMPS 4.10.0
          IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(1)  = abs(CNTL(3))
          ELSE IF  (CNTL3 .GT. ZERO) THEN
            id%DKEEP(1)  = CNTL3*ANORMINF
          ELSE !  (CNTL(3) .EQ. ZERO) THEN
c          id%DKEEP(1)  = NPIV_CRITICAL_PATH*EPS*ANORMINF
           CALL MUMPS_NPIV_CRITICAL_PATH(
     &         N, KEEP(28), id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     &         id%NA(1), id%LNA, id%NE_STEPS(1), NPIV_CRITICAL_PATH )
           id%DKEEP(1)  = sqrt(dble(NPIV_CRITICAL_PATH))*EPS*ANORMINF
          ENDIF
        ENDIF ! fin rank revealing        
        IF ((KEEP(110).NE.0).AND.(PROKG)) THEN
           WRITE(MPG, '(A,I16)')
     &    ' ICNTL(24) null pivot rows detection        =',KEEP(110)   
           WRITE(MPG,'(A,1PD16.4)')
     &    ' ...Zero pivot detection threshold          =',id%DKEEP(1)
        ENDIF 
        IF (CNTL5.GT.ZERO) THEN
            id%DKEEP(2) = CNTL5 * ANORMINF
            IF (PROKG) WRITE(MPG,'(A,1PD10.3)')
     &    ' ...Fixation for null pivots               =',id%DKEEP(2)
         ELSE
            IF (PROKG) WRITE(MPG,*) '...Infinite fixation '
            IF (id%KEEP(50).EQ.0) THEN
C             Unsym
            ! the user let us choose a fixation. set in NEGATIVE
            ! to detect during facto when to set row to zero !
             id%DKEEP(2) = -max(1.0D10*ANORMINF, 
     &                sqrt(huge(ANORMINF))/1.0D8)
            ELSE
C             Sym
            id%DKEEP(2) = ZERO
            ENDIF
         ENDIF
       ENDIF ! fin null pivot detection.
C     Find id of root node if RR is on 
      IF (KEEP(53).NE.0) THEN
        ID_ROOT =MUMPS_PROCNODE(id%PROCNODE_STEPS(id%STEP(KEEP(20))),
     &                          id%KEEP(199))
        IF ( KEEP( 46 )  .NE. 1 ) THEN
          ID_ROOT = ID_ROOT + 1
        END IF
      ENDIF
C Second pass:  set parameters for null pivot detection
C Allocate PIVNUL_LIST in case of null pivot detection
      LPN_LIST = 1
      IF ( associated( id%PIVNUL_LIST) ) DEALLOCATE(id%PIVNUL_LIST)
      IF(KEEP(110) .EQ. 1) THEN
         LPN_LIST = N
      ENDIF
      IF (KEEP(19).NE.0 .AND.
     &   (ID_ROOT.EQ.id%MYID .OR. id%MYID.EQ.MASTER)) THEN
         LPN_LIST = N
      ENDIF
      ALLOCATE( id%PIVNUL_LIST(LPN_LIST),stat = IERR )
      IF ( IERR .GT. 0 ) THEN
        id%INFO(1)=-13
        id%INFO(2)=LPN_LIST
      END IF
      id%PIVNUL_LIST(1:LPN_LIST) = 0
      KEEP(109) = 0
C end set parameter for null pivot detection
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 517
C   --------------------------------------------------------------
C   STATIC PIVOTING
C     -- Static pivoting only when RR and Null pivot detection OFF
C   --------------------------------------------------------------
      IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
        IF (id%MYID .EQ. MASTER) CNTL4 = id%CNTL(4)
        CALL MPI_BCAST( CNTL4, 1, MPI_DOUBLE_PRECISION,
     &                MASTER, id%COMM, IERR )
C 
        IF ( CNTL4 .GE. ZERO ) THEN
         KEEP(97) = 1
         IF ( CNTL4 .EQ. ZERO ) THEN
C           -- set seuil to sqrt(eps)*||A||
            IF(ANORMINF .EQ. ZERO) THEN
               CALL DMUMPS_ANORMINF(  id , ANORMINF, LSCAL )
            ENDIF
            SEUIL = sqrt(EPS) * ANORMINF
         ELSE
            SEUIL = CNTL4
         ENDIF
         SEUIL_LDLT_NIV2 = SEUIL
C
        ELSE 
         SEUIL = ZERO
        ENDIF
      ENDIF
C     set number of tiny pivots / 2x2 pivots in types 1 /
C     2x2 pivots in types 2, to zero. This is because the
C     user can call the factorization step several times.
      KEEP(98)  = 0
      KEEP(103) = 0
      KEEP(105) = 0
      MAXS      = 1_8
*
*     Start allocations
*     *****************
*
C
C  The slaves can now perform the factorization
C
C
C  Allocate id%S on all nodes
C  or point to user provided data WK_USER when LWK_USER>0
C  =======================
C
C     Compute BLR_STRAT and a first estimation 
C     of MAXS, the size of id%S
      CALL  DMUMPS_SET_BLRSTRAT_AND_MAXS_K8 (
     &           MAXS_BASE8, MAXS_BASE_RELAXED8,
     &           BLR_STRAT,
     &           id%KEEP(1), id%KEEP8(1))
C   
      MAXS = MAXS_BASE_RELAXED8
      IF (WK_USER_PROVIDED) THEN
C       -- Set MAXS to size of WK_USER_
        MAXS = id%KEEP8(24)
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 517
      ENDIF
C
      id%KEEP8(75) = huge(id%KEEP8(75))
      id%KEEP8(76) = huge(id%KEEP8(76))
      IF ((.NOT.WK_USER_PROVIDED).AND.(I_AM_SLAVE)) THEN
C
        IF (id%KEEP8(4) .NE. 0_8) THEN
C        -------------------------
C        WE TRY TO USE MEM_ALLOWED (KEEP8(4)/1D6)
C        -------------------------
C        Set MAXS given BLR_STRAT, KEEP(201) and MAXS_BASE_RELAXED8
         CALL DMUMPS_MEM_ALLOWED_SET_MAXS (
     &           MAXS,
     &           BLR_STRAT, id%KEEP(201), MAXS_BASE_RELAXED8,
     &           id%KEEP(1), id%KEEP8(1), id%MYID, id%N, id%NELT, 
     &           id%NA(1), id%LNA, id%NSLAVES, 
     &           KEEP464COPY,  KEEP465COPY,
     &           id%INFO(1), id%INFO(2)
     &           )
        ENDIF ! MEM_ALLOWED
C
      ENDIF ! (.NOT.WK_USER_PROVIDED).AND.(I_AM_SLAVE)) THEN
C
      IF (I_AM_SLAVE) THEN
      ENDIF ! I_AM_SLAVE)
C
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 517
      ENDIF
      CALL MUMPS_SETI8TOI4(MAXS, id%INFO(39))
      CALL DMUMPS_AVGMAX_STAT8(PROKG, MPG, MAXS, id%NSLAVES,
     & PRINT_MAXAVG,
     & id%COMM, " Effective size of S     (based on INFO(39))=   ")
C
      IF ( I_AM_SLAVE ) THEN
C       ------------------
C       Dynamic scheduling
C       ------------------
        CALL DMUMPS_LOAD_SET_INICOST( dble(id%COST_SUBTREES),
     &        KEEP(64), id%DKEEP(15), KEEP(375), MAXS )
        K28=KEEP(28)
        MEMORY_MD_ARG = min(int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8 ),
C       Restrict freedom from dynamic scheduler when
C       MEM_ALLOWED=ICNTL(23) is small (case where KEEP8(4)-MAXS_BASE8
C       is negative after call to DMUMPS_MAX_MEM)
     &                      max(0_8, MAXS-MAXS_BASE8))
        CALL DMUMPS_LOAD_INIT( id, MEMORY_MD_ARG, MAXS )
C
C       Out-Of-Core (OOC) issues. Case where we ran one factorization OOC
C       and the second one is in-core: we try to free OOC
C       related data from previous factorization.
C
        CALL DMUMPS_CLEAN_OOC_DATA(id, IERR)
        IF (IERR < 0) THEN
          id%INFO(1) = -90
          id%INFO(2) = 0
          GOTO 112
        ENDIF
        IF (KEEP(201) .GT. 0) THEN
C          -------------------
C          OOC initializations
C          -------------------
           IF (KEEP(201).EQ.1 !PANEL Version
     &         .AND.KEEP(50).EQ.0 ! Unsymmetric
     &         .AND.KEEP(251).NE.2 ! Store L to disk
     &         ) THEN
              id%OOC_NB_FILE_TYPE=2 ! declared in MUMPS_OOC_COMMON
           ELSE
              id%OOC_NB_FILE_TYPE=1 ! declared in MUMPS_OOC_COMMON
           ENDIF
C          ------------------------------
C          Dimension IO buffer, KEEP(100)
C          ------------------------------
           IF (KEEP(205) .GT. 0) THEN
             KEEP(100) = KEEP(205)
           ELSE
             IF (KEEP(201).EQ.1) THEN ! PANEL version
               I8TMP = int(id%OOC_NB_FILE_TYPE,8) *
     &               2_8 * int(KEEP(226),8)
             ELSE
               I8TMP = 2_8 * id%KEEP8(119)
             ENDIF
             I8TMP = I8TMP +  int(max(KEEP(12),0),8) *
     &               (I8TMP/100_8+1_8)
C            we want to avoid too large IO buffers.
C            12M corresponds to 100Mbytes given to buffers.
             I8TMP = min(I8TMP, 12000000_8)
             KEEP(100)=int(I8TMP)
           ENDIF
           IF (KEEP(201).EQ.1) THEN
C            Panel version. Force the use of a buffer. 
             IF ( KEEP(99) < 3 ) THEN
               KEEP(99) = KEEP(99) + 3
             ENDIF
           ENDIF
C          --------------------------
C          Reset KEEP(100) to 0 if no
C          buffer is used for OOC.
C          --------------------------
           IF (KEEP(99) .LT.3) KEEP(100)=0
           IF((dble(KEEP(100))*dble(KEEP(35))/dble(2)).GT.
     &          (dble(1999999999)))THEN
             IF (PROKG) THEN
               WRITE(MPG,*)id%MYID,': Warning: DIM_BUF_IO might be
     &  too big for Filesystem'
             ENDIF
           ENDIF
           ALLOCATE (id%OOC_INODE_SEQUENCE(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_INODE_SEQUENCE)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_TOTAL_NB_NODES(id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE
              NULLIFY(id%OOC_TOTAL_NB_NODES)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_SIZE_OF_BLOCK(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_SIZE_OF_BLOCK)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_VADDR(KEEP(28),id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_VADDR)
              GOTO 112
           ENDIF
        ENDIF
      ENDIF
 112  CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) < 0) THEN
C       LOAD_END must be done but not OOC_END_FACTO
        GOTO 513
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (KEEP(201) .GT. 0) THEN
           IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
             CALL DMUMPS_OOC_INIT_FACTO(id,MAXS)
          ELSE
             WRITE(*,*) "Internal error in DMUMPS_FAC_DRIVER"
             CALL MUMPS_ABORT()
           ENDIF
           IF(id%INFO(1).LT.0)THEN
              GOTO 111
           ENDIF
        ENDIF
C       First increment corresponds to the number of
C       floating-point operations for subtrees allocated
C       to the local processor.
        CALL DMUMPS_LOAD_UPDATE(0,.FALSE.,dble(id%COST_SUBTREES),
     &          id%KEEP(1),id%KEEP8(1))
        IF (id%INFO(1).LT.0) GOTO 111
      END IF
C     -----------------------
C     Manage main workarray S
C     -----------------------
      EARLYT3ROOTINS = KEEP(200) .EQ.0
#if defined (LARGEMATRICES)
      IF ( id%MYID .ne. MASTER ) THEN
#endif
      IF (.NOT.WK_USER_PROVIDED) THEN
        IF ( EARLYT3ROOTINS ) THEN
C         Standard allocation strategy
          ALLOCATE (id%S(MAXS),stat=IERR)
          id%KEEP8(23) = MAXS
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -13
            CALL MUMPS_SET_IERROR(MAXS, id%INFO(2))
C           On some platforms (IBM for example), an
C           allocation failure returns a non-null pointer.
C           Therefore we nullify S
            NULLIFY(id%S)
            id%KEEP8(23)=0_8
          ENDIF
        ENDIF
      ELSE
       id%S => id%WK_USER(1:id%KEEP8(24))
       id%KEEP8(23) = 0_8
      ENDIF
#if defined (LARGEMATRICES)
      END IF
#endif
C
 111  CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 514
C     --------------------------
C     Initialization of modules
C     related to data management
C     --------------------------
      NB_ACTIVE_FRONTS_ESTIM = 3
      IF (I_AM_SLAVE) THEN
C
        CALL MUMPS_FDM_INIT('A',NB_ACTIVE_FRONTS_ESTIM, id%INFO)
C
        IF ( (KEEP(486).EQ.2) 
     &     .OR. ((KEEP(489).NE.0).AND.(KEEP(400).GT.1))
     &    ) THEN
C         In case of LRSOLVE or CompressCB, 
C         initialize nb of handlers to nb of BLR
C         nodes estimated at analysis
          NB_FRONTS_F_ESTIM = KEEP(470)
        ELSE
          IF (KEEP(489).NE.0) THEN
C          Compress CB and no L0 OMP (or 1 thread under L0): 
C          NB_ACTIVE_FRONTS_ESTIM is too small, 
C          to limit nb of reallocations make it twice larger
           NB_FRONTS_F_ESTIM = 2*NB_ACTIVE_FRONTS_ESTIM
          ELSE
           NB_FRONTS_F_ESTIM = NB_ACTIVE_FRONTS_ESTIM
          ENDIF
        ENDIF
        CALL MUMPS_FDM_INIT('F',NB_FRONTS_F_ESTIM, id%INFO )
        IF (id%INFO(1) .LT. 0 ) GOTO 114
#if ! defined(NO_FDM_DESCBAND)
C         Storage of DESCBAND information
          CALL MUMPS_FDBD_INIT( NB_ACTIVE_FRONTS_ESTIM, id%INFO )
#endif
#if ! defined(NO_FDM_MAPROW)
C         Storage of MAPROW and ROOT2SON information
          CALL MUMPS_FMRD_INIT( NB_ACTIVE_FRONTS_ESTIM, id%INFO )
#endif
        CALL DMUMPS_BLR_INIT_MODULE( NB_FRONTS_F_ESTIM, id%INFO )
 114    CONTINUE
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                       id%COMM, id%MYID )
C     GOTO 500: one of the above module initializations failed
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C
C  Allocate space for matrix in arrowhead 
C  ======================================
C
C  CASE 1 : Matrix is assembled
C  CASE 2 : Matrix is elemental
C
      IF ( KEEP(55) .eq. 0 ) THEN
C       ------------------------------------
C       Space has been allocated already for
C       the integer part during analysis
C       Only slaves need the arrowheads.
C       ------------------------------------
        IF (associated( id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
        IF ( I_AM_SLAVE .and. id%KEEP8(26) .ne. 0_8 ) THEN
          ALLOCATE( id%DBLARR( id%KEEP8(26) ), stat = IERR )
        ELSE
          ALLOCATE( id%DBLARR( 1 ), stat =IERR )
        END IF
        IF ( IERR .NE. 0 ) THEN
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for DBLARR(',id%KEEP8(26),')'
          ENDIF
          id%INFO(1)=-13
          CALL MUMPS_SET_IERROR(id%KEEP8(26), id%INFO(2))
          NULLIFY(id%DBLARR)
          GOTO 100
        END IF
      ELSE
C        ----------------------------------------
C        Allocate variable lists. Systematically.
C        ----------------------------------------
         IF ( associated( id%INTARR ) ) THEN
           DEALLOCATE( id%INTARR )
           NULLIFY( id%INTARR )
         END IF
         IF ( I_AM_SLAVE .and. id%KEEP8(27) .ne. 0_8 ) THEN
           ALLOCATE( id%INTARR( id%KEEP8(27) ), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             CALL MUMPS_SET_IERROR(id%KEEP8(27), id%INFO(2))
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         ELSE
           ALLOCATE( id%INTARR(1),stat=allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         END IF
C        -----------------------------
C        Allocate real values.
C        On master, if hybrid host and
C        no scaling, avoid the copy.
C        -----------------------------
         IF (associated( id%DBLARR)) THEN
           DEALLOCATE(id%DBLARR)
           NULLIFY(id%DBLARR)
         ENDIF
         IF ( I_AM_SLAVE ) THEN
           IF (      id%MYID_NODES .eq. MASTER
     &       .AND.   KEEP(46)   .eq. 1
     &       .AND.   KEEP(52)   .eq. 0 ) THEN
C            --------------------------
C            Simple pointer association
C            --------------------------
             id%DBLARR => id%A_ELT
           ELSE
C            ----------
C            Allocation
C            ----------
             IF ( id%KEEP8(26) .ne. 0_8 ) THEN
               ALLOCATE( id%DBLARR( id%KEEP8(26) ), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 CALL MUMPS_SET_IERROR(id%KEEP8(26), id%INFO(2))
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF
             ELSE
               ALLOCATE( id%DBLARR(1), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = 1
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF 
             END IF
           END IF
         ELSE
           ALLOCATE( id%DBLARR(1), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%DBLARR)
             GOTO 100
           END IF
         END IF
      END IF
C     -----------------
C     Also prepare some
C     data for the root
C     -----------------
      IF ( KEEP(38).NE.0 .AND.  I_AM_SLAVE ) THEN
         CALL DMUMPS_INIT_ROOT_FAC( id%N,
     &   id%root, id%FILS(1), KEEP(38), id%KEEP(1), id%INFO(1) )
      END IF
C
C
 100  CONTINUE
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C       -----------------------------------
C
C       DISTRIBUTION OF THE ORIGINAL MATRIX
C
C       -----------------------------------
C
C     TIMINGS: computed (and printed) on the host
C     Next line: global time for distrib(arrowheads,elts)
C     on the host. Synchronization has been performed.
      IF (id%MYID.EQ.MASTER) CALL MUMPS_SECDEB(TIME)
C     -------------------------------------------
C     S_PTR_ARG / MAXS_ARG will be used for id%S
C     argument to arrowhead/element distribution
C     routines: if id%S is not allocated, we pass
C     S_DUMMY_ARG instead, which is not accessed.
C     -------------------------------------------
      IF (EARLYT3ROOTINS) THEN
            S_PTR_ARG => id%S
            MAXS_ARG = MAXS
      ELSE
            S_PTR_ARG => S_DUMMY_ARG
            MAXS_ARG = 1
      ENDIF
C
      IF ( KEEP( 55 ) .eq. 0 ) THEN
C       ----------------------------
C       Original matrix is assembled
C       Arrowhead format to be used.
C       ----------------------------
C       KEEP8(26) and KEEP8(27) hold the number of entries for real/integer
C       for the matrix in arrowhead format. They have been set by the
C       analysis phase (DMUMPS_ANA_F and DMUMPS_ANA_G)
C
C       ------------------------------------------------------------------
C       Blocking is used for sending arrowhead records (I,J,VAL)
C              buffer(1) is used to store number of bytes already packed
C              buffer(2) number of records already packed
C       KEEP(39) : Number of records (blocking factor)
C       ------------------------------------------------------------------
C
C     ---------------------------------------------
C     In case of parallel root compute minimum
C     size of workspace to receive arrowheads
C     of root node. Will be used to check that
C     MAXS is large enough for arrowheads (case
C     of EARLYT3ROOTINS (KEEP(200)=0); if .NOT.
C     EARLYT3ROOTINS (KEEP(200)=1), root will
C     be assembled into id%S later and size of
C     id%S will be checked later)
C     ---------------------------------------------
      IF (EARLYT3ROOTINS .AND. KEEP(38).NE.0 .AND.
     &    KEEP(60) .EQ.0 .AND. I_AM_SLAVE) THEN
        LWK = int(numroc( id%root%ROOT_SIZE, id%root%MBLOCK,
     &             id%root%MYROW, 0, id%root%NPROW ),8)
        LWK = max( 1_8, LWK )
        LWK = LWK*
     &        int(numroc( id%root%ROOT_SIZE, id%root%NBLOCK,
     &        id%root%MYCOL, 0, id%root%NPCOL ),8)
        LWK = max( 1_8, LWK )
      ELSE
        LWK = 1_8
      ENDIF
C     MAXS must be at least 1, and in case of
C     parallel root, large enough to receive
C     arrowheads of root.
      IF (MAXS .LT. int(LWK,8)) THEN
           id%INFO(1) = -9
           CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
      IF ( KEEP(54) .eq. 0 ) THEN
C       ================================================
C       FIRST CASE : MATRIX IS NOT INITIALLY DISTRIBUTED
C       ================================================
C       A small integer workspace is needed to
C       send the arrowheads.
        IF ( id%MYID .eq. MASTER ) THEN
          ALLOCATE(IWK(id%N), stat=allocok)
          IF ( allocok .NE. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=id%N
          END IF
#if defined(LARGEMATRICES)
          ALLOCATE (WK(LWK),stat=IERR)
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -13
            CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
            write(6,*) ' PB1 ALLOC LARGEMAT'
          ENDIF
#endif
        ENDIF
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
        IF ( id%INFO(1).LT.0 ) GOTO 500
        IF ( id%MYID .eq. MASTER ) THEN
C
C         --------------------------------
C         MASTER sends arowheads using the
C         global communicator with ranks
C         also in global communicator
C         IWK is used as temporary
C         workspace of size N.
C         --------------------------------
           IF ( .not. associated( id%INTARR ) ) THEN
              ALLOCATE( id%INTARR( 1 ),stat=IERR)
              IF ( IERR .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = 1
                 NULLIFY(id%INTARR)
                 write(6,*) ' PB2 ALLOC INTARR'
                 CALL MUMPS_ABORT()
              ENDIF
           ENDIF
          NBRECORDS = KEEP(39)
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(id%KEEP8(28))
          ENDIF
#if defined(LARGEMATRICES)
          CALL DMUMPS_FACTO_SEND_ARROWHEADS(id%N, id%KEEP8(28), id%A(1),
     &      id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &      LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &      id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &      NBRECORDS,
     &      LP, id%COMM, id%root, KEEP,id%KEEP8,
     &      id%FILS(1), IWK(1), ! workspace of size N
     &
     &      id%INTARR(1), id%KEEP8(27), id%DBLARR(1), id%KEEP8(26),
     &      id%PTRAR(1), id%PTRAR(id%N+1),
     &      id%FRERE_STEPS(1), id%STEP(1), WK(1), LWK,
     &      id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &      id%CANDIDATES(1,1)) 
C         write(6,*) '!!! A,IRN,JCN are freed during factorization '
          DEALLOCATE (id%A)
          NULLIFY(id%A)
          DEALLOCATE (id%IRN)
          NULLIFY (id%IRN)
          DEALLOCATE (id%JCN)
          NULLIFY (id%JCN)
          IF (.NOT.WK_USER_PROVIDED) THEN
            IF (EARLYT3ROOTINS) THEN
              ALLOCATE (id%S(MAXS),stat=IERR)
              id%KEEP8(23) = MAXS
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                id%INFO(2) = MAXS
                NULLIFY(id%S)
                id%KEEP8(23)=0_8
                write(6,*) ' PB2 ALLOC LARGEMAT',MAXS
                CALL MUMPS_ABORT()
               ENDIF
             ENDIF
            ENDIF
          ELSE
            id%S => id%WK_USER(1:id%KEEP8(24))
          ENDIF
          IF (EARLYT3ROOTINS) THEN
            id%S(MAXS-LWK+1_8:MAXS) = WK(1_8:LWK)
          ENDIF
          DEALLOCATE (WK)
#else
          CALL DMUMPS_FACTO_SEND_ARROWHEADS(id%N, id%KEEP8(28), id%A(1),
     &    id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &    LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &    id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &    NBRECORDS,
     &    LP, id%COMM, id%root, KEEP(1),id%KEEP8(1),
     &    id%FILS(1), IWK(1),
     &
     &    id%INTARR(1), id%KEEP8(27), id%DBLARR(1), id%KEEP8(26),
     &    id%PTRAR(1), id%PTRAR(id%N+1),
     &    id%FRERE_STEPS(1), id%STEP(1), S_PTR_ARG(1), MAXS_ARG,
     &    id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &    id%CANDIDATES(1,1) ) 
#endif
          DEALLOCATE(IWK)
        ELSE
          NBRECORDS = KEEP(39)
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(id%KEEP8(28))
          ENDIF
          CALL DMUMPS_FACTO_RECV_ARROWHD2( id%N,
     &       id%DBLARR(1), id%KEEP8(26),
     &       id%INTARR(1), id%KEEP8(27),
     &       id%PTRAR( 1 ),
     &       id%PTRAR(id%N+1),
     &       KEEP( 1 ), id%KEEP8(1), id%MYID, id%COMM,
     &       NBRECORDS,
     &
     &       S_PTR_ARG(1), MAXS_ARG,
     &       id%root,
     &       id%PROCNODE_STEPS(1), id%NSLAVES,
     &       id%SYM_PERM(1), id%FRERE_STEPS(1), id%STEP(1),
     &       id%INFO(1), id%INFO(2) )
        ENDIF
      ELSE
C
C     =============================================
C     SECOND CASE : MATRIX IS INITIALLY DISTRIBUTED
C     =============================================
C     Timing on master.
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECDEB(TIME)
      END IF
      IF ( I_AM_SLAVE ) THEN
C       ---------------------------------------------------
C       In order to have possibly IRN_loc/JCN_loc/A_loc
C       of size 0, avoid to pass them inside REDISTRIBUTION
C       and pass id instead
C       NZ_locMAX8 gives as a maximum buffer size (send/recv) used
C        an upper bound to limit buffers on small matrices
C       ---------------------------------------------------
       CALL MPI_ALLREDUCE(id%KEEP8(29), NZ_locMAX8, 1, MPI_INTEGER8,
     &                   MPI_MAX, id%COMM_NODES, IERR)
       NBRECORDS = KEEP(39)
       IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(NZ_locMAX8)
       ENDIF
        CALL DMUMPS_REDISTRIBUTION( id%N,
     &  id%KEEP8(29),
     &  id,
     &  id%DBLARR(1), id%KEEP8(26), id%INTARR(1),
     &  id%KEEP8(27), id%PTRAR(1), id%PTRAR(id%N+1),
     &  KEEP(1), id%KEEP8(1), id%MYID_NODES,
     &  id%COMM_NODES, NBRECORDS,
     &  S_PTR_ARG(1), MAXS_ARG, id%root, id%PROCNODE_STEPS(1),
     &  id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     &  id%ICNTL(1), id%INFO(1), NSEND8, NLOCAL8,
     &  id%ISTEP_TO_INIV2(1),
     &  id%CANDIDATES(1,1) )
        IF ( ( KEEP(52).EQ.7 ).OR. (KEEP(52).EQ.8) ) THEN
C         -------------------------------------------------
C         In that case, scaling arrays have been allocated
C         on all processors. They were useful for matrix
C         distribution. But we now really only need them
C         on the host. In case of distributed solution, we
C         will have to broadcast either ROWSCA or COLSCA
C         (depending on MTYPE) but this is done later.
C
C         In other words, on exit from the factorization,
C         we want to have scaling arrays available only
C         on the host.
C         -------------------------------------------------
          IF ( id%MYID > 0 ) THEN
            IF (associated(id%ROWSCA)) THEN
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
            IF (associated(id%COLSCA)) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ENDIF
          ENDIF
        ENDIF
#if defined(LARGEMATRICES)
C      deallocate id%IRN_loc, id%JCN(loc) to free extra space
C      Note that in this case IRN_loc cannot be used
C      anymore during the solve phase for IR and Error analysis.
         IF (associated(id%IRN_loc)) THEN
            DEALLOCATE(id%IRN_loc)
            NULLIFY(id%IRN_loc)
         ENDIF
         IF (associated(id%JCN_loc)) THEN
            DEALLOCATE(id%JCN_loc)
            NULLIFY(id%JCN_loc)
         ENDIF
         IF (associated(id%A_loc)) THEN
            DEALLOCATE(id%A_loc)
            NULLIFY(id%A_loc)
         ENDIF
       write(6,*) ' Warning :', 
     &        ' id%A_loc, IRN_loc, JCN_loc deallocated !!! '
#endif
      IF (PROK) THEN
        WRITE(MP,120) NLOCAL8, NSEND8
      END IF
      END IF
      IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
C       ------------------------------
C       The host is not working -> had
C       no data from initial matrix
C       ------------------------------
        NSEND8  = 0_8
        NLOCAL8 = 0_8
      END IF
C     --------------------------
C     Put into some info/infog ?
C     --------------------------
      CALL MPI_REDUCE( NSEND8, NSEND_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( NLOCAL8, NLOCAL_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
        WRITE(MPG,125) NLOCAL_TOT8, NSEND_TOT8
      END IF
C
C     -------------------------
C     Check for possible errors
C     -------------------------
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C
      ENDIF
      ELSE
C       -------------------
C       Matrix is elemental,
C       provided on the
C       master only
C       -------------------
        IF ( id%MYID.eq.MASTER)
     &   CALL DMUMPS_MAXELT_SIZE( id%ELTPTR(1),
     &                        id%NELT,
     &                        MAXELT_SIZE )
C
C         Perform the distribution of the elements.
C         A this point,
C           PTRAIW/PTRARW have been computed.
C           INTARR/DBLARR have been allocated
C           ELTPROC gives the mapping of elements
C
        CALL DMUMPS_ELT_DISTRIB( id%N, id%NELT, id%KEEP8(30),
     &     id%COMM, id%MYID,
     &     id%NSLAVES, id%PTRAR(1),
     &     id%PTRAR(id%NELT+2),
     &     id%INTARR(1), id%DBLARR(1), id%KEEP8(27), id%KEEP8(26),
     &     id%KEEP(1), id%KEEP8(1), MAXELT_SIZE,
     &     id%FRTPTR(1), id%FRTELT(1),
     &     S_PTR_ARG(1), MAXS_ARG, id%FILS(1),
     &     id, id%root )
C       ----------------
C       Broadcast errors
C       ----------------
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
      END IF ! Element entry
C     ------------------------
C     Time the redistribution:
C     ------------------------
      IF ( id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(93) = TIME
        IF (PROKG) WRITE(MPG,160) id%DKEEP(93)
      END IF
C
C     TIMINGS:
C     Next line: elapsed time for factorization
      IF (id%MYID.EQ.MASTER)  CALL MUMPS_SECDEB(TIME)
C
C  Allocate buffers on the workers
C  ===============================
C
      IF ( I_AM_SLAVE )  THEN
        CALL DMUMPS_BUF_INI_MYID(id%MYID_NODES)
C
C  Some buffers are required to pack/unpack data and for
C  receiving MPI messages.
C  For packing/unpacking : the buffer must be large
C  enough to send several messages while receives might not
C  be posted yet.
C  It is assumed that the size of an integer is held in KEEP(34)
C  while the size of a complex is held in KEEP(35).
C  BUFR and LBUFR are declared of type integer, since byte is not
C  a standard datatype.
C  We now use KEEP(43) or KEEP(379) and KEEP(44) or KEEP(380)
C  as estimated at analysis to allocate appropriate buffer sizes
C
C  Reception buffer
C  ----------------
        IF (KEEP(486).NE.0) THEN
          DMUMPS_LBUFR_BYTES8 = int(KEEP( 380 ),8) * int(KEEP( 35 ),8)
        ELSE
          DMUMPS_LBUFR_BYTES8 = int(KEEP( 44 ),8) * int(KEEP( 35 ),8)
        ENDIF
C       ---------------------------------------
C       Ensure a reasonable minimal buffer size
C       ---------------------------------------
        DMUMPS_LBUFR_BYTES8 = max( DMUMPS_LBUFR_BYTES8,
     &                      100000_8 )
C
C  If there is pivoting, size of the message might still increase.
C  We use a relaxation (so called PERLU) to increase the estimate.
C
C  Note: PERLU is a global estimate for pivoting. 
C  It may happen that one large contribution block size is increased by more than that.
C  This is why we use an extra factor 2 relaxation coefficient for the relaxation of
C  the reception buffer in the case where pivoting is allowed.
C  A more dynamic strategy could be applied: if message to
C  be received is larger than expected, reallocate a larger
C  buffer. (But this won't work with IRECV.)
C  Finally, one may want (as we are currently doing it for moste messages)
C  to cut large messages into a series of smaller ones.
C
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU = 2
        ELSE
          MIN_PERLU = 0
        ENDIF
C
        DMUMPS_LBUFR_BYTES8 = DMUMPS_LBUFR_BYTES8
     &        + int( 2.0D0 * dble(max(PERLU,MIN_PERLU))*
     &        dble(DMUMPS_LBUFR_BYTES8)/100D0, 8)
        DMUMPS_LBUFR_BYTES8 = min(DMUMPS_LBUFR_BYTES8,
     &                            int(huge (KEEP(44))-100,8))
        DMUMPS_LBUFR_BYTES  = int( DMUMPS_LBUFR_BYTES8 )
        IF (KEEP(48)==5) THEN
C          Since the buffer is going to be allocated, use
C          it as the constraint for memory/granularity
C          in hybrid scheduler
C
           id%KEEP8(21) = id%KEEP8(22) +
     &        int( dble(max(PERLU,MIN_PERLU))*
     &        dble(id%KEEP8(22))/100D0,8)
        ENDIF
C
C  Now estimate the size for the buffer for asynchronous
C  sends of contribution blocks (so called CB). We want to be able to send at
C  least KEEP(213)/100 (two in general) messages at the
C  same time.
C
C   Send buffer
C   -----------
        IF (KEEP(486).NE.0) THEN
         DMUMPS_LBUF8 = int( dble(KEEP(213)) / 100.0D0 *
     &                      dble(KEEP(379)) * dble(KEEP(35)), 8  )
        ELSE
         DMUMPS_LBUF8 = int( dble(KEEP(213)) / 100.0D0 *
     &                      dble(KEEP(43)) * dble(KEEP(35)), 8  )
        ENDIF
        DMUMPS_LBUF8 = max( DMUMPS_LBUF8, 100000_8 )
        DMUMPS_LBUF8 = DMUMPS_LBUF8
     &                 + int( 2.0D0 * dble(max(PERLU,MIN_PERLU))*
     &                   dble(DMUMPS_LBUF8)/100D0, 8)
C       Make DMUMPS_LBUF8 small enough to be stored in a standard integer
        DMUMPS_LBUF8 = min(DMUMPS_LBUF8, int(huge (KEEP(43))-100,8))
C
C       No reason to have send buffer smaller than receive buffer.
C       This should never occur with the formulas above but just
C       in case:
        DMUMPS_LBUF8 = max(DMUMPS_LBUF8, DMUMPS_LBUFR_BYTES8+3*KEEP(34))
        DMUMPS_LBUF  = int(DMUMPS_LBUF8)
        IF(id%KEEP(48).EQ.4)THEN
           DMUMPS_LBUFR_BYTES=DMUMPS_LBUFR_BYTES*5
           DMUMPS_LBUF=DMUMPS_LBUF*5
        ENDIF
C
C  Estimate size of buffer for small messages 
C  Each node can send ( NSLAVES - 1 ) messages to (NSLAVES-1) nodes
C
C  KEEP(56) is the number of nodes of level II.
C  Messages will be sent for the symmetric case
C  for synchronisation issues.
C
C  We take an upperbound
C
        DMUMPS_LBUF_INT = ( KEEP(56) + id%NSLAVES * id%NSLAVES ) * 5
     &               * KEEP(34)
        IF ( KEEP( 38 ) .NE. 0 ) THEN
C
C
          KKKK = MUMPS_PROCNODE( id%PROCNODE_STEPS(id%STEP(KEEP(38))),
     &                           id%KEEP(199) )
          IF ( KKKK .EQ. id%MYID_NODES ) THEN
             DMUMPS_LBUF_INT = DMUMPS_LBUF_INT + 4 * KEEP(34) *
     &         ( id%NSLAVES + id%NE_STEPS(id%STEP(KEEP(38)))
     &      + min(KEEP(56), id%NE_STEPS(id%STEP(KEEP(38)))) * id%NSLAVES
     &         )
          END IF
        END IF
C       At this point, DMUMPS_LBUFR_BYTES, DMUMPS_LBUF
C       and DMUMPS_LBUF_INT have been computed (all
C       are in numbers of bytes).
        IF ( PROK ) THEN
          WRITE( MP, 9999 ) DMUMPS_LBUFR_BYTES,
     &                      DMUMPS_LBUF, DMUMPS_LBUF_INT
        END IF
 9999   FORMAT( /,' Allocated buffers',/,' ------------------',/,
     &  ' Size of reception buffer in bytes ...... = ', I10,
     &  /,
     &  ' Size of async. emission buffer (bytes).. = ', I10,/,
     &  ' Small emission buffer (bytes) .......... = ', I10)
C       --------------------------
C       Allocate small send buffer
C       required for DMUMPS_FAC_B
C       --------------------------
        CALL DMUMPS_BUF_ALLOC_SMALL_BUF( DMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)= -13
C         convert to size in integer  id%INFO(2)= DMUMPS_LBUF_INT
          id%INFO(2)= (DMUMPS_LBUF_INT+KEEP(34)-1)/KEEP(34)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &     ':Allocation error in DMUMPS_BUF_ALLOC_SMALL_BUF'
     &     ,id%INFO(2)
          ENDIF
          GO TO 110
        END IF
C
C       --------------------------------------
C       Allocate reception buffer on all procs
C       This is done now.
C       --------------------------------------
        DMUMPS_LBUFR = (DMUMPS_LBUFR_BYTES+KEEP(34)-1)/KEEP(34)
        ALLOCATE( BUFR( DMUMPS_LBUFR ),stat=IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = DMUMPS_LBUFR
          IF (LPOK) THEN
            WRITE(LP,*)
     &     ': Allocation error for BUFR(', DMUMPS_LBUFR,
     &     ') on MPI process',id%MYID
          ENDIF
          GO TO 110
        END IF
C       -----------------------------------------
C       Estimate MAXIS. IS will be allocated in
C       DMUMPS_FAC_B. It will contain factors and
C       contribution blocks integer information
C       -----------------------------------------
C       Relax integer workspace based on PERLU
        PERLU          = KEEP( 12 )
        IF (KEEP(201).GT.0) THEN
C         OOC panel or non panel (note that
C         KEEP(15)=KEEP(225) if non panel)
          MAXIS_ESTIM   = KEEP(225)
        ELSE
C         In-core or reals for factors not stored
          MAXIS_ESTIM   = KEEP(15)
        ENDIF
        MAXIS = max( 1,
     &       MAXIS_ESTIM + 3 * max(PERLU,10) * 
     &          ( MAXIS_ESTIM / 100 + 1 )
     &  )
C       ----------------------------
C       Allocate PTLUST_S and PTRFAC
C       They will be used to access
C       factors in the solve phase.
C       ----------------------------
        ALLOCATE( id%PTLUST_S( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=id%KEEP(28)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for id%PTLUST_S(', id%KEEP(28),')'
          ENDIF
          NULLIFY(id%PTLUST_S)
          GOTO 110
        END IF
        ALLOCATE( id%PTRFAC( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=id%KEEP(28)
          NULLIFY(id%PTRFAC)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for id%PTRFAC(', id%KEEP(28),')'
          ENDIF
          GOTO 110
        END IF
C       -----------------------------
C       Reserve temporary workspace :
C       IPOOL, PTRWB, ITLOC, PTRIST
C       PTRWB will be subdivided again
C       in routine DMUMPS_FAC_B
C       -----------------------------
        PTRIST = 1
        PTRWB  = PTRIST + id%KEEP(28)
        ITLOC  = PTRWB  + 2 * id%KEEP(28)
C Fwd in facto: ITLOC of size id%N + id%KEEP(253)
        IPOOL  = ITLOC  + id%N + id%KEEP(253)
C
C       --------------------------------
C       NA(1) is an upperbound for LPOOL
C       --------------------------------
C       Structure of the pool:
C     ____________________________________________________
C    | Subtrees   |         | Top nodes           | 1 2 3 |
C     ----------------------------------------------------
        LPOOL = MUMPS_GET_POOL_LENGTH(id%NA(1), id%KEEP(1),id%KEEP8(1))
        ALLOCATE( IWK(  IPOOL + LPOOL - 1 ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=IPOOL + LPOOL - 1
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for IWK(',IPOOL+LPOOL-1,')'
          ENDIF
          GOTO 110
        END IF
        ALLOCATE(IWK8( 2 * id%KEEP(28)), stat = IERR)
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=2 * id%KEEP(28)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for IWKB(', 2*id%KEEP(28),')'
          ENDIF
          GOTO 110
        END IF
C
C  Return to SPMD
C
      ENDIF
C
 110  CONTINUE
C     ----------------
C     Broadcast errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C
      IF ( I_AM_SLAVE )  THEN
C       Store size of receive buffers in DMUMPS_LBUF module
        CALL DMUMPS_BUF_DIST_IRECV_SIZE( DMUMPS_LBUFR_BYTES )
        IF (PROK) THEN
          WRITE( MP, 170 ) MAXS, MAXIS, id%KEEP8(12), KEEP(15),
     &    id%KEEP8(26), id%KEEP8(27), id%KEEP8(11), KEEP(26), KEEP(27)
        ENDIF
      END IF
C     ===============================================================
C     Before calling the main driver, DMUMPS_FAC_B,
C     some statistics should be initialized to 0,
C     even on the host node because they will be
C     used in REDUCE operations afterwards.
C     --------------------------------------------
C     Size of factors written. It will be set to POSFAC in
C     IC, otherwise we accumulate written factors in it.
      id%KEEP8(31)= 0_8
C     Size of factors under L0 will be returned
C     in id%KEEP8(64), not included in KEEP8(31))
C     Number of entries in factors
      id%KEEP8(10) = 0_8
C     KEEP8(8) will hold the volume of extra copies due to
C              in-place stacking in fac_mem_stack.F
      id%KEEP8(8)=0_8
      id%INFO(9:14)=0
      RINFO(2:3)=ZERO
      IF ( I_AM_SLAVE ) THEN
C       ------------------------------------
C       Call effective factorization routine
C       ------------------------------------
        IF ( KEEP(55) .eq. 0 ) THEN
          LDPTRAR = id%N
        ELSE
          LDPTRAR = id%NELT + 1
        END IF
        IF ( id%KEEP(55) .NE. 0 ) THEN
          NELT_arg = id%NELT
        ELSE
C         ------------------------------
C         Use size 1 to avoid complaints
C         when using check bound options
C         ------------------------------
          NELT_arg = 1
        END IF
      ENDIF
C     Compute DKEEP(17) 
      AVG_FLOPS    =  RINFOG(1)/(dble(id%NSLAVES))
      id%DKEEP(17) = max ( id%DKEEP(18), AVG_FLOPS/dble(50) )
     &               
      IF (PROK.AND.id%MYID.EQ.MASTER) THEN
        IF (id%NSLAVES.LE.1) THEN
         WRITE(MPG,'(/A,A,1PD10.3)')
     &' Start factorization with total', 
     &' estimated flops (RINFOG(1))                         = ', 
     &   RINFOG(1)
        ELSE
         WRITE(MP,'(/A,A,1PD10.3,A,1PD10.3)')
     &' Start factorization with total', 
     &' estimated flops RINFOG(1) / Average per MPI proc    = ', 
     &      RINFOG(1), ' / ', AVG_FLOPS
        ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
C       IS/S pointers passed to DMUMPS_FAC_B with
C       implicit interface through intermediate
C       structure S_IS_POINTERS. IS will be allocated
C       during DMUMPS_FAC_B.
        S_IS_POINTERS%IW => id%IS; NULLIFY(id%IS)
        S_IS_POINTERS%A  => id%S ; NULLIFY(id%S)
        CALL DMUMPS_FAC_B(id%N,S_IS_POINTERS,MAXS,MAXIS,id%SYM_PERM(1),
     &  id%NA(1),id%LNA,id%NE_STEPS(1),id%ND_STEPS(1), id%FILS(1),
     &  id%STEP(1),id%FRERE_STEPS(1),id%DAD_STEPS(1),id%CANDIDATES(1,1),
     &  id%ISTEP_TO_INIV2(1),id%TAB_POS_IN_PERE(1,1), id%PTRAR(1),
     &  LDPTRAR,IWK(PTRIST),id%PTLUST_S(1),id%PTRFAC(1),IWK(PTRWB),IWK8,
     &  IWK(ITLOC),RHS_MUMPS(1),IWK(IPOOL),LPOOL,CNTL1,ICNTL(1),
     &  id%INFO(1), RINFO(1),KEEP(1),id%KEEP8(1),id%PROCNODE_STEPS(1),
     &  id%NSLAVES,id%COMM_NODES,id%MYID,id%MYID_NODES,BUFR,DMUMPS_LBUFR
     &  , DMUMPS_LBUFR_BYTES, DMUMPS_LBUF, id%INTARR(1),id%DBLARR(1),
     &  id%root, NELT_arg, id%FRTPTR(1), id%FRTELT(1),id%COMM_LOAD,
     &  id%ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2, id%MEM_DIST(0),
     &  id%DKEEP(1), id%PIVNUL_LIST(1), LPN_LIST, id%LRGROUPS(1)
     &     )
        id%IS => S_IS_POINTERS%IW; NULLIFY(S_IS_POINTERS%IW)
        id%S  => S_IS_POINTERS%A ; NULLIFY(S_IS_POINTERS%A)
C
C       ------------------------------
C       Deallocate temporary workspace
C       ------------------------------
        DEALLOCATE( IWK  )
        DEALLOCATE( IWK8 )
      ENDIF
C     ---------------------------------
C     Free some workspace corresponding
C     to the original matrix in
C     arrowhead or elemental format.
C                  -----
C     Note : INTARR was not allocated
C     during factorization in the case
C     of an assembled matrix.
C     ---------------------------------
        IF ( KEEP(55) .eq. 0 ) THEN
C
C         ----------------
C         Assembled matrix
C         ----------------
          IF (associated( id%DBLARR)) THEN
            DEALLOCATE(id%DBLARR)
            NULLIFY(id%DBLARR)
          ENDIF
C
        ELSE
C
C         ----------------
C         Elemental matrix
C         ----------------
          IF (associated(id%INTARR)) THEN
            DEALLOCATE( id%INTARR)
            NULLIFY( id%INTARR )
          ENDIF
C         ------------------------------------
C         For the master from an hybrid host
C         execution without scaling, then real
C         values have not been copied !
C         -------------------------------------
          IF (      id%MYID_NODES .eq. MASTER
     &      .AND.   KEEP(46)   .eq. 1
     &      .AND.   KEEP(52)   .eq. 0 ) THEN
            NULLIFY( id%DBLARR )
          ELSE
            IF (associated( id%DBLARR)) THEN
              DEALLOCATE(id%DBLARR)
              NULLIFY(id%DBLARR)
            ENDIF
          END IF
        END IF
C     Memroy statistics
C     -----------------------------------
C     If QR (Keep(19)) is not zero, and if
C     the host does not have the information
C     (ie is not slave), send information
C     computed on the slaves during facto
C     to the host.
C     -----------------------------------
      IF ( KEEP(19) .NE. 0 ) THEN
        IF ( KEEP(46) .NE. 1 ) THEN
C         Host was not working during facto_root
C         Send him the information
          IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( KEEP(17), 1, MPI_INTEGER, 1, DEFIC_TAG,
     &                   id%COMM, STATUS, IERR )
          ELSE IF ( id%MYID .EQ. 1 ) THEN
            CALL MPI_SEND( KEEP(17), 1, MPI_INTEGER, 0, DEFIC_TAG,
     &                   id%COMM, IERR )
          END IF
        END IF
      END IF
C     --------------------------------
C     Deallocate communication buffers
C     They will be reallocated
C     in the solve.
C     --------------------------------
      IF (allocated(BUFR)) DEALLOCATE(BUFR)
      CALL DMUMPS_BUF_DEALL_SMALL_BUF( IERR )
C//PIV
      IF (KEEP(219).NE.0) THEN
      CALL DMUMPS_BUF_DEALL_MAX_ARRAY()
      ENDIF
C
C     Check for errors.
C     After DMUMPS_FAC_B every slave is aware of an error. 
C     If master is included in computations, the call below should
C     not be necessary.
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
C
      CALL DMUMPS_EXTRACT_SCHUR_REDRHS(id)
      IF (KEEP(201) .GT. 0) THEN
         IF ((KEEP(201).EQ.1) .OR. (KEEP(201).EQ.2)) THEN
            IF ( I_AM_SLAVE ) THEN
               CALL DMUMPS_OOC_CLEAN_PENDING(IERR)
               IF(IERR.LT.0)THEN
                  id%INFO(1)=IERR
                  id%INFO(2)=0
               ENDIF
            ENDIF
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
C             We want to collect statistics even in case of 
C             error to understand if it is due to numerical
C             issues
CC            IF ( id%INFO(1) < 0 ) GOTO 500
         END IF
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(94)=TIME
      ENDIF
C     =====================================================================
C     COMPUTE MEMORY ALLOCATED BY MUMPS, INFO(16)
C     ---------------------------------------------
      MEM_EFF_ALLOCATED = .TRUE.
      CALL DMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     BLR_STRAT, .TRUE., TOTAL_BYTES,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .FALSE.   ! UNDER_L0_OMP
     &     )
       IF (id%KEEP8(24).NE.0) THEN
C         WK_USER is not part of memory allocated by MUMPS
C         and is not counted, id%KEEP8(23) should be zero
          id%INFO(16) = TOTAL_MBYTES
       ELSE
C         Note that even for the case of ICNTL(23)>0 
C         we report here the memory effectively allocated
C         that can be smaller than ICNTL(23) !
          id%INFO(16) = TOTAL_MBYTES
       ENDIF
C       ----------------------------------------------------
C       Centralize memory statistics on the host
C       id%INFOG(18) = size of mem in Mbytes for facto,
C                   for the processor using largest memory
C       id%INFOG(19) = size of mem in Mbytes for facto,
C                   sum over all processors
C       ----------------------------------------------------
        CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                           id%INFO(16), id%INFOG(18), IRANK )
        CALL DMUMPS_PRINT_ALLOCATED_MEM( PROK, PROKG, PRINT_MAXAVG,
     &       MP, MPG, id%INFO(16), id%INFOG(18), id%INFOG(19),
     &       id%NSLAVES, IRANK,
     &       id%KEEP(1) )
C FIXME Check if WK_USER used and indicate, total space to WK_USER
      IF (PROK ) THEN
          WRITE(MP,'(A,I12) ')
     & ' ** Eff. min. Space MBYTES for facto             (INFO(16)):',
     &                TOTAL_MBYTES
      ENDIF
C     ========================(INFO(16) RELATED)======================
C     ---------------------------------------
C     COMPUTE EFFECTIVE MEMORY USED INFO(22)
C     ---------------------------------------
      PERLU_ON = .TRUE.
      MEM_EFF_ALLOCATED = .FALSE.
      CALL DMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     BLR_STRAT, PERLU_ON, TOTAL_BYTES,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .FALSE.   ! UNDER_L0_OMP
     &     )
C     -- TOTAL_BYTES and TOTAL_MBYTES includes both static 
C     -- (MAXS) and BLR structures computed as the SUM of the PEAKS 
C     -- (KEEP8(67) + KEEP8(70))
      id%KEEP8(7) = TOTAL_BYTES
C     -- INFO(22) holds the effective space (in Mbytes) used by MUMPS
C     -- (it includes part of WK_USER used if provided by user)
      id%INFO(22) = TOTAL_MBYTES
C     ----------------------------------------------------
C     Centralize memory statistics on the host
C       INFOG(21) = size of effective mem (Mbytes) for facto,
C                   for the processor using largest memory
C       INFOG(22) = size of effective mem (Mbytes) for facto,
C                   sum over all processors
C     ----------------------------------------------------
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                    id%INFO(22), id%INFOG(21), IRANK )
      IF ( PROKG ) THEN
       IF (PRINT_MAXAVG) THEN
        WRITE( MPG,'(A,I12) ')
     &  ' ** Memory effectively used, max in  Mbytes     (INFOG(21)):',
     &  id%INFOG(21)
       ENDIF
       WRITE( MPG,'(A,I12) ')
     &  ' ** Memory effectively used, total in Mbytes    (INFOG(22)):',
     &  id%INFOG(22)
      END IF
C
      IF (I_AM_SLAVE) THEN
       K67 = id%KEEP8(67)
       K68 = id%KEEP8(68)
       K70 = id%KEEP8(70)
       K74 = id%KEEP8(74)
       K75 = id%KEEP8(75)
      ELSE
       K67 = 0_8
       K68 = 0_8
       K70 = 0_8
       K74 = 0_8
       K75 = 0_8
      ENDIF
C     -- Save the number of entries effectively used
C        in main working array S
      CALL MUMPS_SETI8TOI4(K67,id%INFO(21))
C
C
      IF ( PROKG ) THEN
          IF (id%INFO(1) .GE.0) THEN
            WRITE(MPG,180) id%DKEEP(94)
          ELSE
            WRITE(MPG,185) id%DKEEP(94)
          ENDIF
      ENDIF
C
C  Sum RINFO(2) : total number of flops for assemblies
C  Sum RINFO(3) : total number of flops for eliminations
C     Initialize RINFO(4) in case BLR was not activated
      RINFO(4) = RINFO(3)
C 
C  Should work even if the master does some work
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
C     Reduce needed to dimension small working array
C     on all procs during DMUMPS_GATHER_SOLUTION
      KEEP(247) = 0
      CALL MPI_REDUCE( KEEP(246), KEEP(247), 1, MPI_INTEGER, 
     &                 MPI_MAX, MASTER, id%COMM, IERR)
C
C     Reduce compression times: get max compression times
      CALL MPI_REDUCE( id%DKEEP(97), id%DKEEP(98), 1,
     &     MPI_DOUBLE_PRECISION,
     &     MPI_MAX, MASTER, id%COMM, IERR)
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MUMPS_REDUCEI8( id%KEEP8(31)+id%KEEP8(64),id%KEEP8(6),
     &                     MPI_SUM, MASTER, id%COMM )
C
      IF (id%MYID.EQ.0) THEN
C      In MegaBytes
       RINFOG(16) = dble(id%KEEP8(6)*int(KEEP(35),8))/dble(1D6)
       IF (KEEP(201).LE.0) THEN
        RINFOG(16) = ZERO
       ENDIF
      ENDIF
      CALL MUMPS_REDUCEI8( id%KEEP8(48),id%KEEP8(148), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(id%KEEP8(148), INFOG(9))
C
      CALL MPI_REDUCE( int(id%INFO(10),8), id%KEEP8(128),
     &                 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SETI8TOI4(id%KEEP8(128), id%INFOG(10))
      ENDIF
C     Use MPI_MAX for this one to get largest front size
      CALL MPI_ALLREDUCE( id%INFO(11), INFOG(11), 1, MPI_INTEGER,
     &                 MPI_MAX, id%COMM, IERR)
C     make maximum effective frontal size available on all procs
C     for solve phase
C     (Note that INFO(11) includes root size on root master)
      KEEP(133) = INFOG(11)
      CALL MPI_REDUCE( id%INFO(12), INFOG(12), 3, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( KEEP(103), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(229) = INFOG(25)
      CALL MPI_REDUCE( KEEP(105), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(230) = INFOG(25)
C
      id%INFO(25) = KEEP(98)
      CALL MPI_ALLREDUCE( id%INFO(25), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, id%COMM, IERR)
C     Extra copies due to in-place stacking
      CALL MUMPS_REDUCEI8( id%KEEP8(8), id%KEEP8(108), MPI_SUM,
     &                     MASTER, id%COMM )
C     Entries in factors
      CALL MUMPS_SETI8TOI4(id%KEEP8(10), id%INFO(27))
      CALL MUMPS_REDUCEI8( id%KEEP8(10),id%KEEP8(110), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(id%KEEP8(110), INFOG(29))
C     Initialize INFO(28)/INFOG(35) in case BLR not activated
      id%INFO(28)  = id%INFO(27)
      INFOG(35)    = INFOG(29)
C     ==============================
C     LOW-RANK
C     ==============================
      IF ( KEEP(486) .NE. 0 ) THEN  !LR is activated
C Compute and Save local amount of flops in case of BLR
            RINFO(4) = dble(FLOP_FRFRONTS + FLOP_FACTO_FR - FLOP_LRGAIN
     &                 + FLOP_COMPRESS  + FLOP_FRFRONTS)
C
C Compute and Save local number of entries in compressed factors
C 
            ITMP8 =  id%KEEP8(10) - int(MRY_LU_LRGAIN,8)
            CALL MUMPS_SETI8TOI4( ITMP8, id%INFO(28))
C
            CALL MPI_REDUCE( MRY_LU_LRGAIN, TMP_MRY_LU_LRGAIN
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( MRY_LU_FR, TMP_MRY_LU_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( MRY_CB_FR, TMP_MRY_CB_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( MRY_CB_LRGAIN, TMP_MRY_CB_LRGAIN
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_LRGAIN, TMP_FLOP_LRGAIN
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_TRSM_FR, TMP_FLOP_TRSM_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_TRSM_LR, TMP_FLOP_TRSM_LR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_UPDATE_FR, TMP_FLOP_UPDATE_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_UPDATE_LR, TMP_FLOP_UPDATE_LR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_FRSWAP_COMPRESS,
     &                       TMP_FLOP_FRSWAP_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_MIDBLK_COMPRESS,
     &                       TMP_FLOP_MIDBLK_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_UPDATE_LRLR3, TMP_FLOP_UPDATE_LRLR3
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(FLOP_ACCUM_COMPRESS, TMP_FLOP_ACCUM_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_TRSM, TMP_FLOP_TRSM
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_PANEL, TMP_FLOP_PANEL
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_FRFRONTS, TMP_FLOP_FRFRONTS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_COMPRESS, TMP_FLOP_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_DECOMPRESS, TMP_FLOP_DECOMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_CB_COMPRESS, TMP_FLOP_CB_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_CB_DECOMPRESS,TMP_FLOP_CB_DECOMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_FACTO_FR, TMP_FLOP_FACTO_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_SOLFWD_FR, TMP_FLOP_SOLFWD_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_SOLFWD_LR, TMP_FLOP_SOLFWD_LR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( CNT_NODES,TMP_CNT_NODES
     &                      , 1, MPI_INTEGER,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            IF (id%NPROCS.GT.1) THEN
              FLOP_FACTO_LR = FLOP_FACTO_FR - FLOP_LRGAIN
     &                          + FLOP_COMPRESS + FLOP_FRFRONTS
              CALL MPI_REDUCE( FLOP_FACTO_LR, AVG_FLOP_FACTO_LR
     &                        , 1, MPI_DOUBLE_PRECISION,
     &                         MPI_SUM, MASTER, id%COMM, IERR)
              IF (id%MYID.EQ.MASTER) THEN
                AVG_FLOP_FACTO_LR = AVG_FLOP_FACTO_LR/id%NPROCS
              ENDIF
              CALL MPI_REDUCE( FLOP_FACTO_LR, MIN_FLOP_FACTO_LR
     &                        , 1, MPI_DOUBLE_PRECISION,
     &                         MPI_MIN, MASTER, id%COMM, IERR)
              CALL MPI_REDUCE( FLOP_FACTO_LR, MAX_FLOP_FACTO_LR
     &                        , 1, MPI_DOUBLE_PRECISION,
     &                         MPI_MAX, MASTER, id%COMM, IERR)
            ENDIF ! NPROCS > 1
            CALL MPI_REDUCE( TIME_UPDATE, TMP_TIME_UPDATE
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_LRLR1, TMP_TIME_UPDATE_LRLR1
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_LRLR2, TMP_TIME_UPDATE_LRLR2
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_LRLR3, TMP_TIME_UPDATE_LRLR3
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_FRLR, TMP_TIME_UPDATE_FRLR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_FRFR, TMP_TIME_UPDATE_FRFR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DIAGCOPY, TMP_TIME_DIAGCOPY
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_COMPRESS,TMP_TIME_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_MIDBLK_COMPRESS, 
     &                       TMP_TIME_MIDBLK_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FRSWAP_COMPRESS,
     &                       TMP_TIME_FRSWAP_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_CB_COMPRESS, TMP_TIME_CB_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION, 
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DECOMP, TMP_TIME_DECOMP
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DECOMP_UCFS, TMP_TIME_DECOMP_UCFS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DECOMP_ASM1, TMP_TIME_DECOMP_ASM1
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_DECOMP_LOCASM2, TMP_TIME_DECOMP_LOCASM2
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_DECOMP_MAPLIG1, TMP_TIME_DECOMP_MAPLIG1
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DECOMP_ASMS2S, TMP_TIME_DECOMP_ASMS2S
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DECOMP_ASMS2M, TMP_TIME_DECOMP_ASMS2M
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_PANEL, TMP_TIME_PANEL
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FAC_I, TMP_TIME_FAC_I
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FAC_MQ, TMP_TIME_FAC_MQ
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FAC_SQ, TMP_TIME_FAC_SQ
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_LRTRSM, TMP_TIME_LRTRSM
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FRTRSM, TMP_TIME_FRTRSM
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FRFRONTS, TMP_TIME_FRFRONTS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_LR_MODULE, TMP_TIME_LR_MODULE
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
        IF (id%MYID.EQ.MASTER) THEN
         IF (id%NPROCS.GT.1) THEN
C rename the stat variable so that COMPUTE_GLOBAL_GAINS can work for any
C number of procs
            MRY_LU_FR            = TMP_MRY_LU_FR
            MRY_LU_LRGAIN        = TMP_MRY_LU_LRGAIN
            MRY_CB_FR            = TMP_MRY_CB_FR
            MRY_CB_LRGAIN        = TMP_MRY_CB_LRGAIN
            FLOP_LRGAIN          = TMP_FLOP_LRGAIN
            FLOP_PANEL           = TMP_FLOP_PANEL
            FLOP_TRSM            = TMP_FLOP_TRSM
            FLOP_TRSM_FR         = TMP_FLOP_TRSM_FR
            FLOP_TRSM_LR         = TMP_FLOP_TRSM_LR
            FLOP_UPDATE_FR       = TMP_FLOP_UPDATE_FR
            FLOP_UPDATE_LR       = TMP_FLOP_UPDATE_LR
            FLOP_UPDATE_LRLR3    = TMP_FLOP_UPDATE_LRLR3
            FLOP_COMPRESS        = TMP_FLOP_COMPRESS
            FLOP_MIDBLK_COMPRESS = TMP_FLOP_MIDBLK_COMPRESS
            FLOP_FRSWAP_COMPRESS = TMP_FLOP_FRSWAP_COMPRESS
            FLOP_ACCUM_COMPRESS  = TMP_FLOP_ACCUM_COMPRESS
            FLOP_CB_COMPRESS     = TMP_FLOP_CB_COMPRESS
            FLOP_DECOMPRESS      = TMP_FLOP_DECOMPRESS
            FLOP_CB_DECOMPRESS   = TMP_FLOP_CB_DECOMPRESS
            FLOP_FRFRONTS        = TMP_FLOP_FRFRONTS
            FLOP_SOLFWD_FR       = TMP_FLOP_SOLFWD_FR
            FLOP_SOLFWD_LR       = TMP_FLOP_SOLFWD_LR
            FLOP_FACTO_FR        = TMP_FLOP_FACTO_FR
            CNT_NODES            = TMP_CNT_NODES
            TIME_UPDATE          = TMP_TIME_UPDATE         /id%NPROCS
            TIME_UPDATE_LRLR1    = TMP_TIME_UPDATE_LRLR1   /id%NPROCS
            TIME_UPDATE_LRLR2    = TMP_TIME_UPDATE_LRLR2   /id%NPROCS
            TIME_UPDATE_LRLR3    = TMP_TIME_UPDATE_LRLR3   /id%NPROCS
            TIME_UPDATE_FRLR     = TMP_TIME_UPDATE_FRLR    /id%NPROCS
            TIME_UPDATE_FRFR     = TMP_TIME_UPDATE_FRFR    /id%NPROCS
            TIME_COMPRESS        = TMP_TIME_COMPRESS       /id%NPROCS
            TIME_MIDBLK_COMPRESS = TMP_TIME_MIDBLK_COMPRESS/id%NPROCS
            TIME_FRSWAP_COMPRESS = TMP_TIME_FRSWAP_COMPRESS/id%NPROCS
            TIME_DIAGCOPY        = TMP_TIME_DIAGCOPY       /id%NPROCS
            TIME_CB_COMPRESS     = TMP_TIME_CB_COMPRESS    /id%NPROCS
            TIME_PANEL           = TMP_TIME_PANEL          /id%NPROCS
            TIME_FAC_I           = TMP_TIME_FAC_I          /id%NPROCS
            TIME_FAC_MQ          = TMP_TIME_FAC_MQ         /id%NPROCS
            TIME_FAC_SQ          = TMP_TIME_FAC_SQ         /id%NPROCS
            TIME_LRTRSM          = TMP_TIME_LRTRSM         /id%NPROCS
            TIME_FRTRSM          = TMP_TIME_FRTRSM         /id%NPROCS
            TIME_FRFRONTS        = TMP_TIME_FRFRONTS       /id%NPROCS
            TIME_LR_MODULE       = TMP_TIME_LR_MODULE      /id%NPROCS
            TIME_DECOMP          = TMP_TIME_DECOMP         /id%NPROCS
            TIME_DECOMP_UCFS     = TMP_TIME_DECOMP_UCFS    /id%NPROCS
            TIME_DECOMP_ASM1     = TMP_TIME_DECOMP_ASM1    /id%NPROCS
            TIME_DECOMP_LOCASM2  = TMP_TIME_DECOMP_LOCASM2 /id%NPROCS
            TIME_DECOMP_MAPLIG1  = TMP_TIME_DECOMP_MAPLIG1 /id%NPROCS
            TIME_DECOMP_ASMS2S   = TMP_TIME_DECOMP_ASMS2S  /id%NPROCS
            TIME_DECOMP_ASMS2M   = TMP_TIME_DECOMP_ASMS2M  /id%NPROCS
         ENDIF
         CALL COMPUTE_GLOBAL_GAINS(id%KEEP8(110),id%RINFOG(3),
     &        id%KEEP8(49), PROKG, MPG)
C         Number of entries in factor  INFOG(35) in 
C         compressed form is updated as long as 
C         BLR is activated,  this independently of the 
C         fact that factors are saved in LR.
          CALL MUMPS_SETI8TOI4(id%KEEP8(49),  id%INFOG(35))
          FRONTWISE = 0
C         WRITE gains also compute stats stored in DKEEP array
          IF (LPOK) THEN
            IF (CNTL(7) < 0.0D0) THEN
C           Warning : using negative values is an experimental and 
C            non recommended setting.
             WRITE(LP,'(/A/,A/,A/,A,A)') 
     &  ' WARNING in BLR input setting',
     &  '          CNTL(7) < 0 is experimental: ',
     &  '          RRQR precision = |CNTL(7| x ||A_pre||, ',
     &  '          where A_pre is the preprocessed matrix as defined',
     &  ' in the Users guide '
            ENDIF
          ENDIF
          CALL SAVEandWRITE_GAINS(FRONTWISE,
     &                KEEP(489), id%DKEEP, N,  id%ICNTL(36),
     &                KEEP(487), KEEP(488), KEEP(490),
     &                KEEP(491), KEEP(50), KEEP(486), KEEP(472),
     &                KEEP(475), KEEP(478), KEEP(480), KEEP(481), 
     &                KEEP(483), KEEP(484), 
     &                id%KEEP8(110), id%KEEP8(49),
     &                KEEP(28), id%NPROCS, MPG, PROKG)
C           flops when BLR activated
            RINFOG(14) = id%DKEEP(56)
          ELSE
            RINFOG(14) = 0.0D00
          ENDIF
      ENDIF
C     ==============================
C     NULL PIVOTS AND RANK-REVEALING
C     ==============================
      IF(KEEP(110) .EQ. 1) THEN
C        -- make available to users the local number of null pivots detected 
C        -- with ICNTL(24) = 1.
         id%INFO(18) = KEEP(109)
         CALL MPI_ALLREDUCE( KEEP(109), KEEP(112), 1, MPI_INTEGER,
     &        MPI_SUM, id%COMM, IERR)
      ELSE
         id%INFO(18)  = 0
         KEEP(109) = 0
         KEEP(112) = 0
      ENDIF
      IF (id%MYID.EQ.MASTER) THEN
C      INFOG(28) deficiency resulting from ICNTL(24) and ICNTL(56).
       INFOG(28)=KEEP(112)+KEEP(17)
      ENDIF
C     ========================================
C     We now provide to the host the part of
C     PIVNUL_LIST resulting from the processing
C     of the root node and we update id%INFO(18)
C     on the processor holding the root to
C     include null pivots relative to the root
C     ========================================
      IF (KEEP(17) .NE. 0) THEN
        IF (id%MYID .EQ. ID_ROOT) THEN
C         Include in id%INFO(18) null pivots resulting
C         from deficiency on the root. In this way,
C         the sum of all id%INFO(18) is equal to INFOG(28).
          id%INFO(18)=id%INFO(18)+KEEP(17)
        ENDIF
        IF (ID_ROOT .EQ. MASTER) THEN
          IF (id%MYID.EQ.MASTER) THEN
C           --------------------------------------------------
C           Null pivots of root have been stored in
C           PIVNUL_LIST(KEEP(109)+1:KEEP(109)+KEEP(17).
C           Shift them at the end of the list because:
C           * this is what we need to build the null space
C           * we would otherwise overwrite them on the host
C             when gathering null pivots from other processors
C           --------------------------------------------------
            DO I=1, KEEP(17)
              id%PIVNUL_LIST(KEEP(112)+I)=id%PIVNUL_LIST(KEEP(109)+I)
            ENDDO
          ENDIF
        ELSE
C         ---------------------------------
C         Null pivots of root must be sent
C         from the processor responsible of
C         the root to the host (or MASTER).
C         ---------------------------------
          IF (id%MYID .EQ. ID_ROOT) THEN
            CALL MPI_SEND(id%PIVNUL_LIST(KEEP(109)+1), KEEP(17),
     &                    MPI_INTEGER, MASTER, ZERO_PIV,
     &                    id%COMM, IERR)
          ELSE IF (id%MYID .EQ. MASTER) THEN
            CALL MPI_RECV(id%PIVNUL_LIST(KEEP(112)+1), KEEP(17),
     &                    MPI_INTEGER, ID_ROOT, ZERO_PIV,
     &                    id%COMM, STATUS, IERR )
          ENDIF
        ENDIF
      ENDIF
C     ===========================
C     gather zero pivots indices
C     on the host node
C     ===========================
C     In case of non working host, the following code also
C     works considering that KEEP(109) is equal to 0 on
C     the non-working host
      IF(KEEP(110) .EQ. 1) THEN
         ALLOCATE(ITMP2(id%NPROCS),stat = IERR )  ! deallocated in 490
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=id%NPROCS
         END IF
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 490
         CALL MPI_GATHER ( KEEP(109),1, MPI_INTEGER, 
     &        ITMP2(1), 1, MPI_INTEGER, 
     &        MASTER, id%COMM, IERR)
         IF(id%MYID .EQ. MASTER) THEN
            POSBUF = ITMP2(1)+1
C           First null pivot of master is in
C           position 1 of global list
            KEEP(220)=1
            DO I = 1,id%NPROCS-1
               CALL MPI_RECV(id%PIVNUL_LIST(POSBUF), ITMP2(I+1), 
     &              MPI_INTEGER,I, 
     &              ZERO_PIV, id%COMM, STATUS, IERR)
C              Send position POSBUF of first null pivot of proc I
C              in global list. Will allow to quickly identify during
C              the solve step if one is concerned by a global position
C              K, 0 <= K <= INFOG(28).
               CALL MPI_SEND(POSBUF, 1, MPI_INTEGER, I, ZERO_PIV,
     &              id%COMM, IERR)
               POSBUF = POSBUF + ITMP2(I+1)
            ENDDO
         ELSE
            CALL MPI_SEND( id%PIVNUL_LIST(1), KEEP(109), MPI_INTEGER,
     &           MASTER,ZERO_PIV, id%COMM, IERR)
            CALL MPI_RECV( KEEP(220), 1, MPI_INTEGER, MASTER, ZERO_PIV,
     &           id%COMM, STATUS, IERR )
         ENDIF
      ENDIF
C     =====================================
C     Statistics relative to min/max pivots
C     =====================================
      CALL MPI_REDUCE( id%DKEEP(19), RINFOG(19), 1, 
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_MIN, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( id%DKEEP(20), RINFOG(20), 1, 
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_MIN, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( id%DKEEP(21), RINFOG(21), 1, 
     &                 MPI_DOUBLE_PRECISION,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
C     =========================================
C     Centralized number of swaps for pivoting
C     =========================================
      CALL MPI_REDUCE( id%KEEP8(80), ITEMP8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF (id%MYID .EQ. MASTER) THEN
        CALL MUMPS_SETI8TOI4(ITEMP8,id%INFOG(48))
      ENDIF
C     ==========================================
C     Centralized largest increase of panel size
C     ==========================================
      CALL MPI_REDUCE( id%KEEP(425), id%INFOG(49), 1, MPI_INTEGER,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
C     =====================================
C     Statistics concerning the determinant
C     =====================================
C
C     1/ on the host better take into account null pivots if scaling:
C
C     Since null pivots are excluded from the computation
C     of the determinant, we also exclude the corresponding
C     scaling entries. Since those entries have already been
C     taken into account before the factorization, we multiply
C     the determinant on the host by the scaling values corresponding
C     to pivots in PIVNUL_LIST.
      IF (id%MYID.EQ.MASTER .AND. LSCAL. AND. KEEP(258).NE.0) THEN
        DO I = 1, id%INFOG(28)
          CALL DMUMPS_UPDATEDETER(id%ROWSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
          CALL DMUMPS_UPDATEDETER(id%COLSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
        ENDDO
      ENDIF
C
C     2/ Swap signs depending on pivoting on each proc
C
      IF (KEEP(258).NE.0) THEN
C       Return the determinant in INFOG(34) and RINFOG(12/13)
C       In case of real arithmetic, initialize
C       RINFOG(13) to 0 (no imaginary part and
C       not touched by DMUMPS_DETER_REDUCTION)
        RINFOG(13)=0.0D0
        IF (KEEP(260).EQ.-1) THEN ! Local to each processor
          id%DKEEP(6)=-id%DKEEP(6)
        ENDIF
C
C       3/ Perform a reduction
C
        CALL DMUMPS_DETER_REDUCTION(
     &           id%COMM, id%DKEEP(6), KEEP(259),
     &           RINFOG(12), INFOG(34), id%NPROCS)
C
C       4/ Swap sign if needed
C
        IF (id%KEEP(50).EQ.0 .AND. id%MYID.EQ. MASTER) THEN
C         Modify sign of determinant according
C         to unsymmetric permutation (max-trans
C         of max-weighted matching)
          IF (id%KEEP(23).NE.0) THEN
            CALL DMUMPS_DETER_SIGN_PERM(
     &           RINFOG(12), id%N,
C           id%STEP: used as workspace of size N still
C                    allocated on master; restored on exit
     &           id%STEP(1),
     &           id%UNS_PERM(1) )
C           Remark that RINFOG(12/13) are modified only
C           on the host but will be broadcast on exit
C           from MUMPS (see DMUMPS_DRIVER)
          ENDIF
        ENDIF
      ENDIF
 490  IF (allocated(ITMP2)) DEALLOCATE(ITMP2)
      IF ( PROKG ) THEN
C     -----------------------------
C     PRINT STATISTICS  (on master)
C     -----------------------------
          WRITE(MPG,99984) RINFOG(2),RINFOG(3),KEEP(52), 
     &                    id%KEEP8(148),
     &                    id%KEEP8(128), INFOG(11), id%KEEP8(110)
          IF (id%KEEP(50) == 1 .OR. id%KEEP(50) == 2) THEN
            ! negative pivots
            WRITE(MPG, 99987) INFOG(12)
          END IF
          IF (id%KEEP(50) == 0) THEN
            ! off diag pivots
            WRITE(MPG, 99985) INFOG(12)
          END IF
          IF (id%KEEP(50) .NE. 1) THEN
            ! delayed pivots
            WRITE(MPG, 99982) INFOG(13)
          END IF
          IF (KEEP(97) .NE. 0) THEN
            ! tiny pivots
            WRITE(MPG, 99986) INFOG(25)
          ENDIF
          IF (id%KEEP(50) == 2) THEN
            !number of 2x2 pivots in type 1 nodes
             WRITE(MPG, 99988) KEEP(229)
            !number of 2x2 pivots in type 2 nodes
             WRITE(MPG, 99989) KEEP(230)
          ENDIF
          !number of zero pivots
          IF (KEEP(110) .NE.0) THEN
              WRITE(MPG, 99991) KEEP(112)
          ENDIF
          !Deficiency on root
          IF ( KEEP(19) .ne. 0 )
c         IF ( KEEP(17) .ne. 0 )
     &    WRITE(MPG, 99983) KEEP(17)
          !Total deficiency
          IF (KEEP(110).NE.0.OR.KEEP(19).NE.0)
c          IF (KEEP(110).NE.0.OR.KEEP(17).NE.0)
     &    WRITE(MPG, 99992) KEEP(17)+KEEP(112)
          ! Memory compress
          WRITE(MPG, 99981) INFOG(14)
          ! Extra copies due to ip stack in unsym case
          ! in core case (or OLD_OOC_PANEL)
          IF (id%KEEP8(108) .GT. 0_8) THEN
            WRITE(MPG, 99980) id%KEEP8(108)
          ENDIF
          IF  ((KEEP(60).NE.0) .AND. INFOG(25).GT.0) THEN
          !  Schur on and tiny pivots set in last level 
          ! before the Schur if KEEP(114)=0
           WRITE(MPG, '(A)') 
     & " ** Warning Static pivoting was necessary"
           WRITE(MPG, '(A)') 
     & " ** to factor interior variables with Schur ON"
          ENDIF
          IF (KEEP(258).NE.0) THEN
            WRITE(MPG,99978) RINFOG(12)
            WRITE(MPG,99977) INFOG(34)
          ENDIF
      END IF
* ==========================================
*
*  End of Factorization Phase
*
* ==========================================
C
C  Goto 500 is done when
C  LOAD_INIT
C  OOC_INIT_FACTO
C  MUMPS_FDM_INIT
#if ! defined(NO_FDM_DESCBAND)
C  MUMPS_FDBD_INIT
#endif
#if ! defined(NO_FDM_MAPROW)
C  MUMPS_FMRD_INIT
#endif
C  are all called.
C
 500  CONTINUE
C     Redo free DBLARR (as in end_driver.F)
C     in case an error occurred after allocating
C     DBLARR and before freeing it above.
      IF (id%KEEP(46).EQ.1 .AND.
     &    id%KEEP(55).NE.0 .AND.
     &    id%MYID.EQ.MASTER .AND.
     &    id%KEEP(52) .EQ. 0) THEN
        NULLIFY(id%DBLARR)
      ELSE
        IF (associated(id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
      ENDIF
#if ! defined(NO_FDM_DESCBAND)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDBD_END(id%INFO(1))  ! INFO(1): input only
      ENDIF
#endif
#if ! defined(NO_FDM_MAPROW)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FMRD_END(id%INFO(1))  ! INFO(1): input only
      ENDIF
#endif
      IF (I_AM_SLAVE) THEN
C       Terminate BLR module except if it is still needed for solve
        IF ( 
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND. id%INFO(1).GE.0
     &     ) THEN
C         Store pointer to BLR_ARRAY in MUMPS structure
C         (requires successful factorization otherwise module is freed)
          CALL DMUMPS_BLR_MOD_TO_STRUC(id%BLRARRAY_ENCODING)
        ELSE
C         INFO(1) positive or negative
          CALL DMUMPS_BLR_END_MODULE(id%INFO(1), id%KEEP8)
        ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDM_END('A')
C       Terminate BLR module except if it is still needed for solve
        IF ( 
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND. id%INFO(1).GE.0
     &     ) THEN
           CALL MUMPS_FDM_MOD_TO_STRUC('F', id%FDM_F_ENCODING,
     &        id%INFO(1))
           IF (.NOT. associated(id%FDM_F_ENCODING)) THEN
             WRITE(*,*) "Internal error 2 in DMUMPS_FAC_DRIVER"
           ENDIF
        ELSE
           CALL MUMPS_FDM_END('F')
        ENDIF
      ENDIF
C
C  Goto 514 is done when an
C  error occurred in MUMPS_FDM_INIT
C  or (after FDM_INIT but before
C  OOC_INIT)
C
 514  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
            CALL DMUMPS_OOC_END_FACTO(id,IERR)
            IF (id%ASSOCIATED_OOC_FILES) THEN
              id%ASSOCIATED_OOC_FILES = .FALSE.
            ENDIF
            IF (IERR.LT.0 .AND. id%INFO(1) .GE. 0) id%INFO(1) = IERR
         ENDIF
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always freed when WK_USER provided
            NULLIFY(id%S)
         ELSE IF (KEEP(201).NE.0) THEN
C           ----------------------------------------
C           In OOC or if KEEP(201).EQ.-1 we always
C           free S at end of factorization. As id%S
C           may be unassociated in case of error
C           during or before the allocation of id%S,
C           we only free S when it was associated.
C           ----------------------------------------
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   ! in all cases
            id%KEEP8(23)=0_8
         ENDIF
      ELSE  ! host not working
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always freed when WK_USER provided
            NULLIFY(id%S)
         ELSE
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   ! in all cases
            id%KEEP8(23)=0_8
         END IF
      END IF
C
C     Goto 513 is done in case of error where LOAD_INIT was
C     called but not OOC_INIT_FACTO.
 513  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         CALL DMUMPS_LOAD_END( id%INFO(1), id%NSLAVES, IERR )
         IF (IERR.LT.0 .AND. id%INFO(1) .GE. 0) id%INFO(1) = IERR
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
C
C     Goto 517 is done when an error occurs when GPU initialization
C     has been performed but not LOAD_INIT or OOC_INIT_FACTO
C
 517  CONTINUE
C
C     Goto 530 is done when an error occurs before
C     the calls to GPU_INIT, LOAD_INIT and OOC_INIT_FACTO
 530  CONTINUE
C  Fwd in facto: free RHS_MUMPS in case
C  it was allocated.
      IF (RHS_MUMPS_ALLOCATED) DEALLOCATE(RHS_MUMPS)
      NULLIFY(RHS_MUMPS)
C
      id%KEEP8(26) = KEEP826_SAVE
      RETURN
 120  FORMAT(/' Local redistrib: data local/sent           =',I16,I16)
 125  FORMAT(/' Redistrib: total data local/sent           =',I16,I16)
 130  FORMAT(//'****** FACTORIZATION STEP ********'/)
 160  FORMAT(
     & /' Elapsed time to reformat/distribute matrix =',F12.4)
 166  FORMAT(' Max difference from 1 after scaling the entries',
     &       ' for ONE-NORM (option 7/8)   =',D9.2)
 170  FORMAT(' STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Size of internal working array S           =',I16/
     &        ' Size of internal working array IS          =',I16/
     &        ' Minimum (ICNTL(14)=0) size of S            =',I16/
     &        ' Minimum (ICNTL(14)=0) size of IS           =',I16/
     &        ' Real space for original matrix             =',I16/
     &        ' Integer space for original matrix          =',I16/
     &        ' INFO(3) Real space for factors (estimated) =',I16/ 
     &        ' INFO(4) Integer space for factors (estim.) =',I16/
     &        ' Maximum frontal size (estimated)           =',I16)
 172  FORMAT(' GLOBAL STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Number of working processes                =',I16/
     &        ' ICNTL(22) Out-of-core option               =',I16/
     &        ' ICNTL(35) BLR activation (eff. choice)     =',I16/
     &        ' ICNTL(14) Memory relaxation                =',I16/
     &        ' INFOG(3) Real space for factors (estimated)=',I16/
     &        ' INFOG(4) Integer space for factors (estim.)=',I16/
     &        ' Maximum frontal size (estimated)           =',I16/
     &        ' Number of nodes in the tree                =',I16/
     &        ' Memory allowed (MB -- 0: N/A )             =',I16/
     &        ' Memory provided by user, sum of LWK_USER   =',I16/
     &        ' Effective threshold for pivoting, CNTL(1)  =',D16.4)
 173  FORMAT( ' Perform forward during facto, NRHS         =',I16)
 174  FORMAT( ' KEEP(268) Relaxed pivoting effective value =',I16)
 180  FORMAT(/' Elapsed time for factorization             =',F12.4)
 185  FORMAT(/' Elapsed time for (failed) factorization    =',F12.4)
99977 FORMAT( ' INFOG(34)  Determinant (base 2 exponent)   =',I16)
99978 FORMAT( ' RINFOG(12) Determinant (real part)         =',F16.8)
99980 FORMAT( ' Extra copies due to In-Place stacking      =',I16)
99981 FORMAT( ' INFOG(14)  Number of memory compress       =',I16)
99982 FORMAT( ' INFOG(13)  Number of delayed pivots        =',I16)
99983 FORMAT( ' Nb of singularities detected by ICNTL(56)  =',I16)
99991 FORMAT( ' Nb of null pivots detected by ICNTL(24)    =',I16)
99992 FORMAT( ' INFOG(28)  Estimated deficiency            =',I16)
99984 FORMAT(/'Leaving factorization with ...'/
     &        ' RINFOG(2)  Operations in node assembly     =',1PD10.3/
     &        ' ------(3)  Operations in node elimination  =',1PD10.3/
     &        ' ICNTL (8)  Scaling effectively used        =',I16/
     &        ' INFOG (9)  Real space for factors          =',I16/
     &        ' INFOG(10)  Integer space for factors       =',I16/
     &        ' INFOG(11)  Maximum front size              =',I16/
     &        ' INFOG(29)  Number of entries in factors    =',I16)
99985 FORMAT( ' INFOG(12)  Number of off diagonal pivots   =',I16)
99986 FORMAT( ' INFOG(25)  Number of tiny pivots(static)   =',I16)
99987 FORMAT( ' INFOG(12)  Number of negative pivots       =',I16)
99988 FORMAT( ' Number of 2x2 pivots in type 1 nodes       =',I16)
99989 FORMAT( ' Number of 2x2 pivots in type 2 nodes       =',I16)
      END SUBROUTINE DMUMPS_FAC_DRIVER
C
      SUBROUTINE DMUMPS_PRINT_ALLOCATED_MEM( PROK, PROKG, PRINT_MAXAVG,
     &       MP, MPG, INFO16, INFOG18, INFOG19, NSLAVES, IRANK, KEEP )
      IMPLICIT NONE
C
C  Purpose:
C  =======
C     Print memory allocated during factorization
C     - called at beginning of factorization in full-rank
C     - called at end of factorization in low-rank (because
C       of dynamic allocations)
C
      LOGICAL, INTENT(IN) :: PROK, PROKG, PRINT_MAXAVG
      INTEGER, INTENT(IN) :: MP, MPG, INFO16, INFOG18, INFOG19
      INTEGER, INTENT(IN) :: IRANK, NSLAVES
      INTEGER, INTENT(IN) :: KEEP(500)
C
      IF ( PROKG ) THEN
        IF (PRINT_MAXAVG) THEN
         WRITE( MPG,'(A,I12) ')
     &   ' ** Memory allocated, max in Mbytes             (INFOG(18)):',
     &   INFOG18
        ENDIF
        WRITE( MPG,'(/A,I12) ')
     &   ' ** Memory allocated, total in Mbytes           (INFOG(19)):',
     &    INFOG19
      END IF
      RETURN
      END SUBROUTINE DMUMPS_PRINT_ALLOCATED_MEM
      SUBROUTINE DMUMPS_AVGMAX_STAT8(PROKG, MPG, VAL, NSLAVES,
     &     PRINT_MAXAVG, COMM, MSG)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL, intent(in) :: PROKG
      INTEGER, intent(in) :: MPG
      INTEGER(8), intent(in) :: VAL
      INTEGER, intent(in) :: NSLAVES
      LOGICAL, intent(in) :: PRINT_MAXAVG
      INTEGER, intent(in) :: COMM
      CHARACTER*48 MSG 
C  Local
      INTEGER(8) MAX_VAL
      INTEGER IERR, MASTER
      DOUBLE PRECISION LOC_VAL, AVG_VAL
      PARAMETER(MASTER=0)
C
      CALL MUMPS_REDUCEI8( VAL, MAX_VAL, MPI_MAX, MASTER, COMM)
      LOC_VAL = dble(VAL)/dble(NSLAVES)
      CALL MPI_REDUCE( LOC_VAL, AVG_VAL, 1, MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MASTER, COMM, IERR )
      IF (PROKG) THEN
        IF (PRINT_MAXAVG) THEN
          WRITE(MPG,100) " Average", MSG, int(AVG_VAL,8)
        ELSE
          WRITE(MPG,110)  MSG, MAX_VAL
        ENDIF
      ENDIF
      RETURN
 100  FORMAT(A8,A48,I18)
 110  FORMAT(A48,I18)
      END SUBROUTINE DMUMPS_AVGMAX_STAT8
C
      SUBROUTINE DMUMPS_EXTRACT_SCHUR_REDRHS(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C     Extract the Schur and possibly also the reduced right-hand side
C     (if Fwd in facto) from the processor working on Schur and copy
C     it into the user datastructures id%SCHUR and id%REDRHS on the host.
C     This routine assumes that the integer list of the Schur has not
C     been permuted and still corresponds to LISTVAR_SCHUR.
C
C     If the Schur is centralized, the master of the Schur holds the
C     Schur and possibly also the reduced right-hand side.
C     If the Schur is distribued (already built in user's datastructure),
C     then the master of the Schur may hold the reduced right-hand side,
C     in which case it is available in root%RHS_CNTR_MASTER_ROOT.
C     
      TYPE(DMUMPS_STRUC) :: id
C
C  Local variables
C  ===============
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
      INTEGER :: ID_SCHUR, SIZE_SCHUR, LD_SCHUR, IB, BL4
      INTEGER :: ROW_LENGTH, I
      INTEGER(8) :: SURFSCHUR8, BL8, SHIFT8
      INTEGER(8) :: ISCHUR_SRC, ISCHUR_DEST, ISCHUR_SYM, ISCHUR_UNS
C
C  External functions
C  ==================
C
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
C     Quick return in case factorization did not terminate correctly
      IF (id%INFO(1) .LT. 0) RETURN
C     Quick return if Schur option off
      IF (id%KEEP(60) .EQ. 0) RETURN
C     Get Schur id
      ID_SCHUR =MUMPS_PROCNODE(
     &    id%PROCNODE_STEPS(id%STEP(max(id%KEEP(20),id%KEEP(38)))),
     &    id%KEEP(199))
      IF ( id%KEEP( 46 )  .NE. 1 ) THEN
        ID_SCHUR = ID_SCHUR + 1
      END IF
C     Get size of Schur
      IF (id%MYID.EQ.ID_SCHUR) THEN
        IF (id%KEEP(60).EQ.1) THEN
C         Sequential Schur
          LD_SCHUR =
     &    id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+2+id%KEEP(IXSZ))
          SIZE_SCHUR = LD_SCHUR - id%KEEP(253)
        ELSE
C         Parallel Schur
          LD_SCHUR   = -999999 ! not used
          SIZE_SCHUR = id%root%TOT_ROOT_SIZE
        ENDIF
      ELSE IF (id%MYID .EQ. MASTER) THEN
        SIZE_SCHUR = id%KEEP(116)
        LD_SCHUR = -44444 ! Not used
      ELSE
C       Proc is not concerned with Schur, return
        RETURN
      ENDIF
      SURFSCHUR8 = int(SIZE_SCHUR,8)*int(SIZE_SCHUR,8)
C     =================================
C     Case of parallel Schur: if REDRHS
C     was requested, obtain it directly
C     from id%root%RHS_CNTR_MASTER_ROOT
C     =================================
      IF (id%KEEP(60) .GT. 1) THEN
        IF (id%KEEP(221).EQ.1 .AND. id%KEEP(252).GT.0) THEN
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
              CALL dcopy(SIZE_SCHUR,
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1), 1,
     &             id%REDRHS((I-1)*id%LREDRHS+1), 1)
            ELSE
              IF (id%MYID.EQ.ID_SCHUR) THEN
C               Send
                CALL MPI_SEND(
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1),
     &             SIZE_SCHUR,
     &             MPI_DOUBLE_PRECISION,
     &             MASTER, TAG_SCHUR,
     &             id%COMM, IERR )
              ELSE ! MYID.EQ.MASTER
C               Receive
                CALL MPI_RECV( id%REDRHS((I-1)*id%LREDRHS+1),
     &             SIZE_SCHUR,
     &             MPI_DOUBLE_PRECISION, ID_SCHUR, TAG_SCHUR,
     &             id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
          ENDDO
C         ------------------------------
C         In case of parallel Schur, we
C         free root%RHS_CNTR_MASTER_ROOT
C         ------------------------------
          IF (id%MYID.EQ.ID_SCHUR) THEN
            DEALLOCATE(id%root%RHS_CNTR_MASTER_ROOT)
            NULLIFY   (id%root%RHS_CNTR_MASTER_ROOT)
          ENDIF
        ENDIF
C       return because this is all we need to do
C       in case of parallel Schur complement
        RETURN
      ENDIF
C     ============================
C     Centralized Schur complement
C     ============================
C     PTRAST has been freed at the moment of calling this
C     routine. Schur is available through
C     PTRFAC(IW( PTLUST_S( STEP(KEEP(20)) ) + 4 +KEEP(IXSZ) ))
      IF (id%KEEP(252).EQ.0) THEN
C       CASE 1 (ORIGINAL CODE):
C       Schur is contiguous on ID_SCHUR
        IF ( ID_SCHUR .EQ. MASTER ) THEN ! Necessarily equals id%MYID
C         ---------------------
C         Copy Schur complement
C         ---------------------
          CALL DMUMPS_COPYI8SIZE( SURFSCHUR8,
     &      id%S(id%PTRFAC(id%STEP(id%KEEP(20)))),
     &      id%SCHUR(1) )
        ELSE
C         -----------------------------------------
C         The processor responsible of the Schur
C         complement sends it to the host processor
C         -----------------------------------------
          BL8=int(huge(BL4)/id%KEEP(35)/10,8)
          DO IB=1, int((SURFSCHUR8+BL8-1_8) / BL8)
            SHIFT8 = int(IB-1,8) * BL8                ! Where to send
            BL4    = int(min(BL8,SURFSCHUR8-SHIFT8)) ! Size of block
            IF ( id%MYID .eq. ID_SCHUR ) THEN
C             Send Schur complement
              CALL MPI_SEND( id%S( SHIFT8 +
     &          id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ)))),
     &          BL4,
     &          MPI_DOUBLE_PRECISION,
     &          MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
            ELSE IF ( id%MYID .eq. MASTER ) THEN
C             Receive Schur complement
              CALL MPI_RECV( id%SCHUR(1_8 + SHIFT8),
     &                     BL4,
     &                     MPI_DOUBLE_PRECISION, ID_SCHUR, TAG_SCHUR,
     &                     id%COMM, STATUS, IERR )
            END IF
          ENDDO
        END IF
      ELSE
C       CASE 2 (Fwd in facto): Schur is not contiguous on ID_SCHUR,
C       process it row by row.
C
C       2.1: We first centralize Schur complement into id%SCHUR
        ISCHUR_SRC = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &               +4+id%KEEP(IXSZ)))
        ISCHUR_DEST= 1_8
        DO I=1, SIZE_SCHUR
          ROW_LENGTH = SIZE_SCHUR
          IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
            CALL dcopy(ROW_LENGTH, id%S(ISCHUR_SRC), 1,
     &                 id%SCHUR(ISCHUR_DEST),1)
          ELSE
            IF (id%MYID.EQ.ID_SCHUR) THEN
C             Send
              CALL MPI_SEND( id%S(ISCHUR_SRC), ROW_LENGTH,
     &        MPI_DOUBLE_PRECISION,
     &        MASTER, TAG_SCHUR,
     &        id%COMM, IERR )
            ELSE
C             Recv
              CALL MPI_RECV( id%SCHUR(ISCHUR_DEST),
     &                   ROW_LENGTH,
     &                   MPI_DOUBLE_PRECISION, ID_SCHUR, TAG_SCHUR,
     &                   id%COMM, STATUS, IERR )
            ENDIF
          ENDIF
          ISCHUR_SRC = ISCHUR_SRC+int(LD_SCHUR,8)
          ISCHUR_DEST= ISCHUR_DEST+int(SIZE_SCHUR,8)
        ENDDO
C       2.2: Get REDRHS on host
C       2.2.1: Symmetric => REDRHS is available in last KEEP(253)
C              rows of Schur structure on ID_SCHUR
C       2.2.2: Unsymmetric => REDRHS corresponds to last KEEP(253)
C              columns. However it must be transposed.
        IF (id%KEEP(221).EQ.1) THEN ! Implies Fwd in facto
          ISCHUR_SYM = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8) *
     &                    int(LD_SCHUR,8)
          ISCHUR_UNS =
     &                 id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8)
          ISCHUR_DEST = 1_8
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR .EQ. MASTER) THEN ! necessarily = id%MYID
              IF (id%KEEP(50) .EQ. 0) THEN
                CALL dcopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ELSE
                CALL dcopy(SIZE_SCHUR, id%S(ISCHUR_SYM), 1,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ENDIF
            ELSE
              IF (id%MYID .NE. MASTER) THEN
                IF (id%KEEP(50) .EQ. 0) THEN
C                 Use id%S(ISCHUR_SYM) as temporary contig. workspace
C                 of size SIZE_SCHUR. 
                  CALL dcopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &            id%S(ISCHUR_SYM), 1)
                ENDIF
                CALL MPI_SEND(id%S(ISCHUR_SYM), SIZE_SCHUR,
     &          MPI_DOUBLE_PRECISION, MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
              ELSE
                CALL MPI_RECV(id%REDRHS(ISCHUR_DEST),
     &          SIZE_SCHUR, MPI_DOUBLE_PRECISION, ID_SCHUR, TAG_SCHUR,
     &          id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
            IF (id%KEEP(50).EQ.0) THEN
              ISCHUR_UNS = ISCHUR_UNS + int(LD_SCHUR,8)
            ELSE
              ISCHUR_SYM = ISCHUR_SYM + int(LD_SCHUR,8)
            ENDIF
            ISCHUR_DEST = ISCHUR_DEST + int(id%LREDRHS,8)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_EXTRACT_SCHUR_REDRHS
