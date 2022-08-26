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
      MODULE DMUMPS_LOAD
      implicit none
      PUBLIC :: DMUMPS_LOAD_SET_INICOST, DMUMPS_LOAD_INIT,
     &  DMUMPS_LOAD_SET_SLAVES, DMUMPS_LOAD_UPDATE,
     &  DMUMPS_LOAD_END, DMUMPS_LOAD_PROCESS_MESSAGE, 
     &  DMUMPS_LOAD_LESS, DMUMPS_LOAD_LESS_CAND,
     &  DMUMPS_LOAD_SET_SLAVES_CAND, DMUMPS_LOAD_MASTER_2_ALL,
     &  DMUMPS_LOAD_RECV_MSGS, DMUMPS_LOAD_MEM_UPDATE,
     &  DMUMPS_LOAD_SET_PARTITION,
     &  DMUMPS_SPLIT_PREP_PARTITION, DMUMPS_SPLIT_POST_PARTITION,
     &  DMUMPS_SPLIT_PROPAGATE_PARTI, DMUMPS_LOAD_POOL_UPD_NEW_POOL,
     &  DMUMPS_LOAD_SBTR_UPD_NEW_POOL, DMUMPS_LOAD_POOL_CHECK_MEM,
     &  DMUMPS_LOAD_SET_SBTR_MEM,
     &  DMUMPS_REMOVE_NODE, DMUMPS_UPPER_PREDICT
     &  ,DMUMPS_LOAD_SEND_MD_INFO,
     &  DMUMPS_LOAD_CLEAN_MEMINFO_POOL, DMUMPS_LOAD_COMP_MAXMEM_POOL,
     &  DMUMPS_LOAD_CHK_MEMCST_POOL, DMUMPS_CHECK_SBTR_COST,
     &  DMUMPS_FIND_BEST_NODE_FOR_MEM,
     &  DMUMPS_LOAD_INIT_SBTR_STRUCT
      DOUBLE PRECISION, DIMENSION(:),
     &       ALLOCATABLE, SAVE, PRIVATE :: LOAD_FLOPS
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PRIVATE :: BUF_LOAD_RECV
      INTEGER, SAVE, PRIVATE :: LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES
      INTEGER, SAVE, PRIVATE :: K50, K69, K35
      INTEGER(8), SAVE, PRIVATE :: MAX_SURF_MASTER
      LOGICAL, SAVE, PRIVATE :: BDC_MEM, BDC_POOL, BDC_SBTR, 
     &     BDC_POOL_MNG,
     &     BDC_M2_MEM,BDC_M2_FLOPS,BDC_MD,REMOVE_NODE_FLAG,
     &     REMOVE_NODE_FLAG_MEM
      DOUBLE PRECISION, SAVE, PRIVATE :: REMOVE_NODE_COST,
     &     REMOVE_NODE_COST_MEM
      INTEGER, SAVE, PRIVATE :: SBTR_WHICH_M
      DOUBLE PRECISION, DIMENSION(:),
     &       ALLOCATABLE, TARGET, SAVE, PRIVATE :: WLOAD
      DOUBLE PRECISION, SAVE, PRIVATE :: DELTA_LOAD, DELTA_MEM
      LOGICAL, SAVE, PRIVATE :: IS_MUMPS_LOAD_ENABLED
      PUBLIC:: MUMPS_LOAD_ENABLE, MUMPS_LOAD_DISABLE
      INTEGER(8), SAVE, PRIVATE :: CHECK_MEM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, TARGET, PRIVATE :: 
     &          IDWLOAD
      DOUBLE PRECISION, SAVE, PRIVATE :: COST_SUBTREE
      DOUBLE PRECISION, SAVE, PRIVATE :: ALPHA
      DOUBLE PRECISION, SAVE, PRIVATE :: BETA
      INTEGER, SAVE, PRIVATE :: MYID, NPROCS, COMM_LD
      INTEGER, SAVE, PRIVATE :: COMM_NODES
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE, 
     &           PRIVATE :: POOL_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE, 
     &           SAVE :: SBTR_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, 
     &           PRIVATE, SAVE :: SBTR_CUR
      INTEGER, DIMENSION(:), ALLOCATABLE, 
     &           PRIVATE, SAVE :: NB_SON
      DOUBLE PRECISION, 
     &           PRIVATE, SAVE :: SBTR_CUR_LOCAL
      DOUBLE PRECISION, 
     &           PRIVATE, SAVE :: PEAK_SBTR_CUR_LOCAL
      DOUBLE PRECISION, 
     &           PRIVATE, SAVE :: MAX_PEAK_STK
      DOUBLE PRECISION, SAVE, 
     &           PRIVATE :: POOL_LAST_COST_SENT
      DOUBLE PRECISION, SAVE, 
     &           PRIVATE :: MIN_DIFF
      INTEGER, SAVE :: POS_ID,POS_MEM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: CB_COST_ID
      INTEGER(8), DIMENSION(:), ALLOCATABLE, SAVE
     &           :: CB_COST_MEM
      PUBLIC :: CB_COST_ID, CB_COST_MEM,POS_MEM,POS_ID
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: LU_USAGE
      INTEGER(8), DIMENSION(:), ALLOCATABLE, SAVE,
     &        PRIVATE::MD_MEM, TAB_MAXS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE ::MEM_SUBTREE
      INTEGER  :: NB_SUBTREES,NIV1_FLAG
      INTEGER, PRIVATE  :: INDICE_SBTR,INDICE_SBTR_ARRAY
      INTEGER :: POOL_NIV2_SIZE
      INTEGER,SAVE :: INSIDE_SUBTREE
      PUBLIC :: NB_SUBTREES,MEM_SUBTREE,INSIDE_SUBTREE,NIV1_FLAG
      DOUBLE PRECISION, SAVE, PRIVATE :: DM_SUMLU,
     &                   DM_THRES_MEM
      DOUBLE PRECISION, DIMENSION(:),
     &   ALLOCATABLE, SAVE , PRIVATE:: DM_MEM
      INTEGER, SAVE, PRIVATE :: POOL_SIZE,ID_MAX_M2
      DOUBLE PRECISION, SAVE, PRIVATE :: MAX_M2,TMP_M2
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE, PRIVATE:: POOL_NIV2
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,SAVE, 
     &      PRIVATE :: POOL_NIV2_COST, NIV2
      DOUBLE PRECISION, SAVE, PRIVATE  ::      CHK_LD
      INTEGER, DIMENSION(:),POINTER, SAVE, PRIVATE  :: 
     &         PROCNODE_LOAD, STEP_TO_NIV2_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE, PRIVATE  :: KEEP_LOAD
      INTEGER, SAVE, PRIVATE :: N_LOAD
      INTEGER(8), DIMENSION(:), POINTER, SAVE, PRIVATE:: KEEP8_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE :: 
     &         FILS_LOAD, STEP_LOAD,
     &         FRERE_LOAD, ND_LOAD,
     &         NE_LOAD,DAD_LOAD
      INTEGER, DIMENSION(:,:),POINTER, SAVE, PRIVATE :: CAND_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE, 
     &         PRIVATE :: MY_FIRST_LEAF,MY_NB_LEAF, MY_ROOT_SBTR
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE, 
     &         PRIVATE ::SBTR_FIRST_POS_IN_POOL
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,SAVE, 
     &         PRIVATE ::SBTR_PEAK_ARRAY,
     &     SBTR_CUR_ARRAY
      DOUBLE PRECISION,DIMENSION(:),POINTER, SAVE :: COST_TRAV
      INTEGER, DIMENSION(:),POINTER, SAVE :: DEPTH_FIRST_LOAD,
     &     DEPTH_FIRST_SEQ_LOAD,SBTR_ID_LOAD
      PUBLIC :: DEPTH_FIRST_LOAD,COST_TRAV, FILS_LOAD,STEP_LOAD,
     &     FRERE_LOAD, ND_LOAD,NE_LOAD,DAD_LOAD,
     &     DEPTH_FIRST_SEQ_LOAD,SBTR_ID_LOAD
      INTEGER, SAVE     :: ROOT_CURRENT_SUBTREE,CURRENT_BEST,
     &     SECOND_CURRENT_BEST
      PUBLIC :: ROOT_CURRENT_SUBTREE,CURRENT_BEST,
     &     SECOND_CURRENT_BEST
      CONTAINS
      SUBROUTINE MUMPS_LOAD_ENABLE()
      IMPLICIT NONE
      IS_MUMPS_LOAD_ENABLED = .TRUE.
      RETURN
      END SUBROUTINE MUMPS_LOAD_ENABLE
      SUBROUTINE MUMPS_LOAD_DISABLE()
      IMPLICIT NONE
      IS_MUMPS_LOAD_ENABLED = .FALSE.
      RETURN
      END SUBROUTINE MUMPS_LOAD_DISABLE
      SUBROUTINE DMUMPS_LOAD_SET_INICOST( COST_SUBTREE_ARG, K64, DK15,
     &     K375, MAXS )
      IMPLICIT NONE
      DOUBLE PRECISION COST_SUBTREE_ARG
      INTEGER, INTENT(IN) :: K64, K375
      DOUBLE PRECISION, INTENT(IN) :: DK15
      INTEGER(8)::MAXS
      DOUBLE PRECISION :: T64, T66
      LOGICAL :: AVOID_LOAD_MESSAGES
      T64 = max ( dble(K64), dble(1) )
      T64 = min ( T64, dble(1000)  )
      T66 = max (dble(DK15), dble(100))
      MIN_DIFF     =  ( T64 / dble(1000)  )* 
     &                  T66 * dble(1000000)
      DM_THRES_MEM = dble(MAXS/300_8)
      COST_SUBTREE = COST_SUBTREE_ARG
      AVOID_LOAD_MESSAGES = .FALSE.
      IF (K375.EQ.1) THEN
        AVOID_LOAD_MESSAGES = .TRUE.
      ENDIF
      IF (AVOID_LOAD_MESSAGES) THEN
        MIN_DIFF = MIN_DIFF * 1000.D0
        DM_THRES_MEM = DM_THRES_MEM * 1000_8
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_SET_INICOST
      SUBROUTINE DMUMPS_SPLIT_PREP_PARTITION ( 
     &      INODE, STEP, N, SLAVEF,
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      CAND, ICNTL, COPY_CAND,
     &      NBSPLIT, NUMORG_SPLIT, SLAVES_LIST,
     &      SIZE_SLAVES_LIST 
     &                                    )
      IMPLICIT NONE
       INTEGER, intent(in) :: INODE, N, SIZE_SLAVES_LIST, SLAVEF, 
     &                        KEEP(500)
       INTEGER, intent(in) :: STEP(N), DAD (KEEP(28)), ICNTL(60),
     &                        PROCNODE_STEPS(KEEP(28)), CAND(SLAVEF+1),
     &                        FILS(N)
       INTEGER, intent(out)   :: NBSPLIT, NUMORG_SPLIT
       INTEGER, intent(inout) :: SLAVES_LIST(SIZE_SLAVES_LIST),
     &                           COPY_CAND(SLAVEF+1)
       INTEGER :: IN, LP, II
       INTEGER  MUMPS_TYPESPLIT
       EXTERNAL MUMPS_TYPESPLIT
       LP = ICNTL(1) 
       IN = INODE
       NBSPLIT = 0
       NUMORG_SPLIT = 0    
       DO WHILE 
     &      (
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),KEEP(199))
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),KEEP(199))
     &           .EQ.6  
     &        ) 
     &      )  
           NBSPLIT = NBSPLIT + 1
           IN = DAD(STEP(IN))
           II = IN
           DO WHILE (II.GT.0)
             NUMORG_SPLIT = NUMORG_SPLIT + 1
             II = FILS(II)
           ENDDO
       END DO
      SLAVES_LIST(1:NBSPLIT) = CAND(1:NBSPLIT)
      COPY_CAND(1:SIZE_SLAVES_LIST-NBSPLIT) = 
     &                   CAND(1+NBSPLIT:SIZE_SLAVES_LIST)
      COPY_CAND(SIZE_SLAVES_LIST-NBSPLIT+1:SLAVEF) = -1
      COPY_CAND(SLAVEF+1) = SIZE_SLAVES_LIST-NBSPLIT
      RETURN
      END SUBROUTINE DMUMPS_SPLIT_PREP_PARTITION
      SUBROUTINE DMUMPS_SPLIT_POST_PARTITION ( 
     &      INODE, STEP, N, SLAVEF, NBSPLIT, NCB, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS, ICNTL,
     &      TAB_POS, NSLAVES_NODE
     &                                    )
      IMPLICIT NONE
       INTEGER, intent(in) :: INODE, N, SLAVEF, NCB, 
     &                        KEEP(500), NBSPLIT
       INTEGER, intent(in) :: STEP(N), DAD (KEEP(28)), ICNTL(60),
     &                        PROCNODE_STEPS(KEEP(28)), 
     &                        FILS(N)
       INTEGER, intent(inout) :: TAB_POS ( SLAVEF+2 ), NSLAVES_NODE
       INTEGER :: IN, LP, II, NUMORG, NBSPLIT_LOC, I
       INTEGER  MUMPS_TYPESPLIT
       EXTERNAL MUMPS_TYPESPLIT
       DO I= NSLAVES_NODE+1, 1, -1
          TAB_POS(I+NBSPLIT) = TAB_POS(I) 
       END DO
       LP = ICNTL(1) 
       IN = INODE
       NBSPLIT_LOC = 0
       NUMORG = 0
       TAB_POS(1) = 1
       DO WHILE 
     &      (
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),KEEP(199))
     &           .EQ.5 
     &        ) 
     &        .OR.
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IN)))),KEEP(199))
     &           .EQ.6  
     &        ) 
     &      )  
           NBSPLIT_LOC = NBSPLIT_LOC + 1
           IN = DAD(STEP(IN))
           II = IN
           DO WHILE (II.GT.0)
             NUMORG = NUMORG + 1
             II = FILS(II)
           ENDDO
           TAB_POS(NBSPLIT_LOC+1) = NUMORG + 1
       END DO
       DO I = NBSPLIT+2, NBSPLIT+NSLAVES_NODE+1
         TAB_POS(I) = TAB_POS(I) + NUMORG
       ENDDO
      NSLAVES_NODE = NSLAVES_NODE + NBSPLIT
      TAB_POS (NSLAVES_NODE+2:SLAVEF+1) = -9999
      TAB_POS ( SLAVEF+2 ) =  NSLAVES_NODE
      RETURN
      END SUBROUTINE DMUMPS_SPLIT_POST_PARTITION
      SUBROUTINE DMUMPS_SPLIT_PROPAGATE_PARTI (
     &      INODE, TYPESPLIT, IFSON, 
     &      CAND, SIZE_CAND,
     &      SON_SLAVE_LIST, NSLSON,
     &      STEP, N, SLAVEF, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS, ICNTL,
     &      ISTEP_TO_INIV2, INIV2,
     &      TAB_POS_IN_PERE, NSLAVES_NODE,
     &      SLAVES_LIST, SIZE_SLAVES_LIST
     &                                    )
      IMPLICIT NONE
       INTEGER, intent(in) :: INODE, TYPESPLIT, IFSON, N, SLAVEF, 
     &                        KEEP(500), 
     &                        NSLSON, SIZE_SLAVES_LIST, SIZE_CAND
       INTEGER, intent(in) :: STEP(N), DAD (KEEP(28)), ICNTL(60),
     &                        PROCNODE_STEPS(KEEP(28)), 
     &                        FILS(N), INIV2,
     &                        SON_SLAVE_LIST (NSLSON),
     &                        ISTEP_TO_INIV2(KEEP(71)),
     &                        CAND(SIZE_CAND)
       INTEGER, intent(out)   ::  NSLAVES_NODE
       INTEGER, intent(inout) ::
     &                   TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
       INTEGER, intent(out)   :: SLAVES_LIST (SIZE_SLAVES_LIST)
       INTEGER :: IN, LP, I, NSLAVES_SONS,
     &            INIV2_FILS, ISHIFT
       LP = ICNTL(1) 
       IN = INODE
      INIV2_FILS = ISTEP_TO_INIV2( STEP( IFSON ))
      NSLAVES_SONS = TAB_POS_IN_PERE (SLAVEF+2, INIV2_FILS)
      TAB_POS_IN_PERE (1,INIV2) = 1
      ISHIFT  = TAB_POS_IN_PERE (2, INIV2_FILS) -1 
      DO I = 2, NSLAVES_SONS
         TAB_POS_IN_PERE (I,INIV2) = 
     &            TAB_POS_IN_PERE (I+1,INIV2_FILS) - ISHIFT
         SLAVES_LIST(I-1) =  SON_SLAVE_LIST (I)
      END DO
      TAB_POS_IN_PERE(NSLAVES_SONS+1:SLAVEF+1,INIV2) = -9999
      NSLAVES_NODE = NSLAVES_SONS - 1
      TAB_POS_IN_PERE (SLAVEF+2, INIV2) = NSLAVES_NODE
      RETURN
      END SUBROUTINE DMUMPS_SPLIT_PROPAGATE_PARTI 
      SUBROUTINE DMUMPS_LOAD_SET_PARTITION(
     &  NCBSON_MAX, SLAVEF,
     &  KEEP,KEEP8,ICNTL,
     &  CAND_OF_NODE,
     &  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,INODE
     &)
       IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: ICNTL(60)
      INTEGER, intent(in) :: SLAVEF, NFRONT
      INTEGER, intent (inout) ::NCB   
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(in) :: NCBSON_MAX
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER i
      INTEGER LP,MP
      INTEGER(8) DUMMY1
      INTEGER DUMMY2
      INTEGER TMP_ARRAY(2)
      LP=ICNTL(4)
      MP=ICNTL(2)
      IF ( KEEP(48) == 0 .OR. KEEP(48) .EQ. 3 ) THEN
         CALL DMUMPS_LOAD_PARTI_REGULAR(
     &        SLAVEF,
     &        KEEP,KEEP8,
     &        CAND_OF_NODE,
     &        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      ELSE IF ( KEEP(48) == 4 ) THEN
         CALL DMUMPS_SET_PARTI_ACTV_MEM(
     &        SLAVEF,
     &        KEEP,KEEP8,
     &        CAND_OF_NODE,
     &        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID)
         DO i=1,NSLAVES_NODE
            IF(TAB_POS(i+1)-TAB_POS(i).LE.0)THEN
               WRITE(*,*)'probleme de partition dans 
     &DMUMPS_LOAD_SET_PARTI_ACTV_MEM'
               CALL MUMPS_ABORT()
            ENDIF
         ENDDO
      ELSE IF ( KEEP(48) == 5 ) THEN
         IF (KEEP(375).EQ.1) THEN 
           GOTO 458
         ENDIF
         CALL DMUMPS_SET_PARTI_FLOP_IRR(
     &        NCBSON_MAX,
     &        SLAVEF,
     &        KEEP,KEEP8,
     &        CAND_OF_NODE,
     &        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID,INODE,
     &        MP,LP)
         DO i=1,NSLAVES_NODE
            IF(TAB_POS(i+1)-TAB_POS(i).LE.0)THEN
               WRITE(*,*)'problem with partition in
     &DMUMPS_SET_PARTI_FLOP_IRR'
               CALL MUMPS_ABORT()
            ENDIF
         ENDDO
         GOTO 457
 458     CONTINUE
         IF ( KEEP(375).EQ.1 )THEN 
           TMP_ARRAY(1)=0
           TMP_ARRAY(2)=0
         ENDIF
         CALL DMUMPS_SET_PARTI_REGULAR(
     &        SLAVEF,
     &        KEEP,KEEP8,
     &        CAND_OF_NODE,
     &        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID,INODE,
     &        TAB_MAXS,TMP_ARRAY,DUMMY1,DUMMY2
     &        )
      ELSE
        WRITE(*,*) "Strategy 6 not implemented"
        CALL MUMPS_ABORT()
      ENDIF
 457  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_LOAD_SET_PARTITION
      SUBROUTINE DMUMPS_LOAD_PARTI_REGULAR(
     &  SLAVEF,
     &  KEEP,KEEP8,
     &  CAND_OF_NODE,
     &  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER ITEMP, NMB_OF_CAND, NSLAVES_LESS
      DOUBLE PRECISION MSG_SIZE
      LOGICAL FORCE_CAND
      INTEGER  MUMPS_REG_GET_NSLAVES
      EXTERNAL MUMPS_REG_GET_NSLAVES
      IF ( KEEP(48) == 0 .AND. KEEP(50) .NE. 0) THEN
      write(*,*) "Internal error 2 in DMUMPS_LOAD_PARTI_REGULAR."
      CALL MUMPS_ABORT()
      END IF
      IF ( KEEP(48) == 3 .AND. KEEP(50) .EQ. 0) THEN
      write(*,*) "Internal error 3 in DMUMPS_LOAD_PARTI_REGULAR."
      CALL MUMPS_ABORT()
      END IF
      MSG_SIZE = dble( NFRONT - NCB ) * dble(NCB)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
      END IF
      IF (FORCE_CAND) THEN
        ITEMP=DMUMPS_LOAD_LESS_CAND
     &       (MEM_DISTRIB,
     &        CAND_OF_NODE,
     &
     &        KEEP(69), SLAVEF, MSG_SIZE,
     &        NMB_OF_CAND )
      ELSE
        ITEMP=DMUMPS_LOAD_LESS(KEEP(69),MEM_DISTRIB,MSG_SIZE)
        NMB_OF_CAND = SLAVEF - 1
      END IF
      NSLAVES_LESS = max(ITEMP,1)
      NSLAVES_NODE = MUMPS_REG_GET_NSLAVES(KEEP8(21), KEEP(48),
     &          KEEP(50),SLAVEF,
     &          NCB, NFRONT, NSLAVES_LESS, NMB_OF_CAND,
     &          KEEP(375), KEEP(119))
      CALL MUMPS_BLOC2_SETPARTITION(
     &            KEEP,KEEP8, SLAVEF,
     &            TAB_POS,
     &            NSLAVES_NODE, NFRONT, NCB
     &             )
      IF (FORCE_CAND) THEN
        CALL DMUMPS_LOAD_SET_SLAVES_CAND(MEM_DISTRIB(0),
     &       CAND_OF_NODE, SLAVEF, NSLAVES_NODE,
     &       SLAVES_LIST)
      ELSE
        CALL DMUMPS_LOAD_SET_SLAVES(MEM_DISTRIB(0),
     &       MSG_SIZE, SLAVES_LIST, NSLAVES_NODE)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_PARTI_REGULAR
      SUBROUTINE DMUMPS_LOAD_INIT( id, MEMORY_MD_ARG, MAXS )
      USE DMUMPS_BUF
      USE DMUMPS_STRUC_DEF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER(8), intent(in) :: MEMORY_MD_ARG
      INTEGER(8), intent(in) :: MAXS
      INTEGER K34_LOC
      INTEGER(8) :: I8SIZE
      INTEGER allocok, IERR, IERR_MPI, i, BUF_LOAD_SIZE
      DOUBLE PRECISION :: MAX_SBTR
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION MEMORY_SENT
      PARAMETER( ZERO=0.0d0 )
      DOUBLE PRECISION SIZE_DBLE(2) 
      INTEGER WHAT
      INTEGER(8) MEMORY_MD, LA
      CALL MUMPS_LOAD_ENABLE()
      STEP_TO_NIV2_LOAD=>id%ISTEP_TO_INIV2
      CAND_LOAD=>id%CANDIDATES
      ND_LOAD=>id%ND_STEPS
      KEEP_LOAD=>id%KEEP
      KEEP8_LOAD=>id%KEEP8
      FILS_LOAD=>id%FILS
      FRERE_LOAD=>id%FRERE_STEPS
      DAD_LOAD=>id%DAD_STEPS
      PROCNODE_LOAD=>id%PROCNODE_STEPS
      STEP_LOAD=>id%STEP
      NE_LOAD=>id%NE_STEPS
      N_LOAD=id%N
      ROOT_CURRENT_SUBTREE=-9999
      MEMORY_MD=MEMORY_MD_ARG
      LA=MAXS
      MAX_SURF_MASTER=id%MAX_SURF_MASTER+
     & (int(id%KEEP(12),8)*int(id%MAX_SURF_MASTER,8)/int(100,8))
      COMM_LD    = id%COMM_LOAD
      COMM_NODES = id%COMM_NODES
      MAX_PEAK_STK = 0.0D0
      K69  = id%KEEP(69)
      IF ( id%KEEP(47) .le. 0 .OR. id%KEEP(47) .gt. 4 ) THEN
        write(*,*) "Internal error 1 in DMUMPS_LOAD_INIT"
        CALL MUMPS_ABORT()
      END IF
      CHK_LD=dble(0)
      BDC_MEM      = ( id%KEEP(47) >= 2 )
      BDC_POOL     = ( id%KEEP(47) >= 3 )
      BDC_SBTR     = ( id%KEEP(47) >= 4 )
      BDC_M2_MEM   = ( ( id%KEEP(80) == 2 .OR. id%KEEP(80) == 3 )
     &             .AND. id%KEEP(47) == 4 )
      BDC_M2_FLOPS   = ( id%KEEP(80) == 1 
     &             .AND. id%KEEP(47) .GE. 1 )
      BDC_MD       = (id%KEEP(86)==1)
      SBTR_WHICH_M       = id%KEEP(90)
      REMOVE_NODE_FLAG=.FALSE.
      REMOVE_NODE_FLAG_MEM=.FALSE.
      REMOVE_NODE_COST_MEM=dble(0)
      REMOVE_NODE_COST=dble(0)
      IF (id%KEEP(80) .LT. 0 .OR. id%KEEP(80)>3) THEN
        WRITE(*,*) "Unimplemented KEEP(80) Strategy"
        CALL MUMPS_ABORT()
      ENDIF
      IF ((id%KEEP(80) == 2 .OR. id%KEEP(80)==3).AND. id%KEEP(47).NE.4)
     &  THEN
        WRITE(*,*) "Internal error 3 in DMUMPS_LOAD_INIT"
        CALL MUMPS_ABORT()
      END IF
      IF (id%KEEP(81) == 1 .AND. id%KEEP(47) < 2) THEN
        WRITE(*,*) "Internal error 2 in DMUMPS_LOAD_INIT"
        CALL MUMPS_ABORT()
      ENDIF
      BDC_POOL_MNG = ((id%KEEP(81) == 1).AND.(id%KEEP(47) >= 2))
      IF(id%KEEP(76).EQ.4)THEN
         DEPTH_FIRST_LOAD=>id%DEPTH_FIRST
      ENDIF
      IF(id%KEEP(76).EQ.5)THEN
         COST_TRAV=>id%COST_TRAV
      ENDIF
      IF(id%KEEP(76).EQ.6)THEN
         DEPTH_FIRST_LOAD=>id%DEPTH_FIRST
         DEPTH_FIRST_SEQ_LOAD=>id%DEPTH_FIRST_SEQ
         SBTR_ID_LOAD=>id%SBTR_ID
      ENDIF
      IF (BDC_M2_MEM.OR.BDC_M2_FLOPS) THEN
         POOL_NIV2_SIZE=max(1,min(id%NBSA+id%KEEP(262),id%NA(1)))
         ALLOCATE(NIV2(id%NSLAVES), NB_SON(id%KEEP(28)),
     &            POOL_NIV2(POOL_NIV2_SIZE),
     &        POOL_NIV2_COST(POOL_NIV2_SIZE),
     &            stat=allocok)
         DO i = 1, id%KEEP(28)
           NB_SON(i)=id%NE_STEPS(i)
         ENDDO
         NIV2=dble(0)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES + id%KEEP(28) + 200
           RETURN
         ENDIF
      ENDIF
      K50      = id%KEEP(50)
      CALL MPI_COMM_RANK( COMM_LD, MYID, IERR_MPI )
      NPROCS = id%NSLAVES
      DM_SUMLU=ZERO
      POOL_SIZE=0
      IF(BDC_MD)THEN
         IF ( allocated(MD_MEM) ) DEALLOCATE(MD_MEM)
         ALLOCATE( MD_MEM( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         IF ( allocated(TAB_MAXS) ) DEALLOCATE(TAB_MAXS)
         ALLOCATE( TAB_MAXS( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         TAB_MAXS=0_8
         IF ( allocated(LU_USAGE) ) DEALLOCATE(LU_USAGE)
         ALLOCATE( LU_USAGE( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         LU_USAGE=dble(0)
         MD_MEM=int(0,8)
      ENDIF
      IF((id%KEEP(81).EQ.2).OR.(id%KEEP(81).EQ.3))THEN
         ALLOCATE(CB_COST_MEM(2*2000*id%NSLAVES), 
     &            stat=allocok)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES 
           RETURN
         ENDIF
         CB_COST_MEM=int(0,8)
         ALLOCATE(CB_COST_ID(2000*3), 
     &            stat=allocok)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES 
           RETURN
         ENDIF
         CB_COST_ID=0
         POS_MEM=1
         POS_ID=1
      ENDIF
      ALLOCATE(FUTURE_NIV2(NPROCS), stat=allocok)
      IF (allocok > 0 ) THEN
         WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      ENDIF
      DO i = 1, NPROCS
        FUTURE_NIV2(i) = id%FUTURE_NIV2(i)
        IF(BDC_MD)THEN
           IF(FUTURE_NIV2(i).EQ.0)THEN
              MD_MEM(i-1)=999999999_8
           ENDIF
        ENDIF
      ENDDO
      DELTA_MEM=ZERO
      DELTA_LOAD=ZERO
      CHECK_MEM=0_8
      IF(BDC_SBTR.OR.BDC_POOL_MNG)THEN
         NB_SUBTREES=id%NBSA_LOCAL
         IF (allocated(MEM_SUBTREE)) DEALLOCATE(MEM_SUBTREE)
         ALLOCATE(MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
         IF (allocok > 0 ) THEN
            WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
            id%INFO(1) = -13
            id%INFO(2) = id%NBSA_LOCAL
            RETURN
         ENDIF
         DO i=1,id%NBSA_LOCAL
            MEM_SUBTREE(i)=id%MEM_SUBTREE(i)
         ENDDO
         MY_FIRST_LEAF=>id%MY_FIRST_LEAF
         MY_NB_LEAF=>id%MY_NB_LEAF
         MY_ROOT_SBTR=>id%MY_ROOT_SBTR
         IF (allocated(SBTR_FIRST_POS_IN_POOL))
     &        DEALLOCATE(SBTR_FIRST_POS_IN_POOL)
         ALLOCATE(SBTR_FIRST_POS_IN_POOL(id%NBSA_LOCAL),stat=allocok)
         IF (allocok > 0 ) THEN
            WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
            id%INFO(1) = -13
            id%INFO(2) = id%NBSA_LOCAL
            RETURN
         ENDIF
         INSIDE_SUBTREE=0
         PEAK_SBTR_CUR_LOCAL = dble(0)
         SBTR_CUR_LOCAL      = dble(0)
         IF (allocated(SBTR_PEAK_ARRAY)) DEALLOCATE(SBTR_PEAK_ARRAY)
         ALLOCATE(SBTR_PEAK_ARRAY(id%NBSA_LOCAL),stat=allocok)
         IF (allocok > 0 ) THEN
            WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
            id%INFO(1) = -13
            id%INFO(2) = id%NBSA_LOCAL
            RETURN
         ENDIF
         SBTR_PEAK_ARRAY=dble(0)
         IF (allocated(SBTR_CUR_ARRAY)) DEALLOCATE(SBTR_CUR_ARRAY)
         ALLOCATE(SBTR_CUR_ARRAY(id%NBSA_LOCAL),stat=allocok)
         IF (allocok > 0 ) THEN
            WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
            id%INFO(1) = -13
            id%INFO(2) = id%NBSA_LOCAL
            RETURN
         ENDIF
         SBTR_CUR_ARRAY=dble(0)
         INDICE_SBTR_ARRAY=1
         NIV1_FLAG=0
         INDICE_SBTR=1
      ENDIF
      IF ( allocated(LOAD_FLOPS) ) DEALLOCATE( LOAD_FLOPS )
      ALLOCATE( LOAD_FLOPS( 0: NPROCS - 1 ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( allocated(WLOAD) ) DEALLOCATE( WLOAD )
      ALLOCATE( WLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( allocated(IDWLOAD) ) DEALLOCATE( IDWLOAD )
      ALLOCATE( IDWLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( BDC_MEM ) THEN
        IF ( allocated(DM_MEM) ) DEALLOCATE( DM_MEM )
        ALLOCATE( DM_MEM( 0:NPROCS-1 ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
      END IF
      IF ( BDC_POOL ) THEN
        IF ( allocated(POOL_MEM) ) DEALLOCATE(POOL_MEM)
        ALLOCATE( POOL_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        POOL_MEM = dble(0) 
        POOL_LAST_COST_SENT = dble(0)
      END IF
      IF ( BDC_SBTR ) THEN
        IF ( allocated(SBTR_MEM) ) DEALLOCATE(SBTR_MEM)
        ALLOCATE( SBTR_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        IF ( allocated(SBTR_CUR) ) DEALLOCATE(SBTR_CUR)
        ALLOCATE( SBTR_CUR(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        SBTR_CUR = dble(0) 
        SBTR_MEM = dble(0) 
      END IF
      K34_LOC=id%KEEP(34)
      CALL MUMPS_SIZE_C(SIZE_DBLE(1),SIZE_DBLE(2),I8SIZE)
      K35  = int(I8SIZE)
      BUF_LOAD_SIZE = K34_LOC * 2 * ( NPROCS - 1 ) +
     &                NPROCS * ( K35 + K34_LOC )
      IF (BDC_MEM) THEN
        BUF_LOAD_SIZE = BUF_LOAD_SIZE + NPROCS * K35
      END IF
      IF (BDC_SBTR)THEN
        BUF_LOAD_SIZE = BUF_LOAD_SIZE + NPROCS * K35
      ENDIF
      LBUF_LOAD_RECV = (BUF_LOAD_SIZE+K34_LOC)/K34_LOC
      LBUF_LOAD_RECV_BYTES = LBUF_LOAD_RECV * K34_LOC
      IF ( allocated(BUF_LOAD_RECV) ) DEALLOCATE(BUF_LOAD_RECV)
      ALLOCATE( BUF_LOAD_RECV( LBUF_LOAD_RECV), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in DMUMPS_LOAD_INIT'
        id%INFO(1) = -13
        id%INFO(2) = LBUF_LOAD_RECV
        RETURN
      ENDIF 
      BUF_LOAD_SIZE = BUF_LOAD_SIZE * 20
      CALL DMUMPS_BUF_ALLOC_LOAD_BUFFER( BUF_LOAD_SIZE, IERR )
      IF ( IERR .LT. 0 ) THEN
         id%INFO(1) = -13
         id%INFO(2) = BUF_LOAD_SIZE
         RETURN
      END IF
      DO i = 0, NPROCS - 1
         LOAD_FLOPS( i ) = ZERO
      END DO
      IF ( BDC_MEM ) THEN
        DO i = 0, NPROCS - 1
          DM_MEM( i )=ZERO
        END DO
      ENDIF
      CALL DMUMPS_INIT_ALPHA_BETA(id%KEEP(69))
      IF(BDC_MD)THEN
         MAX_SBTR=0.0D0
         IF(BDC_SBTR)THEN
            DO i=1,id%NBSA_LOCAL
               MAX_SBTR=max(id%MEM_SUBTREE(i),MAX_SBTR)
            ENDDO
         ENDIF
         MD_MEM(MYID)=MEMORY_MD
         WHAT=8
         CALL DMUMPS_BUF_BROADCAST( WHAT,
     &        COMM_LD, NPROCS,
     &        FUTURE_NIV2,
     &        dble(MEMORY_MD),dble(0) ,MYID, id%KEEP, IERR  )
         WHAT=9
         MEMORY_SENT = dble(LA-MAX_SURF_MASTER)-MAX_SBTR
     &      - max( dble(LA) * dble(3) / dble(100),
     &      dble(2) *
     &      dble(max(id%KEEP(5),id%KEEP(6))) * dble(id%KEEP(127)))
         IF (id%KEEP(12) > 25) THEN
           MEMORY_SENT = MEMORY_SENT -
     &                   dble(id%KEEP(12))*0.2d0*dble(LA)/100.0d0
         ENDIF
         IF (id%KEEP(375).EQ.1) THEN  
           MEMORY_SENT=dble(LA)
         ENDIF
         TAB_MAXS(MYID)=int(MEMORY_SENT,8)
         CALL DMUMPS_BUF_BROADCAST( WHAT,
     &        COMM_LD, NPROCS,
     &        FUTURE_NIV2,
     &        MEMORY_SENT,
     &        dble(0),MYID, id%KEEP, IERR  )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_INIT
      SUBROUTINE DMUMPS_LOAD_UPDATE( CHECK_FLOPS,PROCESS_BANDE,
     &     INC_LOAD, KEEP,KEEP8 )
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      DOUBLE PRECISION INC_LOAD
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      LOGICAL PROCESS_BANDE
      LOGICAL :: EXIT_FLAG
      INTEGER CHECK_FLOPS
      INTEGER IERR
      DOUBLE PRECISION ZERO, SEND_MEM, SEND_LOAD,SBTR_TMP
      PARAMETER( ZERO=0.0d0 )
      INTRINSIC max, abs
      IF (.NOT. IS_MUMPS_LOAD_ENABLED) RETURN
      IF (INC_LOAD == 0.0D0) THEN
         IF(REMOVE_NODE_FLAG)THEN
            REMOVE_NODE_FLAG=.FALSE.
         ENDIF
         RETURN
      ENDIF
      IF((CHECK_FLOPS.NE.0).AND.
     &     (CHECK_FLOPS.NE.1).AND.(CHECK_FLOPS.NE.2))THEN
         WRITE(*,*)MYID,': Bad value for CHECK_FLOPS'
         CALL MUMPS_ABORT()
      ENDIF
      IF(CHECK_FLOPS.EQ.1)THEN
         CHK_LD=CHK_LD+INC_LOAD
      ELSE 
         IF(CHECK_FLOPS.EQ.2)THEN
            RETURN
         ENDIF
      ENDIF
      IF ( PROCESS_BANDE ) THEN
         RETURN                 
      ENDIF
      LOAD_FLOPS( MYID ) = max( LOAD_FLOPS( MYID ) + INC_LOAD, ZERO)
      IF(BDC_M2_FLOPS.AND.REMOVE_NODE_FLAG)THEN
         IF(INC_LOAD.NE.REMOVE_NODE_COST)THEN
            IF(INC_LOAD.GT.REMOVE_NODE_COST)THEN
               DELTA_LOAD = DELTA_LOAD +
     &              (INC_LOAD-REMOVE_NODE_COST)
               GOTO 888
            ELSE
               DELTA_LOAD = DELTA_LOAD -
     &              (REMOVE_NODE_COST-INC_LOAD)
               GOTO 888
            ENDIF
         ENDIF
         GOTO 333
      ENDIF
      DELTA_LOAD = DELTA_LOAD + INC_LOAD
 888  CONTINUE
      IF ( DELTA_LOAD > MIN_DIFF .OR. DELTA_LOAD < -MIN_DIFF) THEN
         SEND_LOAD = DELTA_LOAD
         IF (BDC_MEM) THEN
           SEND_MEM = DELTA_MEM
         ELSE
           SEND_MEM = ZERO
         END IF
         IF(BDC_SBTR)THEN
           SBTR_TMP=SBTR_CUR(MYID)
         ELSE
           SBTR_TMP=dble(0)
         ENDIF
 111     CONTINUE
         CALL DMUMPS_BUF_SEND_UPDATE_LOAD( BDC_SBTR,BDC_MEM,
     &        BDC_MD,COMM_LD, NPROCS,
     &        SEND_LOAD,
     &        SEND_MEM,SBTR_TMP,
     &        DM_SUMLU,
     &        FUTURE_NIV2,
     &        MYID, KEEP, IERR )
         IF ( IERR == -1 )THEN
             CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
             CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
             IF (EXIT_FLAG) THEN
               GOTO 333
             ELSE
               GOTO 111
             ENDIF
         ELSE IF ( IERR .NE.0 ) THEN
             WRITE(*,*) "Internal Error in DMUMPS_LOAD_UPDATE",IERR
             CALL MUMPS_ABORT()
         ENDIF
         DELTA_LOAD = ZERO
         IF (BDC_MEM) DELTA_MEM  = ZERO
      ENDIF
 333  CONTINUE
      IF(REMOVE_NODE_FLAG)THEN
         REMOVE_NODE_FLAG=.FALSE.
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_UPDATE
      SUBROUTINE DMUMPS_LOAD_MEM_UPDATE( SSARBR,
     &           PROCESS_BANDE_ARG, MEM_VALUE, NEW_LU, INC_MEM_ARG,
     &           KEEP,KEEP8,LRLUS)
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER(8), INTENT(IN) :: MEM_VALUE, INC_MEM_ARG, NEW_LU,LRLUS
      LOGICAL, INTENT(IN) :: PROCESS_BANDE_ARG, SSARBR
      INTEGER IERR, KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION ZERO, SEND_MEM, SBTR_TMP
      PARAMETER( ZERO=0.0d0 )
      INTRINSIC max, abs
      INTEGER(8) :: INC_MEM
      LOGICAL PROCESS_BANDE
      LOGICAL :: EXIT_FLAG
      IF (.NOT. IS_MUMPS_LOAD_ENABLED) RETURN
      PROCESS_BANDE=PROCESS_BANDE_ARG
      INC_MEM = INC_MEM_ARG
      IF ( PROCESS_BANDE .AND. NEW_LU .NE. 0_8) THEN
        WRITE(*,*) " Internal Error in DMUMPS_LOAD_MEM_UPDATE."
        WRITE(*,*) " NEW_LU must be zero if called from PROCESS_BANDE"
        CALL MUMPS_ABORT()
      ENDIF
      DM_SUMLU = DM_SUMLU + dble(NEW_LU)
      IF(KEEP_LOAD(201).EQ.0)THEN
         CHECK_MEM = CHECK_MEM + INC_MEM
      ELSE 
         CHECK_MEM = CHECK_MEM + INC_MEM - NEW_LU
      ENDIF
      IF ( MEM_VALUE .NE. CHECK_MEM ) THEN
         WRITE(*,*)MYID,
     &   ':Problem with increments in DMUMPS_LOAD_MEM_UPDATE',
     &   CHECK_MEM, MEM_VALUE, INC_MEM,NEW_LU
         CALL MUMPS_ABORT()
      ENDIF
      IF (PROCESS_BANDE) THEN
         RETURN
      ENDIF
      IF(BDC_POOL_MNG) THEN
         IF(SBTR_WHICH_M.EQ.0)THEN
            IF (SSARBR) SBTR_CUR_LOCAL = SBTR_CUR_LOCAL+
     &                                   dble(INC_MEM-NEW_LU)
         ELSE
            IF (SSARBR) SBTR_CUR_LOCAL = SBTR_CUR_LOCAL+
     &                                   dble(INC_MEM)
         ENDIF
      ENDIF
      IF ( .NOT. BDC_MEM ) THEN
         RETURN
      ENDIF
      IF (BDC_SBTR .AND. SSARBR) THEN
         IF((SBTR_WHICH_M.EQ.0).AND.(KEEP(201).NE.0))THEN
            SBTR_CUR(MYID) = SBTR_CUR(MYID)+dble(INC_MEM-NEW_LU)
         ELSE
            SBTR_CUR(MYID) = SBTR_CUR(MYID)+dble(INC_MEM)
         ENDIF
         SBTR_TMP = SBTR_CUR(MYID)
      ELSE
        SBTR_TMP=dble(0)
      ENDIF
      IF ( NEW_LU > 0_8 ) THEN
        INC_MEM = INC_MEM - NEW_LU
      ENDIF
      DM_MEM( MYID ) = DM_MEM(MYID) + dble(INC_MEM)
      MAX_PEAK_STK=max(MAX_PEAK_STK,DM_MEM(MYID))
      IF(BDC_M2_MEM.AND.REMOVE_NODE_FLAG_MEM)THEN
         IF(dble(INC_MEM).NE.REMOVE_NODE_COST_MEM)THEN
            IF(dble(INC_MEM).GT.REMOVE_NODE_COST_MEM)THEN
               DELTA_MEM = DELTA_MEM +
     &              (dble(INC_MEM)-REMOVE_NODE_COST_MEM)
               GOTO 888               
            ELSE
               DELTA_MEM = DELTA_MEM -
     &              (REMOVE_NODE_COST_MEM-dble(INC_MEM))
               GOTO 888
            ENDIF
         ENDIF
         GOTO 333
      ENDIF
      DELTA_MEM = DELTA_MEM + dble(INC_MEM)
 888  CONTINUE
      IF ((KEEP(48).NE.5).OR.
     &     ((KEEP(48).EQ.5).AND.(abs(DELTA_MEM)
     &      .GE.0.2d0*dble(LRLUS))))THEN
         IF ( abs(DELTA_MEM) > DM_THRES_MEM ) THEN
            SEND_MEM = DELTA_MEM
 111        CONTINUE
            CALL DMUMPS_BUF_SEND_UPDATE_LOAD( 
     &           BDC_SBTR,
     &           BDC_MEM,BDC_MD, COMM_LD,
     &           NPROCS,
     &           DELTA_LOAD,
     &           SEND_MEM,SBTR_TMP,
     &           DM_SUMLU,
     &           FUTURE_NIV2,
     &           MYID, KEEP, IERR )
            IF ( IERR == -1 )THEN
              CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
              CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
              IF (EXIT_FLAG) THEN
                GOTO 333
              ELSE
                GOTO 111
              ENDIF
            ELSE IF ( IERR .NE. 0 ) THEN
              WRITE(*,*) "Internal Error in DMUMPS_LOAD_MEM_UPDATE",IERR
              CALL MUMPS_ABORT()
            ENDIF
            DELTA_LOAD = ZERO
            DELTA_MEM  = ZERO
         ENDIF
      ENDIF
 333  CONTINUE
      IF(REMOVE_NODE_FLAG_MEM)THEN
         REMOVE_NODE_FLAG_MEM=.FALSE.
      ENDIF
      END SUBROUTINE DMUMPS_LOAD_MEM_UPDATE
      INTEGER FUNCTION DMUMPS_LOAD_LESS( K69, MEM_DISTRIB,MSG_SIZE )
      IMPLICIT NONE
      INTEGER i, NLESS, K69 
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      DOUBLE PRECISION LREF
      DOUBLE PRECISION MSG_SIZE
      NLESS = 0
      DO i=1,NPROCS
         IDWLOAD(i) = i - 1
      ENDDO
      WLOAD(1:NPROCS) = LOAD_FLOPS(0:NPROCS-1)
      IF(BDC_M2_FLOPS)THEN
         DO i=1,NPROCS
            WLOAD(i)=WLOAD(i)+NIV2(i)
         ENDDO
      ENDIF
      IF(K69 .gt. 1) THEN
         CALL DMUMPS_ARCHGENWLOAD(MEM_DISTRIB,MSG_SIZE,IDWLOAD,NPROCS)
      ENDIF
      LREF = LOAD_FLOPS(MYID)
      DO i=1, NPROCS
         IF (WLOAD(i).LT.LREF) NLESS=NLESS+1
      ENDDO
      DMUMPS_LOAD_LESS = NLESS
      RETURN
      END FUNCTION DMUMPS_LOAD_LESS
      SUBROUTINE DMUMPS_LOAD_SET_SLAVES(MEM_DISTRIB,MSG_SIZE,DEST,
     &     NSLAVES)
      IMPLICIT NONE
      INTEGER NSLAVES
      INTEGER DEST(NSLAVES)
      INTEGER, DIMENSION(0:NPROCS - 1) :: MEM_DISTRIB
      INTEGER i,J,NBDEST
      DOUBLE PRECISION MSG_SIZE
      IF ( NSLAVES.eq.NPROCS-1 ) THEN
        J = MYID+1
        DO i=1,NSLAVES
           J=J+1
           IF (J.GT.NPROCS) J=1
           DEST(i) = J - 1
        ENDDO
      ELSE
        DO i=1,NPROCS
           IDWLOAD(i) = i - 1
        ENDDO
        CALL MUMPS_SORT_DOUBLES(NPROCS, WLOAD, IDWLOAD)
         NBDEST = 0
         DO i=1, NSLAVES
            J = IDWLOAD(i)
            IF (J.NE.MYID) THEN
               NBDEST = NBDEST+1
               DEST(NBDEST) = J
            ENDIF
         ENDDO
         IF (NBDEST.NE.NSLAVES) THEN
            DEST(NSLAVES) = IDWLOAD(NSLAVES+1)
         ENDIF
         IF(BDC_MD)THEN
            J=NSLAVES+1
            do i=NSLAVES+1,NPROCS
               IF(IDWLOAD(i).NE.MYID)THEN
                  DEST(J)= IDWLOAD(i)
                  J=J+1
               ENDIF
            end do
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_SET_SLAVES
      SUBROUTINE DMUMPS_LOAD_END( INFO1, NSLAVES, IERR )
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: INFO1
      INTEGER, INTENT(IN)  :: NSLAVES
      INTEGER, INTENT(OUT) :: IERR
      INTEGER :: DUMMY_COMMUNICATOR
      IERR=0
      DUMMY_COMMUNICATOR = -999 
      CALL DMUMPS_CLEAN_PENDING( INFO1, KEEP_LOAD(1), BUF_LOAD_RECV(1),
     &     LBUF_LOAD_RECV,
     &     LBUF_LOAD_RECV_BYTES, DUMMY_COMMUNICATOR, COMM_LD,
     &     NSLAVES,
     &     .FALSE.,  
     &     .TRUE.    
     &     )
      DEALLOCATE( LOAD_FLOPS )
      DEALLOCATE( WLOAD )
      DEALLOCATE( IDWLOAD )
      DEALLOCATE(FUTURE_NIV2)
      IF(BDC_MD)THEN
         DEALLOCATE(MD_MEM)
         DEALLOCATE(LU_USAGE)
         DEALLOCATE(TAB_MAXS)
      ENDIF
      IF ( BDC_MEM ) DEALLOCATE( DM_MEM )
      IF ( BDC_POOL) DEALLOCATE( POOL_MEM )
      IF ( BDC_SBTR) THEN
         DEALLOCATE( SBTR_MEM )
         DEALLOCATE( SBTR_CUR )
         DEALLOCATE(SBTR_FIRST_POS_IN_POOL)
         NULLIFY(MY_FIRST_LEAF)
         NULLIFY(MY_NB_LEAF)
         NULLIFY(MY_ROOT_SBTR)
      ENDIF
      IF(KEEP_LOAD(76).EQ.4)THEN
         NULLIFY(DEPTH_FIRST_LOAD)
      ENDIF
      IF(KEEP_LOAD(76).EQ.5)THEN
         NULLIFY(COST_TRAV)
      ENDIF
      IF((KEEP_LOAD(76).EQ.4).OR.(KEEP_LOAD(76).EQ.6))THEN
         NULLIFY(DEPTH_FIRST_LOAD)
         NULLIFY(DEPTH_FIRST_SEQ_LOAD)
         NULLIFY(SBTR_ID_LOAD)
      ENDIF
      IF (BDC_M2_MEM.OR.BDC_M2_FLOPS) THEN
        DEALLOCATE(NB_SON,POOL_NIV2,POOL_NIV2_COST, NIV2)
      END IF
      IF((KEEP_LOAD(81).EQ.2).OR.(KEEP_LOAD(81).EQ.3))THEN
         DEALLOCATE(CB_COST_MEM)
         DEALLOCATE(CB_COST_ID)
      ENDIF
      NULLIFY(ND_LOAD)
      NULLIFY(KEEP_LOAD)
      NULLIFY(KEEP8_LOAD)
      NULLIFY(FILS_LOAD)
      NULLIFY(FRERE_LOAD)
      NULLIFY(PROCNODE_LOAD)
      NULLIFY(STEP_LOAD)
      NULLIFY(NE_LOAD)
      NULLIFY(CAND_LOAD)
      NULLIFY(STEP_TO_NIV2_LOAD)
      NULLIFY(DAD_LOAD)
      IF (BDC_SBTR.OR.BDC_POOL_MNG) THEN
         DEALLOCATE(MEM_SUBTREE)
         DEALLOCATE(SBTR_PEAK_ARRAY)
         DEALLOCATE(SBTR_CUR_ARRAY)
      ENDIF
      CALL DMUMPS_BUF_DEALL_LOAD_BUFFER( IERR )
      DEALLOCATE(BUF_LOAD_RECV)
      RETURN
      END SUBROUTINE DMUMPS_LOAD_END
      RECURSIVE SUBROUTINE DMUMPS_LOAD_RECV_MSGS(COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGTAG, MSGLEN, MSGSOU,COMM
      INTEGER IERR_MPI
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL FLAG
 10   CONTINUE
      CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &     FLAG, STATUS, IERR_MPI )
      IF (FLAG) THEN
        KEEP_LOAD(65)=KEEP_LOAD(65)+1
        KEEP_LOAD(267)=KEEP_LOAD(267)-1
        MSGTAG = STATUS( MPI_TAG )
        MSGSOU = STATUS( MPI_SOURCE )
        IF ( MSGTAG .NE. UPDATE_LOAD) THEN
          write(*,*) "Internal error 1 in DMUMPS_LOAD_RECV_MSGS",
     &    MSGTAG
          CALL MUMPS_ABORT()
        ENDIF
        CALL MPI_GET_COUNT(STATUS, MPI_PACKED, MSGLEN, IERR_MPI)
        IF ( MSGLEN > LBUF_LOAD_RECV_BYTES ) THEN
          write(*,*) "Internal error 2 in DMUMPS_LOAD_RECV_MSGS",
     &    MSGLEN, LBUF_LOAD_RECV_BYTES
          CALL MUMPS_ABORT()
        ENDIF
        CALL MPI_RECV( BUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES,
     &    MPI_PACKED, MSGSOU, MSGTAG, COMM_LD, STATUS, IERR_MPI)
        CALL DMUMPS_LOAD_PROCESS_MESSAGE( MSGSOU, BUF_LOAD_RECV,
     &  LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES )
        GOTO 10
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_RECV_MSGS
      RECURSIVE SUBROUTINE DMUMPS_LOAD_PROCESS_MESSAGE
     &   ( MSGSOU, BUFR, LBUFR, LBUFR_BYTES )
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER MSGSOU, LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INCLUDE 'mpif.h'
      INTEGER POSITION, WHAT, NSLAVES, i
      INTEGER IERR_MPI
      DOUBLE PRECISION LOAD_RECEIVED
      INTEGER INODE_RECEIVED,NCB_RECEIVED
      DOUBLE PRECISION SURF
      INTEGER, POINTER, DIMENSION (:) :: LIST_SLAVES
      DOUBLE PRECISION, POINTER, DIMENSION (:) :: LOAD_INCR
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     WHAT, 1, MPI_INTEGER,
     &     COMM_LD, IERR_MPI )
      IF ( WHAT == 0 ) THEN
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
      LOAD_FLOPS( MSGSOU ) = LOAD_FLOPS(MSGSOU) + LOAD_RECEIVED
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &       LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     &       COMM_LD, IERR_MPI )
          DM_MEM(MSGSOU)  = DM_MEM(MSGSOU) + LOAD_RECEIVED
          MAX_PEAK_STK=max(MAX_PEAK_STK,DM_MEM(MSGSOU))
        END IF
        IF(BDC_SBTR)THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &       LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     &       COMM_LD, IERR_MPI )
          SBTR_CUR(MSGSOU)=LOAD_RECEIVED
        ENDIF
        IF(BDC_MD)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &          LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     &          COMM_LD, IERR_MPI )
           IF(KEEP_LOAD(201).EQ.0)THEN
              LU_USAGE(MSGSOU)=LOAD_RECEIVED
           ENDIF
        ENDIF
      ELSEIF (( WHAT == 1).OR.(WHAT.EQ.19)) THEN
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     NSLAVES, 1, MPI_INTEGER,
     &     COMM_LD, IERR_MPI )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     INODE_RECEIVED, 1, MPI_INTEGER,
     &     COMM_LD, IERR_MPI )
        LIST_SLAVES => IDWLOAD
        LOAD_INCR => WLOAD
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LIST_SLAVES(1), NSLAVES, MPI_INTEGER,
     &     COMM_LD, IERR_MPI)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI)
        DO i = 1, NSLAVES
            LOAD_FLOPS(LIST_SLAVES(i)) =
     &      LOAD_FLOPS(LIST_SLAVES(i)) + 
     &      LOAD_INCR(i)
        END DO
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI)
          DO i = 1, NSLAVES
              DM_MEM(LIST_SLAVES(i)) = DM_MEM(LIST_SLAVES(i)) + 
     &        LOAD_INCR(i)
              MAX_PEAK_STK=max(MAX_PEAK_STK,DM_MEM(LIST_SLAVES(i)))
          END DO
        END IF
        IF(WHAT.EQ.19)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &          LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &          COMM_LD, IERR_MPI)
           CALL DMUMPS_LOAD_CLEAN_MEMINFO_POOL(INODE_RECEIVED)
           CB_COST_ID(POS_ID)=INODE_RECEIVED
           CB_COST_ID(POS_ID+1)=NSLAVES
           CB_COST_ID(POS_ID+2)=POS_MEM
           POS_ID=POS_ID+3
           DO i=1,NSLAVES
              WRITE(*,*)MYID,':',LIST_SLAVES(i),'->',LOAD_INCR(i)
              CB_COST_MEM(POS_MEM)=int(LIST_SLAVES(i),8)
              POS_MEM=POS_MEM+1
              CB_COST_MEM(POS_MEM)=int(LOAD_INCR(i),8)
              POS_MEM=POS_MEM+1
           ENDDO
        ENDIF
        NULLIFY( LIST_SLAVES )
        NULLIFY( LOAD_INCR )
      ELSE IF (WHAT == 2 ) THEN
         IF ( .not. BDC_POOL ) THEN
          WRITE(*,*) "Internal error 2 in DMUMPS_LOAD_PROCESS_MESSAGE"
          CALL MUMPS_ABORT()
        END IF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
        POOL_MEM(MSGSOU)=LOAD_RECEIVED
      ELSE IF ( WHAT == 3 ) THEN
         IF ( .NOT. BDC_SBTR) THEN
          WRITE(*,*) "Internal error 3 in DMUMPS_LOAD_PROCESS_MESSAGE"
          CALL MUMPS_ABORT()
        ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
        SBTR_MEM(MSGSOU)=SBTR_MEM(MSGSOU)+LOAD_RECEIVED
      ELSE IF (WHAT == 4) THEN
        FUTURE_NIV2(MSGSOU+1)=0
        IF(BDC_MD)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &          SURF, 1, MPI_DOUBLE_PRECISION,
     &          COMM_LD, IERR_MPI )
          MD_MEM(MSGSOU)=999999999_8
          TAB_MAXS(MSGSOU)=TAB_MAXS(MSGSOU)+int(SURF,8)
        ENDIF
        IF(BDC_M2_MEM.OR.BDC_M2_FLOPS)THEN
        ENDIF
      ELSE IF (WHAT == 5) THEN
         IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
            WRITE(*,*) "Internal error 7 in DMUMPS_LOAD_PROCESS_MESSAGE"
            CALL MUMPS_ABORT()
         ENDIF
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     INODE_RECEIVED, 1,
     &     MPI_INTEGER,
     &     COMM_LD, IERR_MPI )
         IF(BDC_M2_MEM) THEN
            CALL DMUMPS_PROCESS_NIV2_MEM_MSG(INODE_RECEIVED)
         ELSEIF(BDC_M2_FLOPS) THEN
            CALL DMUMPS_PROCESS_NIV2_FLOPS_MSG(INODE_RECEIVED)
         ENDIF
         IF((KEEP_LOAD(81).EQ.2).OR.(KEEP_LOAD(81).EQ.3))THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           INODE_RECEIVED, 1,
     &           MPI_INTEGER,
     &           COMM_LD, IERR_MPI )   
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              NCB_RECEIVED, 1,
     &              MPI_INTEGER,
     &              COMM_LD, IERR_MPI )   
            IF(
     &          MUMPS_TYPENODE(PROCNODE_LOAD(STEP_LOAD(INODE_RECEIVED)),
     &                         KEEP_LOAD(199)).EQ.1
     &        )THEN
               CB_COST_ID(POS_ID)=INODE_RECEIVED
               CB_COST_ID(POS_ID+1)=1
               CB_COST_ID(POS_ID+2)=POS_MEM
               POS_ID=POS_ID+3
               CB_COST_MEM(POS_MEM)=int(MSGSOU,8)
               POS_MEM=POS_MEM+1
               CB_COST_MEM(POS_MEM)=int(NCB_RECEIVED,8)*
     &              int(NCB_RECEIVED,8)
               POS_MEM=POS_MEM+1
            ENDIF
         ENDIF
      ELSE IF ( WHAT == 6 ) THEN
         IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
            WRITE(*,*) "Internal error 8 in DMUMPS_LOAD_PROCESS_MESSAGE"
            CALL MUMPS_ABORT()
         ENDIF
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
         IF(BDC_M2_MEM) THEN
            NIV2(MSGSOU+1) = LOAD_RECEIVED
         ELSEIF(BDC_M2_FLOPS) THEN
            NIV2(MSGSOU+1) = NIV2(MSGSOU+1) + LOAD_RECEIVED
            IF(NIV2(MSGSOU+1).LT.0.0D0)THEN
               IF(abs(NIV2(MSGSOU+1)) .LE. 1.0D-3) THEN
                  NIV2(MSGSOU+1)=0.0D0
               ELSE
                  WRITE(*,*)'problem with NIV2_FLOPS message',
     &                 NIV2(MSGSOU+1),MSGSOU,LOAD_RECEIVED
                  CALL MUMPS_ABORT()
               ENDIF
            ENDIF
         ENDIF
      ELSEIF(WHAT == 17)THEN
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
         IF(BDC_M2_MEM) THEN
            NIV2(MSGSOU+1) = LOAD_RECEIVED
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           LOAD_RECEIVED, 1,
     &           MPI_DOUBLE_PRECISION,
     &           COMM_LD, IERR_MPI )
            IF(BDC_MD)THEN
               DM_MEM(MYID)=DM_MEM(MYID)+LOAD_RECEIVED
            ELSEIF(BDC_POOL)THEN
               POOL_MEM(MSGSOU)=LOAD_RECEIVED
            ENDIF
         ELSEIF(BDC_M2_FLOPS) THEN
            NIV2(MSGSOU+1) = NIV2(MSGSOU+1) + LOAD_RECEIVED            
            IF(NIV2(MSGSOU+1).LT.0.0D0)THEN
               IF(abs(NIV2(MSGSOU+1)) .LE. 1.0D-3) THEN
                  NIV2(MSGSOU+1)=0.0D0
               ELSE
                  WRITE(*,*)'problem with NIV2_FLOPS message',
     &                 NIV2(MSGSOU+1),MSGSOU,LOAD_RECEIVED
                  CALL MUMPS_ABORT()
               ENDIF
            ENDIF
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &           LOAD_RECEIVED, 1,
     &           MPI_DOUBLE_PRECISION,
     &           COMM_LD, IERR_MPI )
            LOAD_FLOPS( MSGSOU ) = LOAD_FLOPS(MSGSOU) + LOAD_RECEIVED
         ENDIF
      ELSEIF ( WHAT == 7 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 4
     &in DMUMPS_LOAD_PROCESS_MESSAGE'
            CALL MUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     NSLAVES, 1, MPI_INTEGER,
     &     COMM_LD, IERR_MPI )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     INODE_RECEIVED, 1, MPI_INTEGER,
     &     COMM_LD, IERR_MPI )
        LIST_SLAVES => IDWLOAD
        LOAD_INCR => WLOAD
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LIST_SLAVES(1), NSLAVES, MPI_INTEGER,
     &     COMM_LD, IERR_MPI )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_INCR(1), NSLAVES, MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
        DO i = 1, NSLAVES
            MD_MEM(LIST_SLAVES(i)) =
     &      MD_MEM(LIST_SLAVES(i)) + 
     &      int(LOAD_INCR(i),8)
            IF(FUTURE_NIV2(LIST_SLAVES(i)+1).EQ.0)THEN
               MD_MEM(LIST_SLAVES(i))=999999999_8
            ENDIF
        END DO
      ELSEIF ( WHAT == 8 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 5
     &in DMUMPS_LOAD_PROCESS_MESSAGE'
            CALL MUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
        MD_MEM(MSGSOU)=MD_MEM(MSGSOU)+int(LOAD_RECEIVED,8)
        IF(FUTURE_NIV2(MSGSOU+1).EQ.0)THEN
           MD_MEM(MSGSOU)=999999999_8
        ENDIF
      ELSEIF ( WHAT == 9 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 6
     &in DMUMPS_LOAD_PROCESS_MESSAGE'
            CALL MUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &     LOAD_RECEIVED, 1,
     &     MPI_DOUBLE_PRECISION,
     &     COMM_LD, IERR_MPI )
        TAB_MAXS(MSGSOU)=int(LOAD_RECEIVED,8)
      ELSE
          WRITE(*,*) "Internal error 1 in DMUMPS_LOAD_PROCESS_MESSAGE"
          CALL MUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_PROCESS_MESSAGE
      integer function DMUMPS_LOAD_LESS_CAND
     &                 (MEM_DISTRIB,CAND,
     &                  K69,
     &                  SLAVEF,MSG_SIZE,
     &                  NMB_OF_CAND )
      implicit none      
      integer, intent(in) :: K69, SLAVEF
      INTEGER, intent(in) :: CAND(SLAVEF+1)
      INTEGER, DIMENSION(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      INTEGER, intent(out) :: NMB_OF_CAND
      integer i,nless
      DOUBLE PRECISION lref
      DOUBLE PRECISION MSG_SIZE
      nless = 0 
      NMB_OF_CAND=CAND(SLAVEF+1)
      do i=1,NMB_OF_CAND
         WLOAD(i)=LOAD_FLOPS(CAND(i))
         IF(BDC_M2_FLOPS)THEN
            WLOAD(i)=WLOAD(i)+NIV2(CAND(i)+1)
         ENDIF
      end do
      IF(K69 .gt. 1) THEN
         CALL DMUMPS_ARCHGENWLOAD(MEM_DISTRIB,MSG_SIZE,
     &        CAND,NMB_OF_CAND)
      ENDIF
      lref = LOAD_FLOPS(MYID)
      do i=1, NMB_OF_CAND
         if (WLOAD(i).lt.lref) nless=nless+1
      end do 
      DMUMPS_LOAD_LESS_CAND = nless
      return
      end function DMUMPS_LOAD_LESS_CAND
      subroutine DMUMPS_LOAD_SET_SLAVES_CAND
     &           (MEM_DISTRIB,CAND, 
     &
     &            SLAVEF,
     &            nslaves_inode, DEST)
      implicit none
      integer, intent(in) :: nslaves_inode, SLAVEF
      integer, intent(in) :: CAND(SLAVEF+1)
      integer, dimension(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      integer, intent(out) :: DEST(CAND(SLAVEF+1))
      integer i,j,NMB_OF_CAND
      external MUMPS_SORT_DOUBLES
      NMB_OF_CAND = CAND(SLAVEF+1)
      if(nslaves_inode.ge.NPROCS .or.
     &     nslaves_inode.gt.NMB_OF_CAND) then
         write(*,*)'Internal error in DMUMPS_LOAD_SET_SLAVES_CAND',
     &   nslaves_inode, NPROCS, NMB_OF_CAND
         CALL MUMPS_ABORT()
      end if
      if (nslaves_inode.eq.NPROCS-1) then
         j=MYID+1
         do i=1,nslaves_inode
            if(j.ge.NPROCS) j=0
            DEST(i)=j
            j=j+1
         end do
      else
        do i=1,NMB_OF_CAND
               IDWLOAD(i)=i
        end do
        call MUMPS_SORT_DOUBLES(NMB_OF_CAND,
     &       WLOAD(1),IDWLOAD(1) )
        do i=1,nslaves_inode
           DEST(i)= CAND(IDWLOAD(i))
        end do
        IF(BDC_MD)THEN
           do i=nslaves_inode+1,NMB_OF_CAND
              DEST(i)= CAND(IDWLOAD(i))
           end do
        ENDIF
      end if   
      return
      end subroutine DMUMPS_LOAD_SET_SLAVES_CAND
      SUBROUTINE DMUMPS_INIT_ALPHA_BETA(K69)
      IMPLICIT NONE
      INTEGER K69
      IF (K69 .LE. 4) THEN
         ALPHA = 0.0d0
         BETA = 0.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 5) THEN
         ALPHA = 0.5d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 6) THEN
         ALPHA = 0.5d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 7) THEN
         ALPHA = 0.5d0
         BETA = 150000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 8) THEN
         ALPHA = 1.0d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 9) THEN
         ALPHA = 1.0d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 10) THEN
         ALPHA = 1.0d0
         BETA = 150000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 11) THEN
         ALPHA = 1.5d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 12) THEN
         ALPHA = 1.5d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      ALPHA = 1.5d0
      BETA = 150000.0d0
      RETURN
      END SUBROUTINE DMUMPS_INIT_ALPHA_BETA
      SUBROUTINE DMUMPS_ARCHGENWLOAD(MEM_DISTRIB,MSG_SIZE,ARRAY_ADM,LEN)
      IMPLICIT NONE
      INTEGER i,LEN
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      DOUBLE PRECISION MSG_SIZE,FORBIGMSG
      INTEGER ARRAY_ADM(LEN)
      DOUBLE PRECISION MY_LOAD
      FORBIGMSG = 1.0d0
      IF (K69 .lt.2) THEN
         RETURN
      ENDIF
      IF(BDC_M2_FLOPS)THEN
         MY_LOAD=LOAD_FLOPS(MYID)+NIV2(MYID+1)
      ELSE
         MY_LOAD=LOAD_FLOPS(MYID)
      ENDIF
      IF((MSG_SIZE * dble(K35) ) .gt. 3200000.0d0) THEN
         FORBIGMSG = 2.0d0
      ENDIF
      IF (K69 .le. 4) THEN
         DO i = 1,LEN
            IF ((MEM_DISTRIB(ARRAY_ADM(i)) .EQ. 1) .AND.
     &      WLOAD(i) .LT. MY_LOAD ) THEN
               WLOAD(i) = WLOAD(i)/MY_LOAD
            ELSE
              IF ( MEM_DISTRIB(ARRAY_ADM(i)) .NE. 1 ) THEN
                WLOAD(i) = WLOAD(i) *
     &              dble(MEM_DISTRIB(ARRAY_ADM(i)))
     &              * FORBIGMSG
     &              + dble(2)
              ENDIF
            ENDIF
         ENDDO
         RETURN
      ENDIF
      DO i = 1,LEN
         IF ((MEM_DISTRIB(ARRAY_ADM(i)) .EQ. 1) .AND.
     &        WLOAD(i) .LT. MY_LOAD ) THEN
            WLOAD(i) = WLOAD(i) /  MY_LOAD
         ELSE
            IF(MEM_DISTRIB(ARRAY_ADM(i)) .NE. 1) THEN     
               WLOAD(i) = (WLOAD(i) +
     &              ALPHA * MSG_SIZE * dble(K35)  +
     &              BETA) * FORBIGMSG
            ENDIF
         ENDIF
      ENDDO  
      RETURN
      END SUBROUTINE DMUMPS_ARCHGENWLOAD
      SUBROUTINE DMUMPS_LOAD_MASTER_2_ALL(MYID, SLAVEF, COMM,
     &     TAB_POS, NASS, KEEP,KEEP8, LIST_SLAVES, NSLAVES,INODE)
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: MYID, SLAVEF, COMM, NASS, NSLAVES
      INTEGER, INTENT (IN) :: TAB_POS(SLAVEF+2)
      INTEGER, INTENT (IN) :: LIST_SLAVES( NSLAVES )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NCB, NFRONT, NBROWS_SLAVE
      INTEGER i, IERR,WHAT,INODE, allocok
      LOGICAL :: EXIT_FLAG
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MEM_INCREMENT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FLOPS_INCREMENT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CB_BAND
      ALLOCATE(MEM_INCREMENT(NSLAVES), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of MEM_INCREMENT '
     &        //  'in routine DMUMPS_LOAD_MASTER_2_ALL'
         CALL MUMPS_ABORT()
      endif
      ALLOCATE(FLOPS_INCREMENT(NSLAVES), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of FLOPS_INCREMENT '
     &        //    'in routine DMUMPS_LOAD_MASTER_2_ALL'
         CALL MUMPS_ABORT()
      endif
      ALLOCATE(CB_BAND(NSLAVES), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of CB_BAND '
     &        //    'in routine DMUMPS_LOAD_MASTER_2_ALL'
         CALL MUMPS_ABORT()
      endif
      IF((KEEP(81).NE.2).AND.(KEEP(81).NE.3))THEN
         WHAT=1
      ELSE
         WHAT=19
      ENDIF
      FUTURE_NIV2(MYID+1) = FUTURE_NIV2(MYID+1) - 1
      IF ( FUTURE_NIV2(MYID+1) < 0 ) THEN
        WRITE(*,*) "Internal error in DMUMPS_LOAD_MASTER_2_ALL"
        CALL MUMPS_ABORT()
      ENDIF
      IF ( FUTURE_NIV2(MYID + 1) == 0 ) THEN
 112    CONTINUE
        CALL DMUMPS_BUF_SEND_NOT_MSTR(COMM,MYID,SLAVEF,
     &       dble(MAX_SURF_MASTER),KEEP,IERR)
        IF (IERR == -1 ) THEN
          CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
          CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
          IF (EXIT_FLAG) THEN
             GOTO 100 
          ELSE
             GOTO 112
          ENDIF
       ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in DMUMPS_LOAD_MASTER_2_ALL",
     &    IERR
          CALL MUMPS_ABORT()
        ENDIF
      TAB_MAXS(MYID) = TAB_MAXS(MYID) + int(MAX_SURF_MASTER,8)
      ENDIF
      IF ( NSLAVES /= TAB_POS(SLAVEF + 2) ) THEN
        write(*,*) "Error 1 in DMUMPS_LOAD_MASTER_2_ALL",
     &             NSLAVES, TAB_POS(SLAVEF+2)
        CALL MUMPS_ABORT()
      ENDIF
      NCB = TAB_POS(NSLAVES+1) - 1
      NFRONT = NCB + NASS
      DO i = 1, NSLAVES
         NBROWS_SLAVE = TAB_POS(i+1) - TAB_POS(i)
         IF ( KEEP(50) == 0 ) THEN
            FLOPS_INCREMENT( i ) = (dble(NBROWS_SLAVE)*dble( NASS ))+
     &           dble(NBROWS_SLAVE) * dble(NASS) *
     &           dble(2*NFRONT-NASS-1)
         ELSE
            FLOPS_INCREMENT( i ) = dble(NBROWS_SLAVE) * dble(NASS ) *
     &           dble( 2 * ( NASS + TAB_POS(i+1) - 1 ) 
     &           - NBROWS_SLAVE - NASS + 1 )
         ENDIF
         IF ( BDC_MEM ) THEN
            IF ( KEEP(50) == 0 ) THEN
               MEM_INCREMENT( i ) = dble(NBROWS_SLAVE) *
     &              dble(NFRONT)
            ELSE
               MEM_INCREMENT( i ) = dble(NBROWS_SLAVE) *
     &              dble( NASS + TAB_POS(i+1) - 1 )
            END IF
         ENDIF
         IF((KEEP(81).NE.2).AND.(KEEP(81).NE.3))THEN
            CB_BAND(i)=dble(-999999)
         ELSE
            IF ( KEEP(50) == 0 ) THEN
               CB_BAND( i ) = dble(NBROWS_SLAVE) *
     &              dble(NFRONT-NASS)
            ELSE
               CB_BAND( i ) = dble(NBROWS_SLAVE) *
     &              dble(TAB_POS(i+1)-1)
            END IF
         ENDIF
      END DO
      IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
         CB_COST_ID(POS_ID)=INODE
         CB_COST_ID(POS_ID+1)=NSLAVES
         CB_COST_ID(POS_ID+2)=POS_MEM
         POS_ID=POS_ID+3
         DO i=1,NSLAVES
            CB_COST_MEM(POS_MEM)=int(LIST_SLAVES(i),8)
            POS_MEM=POS_MEM+1
            CB_COST_MEM(POS_MEM)=int(CB_BAND(i),8)
            POS_MEM=POS_MEM+1
         ENDDO
      ENDIF
 111  CONTINUE
      CALL DMUMPS_BUF_BCAST_ARRAY(BDC_MEM, COMM, MYID, SLAVEF,
     &     FUTURE_NIV2,
     &     NSLAVES, LIST_SLAVES,INODE,
     &     MEM_INCREMENT,
     &     FLOPS_INCREMENT,CB_BAND, WHAT, KEEP, IERR)
        IF ( IERR == -1 ) THEN
          CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
          CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
          IF (EXIT_FLAG) THEN
             GOTO 100
           ELSE
             GOTO 111
           ENDIF
       ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in DMUMPS_LOAD_MASTER_2_ALL",
     &    IERR
          CALL MUMPS_ABORT()
        ENDIF
      IF (FUTURE_NIV2(MYID+1) .NE. 0) THEN
        DO i = 1, NSLAVES
          LOAD_FLOPS(LIST_SLAVES(i)) = LOAD_FLOPS(LIST_SLAVES(i))
     &       +  FLOPS_INCREMENT(i)
          IF ( BDC_MEM ) THEN
            DM_MEM(LIST_SLAVES(i)) = DM_MEM(LIST_SLAVES(i))
     &       +  MEM_INCREMENT(i)
          END IF
        ENDDO
      ENDIF
 100  CONTINUE
      DEALLOCATE(MEM_INCREMENT,FLOPS_INCREMENT,CB_BAND)
      RETURN
      END SUBROUTINE DMUMPS_LOAD_MASTER_2_ALL
      SUBROUTINE DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &     POOL, LPOOL,
     &     PROCNODE, KEEP,KEEP8, SLAVEF, COMM, MYID, STEP, N,
     &     ND, FILS )
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER LPOOL, SLAVEF, COMM, MYID 
      INTEGER N, KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER POOL( LPOOL ), PROCNODE( KEEP(28) ), STEP( N )
      INTEGER ND( KEEP(28) ), FILS( N )
      INTEGER i, INODE, NELIM, NFR, LEVEL, IERR, WHAT
      DOUBLE PRECISION COST
      LOGICAL :: EXIT_FLAG
      INTEGER NBINSUBTREE,NBTOP,INSUBTREE
      INTEGER MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPENODE
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF(BDC_MD)THEN
         RETURN
      ENDIF
      IF((KEEP(76).EQ.0).OR.(KEEP(76).EQ.2))THEN
         IF(NBTOP.NE.0)THEN
            DO i = LPOOL-NBTOP-2, min(LPOOL-3,LPOOL-NBTOP-2+3)
               INODE = POOL( i )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=dble(0) 
            GOTO 30
         ELSE
            DO i = NBINSUBTREE, max(1,NBINSUBTREE-3), -1
               INODE = POOL( i )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=dble(0) 
            GOTO 30
         ENDIF
      ELSE
         IF(KEEP(76).EQ.1)THEN
            IF(INSUBTREE.EQ.1)THEN
               DO i = NBINSUBTREE, max(1,NBINSUBTREE-3), -1
                  INODE = POOL( i )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=dble(0) 
               GOTO 30
            ELSE
               DO i = LPOOL-NBTOP-2, min(LPOOL-3,LPOOL-NBTOP-2+3)
                  INODE = POOL( i )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=dble(0) 
               GOTO 30
            ENDIF
         ELSE
            WRITE(*,*)
     &      'Internal error: Unknown pool management strategy'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
 20   CONTINUE
        i = INODE
        NELIM = 0
 10     CONTINUE
        IF ( i > 0 ) THEN
          NELIM = NELIM + 1
          i = FILS(i)
          GOTO 10
        ENDIF
        NFR = ND( STEP(INODE) )
        LEVEL = MUMPS_TYPENODE( PROCNODE(STEP(INODE)), KEEP(199) )
        IF (LEVEL .EQ. 1) THEN
          COST = dble( NFR ) * dble( NFR )
        ELSE
          IF ( KEEP(50) == 0 ) THEN
            COST = dble( NFR ) * dble( NELIM )
          ELSE
            COST = dble( NELIM ) * dble( NELIM )
          ENDIF
        ENDIF
 30   CONTINUE
      IF ( abs(POOL_LAST_COST_SENT-COST).GT.DM_THRES_MEM ) THEN
        WHAT = 2
 111    CONTINUE
        CALL DMUMPS_BUF_BROADCAST( WHAT,
     &         COMM, SLAVEF,
     &               FUTURE_NIV2,
     &         COST, dble(0), MYID, KEEP, IERR  )
        POOL_LAST_COST_SENT = COST
        POOL_MEM(MYID)=COST
        IF ( IERR == -1 )THEN
          CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
          CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
          IF (EXIT_FLAG) THEN
          ELSE
             GOTO 111
          ENDIF
       ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in DMUMPS_LOAD_POOL_UPD_NEW_POOL",
     &    IERR
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_LOAD_POOL_UPD_NEW_POOL
      SUBROUTINE DMUMPS_LOAD_SBTR_UPD_NEW_POOL(
     &     OK,INODE,POOL,LPOOL,MYID,SLAVEF,COMM,KEEP,KEEP8)
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER LPOOL,MYID,SLAVEF,COMM,INODE
      INTEGER POOL(LPOOL),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER WHAT,IERR
      LOGICAL OK
      DOUBLE PRECISION COST
      LOGICAL FLAG, EXIT_FLAG
      EXTERNAL MUMPS_ROOTSSARBR,MUMPS_IN_OR_ROOT_SSARBR
      LOGICAL MUMPS_ROOTSSARBR,MUMPS_IN_OR_ROOT_SSARBR
      IF((INODE.LE.0).OR.(INODE.GT.N_LOAD)) THEN
         RETURN
      ENDIF
      IF (.NOT.MUMPS_IN_OR_ROOT_SSARBR(
     &     PROCNODE_LOAD(STEP_LOAD(INODE)), KEEP(199))
     &   ) THEN
         RETURN         
      ENDIF
      IF(MUMPS_ROOTSSARBR(PROCNODE_LOAD(STEP_LOAD(INODE)),
     &         KEEP(199)))THEN
         IF(NE_LOAD(STEP_LOAD(INODE)).EQ.0)THEN
            RETURN
         ENDIF
      ENDIF
      FLAG=.FALSE.
      IF(INDICE_SBTR.LE.NB_SUBTREES)THEN
         IF(INODE.EQ.MY_FIRST_LEAF(INDICE_SBTR))THEN
            FLAG=.TRUE.
         ENDIF
      ENDIF
      IF(FLAG)THEN
         SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY)=MEM_SUBTREE(INDICE_SBTR)
         SBTR_CUR_ARRAY(INDICE_SBTR_ARRAY)=SBTR_CUR(MYID)
         INDICE_SBTR_ARRAY=INDICE_SBTR_ARRAY+1
         WHAT = 3
         IF(dble(MEM_SUBTREE(INDICE_SBTR)).GE.DM_THRES_MEM)THEN
 111        CONTINUE
            CALL DMUMPS_BUF_BROADCAST(
     &           WHAT, COMM, SLAVEF,
     &           FUTURE_NIV2,
     &           dble(MEM_SUBTREE(INDICE_SBTR)), dble(0),
     &           MYID, KEEP, IERR  )
            IF ( IERR == -1 )THEN
               CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
               CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
               IF (EXIT_FLAG) THEN
               ELSE
                  GOTO 111
               ENDIF
            ELSE IF ( IERR .NE. 0 ) THEN
               WRITE(*,*)
     &         "Internal Error 1 in DMUMPS_LOAD_SBTR_UPD_NEW_POOL",
     &         IERR
               CALL MUMPS_ABORT()
            ENDIF
         ENDIF
         SBTR_MEM(MYID)=SBTR_MEM(MYID)+
     &        dble(MEM_SUBTREE(INDICE_SBTR))
         INDICE_SBTR=INDICE_SBTR+1
         IF(INSIDE_SUBTREE.EQ.0)THEN
            INSIDE_SUBTREE=1
         ENDIF
      ELSE 
         IF(INODE.EQ.MY_ROOT_SBTR(INDICE_SBTR-1))THEN
            WHAT = 3
            COST=-SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY-1)
            IF(abs(COST).GE.DM_THRES_MEM)THEN
 112           CONTINUE
               CALL DMUMPS_BUF_BROADCAST(
     &              WHAT, COMM, SLAVEF,
     &              FUTURE_NIV2,
     &              COST, dble(0), MYID, KEEP, IERR  )
               IF ( IERR == -1 )THEN
                  CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
                  CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
                  IF (EXIT_FLAG) THEN
                  ELSE
                     GOTO 112
                  ENDIF
               ELSE IF ( IERR .NE. 0 ) THEN
                  WRITE(*,*)
     &        "Internal Error 3 in DMUMPS_LOAD_SBTR_UPD_NEW_POOL",
     &        IERR
                  CALL MUMPS_ABORT()
               ENDIF
            ENDIF
            INDICE_SBTR_ARRAY=INDICE_SBTR_ARRAY-1            
            SBTR_MEM(MYID)=SBTR_MEM(MYID)-
     &           SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY)
            SBTR_CUR(MYID)=SBTR_CUR_ARRAY(INDICE_SBTR_ARRAY)
            IF(INDICE_SBTR_ARRAY.EQ.1)THEN
               SBTR_CUR(MYID)=dble(0)
               INSIDE_SUBTREE=0
            ENDIF
         ENDIF
         ENDIF 
      RETURN
      END SUBROUTINE DMUMPS_LOAD_SBTR_UPD_NEW_POOL
      SUBROUTINE DMUMPS_SET_PARTI_ACTV_MEM
     &      (SLAVEF,KEEP,KEEP8,PROCS,MEM_DISTRIB,NCB,NFRONT,
     &       NSLAVES_NODE,TAB_POS,
     &       SLAVES_LIST,SIZE_SLAVES_LIST,MYID)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER NUMBER_OF_PROCS,K47, K48, K50
      INTEGER(8) :: K821
      DOUBLE PRECISION DK821
      INTEGER J
      INTEGER KMIN, KMAX
      INTEGER OTHERS,CHOSEN,SMALL_SET,ACC
      DOUBLE PRECISION SOMME,TMP_SUM
      INTEGER AFFECTED
      INTEGER ADDITIONNAL_ROWS,i,X,REF,POS
      INTEGER(8)::TOTAL_MEM
      LOGICAL FORCE_CAND
      DOUBLE PRECISION TEMP(SLAVEF),PEAK
      INTEGER TEMP_ID(SLAVEF),NB_ROWS(SLAVEF)
      EXTERNAL MPI_WTIME
      DOUBLE PRECISION MPI_WTIME
      IF (KEEP8(21) .GT. 0_8) THEN
      write(*,*)MYID,
     & ": Internal Error 1 in DMUMPS_SET_PARTI_ACTV_MEM"
      CALL MUMPS_ABORT()
      ENDIF
      K821=abs(KEEP8(21))
      DK821=dble(K821)
      K50=KEEP(50)
      K48=KEEP(48)
      K47=KEEP(47)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
      END IF
      IF(K48.NE.4)THEN
         WRITE(*,*)'DMUMPS_COMPUTE_PARTI_ACTV_MEM_K821
     &      should be called with KEEP(48) different from 4'
         CALL MUMPS_ABORT()
      ENDIF
         KMIN=1
         KMAX=int(K821/int(NFRONT,8))
         IF(FORCE_CAND)THEN
            DO i=1,PROCS(SLAVEF+1)
               WLOAD(i)=DM_MEM(PROCS(i))
               IDWLOAD(i)=PROCS(i)
            ENDDO
            NUMBER_OF_PROCS=PROCS(SLAVEF+1)
            OTHERS=NUMBER_OF_PROCS
         ELSE
            NUMBER_OF_PROCS=SLAVEF
            WLOAD(1:SLAVEF) = DM_MEM(0:NUMBER_OF_PROCS-1)
            DO i=1,NUMBER_OF_PROCS
               IDWLOAD(i) = i - 1
            ENDDO
            OTHERS=NUMBER_OF_PROCS-1
         ENDIF
         NB_ROWS=0
         CALL MUMPS_SORT_DOUBLES(NUMBER_OF_PROCS, WLOAD, IDWLOAD)
         TOTAL_MEM=int(NCB,8)*int(NFRONT,8)
         SOMME=dble(0)
         J=1
         PEAK=dble(0)
         DO i=1,NUMBER_OF_PROCS
            IF((IDWLOAD(i).NE.MYID))THEN
               PEAK=max(PEAK,WLOAD(i))
               TEMP_ID(J)=IDWLOAD(i)
               TEMP(J)=WLOAD(i)
                IF(BDC_SBTR)THEN
                   TEMP(J)=TEMP(J)+SBTR_MEM(IDWLOAD(i))-
     &                  SBTR_CUR(IDWLOAD(i))
                ENDIF
                IF(BDC_POOL)THEN
                   TEMP(J)=TEMP(J)+POOL_MEM(TEMP_ID(J))
                ENDIF
                IF(BDC_M2_MEM)THEN
                   TEMP(J)=TEMP(J)+NIV2(TEMP_ID(J)+1)
                ENDIF
                J=J+1
            ENDIF
         ENDDO
         NUMBER_OF_PROCS=J-1
         CALL MUMPS_SORT_DOUBLES(NUMBER_OF_PROCS, TEMP, TEMP_ID)
         IF(K50.EQ.0)THEN
           PEAK=max(PEAK,
     &       DM_MEM(MYID)+dble(NFRONT)*dble(NFRONT-NCB))
         ELSE
           PEAK=max(PEAK,
     &       DM_MEM(MYID)+dble(NFRONT-NCB)*dble(NFRONT-NCB))
         ENDIF
         PEAK=max(PEAK,TEMP(OTHERS))
         SOMME=dble(0)
         DO i=1,NUMBER_OF_PROCS
           SOMME=SOMME+TEMP(OTHERS)-TEMP(i)
         ENDDO               
         IF(SOMME.LE.dble(TOTAL_MEM)) THEN
            GOTO 096
         ENDIF
 096     CONTINUE
         SOMME=dble(0)
         DO i=1,OTHERS
            SOMME=SOMME+TEMP(OTHERS)-TEMP(i)
         ENDDO
         IF(dble(TOTAL_MEM).GE.SOMME) THEN
            AFFECTED=0
            CHOSEN=0
            ACC=0
            DO i=1,OTHERS
               IF(K50.EQ.0)THEN
                  IF((TEMP(OTHERS)-TEMP(i)).GT.DK821)THEN
                     TMP_SUM=DK821
                  ELSE
                     TMP_SUM=TEMP(OTHERS)-TEMP(i)
                  ENDIF
                  X=int(TMP_SUM/dble(NFRONT))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
               ENDIF
               IF(K50.NE.0)THEN
                  IF((TEMP(OTHERS)-TEMP(i)).GT.DK821)THEN
                     TMP_SUM=DK821
                  ELSE
                     TMP_SUM=TEMP(OTHERS)-TEMP(i)
                  ENDIF
                  X=int((-dble(NFRONT-NCB+ACC)
     &                 +sqrt(((dble(NFRONT-NCB+ACC)*
     &                 dble(NFRONT-NCB+ACC))+dble(4)*
     &                 (TMP_SUM))))/
     &                 dble(2))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
                  IF(X.LE.0) THEN
                     WRITE(*,*)"Internal Error 2 in
     &                    DMUMPS_SET_PARTI_ACTV_MEM"
                    CALL MUMPS_ABORT()
                  ENDIF
               ENDIF
               NB_ROWS(i)=X
               CHOSEN=CHOSEN+1
               ACC=ACC+X
               IF(NCB-ACC.LT.KMIN) GOTO 111
               IF(NCB.EQ.ACC) GOTO 111
               ENDDO
 111           CONTINUE
               IF((ACC.GT.NCB))THEN
                  X=0
                  DO i=1,OTHERS
                     X=X+NB_ROWS(i)
                  ENDDO
                  WRITE(*,*)'NCB=',NCB,',SOMME=',X
                  WRITE(*,*)MYID,
     &               ": Internal Error 3 in DMUMPS_SET_PARTI_ACTV_MEM"
                  CALL MUMPS_ABORT()
               ENDIF
               IF((NCB.NE.ACC))THEN
                  IF(K50.NE.0)THEN
                     IF(CHOSEN.NE.0)THEN
                        ADDITIONNAL_ROWS=NCB-ACC
                        NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+ADDITIONNAL_ROWS
                     ELSE
                        TMP_SUM=dble(TOTAL_MEM)/dble(NUMBER_OF_PROCS)
                        CHOSEN=0
                        ACC=0
                        DO i=1,OTHERS
                           X=int((-dble(NFRONT-NCB+ACC)
     &                          +sqrt(((dble(NFRONT-NCB+ACC)*
     &                          dble(NFRONT-NCB+ACC))+dble(4)*
     &                          (TMP_SUM))))/
     &                          dble(2))
                           IF((ACC+X).GT.NCB) X=NCB-ACC
                           NB_ROWS(i)=X
                           CHOSEN=CHOSEN+1
                           ACC=ACC+X
                           IF(NCB-ACC.LT.KMIN) GOTO 002
                           IF(NCB.EQ.ACC) GOTO 002
                        ENDDO
 002                    CONTINUE
                        IF(ACC.LT.NCB)THEN
                           NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+(NCB-ACC)
                        ENDIF
                     ENDIF
                     GOTO 333
                  ENDIF
                  ADDITIONNAL_ROWS=NCB-ACC
                  DO i=CHOSEN,1,-1
                     IF(int(dble(ADDITIONNAL_ROWS)/
     &                    dble(i)).NE.0)THEN
                        GOTO 222 
                     ENDIF
                  ENDDO
 222              CONTINUE
                  X=int(dble(ADDITIONNAL_ROWS)/dble(i))
                  DO J=1,i
                     NB_ROWS(J)=NB_ROWS(J)+X
                     ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                  ENDDO 
                  IF(ADDITIONNAL_ROWS.NE.0) THEN
                     NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                  ENDIF    
               ENDIF
 333           CONTINUE
               IF(NB_ROWS(CHOSEN).EQ.0) CHOSEN=CHOSEN-1
               GOTO 889
            ELSE
               DO i=OTHERS,1,-1
                  SOMME=dble(0)
                  DO J=1,i
                     SOMME=SOMME+TEMP(J)
                  ENDDO
                  SOMME=(dble(i)*TEMP(i))-SOMME
                  IF(dble(TOTAL_MEM).GE.SOMME) GOTO 444
               ENDDO
 444           CONTINUE
               REF=i
               DO J=1,i
                  IF(TEMP(J).EQ.TEMP(i)) THEN
                     SMALL_SET=J
                     GOTO 123
                  ENDIF
               ENDDO
 123           CONTINUE
               IF(i.EQ.1)THEN 
                  NB_ROWS(i)=NCB
                  CHOSEN=1
                  GOTO 666
               ENDIF
 323           CONTINUE
               AFFECTED=0
               CHOSEN=0
               ACC=0
               DO i=1,SMALL_SET
                  IF(K50.EQ.0)THEN
                     IF((TEMP(SMALL_SET)-TEMP(i)).GT.DK821)THEN
                        TMP_SUM=DK821
                     ELSE
                        TMP_SUM=TEMP(SMALL_SET)-TEMP(i)
                     ENDIF
                     X=int(TMP_SUM/dble(NFRONT))
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  IF(K50.NE.0)THEN
                     IF((TEMP(SMALL_SET)-TEMP(i)).GT.DK821)THEN
                        TMP_SUM=DK821
                     ELSE
                        TMP_SUM=TEMP(SMALL_SET)-TEMP(i)
                     ENDIF
                      X=int((-dble(NFRONT-NCB+ACC)
     &                  +sqrt(((dble(NFRONT-NCB+ACC)*
     &                  dble(NFRONT-NCB+ACC))+dble(4)*
     &                  (TMP_SUM))))/
     &                  dble(2))
                      IF(X.LT.0)THEN
                        WRITE(*,*)MYID,
     &             ': Internal error 4 in DMUMPS_SET_PARTI_ACTV_MEM'
                        CALL MUMPS_ABORT()
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  NB_ROWS(i)=X
                  ACC=ACC+X
                  CHOSEN=CHOSEN+1
                  IF(NCB-ACC.LT.KMIN) GOTO 888
                  IF(NCB.EQ.ACC) GOTO 888
                  IF(ACC.GT.NCB) THEN
                    WRITE(*,*)MYID,
     &            ': Internal error 5 in DMUMPS_SET_PARTI_ACTV_MEM'
                    CALL MUMPS_ABORT()
                  ENDIF
               ENDDO
 888           CONTINUE
               SOMME=dble(0)
               X=NFRONT-NCB
               IF((ACC.GT.NCB))THEN
                  WRITE(*,*)MYID,
     &           ':Internal error 6 in DMUMPS_SET_PARTI_ACTV_MEM'
                  CALL MUMPS_ABORT()
               ENDIF
               IF((ACC.LT.NCB))THEN
                  IF(K50.NE.0)THEN
                     IF(SMALL_SET.LT.OTHERS)THEN
                       SMALL_SET=REF+1
                       REF=SMALL_SET
                       GOTO 323
                     ELSE
                       NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+NCB-ACC
                       GOTO 666
                     ENDIF
                 ENDIF
                 ADDITIONNAL_ROWS=NCB-ACC
                 i=CHOSEN+1
                 DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                .AND.(i.LE.NUMBER_OF_PROCS))
                    J=1
                    X=int(ADDITIONNAL_ROWS/(i-1))
                    IF((X.EQ.0).AND.(ADDITIONNAL_ROWS.NE.0))THEN
                       DO WHILE ((J.LT.i).AND.(ADDITIONNAL_ROWS.GT.0))
                         NB_ROWS(J)=NB_ROWS(J)+1
                         ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                         J=J+1
                       ENDDO
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          WRITE(*,*)MYID,
     &             ':Internal error 7 in DMUMPS_SET_PARTI_ACTV_MEM'
                          CALL MUMPS_ABORT()
                      ENDIF
                      GOTO 047
                    ENDIF
                    IF((TEMP(1)+dble((NB_ROWS(1)+X)*NFRONT)).LE.
     &                   TEMP(i))THEN
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(J.LT.i))
                          AFFECTED=X
                          IF((AFFECTED+NB_ROWS(J)).GT.
     &                         KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(J)
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                         AFFECTED
                          J=J+1
                       ENDDO
                    ELSE
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(J.LE.i))
                          AFFECTED=int((TEMP(i)-(TEMP(J)+
     &                         (dble(NB_ROWS(J))*dble(NFRONT))))
     &                         /dble(NFRONT))
                          IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(J)
                          ENDIF
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          J=J+1
                       ENDDO
                    ENDIF
                    i=i+1
                 ENDDO
 047             CONTINUE
                 IF((ADDITIONNAL_ROWS.EQ.0).AND.
     &                (i.LT.NUMBER_OF_PROCS))THEN
                    CHOSEN=i-1
                 ELSE
                    CHOSEN=i-2
                 ENDIF
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     &                 (ADDITIONNAL_ROWS.NE.0))THEN
                    DO i=1,CHOSEN
                       NB_ROWS(i)=NB_ROWS(i)+1
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 048
                    ENDDO
 048                CONTINUE
                 ENDIF
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     &                (ADDITIONNAL_ROWS.NE.0))THEN
                    i=CHOSEN+1
                    DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                   .AND.(i.LE.NUMBER_OF_PROCS))
                       J=1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(J.LE.i))
                          AFFECTED=int((TEMP(i)-(TEMP(J)+
     &                         (dble(NB_ROWS(J))*
     &                         dble(NFRONT))))/dble(NFRONT))
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          J=J+1
                       ENDDO
                       i=i+1
                    ENDDO
                    CHOSEN=i-2
                 ENDIF
                 CONTINUE
              ENDIF
 666          CONTINUE
              SOMME=dble(0)
              X=0
              POS=0
              DO i=1,CHOSEN
                 IF(K50.NE.0) THEN
                    IF((TEMP(i)+dble(NB_ROWS(i))
     &                   *dble(X+NB_ROWS(i)+NFRONT-NCB))
     &                   .GT.PEAK)THEN
                       SMALL_SET=SMALL_SET+1
                    ENDIF
                 ENDIF
                 IF(K50.EQ.0) THEN
                    IF((TEMP(i)+dble(NB_ROWS(i))*dble(NFRONT))
     &                   .GT.PEAK)THEN
                       SMALL_SET=SMALL_SET+1
                    ENDIF
                 ENDIF
                 X=X+NB_ROWS(i)
                 SOMME=SOMME+ dble(NB_ROWS(i))
              ENDDO
           ENDIF
 889       CONTINUE
           J=CHOSEN
           X=0
           DO i=J,1,-1
             IF(NB_ROWS(i).EQ.0)THEN
                IF(X.EQ.1)THEN
                  WRITE(*,*)MYID,
     &         ':Internal error 12 in DMUMPS_SET_PARTI_ACTV_MEM'
                  CALL MUMPS_ABORT()
                ENDIF
                CHOSEN=CHOSEN-1
             ELSE 
                  IF(NB_ROWS(i).GT.0)THEN
                    X=1
                 ELSE
                    WRITE(*,*)
     &            'Internal error 13 in DMUMPS_SET_PARTI_ACTV_MEM'
                    CALL MUMPS_ABORT()
                  ENDIF
             ENDIF
          ENDDO
           NSLAVES_NODE=CHOSEN
           TAB_POS(NSLAVES_NODE+1)= NCB+1
           TAB_POS(SLAVEF+2) = CHOSEN
           POS=1
           DO i=1,CHOSEN
              SLAVES_LIST(i)=TEMP_ID(i)
              TAB_POS(i)=POS
              POS=POS+NB_ROWS(i) 
              IF(NB_ROWS(i).LE.0)THEN
                WRITE(*,*)
     &          'Internal error 14 in DMUMPS_SET_PARTI_ACTV_MEM'
                CALL MUMPS_ABORT()
              ENDIF
           ENDDO
           DO i=CHOSEN+1,NUMBER_OF_PROCS
              SLAVES_LIST(i)=TEMP_ID(i)
           ENDDO
           IF(POS.NE.(NCB+1))THEN
              WRITE(*,*)
     &        'Internal error 15 in DMUMPS_SET_PARTI_ACTV_MEM'
              CALL MUMPS_ABORT()
           ENDIF
      END SUBROUTINE DMUMPS_SET_PARTI_ACTV_MEM
      SUBROUTINE DMUMPS_SET_PARTI_FLOP_IRR
     &      (NCBSON_MAX,SLAVEF,KEEP,KEEP8,
     &       PROCS,MEM_DISTRIB,NCB,NFRONT,
     &       NSLAVES_NODE,TAB_POS,
     &       SLAVES_LIST,SIZE_SLAVES_LIST,MYID,INODE,MP,LP)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: NCBSON_MAX
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(in) :: MP,LP
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER NUMBER_OF_PROCS,K47,K48, K50,K83,K69
      INTEGER(8) :: K821
      INTEGER J
      INTEGER KMIN, KMAX
      INTEGER OTHERS,CHOSEN,SMALL_SET,ACC
      DOUBLE PRECISION SOMME,TMP_SUM,DELTA,A,B,C,MASTER_WORK
      INTEGER AFFECTED
      INTEGER ADDITIONNAL_ROWS,i,X,REF,POS,NELIM
      INTEGER(8) X8
      LOGICAL FORCE_CAND,SMP
      DOUBLE PRECISION BANDE_K821
      INTEGER NB_SAT,NB_ZERO
      DOUBLE PRECISION TEMP(SLAVEF),TOTAL_COST, MAX_MEM_ALLOW
      INTEGER TEMP_ID(SLAVEF),NB_ROWS(SLAVEF)
      INTEGER NSLAVES_REF,NCB_FILS
      EXTERNAL MPI_WTIME,MUMPS_GETKMIN
      INTEGER MUMPS_GETKMIN
      INTEGER POS_MIN_LOAD,SIZE_MY_SMP,WHAT,ADDITIONNAL_ROWS_SPECIAL
      LOGICAL HAVE_TYPE1_SON
      DOUBLE PRECISION MIN_LOAD,MAX_LOAD,TEMP_MAX_LOAD
      DOUBLE PRECISION MPI_WTIME
      DOUBLE PRECISION BUF_SIZE,NELIM_MEM_SIZE
      DOUBLE PRECISION MEM_SIZE_STRONG(SLAVEF),MEM_SIZE_WEAK(SLAVEF)
      K821=abs(KEEP8(21))
      TEMP_MAX_LOAD=dble(0)
      K50=KEEP(50)
      K48=KEEP(48)
      K47=KEEP(47)
      K83=KEEP(83)
      K69=0
      NCB_FILS=NCBSON_MAX
      IF(int(NCB_FILS,8)*int(min(NCB,NCB_FILS),8).GT.K821)THEN
         HAVE_TYPE1_SON=.TRUE.
      ELSE
         HAVE_TYPE1_SON=.FALSE.
      ENDIF
      SMP=(K69.NE.0)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
      END IF
      NELIM=NFRONT-NCB
         KMAX=int(K821/int(NCB,8))
         IF(FORCE_CAND)THEN
           DO i=1,PROCS(SLAVEF+1)
              WLOAD(i)=LOAD_FLOPS(PROCS(i))
              IDWLOAD(i)=PROCS(i)
              WLOAD(i)=max(WLOAD(i),0.0d0)
           ENDDO
           NUMBER_OF_PROCS=PROCS(SLAVEF+1)
            OTHERS=NUMBER_OF_PROCS
         ELSE
            NUMBER_OF_PROCS=SLAVEF
            WLOAD(1:SLAVEF) = LOAD_FLOPS(0:NUMBER_OF_PROCS-1)
            DO i=1,NUMBER_OF_PROCS
               IDWLOAD(i) = i - 1
               IF (WLOAD(i) < -0.5d0 ) THEN
                  IF((MP.GT.0).AND.(LP.GE.2))THEN
                     WRITE(MP,*)MYID,': Negative load ',
     &                    WLOAD(i)
                  ENDIF
               ENDIF 
               WLOAD(i)=max(WLOAD(i),0.0d0)
            ENDDO
            OTHERS=NUMBER_OF_PROCS-1
         ENDIF
         KMAX=int(NCB/OTHERS)
         KMIN=MUMPS_GETKMIN(int(NCB,8)*int(KMAX,8),K50,KMAX,NCB)
         NB_ROWS=0
         CALL MUMPS_SORT_DOUBLES(NUMBER_OF_PROCS, WLOAD, IDWLOAD)
         IF(K50.EQ.0)THEN
            TOTAL_COST=dble( NELIM ) * dble ( NCB ) +
     &           dble(NCB) * dble(NELIM)*dble(2*NFRONT-NELIM-1)
         ELSE
            TOTAL_COST=dble(NELIM) * dble ( NCB ) *
     &           dble(NFRONT+1)
         ENDIF
         CALL MUMPS_GET_FLOPS_COST(NFRONT,NELIM,NELIM,K50,
     &        2,MASTER_WORK)
         SOMME=dble(0)
         J=1
         IF(FORCE_CAND.AND.(NUMBER_OF_PROCS.GT.K83))THEN
            MASTER_WORK=dble(KEEP(88))*MASTER_WORK/dble(100)
         ENDIF
         IF(FORCE_CAND.AND.(NUMBER_OF_PROCS.LE.K83))THEN
            MASTER_WORK=dble(KEEP(87))*MASTER_WORK/dble(100)
         ENDIF
         IF(MASTER_WORK.LT.dble(1))THEN
            MASTER_WORK=dble(1)
         ENDIF
         NSLAVES_REF=int(TOTAL_COST/MASTER_WORK)+1
         IF(FORCE_CAND)THEN
            NSLAVES_REF=min(NSLAVES_REF,NUMBER_OF_PROCS)
         ELSE
            NSLAVES_REF=min(NSLAVES_REF,NUMBER_OF_PROCS-1)
         ENDIF
        DO i=1,NUMBER_OF_PROCS
           IF((IDWLOAD(i).NE.MYID))THEN
              TEMP_ID(J)=IDWLOAD(i)
              TEMP(J)=WLOAD(i)
              IF(BDC_M2_FLOPS)THEN
                 TEMP(J)=TEMP(J)+NIV2(TEMP_ID(J)+1)
              ENDIF
              J=J+1
           ENDIF
        ENDDO               
        NUMBER_OF_PROCS=J-1
        CALL MUMPS_SORT_DOUBLES(NUMBER_OF_PROCS, TEMP, TEMP_ID)
        SOMME=dble(0)
        TMP_SUM=dble(0)
        DO i=1,OTHERS
            SOMME=SOMME+TEMP(OTHERS)-TEMP(i)
            TMP_SUM=TMP_SUM+TEMP(i)
        ENDDO
         TMP_SUM=(TMP_SUM/dble(OTHERS))+
     &        (TOTAL_COST/dble(OTHERS))
         SIZE_MY_SMP=OTHERS
         MIN_LOAD=TEMP(1)
         POS_MIN_LOAD=1
         IF(.NOT.SMP) MAX_LOAD=TEMP(OTHERS)
         IF(SMP)THEN
            J=1
            DO i=1,OTHERS
               IF(MEM_DISTRIB(TEMP_ID(i)).EQ.1)THEN
                  IF(TEMP(i).LE.TMP_SUM)THEN
                     WLOAD(J)=TEMP(i)
                     IDWLOAD(J)=TEMP_ID(i)
                     J=J+1
                  ELSE
                  ENDIF
               ENDIF
            ENDDO
            MAX_LOAD=WLOAD(J-1)
            SIZE_MY_SMP=J-1
            DO i=1,OTHERS
               IF((MEM_DISTRIB(TEMP_ID(i)).NE.1).OR.
     &              ((MEM_DISTRIB(TEMP_ID(i)).EQ.1).AND.
     &              (TEMP(i).GE.TMP_SUM)))THEN
                  WLOAD(J)=TEMP(i)
                  IDWLOAD(J)=TEMP_ID(i)
                  J=J+1
               ENDIF
            ENDDO
            TEMP=WLOAD
            TEMP_ID=IDWLOAD
         ENDIF
        IF(BDC_MD)THEN
           BUF_SIZE=dble(K821)
           IF (KEEP(201).EQ.2) THEN
              A=dble(int((dble(KEEP(100))/dble(2))/dble(NELIM)))
              IF(K50.EQ.0)THEN
                 BUF_SIZE=min(BUF_SIZE,A*dble(NCB))
              ELSE
                 BUF_SIZE=min(BUF_SIZE,A*A)
              ENDIF
           ENDIF
           BUF_SIZE=dble(K821)
           DO i=1,NUMBER_OF_PROCS
              A=dble(MD_MEM(TEMP_ID(i)))/
     &             dble(NELIM)
              A=A*dble(NFRONT)
              IF(K50.EQ.0)THEN
                 B=dble(int(dble(NCB)/dble(NUMBER_OF_PROCS))+1)*
     &                dble(NFRONT)
              ELSE
                 WHAT = 5 
                 CALL MUMPS_MAX_SURFCB_NBROWS(WHAT, KEEP,KEEP8, NCB,
     &                NFRONT, min(NCB,OTHERS), J, X8)
                 B=dble(X8)+(dble(J)*dble(NELIM))
              ENDIF
              NELIM_MEM_SIZE=A+B
              MEM_SIZE_WEAK(i)=NELIM_MEM_SIZE
            IF((SBTR_WHICH_M.EQ.0).OR.(.NOT.BDC_SBTR))THEN
               IF(BDC_M2_MEM)THEN
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)
               ELSE
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))
               ENDIF
            ELSE
               IF(BDC_SBTR)THEN
                  IF(BDC_M2_MEM)THEN
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ELSE
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ENDIF
               ENDIF
            ENDIF
            IF(min(MEM_SIZE_STRONG(i),MEM_SIZE_WEAK(i)).LT.dble(0))THEN
                IF(MEM_SIZE_STRONG(i).LT.0.0d0)THEN
                   MEM_SIZE_STRONG(i)=dble(0)
                ELSE
                   MEM_SIZE_WEAK(i)=dble(0)
                ENDIF
             ENDIF
          ENDDO
       ELSE
          BUF_SIZE=dble(K821)
          DO i=1,NUMBER_OF_PROCS
            IF((SBTR_WHICH_M.EQ.0).OR.(.NOT.BDC_SBTR))THEN
               IF(BDC_M2_MEM)THEN
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)
               ELSE
                  MEM_SIZE_STRONG(i)=
     &                 dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                 LU_USAGE(TEMP_ID(i))
               ENDIF
            ELSE
               IF(BDC_SBTR)THEN
                  IF(BDC_M2_MEM)THEN
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-NIV2(TEMP_ID(i)+1)-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ELSE
                     MEM_SIZE_STRONG(i)=
     &                    dble(TAB_MAXS(TEMP_ID(i)))-DM_MEM(TEMP_ID(i))-
     &                    LU_USAGE(TEMP_ID(i))-
     &                    (SBTR_MEM(TEMP_ID(i))-SBTR_CUR(TEMP_ID(i)))
                  ENDIF
               ENDIF
            ENDIF
            MEM_SIZE_STRONG(i)=max(dble(0),MEM_SIZE_STRONG(i))
            MEM_SIZE_WEAK(i)=huge(MEM_SIZE_WEAK(i))
          ENDDO
       ENDIF
       IF((((NUMBER_OF_PROCS.LE.K83).AND.FORCE_CAND).AND.
     &      (TOTAL_COST.GE.SOMME)).OR.
     &      (.NOT.FORCE_CAND).OR.
     &      (((NUMBER_OF_PROCS+1).GT.K83).AND.FORCE_CAND))THEN
               REF=NSLAVES_REF
               SMALL_SET=NSLAVES_REF
               IF(.NOT.SMP)THEN
                  DO i=NSLAVES_REF,1,-1
                     SOMME=dble(0)
                     DO J=1,i
                        SOMME=SOMME+TEMP(J)
                     ENDDO
                     SOMME=(dble(i)*TEMP(i))-SOMME
                     IF(TOTAL_COST.GE.SOMME) GOTO 444
                  ENDDO
 444              CONTINUE
                  REF=i
                  SMALL_SET=REF
                  MAX_LOAD=TEMP(SMALL_SET)
               ELSE
                  X=min(SIZE_MY_SMP,NSLAVES_REF)
 450              CONTINUE
                  SOMME=dble(0)
                  DO J=1,X
                     SOMME=SOMME+(TEMP(X)-TEMP(J))
                  ENDDO
                  IF(SOMME.GT.TOTAL_COST)THEN
                     X=X-1
                     GOTO 450
                  ELSE
                     IF(X.LT.SIZE_MY_SMP) THEN
                        REF=X
                        SMALL_SET=REF
                        MAX_LOAD=TEMP(SMALL_SET)
                     ELSE
                        X=min(SIZE_MY_SMP,NSLAVES_REF)
                        J=X+1
                        MAX_LOAD=TEMP(X)
                        TMP_SUM=MAX_LOAD
                        DO i=X+1,OTHERS
                           IF(TEMP(i).GT.MAX_LOAD)THEN
                              SOMME=SOMME+(dble(i-1)*(TEMP(i)-MAX_LOAD))
                              TMP_SUM=MAX_LOAD
                              MAX_LOAD=TEMP(i)
                           ELSE
                              SOMME=SOMME+(MAX_LOAD-TEMP(i))
                           ENDIF
                           IF(i.EQ.NSLAVES_REF)THEN
                              SMALL_SET=NSLAVES_REF
                              REF=SMALL_SET
                              GOTO 323
                           ENDIF
                           IF(SOMME.GT.TOTAL_COST)THEN
                              REF=i-1
                              SMALL_SET=i-1
                              MAX_LOAD=TMP_SUM
                              GOTO 323
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
 323           CONTINUE
               MAX_LOAD=dble(0)
               DO i=1,SMALL_SET
                  MAX_LOAD=max(MAX_LOAD,TEMP(i))
               ENDDO
               TEMP_MAX_LOAD=MAX_LOAD
               NB_ROWS=0
               TMP_SUM=dble(0)
               CHOSEN=0
               ACC=0
               NB_SAT=0
               NB_ZERO=0
               DO i=1,SMALL_SET
                  IF(K50.EQ.0)THEN
                     X=int(BUF_SIZE/dble(NCB+1))-1
                     BANDE_K821=dble(X)*dble(NFRONT)
                  ELSE
                     A=dble(1)
                     B=dble(ACC+2)
                     C=-BUF_SIZE+dble(ACC+NELIM)
                     DELTA=(B*B)-(dble(4)*A*C)
                     X=int((-B+sqrt(DELTA))/(dble(2)*A))
                     IF(X.GT.NCB-ACC) X=NCB-ACC
                     BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                  ENDIF
                  IF(HAVE_TYPE1_SON)THEN
                     IF(K50.EQ.0)THEN
                        X=int((BUF_SIZE-dble(NFRONT))/dble(NFRONT+1))
                        BANDE_K821=dble(X)*dble(NFRONT)
                     ELSE
                        A=dble(1)
                        B=dble(ACC+2+NELIM)
                        C=-BUF_SIZE+dble(ACC+NELIM)
                        DELTA=(B*B)-(dble(4)*A*C)
                        X=int((-B+sqrt(DELTA))/(dble(2)*A))
                        IF(X.GT.NCB-ACC) X=NCB-ACC
                        BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                     ENDIF
                  ENDIF
                  MAX_MEM_ALLOW=BANDE_K821
                  IF(BDC_MD)THEN
                     MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                    BANDE_K821)
                     MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                  ENDIF
                  IF(K50.EQ.0)THEN
                     KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                     X=int((MAX_LOAD-TEMP(i))/
     &                    (dble(NELIM)*dble(2*NFRONT-NELIM)))
                     IF(X.GE.KMAX)THEN
                        IF(KMAX.GE.KMIN)THEN
                           X=KMAX
                           NB_SAT=NB_SAT+1
                        ELSE
                           X=0
                        ENDIF
                     ELSE
                        IF(X.LT.KMIN)THEN
                           X=0
                        ENDIF                        
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  IF(K50.NE.0)THEN
                        A=dble(1)
                        B=dble(ACC+NELIM)
                        C=dble(-MAX_MEM_ALLOW)
                        DELTA=((B*B)-(dble(4)*A*C))
                        KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                     A=dble(NELIM)
                     B=dble(NELIM)*(dble(NELIM)+dble(2*ACC+1))
                     C=-(MAX_LOAD-TEMP(i))
                     DELTA=(B*B-(dble(4)*A*C))
                     X=int((-B+sqrt(DELTA))/(dble(2)*A))
                     IF(X.LT.0) THEN
                        WRITE(*,*)MYID,
     &    ': Internal error 1 in DMUMPS_SET_PARTI_FLOP_IRR'
                        CALL MUMPS_ABORT()
                     ENDIF
                     IF(X.GE.KMAX)THEN
                        IF(KMAX.GE.KMIN)THEN
                           X=KMAX
                           NB_SAT=NB_SAT+1
                        ELSE
                           X=0
                        ENDIF
                     ELSE
                        IF(X.LT.KMIN)THEN
                           X=0
                        ENDIF
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF                  
                  NB_ROWS(i)=X
                  ACC=ACC+X
                  CHOSEN=CHOSEN+1
                  IF(SMP)THEN
                     IF(MIN_LOAD.GT.TEMP(i))THEN
                        MIN_LOAD=TEMP(i)
                        POS_MIN_LOAD=i
                     ENDIF
                  ENDIF
                  TMP_SUM=MAX_LOAD
                  IF(K50.EQ.0)THEN
                     MAX_LOAD=max(MAX_LOAD,
     &                    (TEMP(i)+(dble(NELIM) *
     &                    dble(NB_ROWS(i)))+
     &                    (dble(NB_ROWS(i))*dble(NELIM)*
     &                    dble(2*NFRONT-NELIM-1))))
                  ELSE
                     MAX_LOAD=max(MAX_LOAD,
     &               TEMP(i)+(dble(NELIM) * dble(NB_ROWS(i)))*
     &                    dble(2*(NELIM+ACC)-NB_ROWS(i)
     &                    -NELIM+1))
                  ENDIF
                  IF(TMP_SUM.LT.MAX_LOAD)THEN
                  ENDIF
                  IF(NCB-ACC.LT.KMIN) GOTO 888
                  IF(NCB.EQ.ACC) GOTO 888
                  IF(ACC.GT.NCB) THEN
                    WRITE(*,*)MYID,
     &      ': Internal error 2 in DMUMPS_SET_PARTI_FLOP_IRR'
                    CALL MUMPS_ABORT()
                  ENDIF
               ENDDO
 888           CONTINUE
               SOMME=dble(0)
               X=NFRONT-NCB
               IF((ACC.GT.NCB))THEN
                  WRITE(*,*)MYID,
     &          ': Internal error 3 in DMUMPS_SET_PARTI_FLOP_IRR'
                  CALL MUMPS_ABORT()
               ENDIF
               IF((ACC.LT.NCB))THEN
                  IF(K50.NE.0)THEN
                     IF(SMALL_SET.LE.OTHERS)THEN
                       IF((NB_SAT.EQ.SMALL_SET).AND.(SMALL_SET.LT.
     &                      NSLAVES_REF))THEN
                          SMALL_SET=REF+1
                          REF=REF+1
                          NB_ROWS=0
                          GOTO 323
                       ENDIF
                       ADDITIONNAL_ROWS_SPECIAL=NCB-ACC
                       DO i=1,SMALL_SET
                          MAX_LOAD=TEMP_MAX_LOAD
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=dble(NELIM)*
     &                         dble(ADDITIONNAL_ROWS)*
     &                         dble(2*NFRONT-ADDITIONNAL_ROWS-NELIM
     &                         +1)
                          SOMME=SOMME/dble(SMALL_SET-NB_SAT)
                          NB_ROWS=0
                          NB_ZERO=0
                          ACC=0
                          CHOSEN=0
                          NB_SAT=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO J=1,SMALL_SET
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                        min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=(dble(NELIM)*dble(NELIM+2*ACC+1))
                             C=-(MAX_LOAD-TEMP(J)+SOMME)
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             X=X+1
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 4 in DMUMPS_SET_PARTI_FLOP_IRR'
                                CALL MUMPS_ABORT()
                             ENDIF
                             IF(X.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   X=KMAX
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   NB_ZERO=NB_ZERO+1
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.min(KMIN,KMAX))THEN
                                   NB_ZERO=NB_ZERO+1
                                   X=0
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(J)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(J))THEN
                                   MIN_LOAD=TEMP(J)
                                   POS_MIN_LOAD=i
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             TMP_SUM=MAX_LOAD
                             TEMP_MAX_LOAD=max(TEMP_MAX_LOAD,
     &                            TEMP(J)+(dble(NELIM) *
     &                            dble(NB_ROWS(J)))*
     &                            dble(2*(NELIM+
     &                            ACC)-NB_ROWS(J)
     &                            -NELIM+1))
                             IF(REF.LE.NUMBER_OF_PROCS-1)THEN
                                IF(TEMP_MAX_LOAD.GT.TEMP(REF+1))THEN
                                   IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                      SMALL_SET=REF+1
                                      REF=REF+1
                                      NB_ROWS=0
                                      GOTO 323
                                   ENDIF
                                ENDIF
                             ENDIF
                             IF(NCB.EQ.ACC) GOTO 666
                          ENDDO
                          IF(NB_SAT.EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                          IF(NB_ZERO.EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                          IF((NB_SAT+NB_ZERO).EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                       ENDDO
 434                   CONTINUE
                       ADDITIONNAL_ROWS=NCB-ACC
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          IF(ADDITIONNAL_ROWS.LT.KMIN)THEN
                             i=CHOSEN
                             J=ACC
 436                         CONTINUE
                             IF(NB_ROWS(i).NE.0)THEN
                                J=J-NB_ROWS(i)
                                A=dble(1)
                                B=dble(J+2)
                                C=-BUF_SIZE+dble(J+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-J) X=NCB-J
                                BANDE_K821=dble(X)*dble(NELIM+J+X)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=dble(1)
                                   B=dble(J+2+NELIM)
                                   C=-BUF_SIZE+dble(J+NELIM)
                                   DELTA=(B*B)-(dble(4)*A*C)
                                   X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                   IF(X.GT.NCB-J) X=NCB-J
                                   BANDE_K821=dble(X)*dble(NELIM+J+X)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K821
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=min(
     &                         min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                                  BANDE_K821)
                                   MAX_MEM_ALLOW=max(dble(0),
     &                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=dble(1)
                                B=dble(J+NELIM)
                                C=dble(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(dble(4)*A*C))
                                KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(NB_ROWS(i).NE.KMAX)THEN
                                   IF(NCB-J.LE.KMAX)THEN
                                      NB_ROWS(i)=+NCB-J
                                      ADDITIONNAL_ROWS=0
                                   ENDIF
                                ENDIF
                                TEMP_MAX_LOAD=max(TEMP_MAX_LOAD,
     &                               TEMP(i)+
     &                               (dble(NELIM) * dble(NB_ROWS(i)))*
     &                               dble(2*(NELIM+
     &                               ACC)-NB_ROWS(i)
     &                               -NELIM+1))
                                IF(REF.LE.NUMBER_OF_PROCS-1)THEN
                                   IF(TEMP_MAX_LOAD.GT.TEMP(REF+1))THEN
                                      IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                         SMALL_SET=REF+1
                                         REF=REF+1
                                         NB_ROWS=0
                                         GOTO 323
                                      ENDIF
                                   ENDIF
                                ENDIF
                             ELSE
                                i=i-1
                                IF(i.NE.0)GOTO 436
                             ENDIF
                             IF(ADDITIONNAL_ROWS.NE.0)THEN
                                i=CHOSEN
                                IF(i.NE.SMALL_SET)THEN
                                   i=i+1
                                   IF(NB_ROWS(i).NE.0)THEN
                                      WRITE(*,*)MYID,
     &    ': Internal error 5 in DMUMPS_SET_PARTI_FLOP_IRR'
                                      CALL MUMPS_ABORT()
                                   ENDIF
                                ENDIF
                                NB_ROWS(i)=NB_ROWS(i)+ADDITIONNAL_ROWS
                                ADDITIONNAL_ROWS=0
                             ENDIF
                             CHOSEN=i
                          ENDIF
                       ENDIF
                       i=CHOSEN+1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                      .AND.(i.LE.NUMBER_OF_PROCS))
                          IF((TEMP(i).LE.MAX_LOAD))THEN
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=dble(NELIM)*dble(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(i))
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   X=KMAX
                                ELSE
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.KMIN)THEN
                                   X=0
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(i)=X
                             ACC=ACC+X
                             ADDITIONNAL_ROWS=NCB-ACC
                          ELSE IF((TEMP(i).GT.MAX_LOAD))THEN
                             MAX_LOAD=TEMP(i)
                             NB_SAT=0
                             ACC=0
                             NB_ROWS=0
                             DO J=1,i
                                A=dble(1)
                                B=dble(ACC+2)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=dble(1)
                                   B=dble(ACC+2+NELIM)
                                   C=-BUF_SIZE+dble(ACC+NELIM)
                                   DELTA=(B*B)-(dble(4)*A*C)
                                   X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                   IF(X.GT.NCB-ACC) X=NCB-ACC
                                   BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K821
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                                  BANDE_K821)
                                   MAX_MEM_ALLOW=max(dble(0),
     &                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=dble(1)
                                B=dble(ACC+NELIM)
                                C=dble(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(dble(4)*A*C))
                                KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                                A=dble(NELIM)
                                B=dble(NELIM)*dble(NELIM+2*ACC+1)
                                C=-(MAX_LOAD-TEMP(J))
                                DELTA=(B*B-(dble(4)*A*C))
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.LT.0) THEN
                                   WRITE(*,*)MYID,
     &    ': Internal error 6 in DMUMPS_SET_PARTI_FLOP_IRR'
                                   CALL MUMPS_ABORT()
                                ENDIF
                                IF(X.GE.KMAX)THEN
                                   IF(KMAX.GE.KMIN)THEN
                                      X=KMAX
                                      NB_SAT=NB_SAT+1
                                   ELSE
                                      X=0
                                   ENDIF
                                ELSE
                                   IF(X.LT.min(KMIN,KMAX))THEN
                                      X=0
                                   ENDIF
                                ENDIF
                                IF((ACC+X).GT.NCB) X=NCB-ACC
                                NB_ROWS(J)=X
                                IF(SMP)THEN
                                   IF(MIN_LOAD.GT.TEMP(J))THEN
                                      MIN_LOAD=TEMP(J)
                                      POS_MIN_LOAD=i
                                   ENDIF
                                ENDIF
                                ACC=ACC+X
                                MAX_LOAD=max(MAX_LOAD,
     &                               TEMP(J)+
     &                               (dble(NELIM)*dble(NB_ROWS(J)))*
     &                               dble(2*(NELIM+
     &                               ACC)-NB_ROWS(J)
     &                               -NELIM+1))
                                IF(NCB.EQ.ACC) GOTO 741
                                IF(NCB-ACC.LT.KMIN) GOTO 210
                             ENDDO
 210                         CONTINUE
                          ENDIF
 741                      CONTINUE
                          i=i+1
                          ADDITIONNAL_ROWS=NCB-ACC
                       ENDDO
                       CHOSEN=i-1
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=dble(NELIM)*dble(ADDITIONNAL_ROWS)*
     &                         dble(2*NFRONT-ADDITIONNAL_ROWS-
     &                         NELIM+1)
                          SOMME=SOMME/dble(NUMBER_OF_PROCS)
                          NB_ROWS=0
                          ACC=0
                          CHOSEN=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO i=1,OTHERS
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=dble(NELIM)*dble(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(i)+SOMME)
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 7 in DMUMPS_SET_PARTI_FLOP_IRR'
                                CALL MUMPS_ABORT()
                             ENDIF
                             IF(X.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   X=KMAX
                                ELSE
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.min(KMIN,KMAX))THEN
                                   X=min(KMAX,KMIN)
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(i)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(i))THEN
                                   MIN_LOAD=TEMP(i)
                                   POS_MIN_LOAD=i
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             IF(NCB.EQ.ACC) GOTO 666
                             IF(NCB-ACC.LT.KMIN) GOTO 488
                          ENDDO
 488                      CONTINUE
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=dble(NELIM)*
     &                         dble(ADDITIONNAL_ROWS)*
     &                         dble(2*NFRONT-ADDITIONNAL_ROWS-
     &                         NELIM+1)
                          SOMME=SOMME/dble(NUMBER_OF_PROCS)
                          NB_ROWS=0
                          ACC=0
                          CHOSEN=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO i=1,OTHERS
                             A=dble(1)
                             B=dble(ACC+2)
                             C=-BUF_SIZE+dble(ACC+NELIM)
                             DELTA=(B*B)-(dble(4)*A*C)
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             IF(HAVE_TYPE1_SON)THEN
                                A=dble(1)
                                B=dble(ACC+2+NELIM)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(BANDE_K821,
     &                               MEM_SIZE_STRONG(i))
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             A=dble(1)
                             B=dble(ACC+NELIM)
                             C=dble(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(dble(4)*A*C))
                             KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                             A=dble(NELIM)
                             B=dble(NELIM)*dble(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(i)+SOMME)
                             DELTA=(B*B-(dble(4)*A*C))
                             X=int((-B+sqrt(DELTA))/(dble(2)*A))
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 8 in DMUMPS_SET_PARTI_FLOP_IRR'
                                CALL MUMPS_ABORT()
                             ENDIF
                             IF(X.GE.KMAX)THEN
                                X=KMAX
                             ELSE
                                IF(X.LT.KMIN)THEN
                                   X=KMIN
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(i)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(i))THEN
                                   MIN_LOAD=TEMP(i)
                                   POS_MIN_LOAD=i
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             IF(NCB.EQ.ACC) GOTO 666
                             IF(NCB-ACC.LT.KMIN) GOTO 477
                          ENDDO
 477                      CONTINUE
                          IF(ACC.NE.NCB)THEN
                             NB_SAT=0
                             ACC=0
                             CHOSEN=0
                             IF(SMP)THEN
                                MIN_LOAD=TEMP(1)
                                POS_MIN_LOAD=1
                             ENDIF
                             DO i=1,OTHERS
                                A=dble(1)
                                B=dble(ACC+2)
                                C=-BUF_SIZE+dble(ACC+NELIM)
                                DELTA=(B*B)-(dble(4)*A*C)
                                X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=dble(1)
                                   B=dble(ACC+2+NELIM)
                                   C=-BUF_SIZE+dble(ACC+NELIM)
                                   DELTA=(B*B)-(dble(4)*A*C)
                                   X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                   IF(X.GT.NCB-ACC) X=NCB-ACC
                                   BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K821
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=min(BANDE_K821,
     &                                  MEM_SIZE_STRONG(i))
                                   MAX_MEM_ALLOW=max(dble(0),
     &                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=dble(1)
                                B=dble(ACC+NELIM)
                                C=dble(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(dble(4)*A*C))
                                KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                                X=KMAX-NB_ROWS(i)
                                IF((ACC+NB_ROWS(i)+X).GT.NCB)
     &                               X=NCB-(ACC+NB_ROWS(i))
                                NB_ROWS(i)=NB_ROWS(i)+X
                                IF((dble(NB_ROWS(i))*
     &                               dble(NB_ROWS(i)+ACC)).EQ.
     &                               BANDE_K821)THEN
                                   NB_SAT=NB_SAT+1
                                ENDIF
                                ACC=ACC+NB_ROWS(i)
                                IF(SMP)THEN
                                   IF(MIN_LOAD.GT.TEMP(i))THEN
                                      MIN_LOAD=TEMP(i)
                                      POS_MIN_LOAD=i
                                   ENDIF
                                ENDIF
                                CHOSEN=CHOSEN+1
                                IF(NCB.EQ.ACC) GOTO 666
                                IF(NCB-ACC.LT.KMIN) GOTO 834
                             ENDDO
 834                         CONTINUE
                          ENDIF
                          IF(ACC.NE.NCB)THEN
                            ADDITIONNAL_ROWS=NCB-ACC
                            SOMME=dble(NELIM)*
     &                           dble(ADDITIONNAL_ROWS)*
     &                           dble(2*NFRONT-ADDITIONNAL_ROWS-
     &                           NELIM+1)
                            SOMME=SOMME/dble(NUMBER_OF_PROCS-NB_SAT)
                            ACC=0
                            DO i=1,CHOSEN
                               A=dble(1)
                               B=dble(ACC+2)
                               C=-BUF_SIZE+dble(ACC+NELIM)
                               DELTA=(B*B)-(dble(4)*A*C)
                               X=int((-B+sqrt(DELTA))/(dble(2)*A))
                               IF(X.GT.NCB-ACC) X=NCB-ACC
                               BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                               IF(HAVE_TYPE1_SON)THEN
                                  A=dble(1)
                                  B=dble(ACC+2+NELIM)
                                  C=-BUF_SIZE+dble(ACC+NELIM)
                                  DELTA=(B*B)-(dble(4)*A*C)
                                  X=int((-B+sqrt(DELTA))/(dble(2)*A))
                                  IF(X.GT.NCB-ACC) X=NCB-ACC
                                  BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                               ENDIF
                               IF((dble(NB_ROWS(i))*
     &                              dble(NB_ROWS(i)+ACC)).EQ.
     &                              BANDE_K821)THEN
                                  GOTO 102
                               ENDIF
                               A=dble(NELIM)
                               B=dble(NELIM)*
     &                              dble(NELIM+2*(ACC+NB_ROWS(i))+1)
                               C=-(SOMME)
                               DELTA=(B*B-(dble(4)*A*C))
                               X=int((-B+sqrt(DELTA))/(dble(2)*A))
                               A=dble(1)
                               B=dble(ACC+NELIM)
                               C=dble(-BANDE_K821)
                               DELTA=((B*B)-(dble(4)*A*C))
                               KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                               IF(X.LT.0) THEN
                                  WRITE(*,*)MYID,
     &    ': Internal error 9 in DMUMPS_SET_PARTI_FLOP_IRR'
                                  CALL MUMPS_ABORT()
                               ENDIF
                               IF((ACC+X+NB_ROWS(i)).GT.NCB)THEN
                                  IF((NCB-ACC).GT.KMAX)THEN
                                     NB_ROWS(i)=KMAX
                                  ELSE
                                     NB_ROWS(i)=NCB-ACC
                                  ENDIF
                               ELSE
                                  IF((NB_ROWS(i)+X).GT.KMAX)THEN
                                     NB_ROWS(i)=KMAX
                                  ELSE
                                     NB_ROWS(i)=NB_ROWS(i)+X
                                  ENDIF
                               ENDIF
 102                           CONTINUE
                               ACC=ACC+NB_ROWS(i)
                               IF(NCB.EQ.ACC) THEN
                                  CHOSEN=i
                                  GOTO 666
                               ENDIF
                               IF(NCB-ACC.LT.KMIN) THEN
                                  CHOSEN=i
                                  GOTO 007
                               ENDIF
                            ENDDO
 007                        CONTINUE
                            DO i=1,CHOSEN
                               NB_ROWS(i)=NB_ROWS(i)+1
                               ACC=ACC+1
                               IF(ACC.EQ.NCB)GOTO 666
                            ENDDO
                            IF(ACC.LT.NCB)THEN
                               IF(SMP)THEN
                                  NB_ROWS(1)=NB_ROWS(1)+NCB-ACC
                               ELSE
                                  NB_ROWS(POS_MIN_LOAD)=
     &                                 NB_ROWS(POS_MIN_LOAD)+NCB-ACC
                               ENDIF
                            ENDIF
                         ENDIF
                         GOTO 666
                     ENDIF
                  ENDIF
                  GOTO 666
                 ENDIF
                 ADDITIONNAL_ROWS=NCB-ACC
                 i=CHOSEN+1
                 IF(NB_SAT.EQ.SMALL_SET) GOTO 777
                 DO i=1,SMALL_SET
                    IDWLOAD(i)=i
                    AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                    BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                    IF(HAVE_TYPE1_SON)THEN
                       AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                      (dble(NFRONT+1)))
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                    ENDIF
                    MAX_MEM_ALLOW=BANDE_K821
                    IF(BDC_MD)THEN
                       MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                      BANDE_K821)
                       MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                    ENDIF
                    WLOAD(i)=MAX_MEM_ALLOW
                 ENDDO
                 CALL MUMPS_SORT_DOUBLES(SMALL_SET, WLOAD, IDWLOAD)
                 NB_ZERO=0
                 IF((NB_SAT.EQ.SMALL_SET).AND.
     &                (SMALL_SET.LT.NSLAVES_REF))THEN
                    SMALL_SET=REF+1
                    REF=REF+1
                    NB_ROWS=0
                    GOTO 323
                 ENDIF
                 IF((NB_SAT.EQ.SMALL_SET).AND.
     &                (SMALL_SET.LE.NUMBER_OF_PROCS))GOTO 777
                 AFFECTED=int(ADDITIONNAL_ROWS/(SMALL_SET-NB_SAT))
                 AFFECTED=max(AFFECTED,1)
                 DO i=1,SMALL_SET
                    KMAX=int(WLOAD(i)/dble(NFRONT))
                    IF(NB_ROWS(IDWLOAD(i)).EQ.KMAX)THEN 
                       GOTO 912
                    ENDIF
                    IF((NB_ROWS(IDWLOAD(i))+min(AFFECTED,
     &                   ADDITIONNAL_ROWS)).GT.KMAX)THEN
                       IF(NB_ROWS(IDWLOAD(i)).GT.KMAX)THEN
                       ENDIF
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                      (KMAX-NB_ROWS(IDWLOAD(i)))
                       NB_ROWS(IDWLOAD(i))=KMAX
                       NB_SAT=NB_SAT+1
                       IF(NB_SAT.EQ.SMALL_SET)THEN
                          IF(SMALL_SET.NE.NSLAVES_REF)THEN
                             SMALL_SET=REF+1
                             REF=REF+1
                             NB_ROWS=0
                             GOTO 323
                          ELSE
                             MAX_LOAD=max(MAX_LOAD,
     &                            (TEMP(IDWLOAD(i))+(dble(NELIM) *
     &                            dble(NB_ROWS(IDWLOAD(i))))+
     &                            (dble(NB_ROWS(IDWLOAD(i)))*
     &                            dble(NELIM))*
     &                            dble(2*NFRONT-NELIM-1)))
                             GOTO 777
                          ENDIF
                       ENDIF
                       AFFECTED=int(ADDITIONNAL_ROWS/(SMALL_SET-NB_SAT))
                       AFFECTED=max(AFFECTED,1)
                    ELSE
                       IF((NB_ROWS(IDWLOAD(i))+min(AFFECTED,
     &                      ADDITIONNAL_ROWS)).GE.KMIN)THEN
                          X=min(AFFECTED,ADDITIONNAL_ROWS)
                          NB_ROWS(IDWLOAD(i))=NB_ROWS(IDWLOAD(i))+
     &                         X
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                       ELSE
                          X=int((MAX_LOAD-TEMP(IDWLOAD(i)))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          IF(X+AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ELSE
                             X=AFFECTED+X
                          ENDIF
                          IF(X.GE.KMIN)THEN
                             NB_ROWS(IDWLOAD(i))=NB_ROWS(IDWLOAD(i))+
     &                            X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            X
                          ELSE
                             NB_ZERO=NB_ZERO+1
                          ENDIF
                       ENDIF
                    ENDIF
 912                CONTINUE
                    MAX_LOAD=max(MAX_LOAD,
     &                   (TEMP(IDWLOAD(i))+(dble(NELIM)*
     &                   dble(NB_ROWS(IDWLOAD(i))))+
     &                   (dble(NB_ROWS(IDWLOAD(i)))*dble(NELIM))*
     &                   dble(2*NFRONT-NELIM-1)))
                    IF(SMALL_SET.LT.NUMBER_OF_PROCS)THEN
                       IF(MAX_LOAD.GT.TEMP(SMALL_SET+1))THEN
                          IF(SMALL_SET.LT.NSLAVES_REF)THEN
                             SMALL_SET=REF+1
                             REF=REF+1
                             NB_ROWS=0
                             GOTO 323
                          ENDIF
                       ENDIF
                    ENDIF
                    IF(SMALL_SET.EQ.NB_SAT)GOTO 777
                    IF(ADDITIONNAL_ROWS.EQ.0)THEN
                       CHOSEN=SMALL_SET
                       GOTO 049
                    ENDIF
                 ENDDO
 777             CONTINUE
                 IF((NB_ZERO.NE.0).AND.(ADDITIONNAL_ROWS.GE.KMIN))THEN
                    J=NB_ZERO
 732                CONTINUE
                    X=int(ADDITIONNAL_ROWS/(J))
                    IF(X.LT.KMIN)THEN
                       J=J-1
                       GOTO 732
                    ENDIF
                    IF(X*J.LT.ADDITIONNAL_ROWS)THEN
                       X=X+1
                    ENDIF
                    DO i=1,SMALL_SET
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                         dble(BANDE_K821))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF(NB_ROWS(i).EQ.0)THEN
                          IF(X.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ENDIF
                          IF(X.GT.KMAX)THEN
                             X=KMAX
                          ENDIF
                          IF(X.GT.KMIN)THEN
                             NB_ROWS(i)=X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                             MAX_LOAD=max(MAX_LOAD,
     &                            (TEMP(i)+(dble(NELIM) *
     &                            dble(NB_ROWS(i)))+
     &                            (dble(NB_ROWS(i))*dble(NELIM))*
     &                            dble(2*NFRONT-NELIM-1)))
                          ENDIF
                       ENDIF
                    ENDDO
                 ENDIF
                 i=CHOSEN+1
                 DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                .AND.(i.LE.NUMBER_OF_PROCS))
                    IF((TEMP(i).LE.MAX_LOAD))THEN
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                         BANDE_K821)
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       AFFECTED=int((MAX_LOAD-TEMP(i))/
     &                      (dble(NELIM)*dble(2*NFRONT-NELIM)))
                       IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                          AFFECTED=ADDITIONNAL_ROWS
                       ENDIF
                       IF(NB_ROWS(i).LT.KMAX)THEN
                          IF((AFFECTED+NB_ROWS(i)).GT.KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(i)
                             NB_SAT=NB_SAT+1
                          ELSE
                             IF((AFFECTED+NB_ROWS(i)).LT.
     &                            KMIN)THEN
                                AFFECTED=0
                             ENDIF
                          ENDIF
                          NB_ROWS(i)=NB_ROWS(i)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                       ENDIF
                    ELSE IF((TEMP(i).GT.MAX_LOAD))THEN
                       IF(NB_SAT.EQ.i-1) GOTO 218
                       X=(ADDITIONNAL_ROWS/(i-1-NB_SAT))
                       ACC=1
                       DO J=1,i-1
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(J)+X))
     &                         +(dble(NB_ROWS(J)+X)*dble(NELIM))*
     &                         dble(2*NFRONT-NELIM-1))
                          IF((TEMP(J)+TMP_SUM).GT.MAX_LOAD)THEN
                             ACC=0
                          ENDIF
                       ENDDO
                       IF(ACC.EQ.1)THEN
                          MAX_LOAD=TEMP(i)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                         .AND.(J.LT.i))
                             AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             IF(HAVE_TYPE1_SON)THEN
                                AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             ENDIF
                             AFFECTED=X
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                                AFFECTED=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                                   AFFECTED=KMAX-NB_ROWS(J)
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   IF((AFFECTED+NB_ROWS(J)).LT.
     &                                  KMIN)THEN
                                      AFFECTED=0
                                   ENDIF
                                ENDIF
                                NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               AFFECTED
                             ENDIF
                             J=J+1
                          ENDDO                          
                       ELSE
                          MAX_LOAD=TEMP(i)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                         .AND.(J.LT.i))
                             AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             IF(HAVE_TYPE1_SON)THEN
                                AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                             ENDIF
                             TMP_SUM=((dble(NELIM)* dble(NB_ROWS(J)))
     &                            +(dble(NB_ROWS(J))*dble(NELIM))*
     &                            dble(2*NFRONT-NELIM-1))
                             X=int((MAX_LOAD-(TEMP(J)+TMP_SUM))/
     &                            (dble(NELIM)*dble(2*NFRONT-NELIM)))
                             IF(X.LT.0)THEN
                                WRITE(*,*)MYID,
     &    ': Internal error 10 in DMUMPS_SET_PARTI_FLOP_IRR'
                                CALL MUMPS_ABORT()
                             ENDIF
                             AFFECTED=X
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                                AFFECTED=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                                   AFFECTED=KMAX-NB_ROWS(J)
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   IF((AFFECTED+NB_ROWS(J)).LT.
     &                                  KMIN)THEN
                                      AFFECTED=0
                                   ENDIF
                                ENDIF
                                NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               AFFECTED
                             ENDIF
                             J=J+1
                          ENDDO
                       ENDIF
                    ENDIF
 218                CONTINUE
                    i=i+1
                 ENDDO
                 CHOSEN=i-1
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     &                 (ADDITIONNAL_ROWS.NE.0))THEN
                    DO i=1,CHOSEN
                       IF(NB_ROWS(i)+1.GE.KMIN)THEN
                          NB_ROWS(i)=NB_ROWS(i)+1
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                       ENDIF
                       MAX_LOAD=max(MAX_LOAD,
     &                      (TEMP(i)+(dble(NELIM) *
     &                      dble(NB_ROWS(i)))+
     &                      (dble(NB_ROWS(i))*dble(NELIM))*
     &                      dble(2*NFRONT-NELIM-1)))
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 048
                    ENDDO
 048                CONTINUE
                 ENDIF
                 IF((ADDITIONNAL_ROWS.NE.0))THEN
                    IF(CHOSEN.LT.NUMBER_OF_PROCS)THEN
                       i=CHOSEN+1
                    ELSE
                       IF(CHOSEN.NE.NUMBER_OF_PROCS)THEN
                          WRITE(*,*)MYID,
     &    ': Internal error 11 in DMUMPS_SET_PARTI_FLOP_IRR'
                          CALL MUMPS_ABORT()
                       ENDIF
                       i=CHOSEN
                    ENDIF
                    DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                   .AND.(i.LE.NUMBER_OF_PROCS))
                       IF(TEMP(i).LE.MAX_LOAD)THEN
                          AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                            BANDE_K821)
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(i)))
     &                         +(dble(NB_ROWS(i))*dble(NELIM))*
     &                         dble(2*NFRONT-NELIM-1))
                          X=int((MAX_LOAD-(TEMP(i)+TMP_SUM))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          AFFECTED=X
                          IF(X.LT.0)THEN
                             WRITE(*,*)MYID,
     &    ': Internal error 12 in DMUMPS_SET_PARTI_FLOP_IRR'
                             CALL MUMPS_ABORT()
                          ENDIF
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          IF(NB_ROWS(i).LT.KMAX)THEN
                             IF((AFFECTED+NB_ROWS(i)).GT.KMAX)THEN
                                AFFECTED=KMAX-NB_ROWS(i)
                             ELSE
                                IF((AFFECTED+NB_ROWS(i)).LT.
     &                               KMIN)THEN
                                   AFFECTED=0
                                ENDIF
                             ENDIF
                             NB_ROWS(i)=NB_ROWS(i)+AFFECTED
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          ENDIF
                          IF(i.NE.NUMBER_OF_PROCS) GOTO 624
                       ELSE IF((TEMP(i).GT.MAX_LOAD))THEN
                          X=int(ADDITIONNAL_ROWS/i-1)
                          X=max(X,1)
                          IF((MAX_LOAD+((dble(NELIM)*
     &                         dble(X))+(dble(
     &                         X)*dble(NELIM))*dble(
     &                         (2*NFRONT-NELIM-1)))).LE.TEMP(i))THEN
                             AFFECTED=X
                             POS=1
                          ELSE
                             POS=0
                          ENDIF
                          MAX_LOAD=TEMP(i)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     &                         .AND.(J.LT.i))
                             X=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(X)*dble(NFRONT)
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(HAVE_TYPE1_SON)THEN
                                X=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(X)*dble(NFRONT)
                             ENDIF
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(POS.EQ.0)THEN
                                TMP_SUM=((dble(NELIM) *
     &                               dble(NB_ROWS(J)))
     &                               +(dble(NB_ROWS(J))*dble(NELIM))*
     &                               dble(2*NFRONT-NELIM-1))
                                X=int((TEMP(i)-(TEMP(J)+TMP_SUM))/
     &                               (dble(NELIM)*dble(2*NFRONT-
     &                               NELIM)))
                             ELSE
                                X=int(TMP_SUM)
                             ENDIF
                             IF(X.GT.ADDITIONNAL_ROWS)THEN
                                X=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((X+NB_ROWS(J)).GT.KMAX)THEN
                                   X=KMAX-NB_ROWS(J)
                                ELSE
                                   IF((NB_ROWS(J)+X).LT.
     &                                  KMIN)THEN
                                     X=0
                                  ENDIF
                               ENDIF
                               NB_ROWS(J)=NB_ROWS(J)+X
                               ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                            ENDIF
                            J=J+1
                         ENDDO
                       ENDIF
 624                   CONTINUE
                       i=i+1
                    ENDDO
                    CHOSEN=i-1
                    IF(ADDITIONNAL_ROWS.NE.0)THEN
                       ACC=0
                       DO i=1,CHOSEN
                          X=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(X)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             X=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(X)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                            BANDE_K821)
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(i)))
     &                    +(dble(NB_ROWS(i))*dble(NELIM))*
     &                    dble(2*NFRONT-NELIM-1))
                          X=int((MAX_LOAD-
     &                         (TEMP(i)+TMP_SUM))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          IF(X.LT.0)THEN
                             WRITE(*,*)MYID,
     &    ': Internal error 13 in DMUMPS_SET_PARTI_FLOP_IRR'
                             CALL MUMPS_ABORT()
                          ENDIF
                          IF(X.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ENDIF
                          IF(NB_ROWS(i).LT.KMAX)THEN
                             IF((X+NB_ROWS(i)).GE.KMAX)THEN
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               (KMAX-NB_ROWS(i))
                                NB_ROWS(i)=KMAX
                             ELSE
                                IF((X+NB_ROWS(i)).GE.
     &                               KMIN)THEN
                                   NB_ROWS(i)=NB_ROWS(i)+X
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                                   ACC=ACC+1
                                ELSE
                                   ACC=ACC+1
                                ENDIF
                             ENDIF
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       IF(CHOSEN.LT.NUMBER_OF_PROCS)THEN
                          CHOSEN=CHOSEN+1
                       ENDIF
                       IF(ACC.EQ.0)THEN
                          ACC=1
                       ENDIF
                       X=int(ADDITIONNAL_ROWS/ACC)
                       X=max(X,1)
                       ACC=0
                       DO i=1,CHOSEN
                          J=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(J)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             J=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(J)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                            BANDE_K821)
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          TMP_SUM=((dble(NELIM) * dble(NB_ROWS(i)))
     &                         +(dble(NB_ROWS(i))*dble(NELIM))*
     &                         dble(2*NFRONT-NELIM-1))
                          J=int((MAX_LOAD-
     &                         (TEMP(i)+TMP_SUM))/
     &                         (dble(NELIM)*dble(2*NFRONT-NELIM)))
                          IF(NB_ROWS(i).LT.KMAX)THEN
                             IF((min(X,J)+NB_ROWS(i)).GE.KMAX)THEN
                                IF((KMAX-NB_ROWS(i)).GT.
     &                               ADDITIONNAL_ROWS)THEN
                                   NB_ROWS(i)=NB_ROWS(i)+
     &                                  ADDITIONNAL_ROWS
                                   ADDITIONNAL_ROWS=0
                                ELSE
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                                  (KMAX-NB_ROWS(i))
                                   NB_ROWS(i)=KMAX
                                ENDIF
                             ELSE
                                IF((min(X,J)+NB_ROWS(i)).GE.
     &                            KMIN)THEN
                                   NB_ROWS(i)=NB_ROWS(i)+min(X,J)
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                                  min(X,J)
                                   ACC=ACC+1
                                ENDIF
                             ENDIF
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       IF(ACC.GT.0)THEN
                          DO i=1,CHOSEN
                             X=int(BUF_SIZE/dble(NCB+1))-1
                             BANDE_K821=dble(X)*dble(NFRONT)
                             IF(HAVE_TYPE1_SON)THEN
                                X=int((BUF_SIZE-dble(NFRONT))/
     &                               dble(NFRONT+1))
                                BANDE_K821=dble(X)*dble(NFRONT)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K821
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)),
     &                               BANDE_K821)
                                MAX_MEM_ALLOW=max(dble(0),
     &                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                             IF(KMAX-NB_ROWS(i).LT.
     &                            ADDITIONNAL_ROWS)THEN
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               (KMAX-NB_ROWS(i))
                                NB_ROWS(i)=KMAX
                             ELSE
                                IF(NB_ROWS(i).EQ.0)THEN
                                   IF(min(KMIN,KMAX).LT.
     &                                  ADDITIONNAL_ROWS)THEN
                                      NB_ROWS(i)=min(KMIN,KMAX)
                                      ADDITIONNAL_ROWS=
     &                                     ADDITIONNAL_ROWS-
     &                                     min(KMIN,KMAX)
                                   ENDIF
                                ELSE
                                   NB_ROWS(i)=NB_ROWS(i)+
     &                                  ADDITIONNAL_ROWS
                                   ADDITIONNAL_ROWS=0
                                ENDIF
                             ENDIF
                             IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                          ENDDO
                       ENDIF
                       DO i=1,CHOSEN
                          IDWLOAD(i)=i
                          AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          ENDIF                          
                          WLOAD(i)=(BANDE_K821-dble(NB_ROWS(i)*NFRONT))
                       ENDDO
                       CALL MUMPS_SORT_DOUBLES(NUMBER_OF_PROCS, WLOAD,
     &                      IDWLOAD)
                       NB_SAT=0
                       DO i=1,CHOSEN
                          X=int(ADDITIONNAL_ROWS/(CHOSEN-NB_SAT))
                          X=max(X,1)
                          AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                          ENDIF
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(BANDE_K821,
     &                            MEM_SIZE_STRONG(i))
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          IF(NB_ROWS(IDWLOAD(i)).LT.KMAX)THEN
                             IF((NB_ROWS(IDWLOAD(i))+X).LT.KMAX)THEN
                                NB_ROWS(IDWLOAD(i))=
     &                               NB_ROWS(IDWLOAD(i))+X
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                             ELSE
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                               (KMAX-NB_ROWS(IDWLOAD(i)))
                                NB_ROWS(IDWLOAD(i))=KMAX
                             ENDIF
                          ENDIF
                          IF(NB_ROWS(IDWLOAD(i)).EQ.KMAX)THEN
                             NB_SAT=NB_SAT+1
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0) GOTO 049
                       ENDDO
                       DO i=1,CHOSEN
                          X=int(BUF_SIZE/dble(NCB+1))-1
                          BANDE_K821=dble(X)*dble(NFRONT)
                          IF(HAVE_TYPE1_SON)THEN
                             X=int((BUF_SIZE-dble(NFRONT))/
     &                            dble(NFRONT+1))
                             BANDE_K821=dble(X)*dble(NFRONT)
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K821
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=min(BANDE_K821,
     &                            MEM_SIZE_STRONG(i))
                             MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                          IF(KMAX-NB_ROWS(i).LT.ADDITIONNAL_ROWS)THEN
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            (KMAX-NB_ROWS(i))
                             NB_ROWS(i)=KMAX
                          ELSE
                             NB_ROWS(i)=NB_ROWS(i)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       X=int(ADDITIONNAL_ROWS/CHOSEN)
                       X=max(X,1)
                       DO i=1,CHOSEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          NB_ROWS(i)=NB_ROWS(i)+X
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                    ENDIF
                 ENDIF
 049             CONTINUE
              ENDIF
 666          CONTINUE
              SOMME=dble(0)
              X=0
              POS=0
              DO i=1,CHOSEN
                 X=X+NB_ROWS(i)
                 SOMME=SOMME+ dble(NB_ROWS(i))
              ENDDO
              GOTO 890
           ELSE IF((KEEP(83).GE.NUMBER_OF_PROCS).AND.FORCE_CAND)THEN
              MAX_LOAD=dble(0)
              DO i=1,OTHERS
                 MAX_LOAD=max(MAX_LOAD,TEMP(i))
              ENDDO
              ACC=0
              CHOSEN=0
              X=1
              DO i=1,OTHERS
              ENDDO
              DO i=2,OTHERS 
                 IF(TEMP(i).EQ.TEMP(1))THEN
                    X=X+1
                 ELSE 
                    GOTO 329
                 ENDIF
              ENDDO
 329          CONTINUE
              TMP_SUM=TOTAL_COST/dble(X)
              TEMP_MAX_LOAD=dble(0)
              DO i=1,OTHERS
                 IF(K50.EQ.0)THEN
                    X=int(BUF_SIZE/dble(NCB+1))-1
                    BANDE_K821=dble(X)*dble(NFRONT)
                 ELSE
                    A=dble(1)
                    B=dble(ACC+2)
                    C=-BUF_SIZE+dble(ACC+NELIM)
                    DELTA=(B*B)-(dble(4)*A*C)
                    X=int((-B+sqrt(DELTA))/(dble(2)*A))
                    IF(X.GT.NCB-ACC) X=NCB-ACC
                    BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                 ENDIF
                 IF(HAVE_TYPE1_SON)THEN
                    IF(K50.EQ.0)THEN
                       X=int((BUF_SIZE-dble(NFRONT))/dble(NFRONT+1))
                       BANDE_K821=dble(X)*dble(NFRONT)
                    ELSE
                       A=dble(1)
                       B=dble(ACC+2+NELIM)
                       C=-BUF_SIZE+dble(ACC+NELIM)
                       DELTA=(B*B)-(dble(4)*A*C)
                       X=int((-B+sqrt(DELTA))/(dble(2)*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                    ENDIF
                 ENDIF
                 MAX_MEM_ALLOW=BANDE_K821
                 IF(BDC_MD)THEN
                    MAX_MEM_ALLOW=min(BANDE_K821,
     &                    min(MEM_SIZE_WEAK(i),MEM_SIZE_STRONG(i)))
                    MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                 ENDIF
                 IF(K50.EQ.0)THEN
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                    IF(TMP_SUM+TEMP(i).GT.MAX_LOAD)THEN
                       SOMME=MAX_LOAD-TEMP(i)
                    ELSE
                       SOMME=TMP_SUM
                    ENDIF
                    X=int(SOMME/
     &                   (dble(NELIM)*dble(2*NFRONT-NELIM)))
                    IF(X.GT.KMAX)THEN
                       X=KMAX
                    ELSE
                       IF(X.LT.KMIN)THEN
                          X=min(KMIN,KMAX)
                       ENDIF
                    ENDIF
                    IF((ACC+X).GT.NCB) X=NCB-ACC
                 ENDIF
                 IF(K50.NE.0)THEN
                       A=dble(1)
                       B=dble(ACC+NELIM)
                       C=dble(-MAX_MEM_ALLOW)
                       DELTA=((B*B)-(dble(4)*A*C))
                       KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                    A=dble(NELIM)
                    B=dble(NELIM)*dble(NELIM+2*ACC+1)
                    IF(TMP_SUM+TEMP(i).GT.MAX_LOAD)THEN
                       C=-(MAX_LOAD-TEMP(i))
                    ELSE
                       C=-TMP_SUM
                    ENDIF
                    DELTA=(B*B-(dble(4)*A*C))
                    X=int((-B+sqrt(DELTA))/(dble(2)*A))
                    IF(X.LT.0) THEN
                       WRITE(*,*)MYID,
     &    ': Internal error 14 in DMUMPS_SET_PARTI_FLOP_IRR'
                       CALL MUMPS_ABORT()
                    ENDIF
                    IF(X.GE.KMAX)THEN
                       IF(KMAX.GT.KMIN)THEN
                          X=KMAX
                       ELSE
                          X=0
                       ENDIF
                    ELSE
                       IF(X.LE.min(KMIN,KMAX))THEN
                          IF(KMAX.LT.KMIN)THEN
                             X=0
                          ELSE
                             X=min(KMIN,KMAX)
                          ENDIF
                       ENDIF
                    ENDIF
                    IF((ACC+X).GT.NCB) X=NCB-ACC
                 ENDIF
                 TEMP_MAX_LOAD=max(TEMP_MAX_LOAD,TEMP(i))
                 NB_ROWS(i)=X
                 CHOSEN=CHOSEN+1
                 ACC=ACC+X
                 IF(ACC.EQ.NCB) GOTO 541
              ENDDO
 541          CONTINUE
              IF(ACC.LT.NCB)THEN
                 IF(K50.EQ.0)THEN
                    ADDITIONNAL_ROWS=NCB-ACC
                    DO J=1,CHOSEN
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(
     &                    min(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     &                         dble(BANDE_K821))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF((NB_ROWS(J)).LT.KMAX)THEN
                          IF(ADDITIONNAL_ROWS.GT.(KMAX-NB_ROWS(J)))THEN
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            (KMAX-NB_ROWS(J))
                             NB_ROWS(J)=KMAX
                          ELSE
                             NB_ROWS(J)=NB_ROWS(J)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ENDIF
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO 
                    X=int(ADDITIONNAL_ROWS/CHOSEN)
                    X=max(X,1)
                    DO J=1,CHOSEN
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(BANDE_K821,
     &                         MEM_SIZE_STRONG(J))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF((NB_ROWS(J)+X).GT.KMAX)THEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                         (KMAX-NB_ROWS(J))
                          NB_ROWS(J)=KMAX
                       ELSE
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          NB_ROWS(J)=NB_ROWS(J)+X
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO 
                    DO i=1,CHOSEN
                       X=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(X)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          X=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(X)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(BANDE_K821,
     &                         MEM_SIZE_STRONG(i))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF(KMAX-NB_ROWS(i).LT.ADDITIONNAL_ROWS)THEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                         (KMAX-NB_ROWS(i))
                          NB_ROWS(i)=KMAX
                       ELSE
                          NB_ROWS(i)=NB_ROWS(i)+ADDITIONNAL_ROWS
                          ADDITIONNAL_ROWS=0
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO
                    DO i=1,NUMBER_OF_PROCS
                       IDWLOAD(i)=i
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF                          
                       WLOAD(i)=(BANDE_K821-(dble(NB_ROWS(i))*
     &                      dble(NFRONT)))
                    ENDDO
                    CALL MUMPS_SORT_DOUBLES(NUMBER_OF_PROCS, WLOAD,
     &                   IDWLOAD)
                    NB_SAT=0
                    DO i=1,CHOSEN
                       X=int(ADDITIONNAL_ROWS/(CHOSEN-NB_SAT))
                       X=max(X,1)
                       AFFECTED=int(BUF_SIZE/dble(NCB+1))-1
                       BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=int((BUF_SIZE-dble(NFRONT))/
     &                         dble(NFRONT+1))
                          BANDE_K821=dble(AFFECTED)*dble(NFRONT)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       KMAX=int(MAX_MEM_ALLOW/dble(NFRONT))
                       IF(NB_ROWS(IDWLOAD(i)).LT.KMAX)THEN
                          IF((NB_ROWS(IDWLOAD(i))+X).LT.KMAX)THEN
                             NB_ROWS(IDWLOAD(i))=
     &                            NB_ROWS(IDWLOAD(i))+X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          ELSE
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     &                            (KMAX-NB_ROWS(IDWLOAD(i)))
                             NB_ROWS(IDWLOAD(i))=KMAX
                          ENDIF
                       ENDIF
                       IF(NB_ROWS(IDWLOAD(i)).EQ.KMAX)THEN
                          NB_SAT=NB_SAT+1
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 889
                    ENDDO
                    GOTO 994
                 ELSE
                    ACC=0
                    CHOSEN=0
                    DO i=1,OTHERS
                       A=dble(1)
                       B=dble(ACC+2)
                       C=-BUF_SIZE+dble(ACC+NELIM)
                       DELTA=(B*B)-(dble(4)*A*C)
                       X=int((-B+sqrt(DELTA))/(dble(2)*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                       IF(HAVE_TYPE1_SON)THEN
                          A=dble(1)
                          B=dble(ACC+2+NELIM)
                          C=-BUF_SIZE+dble(ACC+NELIM)
                          DELTA=(B*B)-(dble(4)*A*C)
                          X=int((-B+sqrt(DELTA))/(dble(2)*A))
                          IF(X.GT.NCB-ACC) X=NCB-ACC
                          BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K821
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=min(BANDE_K821,
     &                         MEM_SIZE_STRONG(i))
                          MAX_MEM_ALLOW=max(dble(0),MAX_MEM_ALLOW)
                       ENDIF
                       A=dble(1)
                       B=dble(ACC+NELIM)
                       C=dble(-MAX_MEM_ALLOW)
                       DELTA=((B*B)-(dble(4)*A*C))
                       KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                       X=KMAX-NB_ROWS(i)
                       IF((ACC+NB_ROWS(i)+X).GT.NCB)
     &                            X=NCB-(ACC+NB_ROWS(i))
                       NB_ROWS(i)=NB_ROWS(i)+X
                       ACC=ACC+NB_ROWS(i)
                       CHOSEN=CHOSEN+1
                       IF(NCB.EQ.ACC) GOTO 889
                    ENDDO
                    ADDITIONNAL_ROWS=NCB-ACC
                 ENDIF
                 ACC=0
                 CHOSEN=0
                 DO i=1,OTHERS
                    A=dble(1)
                    B=dble(ACC+2)
                    C=-BUF_SIZE+dble(ACC+NELIM)
                    DELTA=(B*B)-(dble(4)*A*C)
                    X=int((-B+sqrt(DELTA))/(dble(2)*A))
                    IF(X.GT.NCB-ACC) X=NCB-ACC
                    BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                    IF(HAVE_TYPE1_SON)THEN
                       A=dble(1)
                       B=dble(ACC+2+NELIM)
                       C=-BUF_SIZE+dble(ACC+NELIM)
                       DELTA=(B*B)-(dble(4)*A*C)
                       X=int((-B+sqrt(DELTA))/(dble(2)*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K821=dble(X)*dble(NELIM+ACC+X)
                    ENDIF
                    MAX_MEM_ALLOW=BANDE_K821
                    A=dble(1)
                    B=dble(ACC+NELIM)
                    C=dble(-MAX_MEM_ALLOW)
                    DELTA=((B*B)-(dble(4)*A*C))
                    KMAX=int((-B+sqrt(DELTA))/(dble(2)*A))
                    X=KMAX-NB_ROWS(i)
                    IF((ACC+NB_ROWS(i)+X).GT.NCB)
     &                   X=NCB-(ACC+NB_ROWS(i))
                    NB_ROWS(i)=NB_ROWS(i)+X
                    ACC=ACC+NB_ROWS(i)
                    CHOSEN=CHOSEN+1
                    IF(NCB.EQ.ACC) GOTO 889
                 ENDDO
                 ADDITIONNAL_ROWS=NCB-ACC
 994             CONTINUE
                 X=int(dble(ADDITIONNAL_ROWS)/dble(OTHERS))
                 IF((X*OTHERS).LT.ADDITIONNAL_ROWS)THEN
                    X=X+1
                 ENDIF
                 DO i=1,OTHERS
                    NB_ROWS(i)=NB_ROWS(i)+X
                    ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                    IF(ADDITIONNAL_ROWS.LT.X)X=ADDITIONNAL_ROWS
                 ENDDO
                 CHOSEN=OTHERS
              ENDIF
           ENDIF
 889       CONTINUE
           MAX_LOAD=TEMP_MAX_LOAD
 890       CONTINUE
           J=CHOSEN
           X=0
              DO i=J,1,-1
                 IF(NB_ROWS(i).EQ.0)THEN
                    CHOSEN=CHOSEN-1
                    ELSE 
                       IF(NB_ROWS(i).GT.0)THEN
                          X=1
                       ELSE
                          WRITE(*,*)MYID,
     &    ': Internal error 15 in DMUMPS_SET_PARTI_FLOP_IRR'
                          CALL MUMPS_ABORT()
                       ENDIF
                    ENDIF
                 ENDDO
           NSLAVES_NODE=CHOSEN
           TAB_POS(NSLAVES_NODE+1)= NCB+1
           TAB_POS(SLAVEF+2) = CHOSEN
           POS=1
           X=1
           DO i=1,J
              IF(NB_ROWS(i).NE.0)THEN
                 SLAVES_LIST(X)=TEMP_ID(i)
                 TAB_POS(X)=POS
                 POS=POS+NB_ROWS(i) 
                 IF(NB_ROWS(i).LE.0)THEN
                    WRITE(*,*)MYID,
     &    ': Internal error 16 in DMUMPS_SET_PARTI_FLOP_IRR'
                    CALL MUMPS_ABORT()
                 ENDIF
                 X=X+1
               ENDIF
           ENDDO
           IF(POS.NE.(NCB+1))THEN
              WRITE(*,*)MYID,
     &    ': Internal error 17 in DMUMPS_SET_PARTI_FLOP_IRR',
     &             POS,NCB+1
              CALL MUMPS_ABORT()
           ENDIF
      END SUBROUTINE DMUMPS_SET_PARTI_FLOP_IRR
      SUBROUTINE DMUMPS_LOAD_POOL_CHECK_MEM
     &      (INODE,UPPER,SLAVEF,KEEP,KEEP8,
     &       STEP,POOL,LPOOL,PROCNODE,N)
      IMPLICIT NONE
      INTEGER INODE, LPOOL, SLAVEF, N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(KEEP(28)), POOL(LPOOL), PROCNODE(KEEP(28))
      LOGICAL UPPER
      INTEGER J
      DOUBLE PRECISION MEM_COST
      INTEGER NBINSUBTREE,i,NBTOP
      EXTERNAL DMUMPS_POOL_EMPTY,
     & MUMPS_IN_OR_ROOT_SSARBR
      LOGICAL DMUMPS_POOL_EMPTY,
     & MUMPS_IN_OR_ROOT_SSARBR
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      IF(KEEP(47).LT.2)THEN
         WRITE(*,*)'DMUMPS_LOAD_POOL_CHECK_MEM must
     &        be called with K47>=2'
         CALL MUMPS_ABORT()
      ENDIF        
      IF((INODE.GT.0).AND.(INODE.LE.N))THEN
      MEM_COST=DMUMPS_LOAD_GET_MEM(INODE)
         IF((DM_MEM(MYID)+dble(MEM_COST)+ PEAK_SBTR_CUR_LOCAL-
     &        SBTR_CUR_LOCAL)
     &        .GT.MAX_PEAK_STK)THEN
            DO i=NBTOP-1,1,-1
               INODE = POOL( LPOOL - 2 - i)
               MEM_COST=DMUMPS_LOAD_GET_MEM(INODE)
               IF((INODE.LT.0).OR.(INODE.GT.N)) THEN
                  DO J=i+1,NBTOP,-1
                     POOL(J-1)=POOL(J)
                  ENDDO
                  UPPER=.TRUE.
                  RETURN
               ENDIF
               IF((DM_MEM(MYID)+dble(MEM_COST)+PEAK_SBTR_CUR_LOCAL-
     &              SBTR_CUR_LOCAL).LE.
     &              MAX_PEAK_STK) THEN
                  DO J=i+1,NBTOP,-1
                     POOL(J-1)=POOL(J)
                  ENDDO
                  UPPER=.TRUE.
                  RETURN
               ENDIF
            ENDDO
            IF(NBINSUBTREE.NE.0)THEN
               INODE = POOL( NBINSUBTREE )
               IF(.NOT.MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &              KEEP(199)))THEN
                  WRITE(*,*)
     &        'Internal error 1 in DMUMPS_LOAD_POOL_CHECK_MEM'
                  CALL MUMPS_ABORT()
               ENDIF
               UPPER=.FALSE.
               RETURN
            ENDIF
            INODE=POOL(LPOOL-2-NBTOP)
            UPPER=.TRUE.
            RETURN
         ENDIF
      ENDIF
      UPPER=.TRUE.
      END SUBROUTINE DMUMPS_LOAD_POOL_CHECK_MEM
      SUBROUTINE DMUMPS_LOAD_SET_SBTR_MEM(WHAT)
      IMPLICIT NONE
      LOGICAL WHAT
      IF(.NOT.BDC_POOL_MNG)THEN
         WRITE(*,*)'DMUMPS_LOAD_SET_SBTR_MEM
     &        should be called when K81>0 and K47>2'
      ENDIF
      IF(WHAT)THEN
         PEAK_SBTR_CUR_LOCAL=PEAK_SBTR_CUR_LOCAL+
     &        dble(MEM_SUBTREE(INDICE_SBTR))
         IF(.NOT.BDC_SBTR) INDICE_SBTR=INDICE_SBTR+1
      ELSE
         PEAK_SBTR_CUR_LOCAL=dble(0)
         SBTR_CUR_LOCAL=dble(0)
      ENDIF
      END SUBROUTINE DMUMPS_LOAD_SET_SBTR_MEM
      DOUBLE PRECISION FUNCTION DMUMPS_LOAD_GET_MEM( INODE )
      IMPLICIT NONE
      INTEGER INODE,LEVEL,i,NELIM,NFR
      DOUBLE PRECISION COST
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      i = INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
        NELIM = NELIM + 1
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      NFR = ND_LOAD( STEP_LOAD(INODE) ) + KEEP_LOAD(253)
      LEVEL = MUMPS_TYPENODE( PROCNODE_LOAD(STEP_LOAD(INODE)),
     &                        KEEP_LOAD(199) )
      IF (LEVEL .EQ. 1) THEN
        COST =  dble(NFR) * dble(NFR)
      ELSE
        IF ( K50 == 0 ) THEN
           COST =  dble(NFR) * dble(NELIM)
        ELSE
           COST = dble(NELIM) * dble(NELIM)
        ENDIF
      ENDIF
      DMUMPS_LOAD_GET_MEM=COST
      RETURN
      END FUNCTION DMUMPS_LOAD_GET_MEM
      RECURSIVE SUBROUTINE DMUMPS_NEXT_NODE(FLAG,COST,COMM)
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER COMM,WHAT,IERR
      LOGICAL FLAG, EXIT_FLAG
      DOUBLE PRECISION COST
      DOUBLE PRECISION TO_BE_SENT
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      IF(FLAG)THEN
         WHAT=17 
         IF(BDC_M2_FLOPS)THEN
            TO_BE_SENT=DELTA_LOAD-COST
            DELTA_LOAD=dble(0)
         ELSE IF(BDC_M2_MEM)THEN
            IF(BDC_POOL.AND.(.NOT.BDC_MD))THEN
               TO_BE_SENT=max(TMP_M2,POOL_LAST_COST_SENT)
               POOL_LAST_COST_SENT=TO_BE_SENT
            ELSE IF(BDC_MD)THEN
               DELTA_MEM=DELTA_MEM+TMP_M2
               TO_BE_SENT=DELTA_MEM
            ELSE
               TO_BE_SENT=dble(0)
            ENDIF
         ENDIF
      ELSE
         WHAT=6
         TO_BE_SENT=dble(0)
      ENDIF
 111  CONTINUE
      CALL DMUMPS_BUF_BROADCAST( WHAT,
     &         COMM, NPROCS,
     &               FUTURE_NIV2,
     &         COST, 
     &         TO_BE_SENT, 
     &         MYID, KEEP_LOAD, IERR  )
      IF ( IERR == -1 )THEN
         CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
         CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
         IF (EXIT_FLAG) THEN
            GOTO 100
         ELSE
            GOTO 111
         ENDIF
      ELSE IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) "Internal Error in DMUMPS_LOAD_POOL_UPD_NEW_POOL",
     &   IERR
         CALL MUMPS_ABORT()
      ENDIF
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_NEXT_NODE
      SUBROUTINE DMUMPS_UPPER_PREDICT(INODE,STEP,NSTEPS,PROCNODE,FRERE,
     &     NE,COMM,SLAVEF,MYID,KEEP,KEEP8,N)
      USE DMUMPS_BUF
      IMPLICIT NONE
      INTEGER INODE,NSTEPS,MYID,SLAVEF,COMM,N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER FRERE(NSTEPS),NE(NSTEPS),STEP(N),PROCNODE(NSTEPS)
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR,MUMPS_PROCNODE
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      INTEGER i,NCB,NELIM
      INTEGER MUMPS_PROCNODE
      INTEGER FATHER_NODE,FATHER,WHAT,IERR
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      LOGICAL :: EXIT_FLAG
      IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
         WRITE(*,*)MYID,': Problem in DMUMPS_UPPER_PREDICT'
         CALL MUMPS_ABORT()
      ENDIF
      IF((INODE.LT.0).OR.(INODE.GT.N)) THEN
         RETURN
      ENDIF
      i=INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
         NELIM = NELIM + 1
         i = FILS_LOAD(i)
         GOTO 10
      ENDIF
      NCB=ND_LOAD(STEP_LOAD(INODE))-NELIM + KEEP_LOAD(253)
      WHAT=5
      FATHER_NODE=DAD_LOAD(STEP_LOAD(INODE))
      IF (FATHER_NODE.EQ.0) THEN
         RETURN
      ENDIF
      IF((FRERE(STEP(FATHER_NODE)).EQ.0).AND.
     &     ((FATHER_NODE.EQ.KEEP(38)).OR.
     &     (FATHER_NODE.EQ.KEEP(20))))THEN
         RETURN 
      ENDIF
      IF(MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(FATHER_NODE)),
     &            KEEP(199))) THEN
         RETURN
      ENDIF
      FATHER=MUMPS_PROCNODE(PROCNODE(STEP(FATHER_NODE)),KEEP(199))
      IF(FATHER.EQ.MYID)THEN
        IF(BDC_M2_MEM)THEN
           CALL DMUMPS_PROCESS_NIV2_MEM_MSG(FATHER_NODE)
        ELSEIF(BDC_M2_FLOPS)THEN
           CALL DMUMPS_PROCESS_NIV2_FLOPS_MSG(FATHER_NODE)
        ENDIF
        IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
           IF(MUMPS_TYPENODE(PROCNODE_LOAD(STEP_LOAD(INODE)),
     &          KEEP(199)).EQ.1)THEN
              CB_COST_ID(POS_ID)=INODE
              CB_COST_ID(POS_ID+1)=1
              CB_COST_ID(POS_ID+2)=POS_MEM
              POS_ID=POS_ID+3
              CB_COST_MEM(POS_MEM)=int(MYID,8)
              POS_MEM=POS_MEM+1
              CB_COST_MEM(POS_MEM)=int(NCB,8)*int(NCB,8)
              POS_MEM=POS_MEM+1
           ENDIF
        ENDIF
        GOTO 666
      ENDIF
 111  CONTINUE      
      CALL DMUMPS_BUF_SEND_FILS(WHAT, COMM, NPROCS,
     &     FATHER_NODE,INODE,NCB, KEEP,MYID,
     &     FATHER, IERR)
      IF (IERR == -1 ) THEN
        CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
        CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
        IF (EXIT_FLAG) THEN
           GOTO 666
        ELSE
           GOTO 111
        ENDIF
      ELSE IF ( IERR .NE. 0 ) THEN
        WRITE(*,*) "Internal Error in DMUMPS_UPPER_PREDICT",
     &  IERR
        CALL MUMPS_ABORT()
      ENDIF 
 666  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_UPPER_PREDICT
      SUBROUTINE DMUMPS_REMOVE_NODE(INODE,NUM_CALL)
      IMPLICIT NONE
      DOUBLE PRECISION MAXI
      INTEGER i,J,IND_MAXI
      INTEGER INODE,NUM_CALL
      IF(BDC_M2_MEM)THEN
         IF(((NUM_CALL.EQ.1).AND.(BDC_MD)).OR.
     &       ((NUM_CALL.EQ.2).AND.(.NOT.BDC_MD)))THEN
            RETURN
         ENDIF
      ENDIF
      IF((FRERE_LOAD(STEP_LOAD(INODE)).EQ.0).AND.
     &     ((INODE.EQ.KEEP_LOAD(38)).OR.
     &     (INODE.EQ.KEEP_LOAD(20)))) THEN
         RETURN
      ENDIF
      DO i=POOL_SIZE,1,-1
         IF(POOL_NIV2(i).EQ.INODE) GOTO 666
      ENDDO
         NB_SON(STEP_LOAD(INODE))=-1
      RETURN         
 666  CONTINUE
      IF(BDC_M2_MEM)THEN
         IF(POOL_NIV2_COST(i).EQ.MAX_M2)THEN
            TMP_M2=MAX_M2
            MAXI=dble(0)
            IND_MAXI=-9999
            DO J=POOL_SIZE,1,-1
               IF(J.NE.i) THEN
                  IF(POOL_NIV2_COST(J).GT.MAXI)THEN
                     MAXI=POOL_NIV2_COST(J)
                     IND_MAXI=J
                  ENDIF
               ENDIF
            ENDDO
            MAX_M2=MAXI
            J=IND_MAXI
            REMOVE_NODE_FLAG_MEM=.TRUE.
            REMOVE_NODE_COST_MEM=TMP_M2
            CALL DMUMPS_NEXT_NODE(REMOVE_NODE_FLAG,MAX_M2,COMM_LD)
            NIV2(MYID+1)=MAX_M2
         ENDIF
      ELSEIF(BDC_M2_FLOPS)THEN
         REMOVE_NODE_COST=POOL_NIV2_COST(i)
         REMOVE_NODE_FLAG=.TRUE.
         CALL DMUMPS_NEXT_NODE(REMOVE_NODE_FLAG,
     &        -POOL_NIV2_COST(i),COMM_LD)
         NIV2(MYID+1)=NIV2(MYID+1)-POOL_NIV2_COST(i)
      ENDIF
      DO J=i+1,POOL_SIZE
         POOL_NIV2(J-1)=POOL_NIV2(J)
         POOL_NIV2_COST(J-1)=POOL_NIV2_COST(J)
      ENDDO
      POOL_SIZE=POOL_SIZE-1
      END SUBROUTINE DMUMPS_REMOVE_NODE
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_NIV2_MEM_MSG(INODE)
      IMPLICIT NONE
      INTEGER INODE
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      IF((INODE.EQ.KEEP_LOAD(20)).OR.
     &     (INODE.EQ.KEEP_LOAD(38)))THEN
         RETURN
      ENDIF
      IF(NB_SON(STEP_LOAD(INODE)).EQ.-1)THEN
         RETURN
      ELSE
         IF(NB_SON(STEP_LOAD(INODE)).LT.0)THEN
            WRITE(*,*)
     &        'Internal error 1 in DMUMPS_PROCESS_NIV2_MEM_MSG'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      NB_SON(STEP_LOAD(INODE))=
     &     NB_SON(STEP_LOAD(INODE))-1
      IF(NB_SON(STEP_LOAD(INODE)).EQ.0)THEN
         IF(POOL_SIZE.EQ.POOL_NIV2_SIZE)THEN
            WRITE(*,*)MYID,': Internal Error 2 in
     &DMUMPS_PROCESS_NIV2_MEM_MSG'
            CALL MUMPS_ABORT()
         ENDIF
         POOL_NIV2(POOL_SIZE+1)=INODE
         POOL_NIV2_COST(POOL_SIZE+1)=
     &        DMUMPS_LOAD_GET_MEM(INODE)
         POOL_SIZE=POOL_SIZE+1
         IF(POOL_NIV2_COST(POOL_SIZE).GT.MAX_M2)THEN
            MAX_M2=POOL_NIV2_COST(POOL_SIZE)
            ID_MAX_M2=POOL_NIV2(POOL_SIZE)
            CALL DMUMPS_NEXT_NODE(REMOVE_NODE_FLAG_MEM,MAX_M2,COMM_LD)
            NIV2(1+MYID)=MAX_M2
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_NIV2_MEM_MSG    
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_NIV2_FLOPS_MSG(INODE)
      IMPLICIT NONE
      INTEGER INODE
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      IF((INODE.EQ.KEEP_LOAD(20)).OR.
     &     (INODE.EQ.KEEP_LOAD(38)))THEN
         RETURN
      ENDIF
      IF(NB_SON(STEP_LOAD(INODE)).EQ.-1)THEN
         RETURN
      ELSE
         IF(NB_SON(STEP_LOAD(INODE)).LT.0)THEN
            WRITE(*,*)
     &        'Internal error 1 in DMUMPS_PROCESS_NIV2_FLOPS_MSG'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      NB_SON(STEP_LOAD(INODE))=
     &     NB_SON(STEP_LOAD(INODE))-1
      IF(NB_SON(STEP_LOAD(INODE)).EQ.0)THEN
         IF(POOL_SIZE.EQ.POOL_NIV2_SIZE)THEN
            WRITE(*,*)MYID,': Internal Error 2 in
     &DMUMPS_PROCESS_NIV2_FLOPS_MSG',POOL_NIV2_SIZE,
     &           POOL_SIZE
            CALL MUMPS_ABORT()
         ENDIF
         POOL_NIV2(POOL_SIZE+1)=INODE
         POOL_NIV2_COST(POOL_SIZE+1)=
     &        DMUMPS_LOAD_GET_FLOPS_COST(INODE)
         POOL_SIZE=POOL_SIZE+1
         MAX_M2=POOL_NIV2_COST(POOL_SIZE)
         ID_MAX_M2=POOL_NIV2(POOL_SIZE)
         CALL DMUMPS_NEXT_NODE(REMOVE_NODE_FLAG,
     &           POOL_NIV2_COST(POOL_SIZE),
     &        COMM_LD)
         NIV2(MYID+1)=POOL_NIV2_COST(POOL_SIZE)+NIV2(MYID+1)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_NIV2_FLOPS_MSG
      DOUBLE PRECISION FUNCTION DMUMPS_LOAD_GET_FLOPS_COST(INODE)
      USE MUMPS_FUTURE_NIV2
      INTEGER INODE
      INTEGER NFRONT,NELIM,i,LEVEL
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      DOUBLE PRECISION COST
      i = INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
        NELIM = NELIM + 1
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      NFRONT = ND_LOAD( STEP_LOAD(INODE) ) + KEEP_LOAD(253)
      LEVEL = MUMPS_TYPENODE( PROCNODE_LOAD(STEP_LOAD(INODE)),
     &        KEEP_LOAD(199) )
      COST=dble(0)
      CALL MUMPS_GET_FLOPS_COST(NFRONT,NELIM,NELIM,
     &                          KEEP_LOAD(50),LEVEL,COST)
      DMUMPS_LOAD_GET_FLOPS_COST=COST
      RETURN
      END FUNCTION DMUMPS_LOAD_GET_FLOPS_COST
      INTEGER FUNCTION DMUMPS_LOAD_GET_CB_FREED( INODE )
      IMPLICIT NONE
      INTEGER INODE,NELIM,NFR,SON,IN,i
      INTEGER COST_CB
      COST_CB=0
      i = INODE
 10   CONTINUE
      IF ( i > 0 ) THEN
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      SON=-i
      DO i=1, NE_LOAD(STEP_LOAD(INODE))
         NFR = ND_LOAD( STEP_LOAD(SON) ) + KEEP_LOAD(253)
         IN=SON
         NELIM = 0
 20      CONTINUE
         IF ( IN > 0 ) THEN
            NELIM = NELIM + 1
            IN = FILS_LOAD(IN)
            GOTO 20
         ENDIF
         COST_CB=COST_CB+((NFR-NELIM)*(NFR-NELIM))
         SON=FRERE_LOAD(STEP_LOAD(SON))
      ENDDO
      DMUMPS_LOAD_GET_CB_FREED=COST_CB
      RETURN
      END FUNCTION DMUMPS_LOAD_GET_CB_FREED
      SUBROUTINE DMUMPS_LOAD_SEND_MD_INFO(SLAVEF,NMB_OF_CAND,
     &     LIST_OF_CAND,
     &     TAB_POS, NASS, KEEP,KEEP8, LIST_SLAVES,
     &     NSLAVES,INODE)
      USE DMUMPS_BUF
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, NASS, NSLAVES
      INTEGER, INTENT (IN) :: NMB_OF_CAND
      INTEGER, INTENT (IN) :: LIST_OF_CAND(NMB_OF_CAND)
      INTEGER, INTENT (IN) :: TAB_POS(SLAVEF+2)
      INTEGER, INTENT (IN) :: LIST_SLAVES(NSLAVES)
      INTEGER KEEP(500),INODE
      INTEGER(8) KEEP8(150)
      INTEGER allocok
      DOUBLE PRECISION MEM_COST,FCT_COST
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: DELTA_MD
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPROC2POSINDELTAMD
      INTEGER, DIMENSION(:), ALLOCATABLE :: P_TO_UPDATE
      INTEGER NBROWS_SLAVE,i,WHAT,IERR
      INTEGER :: NP_TO_UPDATE, K
      LOGICAL FORCE_CAND
      LOGICAL :: EXIT_FLAG
      MEM_COST=dble(0)
      FCT_COST=dble(0)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        FORCE_CAND = .FALSE.
      ELSE
        FORCE_CAND = (mod(KEEP(24),2).eq.0)
      END IF
      CALL DMUMPS_LOAD_GET_ESTIM_MEM_COST(INODE,FCT_COST,
     &        MEM_COST,NMB_OF_CAND,NASS)
      ALLOCATE(IPROC2POSINDELTAMD(0:SLAVEF-1),
     & DELTA_MD(min(SLAVEF, NMB_OF_CAND+NSLAVES)),
     & P_TO_UPDATE(min(SLAVEF, NMB_OF_CAND+NSLAVES)),
     & stat=allocok)
      IF (allocok > 0 ) THEN
        WRITE(*,*) "PB ALLOC IN DMUMPS_LOAD_SEND_MD_INFO",
     &  SLAVEF, NMB_OF_CAND, NSLAVES
        CALL MUMPS_ABORT()
      ENDIF
      IPROC2POSINDELTAMD = -99
      NP_TO_UPDATE = 0
      DO i = 1, NSLAVES
        NP_TO_UPDATE = NP_TO_UPDATE + 1
        IPROC2POSINDELTAMD (LIST_SLAVES(i)) = NP_TO_UPDATE
        NBROWS_SLAVE = TAB_POS(i+1) - TAB_POS(i)
        DELTA_MD(NP_TO_UPDATE)=-dble(NBROWS_SLAVE)*
     &           dble(NASS)
        P_TO_UPDATE(NP_TO_UPDATE) = LIST_SLAVES(i)
      ENDDO
      DO i = 1, NMB_OF_CAND
        K = IPROC2POSINDELTAMD(LIST_OF_CAND(i))
        IF ( K > 0 ) THEN
          DELTA_MD(K)=DELTA_MD(K)+FCT_COST
        ELSE
          NP_TO_UPDATE = NP_TO_UPDATE + 1
          IPROC2POSINDELTAMD (LIST_OF_CAND(i)) = NP_TO_UPDATE
          DELTA_MD   (NP_TO_UPDATE) = FCT_COST
          P_TO_UPDATE(NP_TO_UPDATE) = LIST_OF_CAND(i)
        ENDIF
      ENDDO
      WHAT=7
 111  CONTINUE
      CALL DMUMPS_BUF_BCAST_ARRAY(.FALSE., COMM_LD, MYID, SLAVEF,
     &     FUTURE_NIV2,
     &     NP_TO_UPDATE, P_TO_UPDATE,0,
     &     DELTA_MD, 
     &     DELTA_MD, 
     &     DELTA_MD, 
     &     WHAT, KEEP, IERR)
      IF ( IERR == -1 ) THEN
          CALL DMUMPS_LOAD_RECV_MSGS(COMM_LD)
          CALL MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
          IF (EXIT_FLAG) THEN
             GOTO 100
          ELSE
             GOTO 111
          ENDIF
       ELSE IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) "Internal Error 2 in DMUMPS_LOAD_SEND_MD_INFO",
     &   IERR
         CALL MUMPS_ABORT()
      ENDIF
      IF (FUTURE_NIV2(MYID+1) .NE. 0) THEN
        DO i = 1, NP_TO_UPDATE
           MD_MEM(P_TO_UPDATE(i))=MD_MEM(P_TO_UPDATE(i))+
     &          int(DELTA_MD( i ),8)
           IF(FUTURE_NIV2(P_TO_UPDATE(i)+1).EQ.0)THEN
              MD_MEM(P_TO_UPDATE(i))=999999999_8
           ENDIF
        ENDDO
      ENDIF
 100  CONTINUE
      DEALLOCATE(DELTA_MD,P_TO_UPDATE,IPROC2POSINDELTAMD)
      RETURN
      END SUBROUTINE DMUMPS_LOAD_SEND_MD_INFO
      SUBROUTINE DMUMPS_LOAD_GET_ESTIM_MEM_COST(INODE,FCT_COST,
     &     MEM_COST,NSLAVES,NELIM)
      IMPLICIT NONE
      INTEGER INODE,NSLAVES,NFR,NELIM,IN
      DOUBLE PRECISION MEM_COST,FCT_COST
      NFR=ND_LOAD(STEP_LOAD(INODE)) + KEEP_LOAD(253)
      IN = INODE
      FCT_COST=dble(int(dble(NFR-NELIM)/dble(NSLAVES))+1)*
     &     dble(NELIM)
      MEM_COST=dble(int(dble(NFR-NELIM)/dble(NSLAVES))+1)*
     &     dble(NFR)
      END SUBROUTINE DMUMPS_LOAD_GET_ESTIM_MEM_COST
      SUBROUTINE DMUMPS_LOAD_CLEAN_MEMINFO_POOL(INODE)
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER INODE
      INTEGER i,J,SON,NSLAVES_TEMP,POS_TEMP,K
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      IF((INODE.LT.0).OR.(INODE.GT.N_LOAD))THEN
         RETURN
      ENDIF
      IF(POS_ID.GT.1)THEN
         i=INODE
 10      CONTINUE
         IF ( i > 0 ) THEN
            i = FILS_LOAD(i)
            GOTO 10
         ENDIF
         SON=-i
         IF(POS_ID.LT.NE_LOAD(STEP_LOAD(INODE))*3)THEN
            i=1
         ENDIF
         DO i=1, NE_LOAD(STEP_LOAD(INODE))
            J=1
            DO WHILE (J.LT.POS_ID)
               IF(CB_COST_ID(J).EQ.SON)GOTO 295 
               J=J+3
            ENDDO
 295        CONTINUE
            IF(J.GE.POS_ID)THEN
               IF ( MUMPS_PROCNODE(
     &                PROCNODE_LOAD(STEP_LOAD(INODE)),
     &                KEEP_LOAD(199) ) .EQ. MYID ) THEN
                  IF(INODE.EQ.KEEP_LOAD(38))THEN
                     GOTO 666
                  ELSE
                     IF(FUTURE_NIV2(MYID+1).NE.0)THEN
                        WRITE(*,*)MYID,': i did not find ',SON
                        CALL MUMPS_ABORT()
                     ENDIF
                     GOTO 666
                  ENDIF
               ELSE
                  GOTO 666
               ENDIF
            ENDIF
            NSLAVES_TEMP=CB_COST_ID(J+1)       
            POS_TEMP=CB_COST_ID(J+2)
            DO K=J,POS_ID-1
               CB_COST_ID(K)=CB_COST_ID(K+3)
            ENDDO
            K=POS_TEMP
            DO WHILE (K.LE.POS_MEM-1)
               CB_COST_MEM(K)=CB_COST_MEM(K+2*NSLAVES_TEMP)
               K=K+1
            ENDDO
            POS_MEM=POS_MEM-2*NSLAVES_TEMP
            POS_ID=POS_ID-3
            IF((POS_MEM.LT.1).OR.(POS_ID.LT.1))THEN
               WRITE(*,*)MYID,': negative pos_mem or pos_id'
               CALL MUMPS_ABORT()
            ENDIF
 666        CONTINUE
            SON=FRERE_LOAD(STEP_LOAD(SON))
         ENDDO
      ENDIF
      END SUBROUTINE DMUMPS_LOAD_CLEAN_MEMINFO_POOL
      SUBROUTINE DMUMPS_LOAD_CHK_MEMCST_POOL(FLAG)
      IMPLICIT NONE
      LOGICAL FLAG
      INTEGER i
      DOUBLE PRECISION MEM
      FLAG=.FALSE.
      DO i=0,NPROCS-1
         MEM=DM_MEM(i)+LU_USAGE(i)
         IF(BDC_SBTR)THEN
            MEM=MEM+SBTR_MEM(i)-SBTR_CUR(i)
         ENDIF
         IF((MEM/dble(TAB_MAXS(i))).GT.0.8d0)THEN
            FLAG=.TRUE.
            GOTO 666
         ENDIF
      ENDDO
 666  CONTINUE
      END SUBROUTINE DMUMPS_LOAD_CHK_MEMCST_POOL
      SUBROUTINE DMUMPS_CHECK_SBTR_COST(NBINSUBTREE,INSUBTREE,NBTOP,
     &           MIN_COST,SBTR)
      IMPLICIT NONE
      INTEGER NBINSUBTREE,INSUBTREE,NBTOP
      DOUBLE PRECISION MIN_COST
      LOGICAL SBTR
      INTEGER i
      DOUBLE PRECISION TMP_COST,TMP_MIN
      TMP_MIN=huge(TMP_MIN)
      DO i=0,NPROCS-1
         IF(i.NE.MYID)THEN
            IF(BDC_SBTR)THEN
               TMP_MIN=min(TMP_MIN,dble(TAB_MAXS(i))-(DM_MEM(i)+
     &              LU_USAGE(i))-(SBTR_MEM(i)-SBTR_CUR(i)))
            ELSE
               TMP_MIN=min(TMP_MIN,dble(TAB_MAXS(i))-
     &              (DM_MEM(i)+LU_USAGE(i)))
            ENDIF
         ENDIF
      ENDDO
      IF(NBINSUBTREE.GT.0)THEN
         IF(INSUBTREE.EQ.1)THEN
            TMP_COST=dble(TAB_MAXS(MYID))-(DM_MEM(MYID)+
     &           LU_USAGE(MYID))
     &           -(SBTR_MEM(MYID)-SBTR_CUR(MYID))
         ELSE
            SBTR=.FALSE.
            GOTO 777
         ENDIF
      ENDIF
      TMP_MIN=min(TMP_COST,TMP_MIN)
      IF(TMP_MIN.GT.MIN_COST) SBTR=.TRUE.
 777  CONTINUE
      END SUBROUTINE DMUMPS_CHECK_SBTR_COST
      SUBROUTINE DMUMPS_LOAD_COMP_MAXMEM_POOL(INODE,MAX_MEM,PROC)
      USE MUMPS_FUTURE_NIV2
      IMPLICIT NONE
      INTEGER INODE,PROC
      INTEGER i,POS,NSLAVES,SLAVE,NCAND,J,NELIM,NCB,NFRONT,SON,K
      INTEGER allocok
      EXTERNAL MUMPS_TYPENODE
      INTEGER  MUMPS_TYPENODE
      DOUBLE PRECISION MAX_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MEM_ON_PROCS,
     &     RECV_BUF
      LOGICAL, DIMENSION(:), ALLOCATABLE :: CONCERNED
      DOUBLE PRECISION MAX_SENT_MSG
      IF((FRERE_LOAD(STEP_LOAD(INODE)).EQ.0)
     &           .AND.(INODE.EQ.KEEP_LOAD(38)))THEN
         RETURN
      ENDIF
      ALLOCATE( MEM_ON_PROCS(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in DMUMPS_LOAD_COMP_MAXMEM_POOL'
        CALL MUMPS_ABORT()
      ENDIF 
      ALLOCATE( CONCERNED(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in DMUMPS_LOAD_COMP_MAXMEM_POOL'
        CALL MUMPS_ABORT()
      ENDIF 
      ALLOCATE( RECV_BUF(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in DMUMPS_LOAD_COMP_MAXMEM_POOL'
        CALL MUMPS_ABORT()
      ENDIF 
      RECV_BUF=dble(0)
      MAX_SENT_MSG=dble(0)
      i = INODE
      NELIM = 0
 10   CONTINUE
      IF ( i > 0 ) THEN
        NELIM = NELIM + 1
        i = FILS_LOAD(i)
        GOTO 10
      ENDIF
      SON=-i
      NFRONT=ND_LOAD(STEP_LOAD(INODE)) + KEEP_LOAD(253)
      NCB=NFRONT-NELIM
      IF(MUMPS_TYPENODE(PROCNODE_LOAD(STEP_LOAD(INODE)),
     &     KEEP_LOAD(199)).EQ.2)THEN
         NCAND=CAND_LOAD(NPROCS+1, STEP_TO_NIV2_LOAD(STEP_LOAD(INODE)))
      ENDIF
      DO i=0,NPROCS-1
         IF(i.EQ.MYID)THEN
            MEM_ON_PROCS(i)=dble(TAB_MAXS(i))-(DM_MEM(i)+
     &           LU_USAGE(i)+
     &           DMUMPS_LOAD_GET_MEM(INODE))
            IF(BDC_SBTR)THEN
               MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-(SBTR_MEM(i)-SBTR_CUR(i))
            ENDIF
            CONCERNED(i)=.TRUE.
         ELSE
            MEM_ON_PROCS(i)=dble(TAB_MAXS(i))-(DM_MEM(i)+LU_USAGE(i))
            IF(BDC_SBTR)THEN
               MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-(SBTR_MEM(i)-SBTR_CUR(i))
            ENDIF
            IF(BDC_M2_MEM)THEN
               MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-NIV2(i+1)
            ENDIF
         ENDIF
         IF(MUMPS_TYPENODE(PROCNODE_LOAD(STEP_LOAD(INODE)),
     &        KEEP_LOAD(199)).EQ.2)THEN
            IF(BDC_MD.AND.(KEEP_LOAD(48).EQ.5))THEN
               DO J=1,NCAND
                  IF(CAND_LOAD(J, STEP_TO_NIV2_LOAD(STEP_LOAD(INODE)))
     &                 .EQ.i)THEN
                     MEM_ON_PROCS(i)=MEM_ON_PROCS(i)-
     &                 ((dble(NFRONT)*dble(NCB))/dble(NCAND))
                     CONCERNED(i)=.TRUE.
                     GOTO 666
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
 666     CONTINUE
      ENDDO
      DO K=1, NE_LOAD(STEP_LOAD(INODE))
         i=1
         DO WHILE (i.LE.POS_ID)
            IF(CB_COST_ID(i).EQ.SON)GOTO 295
            i=i+3
         ENDDO
 295     CONTINUE
         IF(i.GE.POS_ID)THEN
            IF(FUTURE_NIV2(MYID+1).NE.0)THEN
               WRITE(*,*)MYID,': ',SON,'has not been found
     & in DMUMPS_LOAD_COMP_MAXMEM_POOL'
               CALL MUMPS_ABORT()
            ENDIF
            GOTO 777
         ENDIF
         NSLAVES=CB_COST_ID(i+1)
         POS=CB_COST_ID(i+2)
         DO i=1,NSLAVES
            SLAVE=int(CB_COST_MEM(POS))
            IF(.NOT.CONCERNED(SLAVE))THEN
               MEM_ON_PROCS(SLAVE)=MEM_ON_PROCS(SLAVE)+
     &              dble(CB_COST_MEM(POS+1))
            ENDIF
            DO J=0,NPROCS-1
               IF(CONCERNED(J))THEN
                  IF(SLAVE.NE.J)THEN
                     RECV_BUF(J)=max(RECV_BUF(J),
     &                    dble(CB_COST_MEM(POS+1)))
                  ENDIF
               ENDIF
            ENDDO
            POS=POS+2
         ENDDO
 777     CONTINUE
         SON=FRERE_LOAD(STEP_LOAD(SON))
      ENDDO
      MAX_MEM=huge(MAX_MEM)
      WRITE(*,*)'NPROCS=',NPROCS,MAX_MEM
      DO i=0,NPROCS-1
         IF(MAX_MEM.GT.MEM_ON_PROCS(i))THEN
            PROC=i
         ENDIF
         MAX_MEM=min(MEM_ON_PROCS(i),MAX_MEM)
      ENDDO
      DEALLOCATE(MEM_ON_PROCS)
      DEALLOCATE(CONCERNED)
      DEALLOCATE(RECV_BUF)
      END SUBROUTINE DMUMPS_LOAD_COMP_MAXMEM_POOL
      SUBROUTINE DMUMPS_FIND_BEST_NODE_FOR_MEM(MIN_PROC,POOL,
     &                      LPOOL,INODE)
      IMPLICIT NONE
      INTEGER INODE,LPOOL,MIN_PROC
      INTEGER POOL(LPOOL)
      EXTERNAL MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      INTEGER i,NBTOP,INSUBTREE,NBINSUBTREE,NODE,FATHER,SON,J
      INTEGER SBTR_NB_LEAF,POS,K,allocok,L
      INTEGER, ALLOCATABLE, DIMENSION (:) ::  TMP_SBTR
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF((KEEP_LOAD(47).EQ.4).AND.
     &     ((NBINSUBTREE.NE.0)))THEN
         DO J=INDICE_SBTR,NB_SUBTREES
            NODE=MY_ROOT_SBTR(J) 
            FATHER=DAD_LOAD(STEP_LOAD(NODE))
            i=FATHER
 110        CONTINUE
            IF ( i > 0 ) THEN
               i = FILS_LOAD(i)
               GOTO 110
            ENDIF
            SON=-i
            i=SON
 120        CONTINUE
            IF ( i > 0 ) THEN
               IF( MUMPS_PROCNODE(PROCNODE_LOAD(STEP_LOAD(i)),
     &             KEEP_LOAD(199)) .EQ.  MIN_PROC ) THEN
                  SBTR_NB_LEAF=MY_NB_LEAF(J)
                  POS=SBTR_FIRST_POS_IN_POOL(J)
                  IF(POOL(POS+SBTR_NB_LEAF).NE.MY_FIRST_LEAF(J))THEN
                     WRITE(*,*)MYID,': The first leaf is not ok'
                     CALL MUMPS_ABORT()
                  ENDIF
                  ALLOCATE (TMP_SBTR(SBTR_NB_LEAF), stat=allocok)
                  IF (allocok > 0 ) THEN
                     WRITE(*,*)MYID,': Not enough space
     &                    for allocation'
                     CALL MUMPS_ABORT()
                  ENDIF
                  POS=SBTR_FIRST_POS_IN_POOL(J)
                  DO K=1,SBTR_NB_LEAF
                     TMP_SBTR(K)=POOL(POS+K-1)
                  ENDDO
                  DO K=POS+1,NBINSUBTREE-SBTR_NB_LEAF
                     POOL(K)=POOL(K+SBTR_NB_LEAF)
                  ENDDO
                  POS=1
                  DO K=NBINSUBTREE-SBTR_NB_LEAF+1,NBINSUBTREE
                     POOL(K)=TMP_SBTR(POS)
                     POS=POS+1
                  ENDDO
                  DO K=INDICE_SBTR,J
                     SBTR_FIRST_POS_IN_POOL(K)=SBTR_FIRST_POS_IN_POOL(K)
     &                    -SBTR_FIRST_POS_IN_POOL(J)
                  ENDDO
                  SBTR_FIRST_POS_IN_POOL(J)=NBINSUBTREE-SBTR_NB_LEAF
                  POS=MY_FIRST_LEAF(J)
                  L=MY_NB_LEAF(J)
                  DO K=INDICE_SBTR,J
                     MY_FIRST_LEAF(J)=MY_FIRST_LEAF(J+1)
                     MY_NB_LEAF(J)=MY_NB_LEAF(J+1)
                  ENDDO
                  MY_FIRST_LEAF(INDICE_SBTR)=POS
                  MY_NB_LEAF(INDICE_SBTR)=L
                  INODE=POOL(NBINSUBTREE)
                  DEALLOCATE(TMP_SBTR)
                  RETURN
               ENDIF
               i = FRERE_LOAD(STEP_LOAD(i))
               GOTO 120
            ENDIF           
         ENDDO
      ENDIF
      DO J=NBTOP,1,-1
         NODE=POOL(LPOOL-2-J)
         FATHER=DAD_LOAD(STEP_LOAD(NODE))
         i=FATHER
 11      CONTINUE
         IF ( i > 0 ) THEN
            i = FILS_LOAD(i)
            GOTO 11
         ENDIF
         SON=-i
         i=SON
 12      CONTINUE
         IF ( i > 0 ) THEN
            IF( MUMPS_PROCNODE(PROCNODE_LOAD(STEP_LOAD(i)),
     &          KEEP_LOAD(199)) .EQ. MIN_PROC ) THEN
               INODE=NODE
               RETURN
            ENDIF
            i = FRERE_LOAD(STEP_LOAD(i))
            GOTO 12
         ENDIF
      ENDDO
      END SUBROUTINE DMUMPS_FIND_BEST_NODE_FOR_MEM
      SUBROUTINE DMUMPS_LOAD_INIT_SBTR_STRUCT(POOL, LPOOL,KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER LPOOL,POOL(LPOOL),KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER i,POS
      EXTERNAL MUMPS_ROOTSSARBR
      LOGICAL MUMPS_ROOTSSARBR
      IF(.NOT.BDC_SBTR) RETURN
      POS=0
      DO i=NB_SUBTREES,1,-1
         DO WHILE(MUMPS_ROOTSSARBR(
     &            PROCNODE_LOAD(STEP_LOAD(POOL(POS+1))),
     &            KEEP(199)))
            POS=POS+1
         ENDDO
         SBTR_FIRST_POS_IN_POOL(i)=POS+1
         POS=POS+MY_NB_LEAF(i)
      ENDDO
      END SUBROUTINE DMUMPS_LOAD_INIT_SBTR_STRUCT
      END MODULE DMUMPS_LOAD 
      SUBROUTINE DMUMPS_SET_PARTI_REGULAR(
     &     SLAVEF,
     &     KEEP,KEEP8,
     &     PROCS,
     &     MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &     TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID,INODE,
     &     TAB_MAXS_ARG,SUP_PROC_ARG,MAX_SURF,NB_ROW_MAX
     &     )
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER(8), intent(in) :: TAB_MAXS_ARG(0:SLAVEF-1)
      INTEGER, intent(in) :: SUP_PROC_ARG(2)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE,NB_ROW_MAX
      INTEGER(8), intent(out):: MAX_SURF
      LOGICAL :: FORCE_LDLTRegular_NIV2
      INTEGER NSLAVES,ACC
      INTEGER i,J,NELIM,NB_SUP,K50,NB_ROWS(PROCS(SLAVEF+1))
      INTEGER TMP_NROW,X,K
      LOGICAL SUP,MEM_CSTR
      DOUBLE PRECISION MAX_LOAD,TOTAL_LOAD,VAR,TMP,A,B,C,DELTA,
     &     LOAD_CORR
      INTEGER IDWLOAD(SLAVEF)
      INTEGER(8) MEM_CONSTRAINT(2)
      K50=KEEP(50)
      FORCE_LDLTRegular_NIV2 = .FALSE.
      MAX_SURF=0
      NB_ROW_MAX=0
      NELIM=NFRONT-NCB
      NB_SUP=0
      TOTAL_LOAD=0.0D0
      SUP=.FALSE.
      IF(SUP_PROC_ARG(1).NE.
     &     0)THEN
         MEM_CONSTRAINT(1)=TAB_MAXS_ARG(PROCS(1))
         TOTAL_LOAD=TOTAL_LOAD+dble(SUP_PROC_ARG(1))/100.0D0
         NB_SUP=NB_SUP+1
      ENDIF
      IF(SUP_PROC_ARG(2).NE.
     &     0)THEN
         MEM_CONSTRAINT(2)=TAB_MAXS_ARG(PROCS(PROCS(SLAVEF+1)))
         TOTAL_LOAD=TOTAL_LOAD+dble(SUP_PROC_ARG(2))/100.0D0
         NB_SUP=NB_SUP+1
      ENDIF
      TOTAL_LOAD=TOTAL_LOAD+(PROCS(SLAVEF+1)-NB_SUP)
      IF(K50.EQ.0)THEN
         MAX_LOAD=dble( NELIM ) * dble ( NCB ) +
     *        dble(NCB) * dble(NELIM)*dble(2*NFRONT-NELIM-1)
      ELSE
         MAX_LOAD=dble(NELIM) * dble ( NCB ) *
     *        dble(NFRONT+1)
      ENDIF
      TMP=min(MAX_LOAD,MAX_LOAD/TOTAL_LOAD)
      J=1
      DO i=1,PROCS(SLAVEF+1)
         IF((NB_SUP.GT.0).AND.(i.EQ.1))THEN
            CYCLE
         ELSEIF((NB_SUP.EQ.2).AND.(i.EQ.PROCS(SLAVEF+1)))THEN
            CYCLE
         ENDIF
         IDWLOAD(J)=PROCS(i)
         J=J+1
      ENDDO
      DO i=1,NB_SUP
         IF(i.EQ.1)THEN
            IDWLOAD(J)=PROCS(1)
         ELSE
            IDWLOAD(J)=PROCS(PROCS(SLAVEF+1))
         ENDIF
         J=J+1
      ENDDO
      IF ((K50.EQ.0).OR.FORCE_LDLTRegular_NIV2) THEN
         ACC=0
         J=PROCS(SLAVEF+1)-NB_SUP+1
         DO i=1,NB_SUP
            VAR=dble(SUP_PROC_ARG(i))/100.0D0
            TMP_NROW=int(dble(MEM_CONSTRAINT(i))/dble(NFRONT))
            NB_ROWS(J)=int(max((VAR*dble(TMP))/
     &           (dble(NELIM)*dble(2*NFRONT-NELIM)),
     &           dble(1)))
            IF(NB_ROWS(J).GT.TMP_NROW)THEN
               NB_ROWS(J)=TMP_NROW
            ENDIF
            IF(NCB-ACC.LT.NB_ROWS(J)) THEN 
               NB_ROWS(J)=NCB-ACC
               ACC=NCB
               EXIT
            ENDIF
            ACC=ACC+NB_ROWS(J)
            J=J+1
         ENDDO
         IF(ACC.EQ.NCB)THEN
            GOTO 777
         ENDIF
         DO i=1,PROCS(SLAVEF+1)-NB_SUP
            VAR=1.0D0
            TMP_NROW=int((dble(TAB_MAXS_ARG(IDWLOAD(i))))/dble(NFRONT))
            NB_ROWS(i)=int((dble(VAR)*dble(TMP))/
     &           (dble(NELIM)*dble(2*NFRONT-NELIM)))
            IF(NB_ROWS(i).GT.TMP_NROW)THEN
               NB_ROWS(i)=TMP_NROW
            ENDIF
            IF(NCB-ACC.LT.NB_ROWS(i)) THEN 
               NB_ROWS(i)=NCB-ACC
               ACC=NCB
               EXIT
            ENDIF
            ACC=ACC+NB_ROWS(i)
         ENDDO
         IF(ACC.NE.NCB)THEN
            IF(PROCS(SLAVEF+1).EQ.NB_SUP)THEN
               TMP_NROW=(NCB-ACC)/PROCS(SLAVEF+1)+1
               DO i=1,PROCS(SLAVEF+1)
                  NB_ROWS(i)=NB_ROWS(i)+TMP_NROW
                  IF(ACC+TMP_NROW.GT.NCB)THEN
                     NB_ROWS(i)=NB_ROWS(i)-TMP_NROW+NCB-ACC
                     ACC=NCB
                     EXIT
                  ENDIF
                  ACC=ACC+TMP_NROW
               ENDDO
            ELSE
               TMP_NROW=(NCB-ACC)/(PROCS(SLAVEF+1)-NB_SUP)+1
               DO i=1,PROCS(SLAVEF+1)-NB_SUP
                  NB_ROWS(i)=NB_ROWS(i)+TMP_NROW
                  ACC=ACC+TMP_NROW
                  IF(ACC.GT.NCB) THEN 
                     NB_ROWS(i)=NB_ROWS(i)-TMP_NROW+
     &                    (NCB-(ACC-TMP_NROW))
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSE
         ACC=0
         i=PROCS(SLAVEF+1)-NB_SUP+1
         X=NCB
         LOAD_CORR=0.0D0
         MEM_CSTR=.FALSE.
         DO J=1,NB_SUP
            VAR=DBLE(SUP_PROC_ARG(J))/DBLE(100)
            A=1.0D0
            B=dble(X+NELIM)
            C=-dble(max(MEM_CONSTRAINT(J),0_8))
            DELTA=((B*B)-(4*A*C))
            TMP_NROW=int((-B+sqrt(DELTA))/(2*A))
            A=dble(-NELIM)
            B=dble(NELIM)*(dble(-NELIM)+dble(2*(X+NELIM)+1))
            C=-(VAR*TMP)
            DELTA=(B*B-(4*A*C))
            NB_ROWS(i)=int((-B+sqrt(DELTA))/(2*A))
            IF(NB_ROWS(i).GT.TMP_NROW)THEN
               NB_ROWS(i)=TMP_NROW
               MEM_CSTR=.TRUE.               
            ENDIF
            IF(ACC+NB_ROWS(i).GT.NCB)THEN
               NB_ROWS(i)=NCB-ACC
               ACC=NCB
               X=0
               EXIT
            ENDIF
            X=X-NB_ROWS(i)
            ACC=ACC+NB_ROWS(i) 
            LOAD_CORR=LOAD_CORR+(dble(NELIM) * dble (NB_ROWS(i)) *
     *           dble(2*(X+NELIM) - NELIM - NB_ROWS(i) + 1))
            i=i+1
         ENDDO
         IF(ACC.EQ.NCB)THEN
            GOTO 777
         ENDIF
            IF((PROCS(SLAVEF+1).NE.NB_SUP).AND.MEM_CSTR)THEN
               TMP=(MAX_LOAD-LOAD_CORR)/(PROCS(SLAVEF+1)-NB_SUP)
            ENDIF
         X=ACC
         ACC=0
         DO i=1,PROCS(SLAVEF+1)-NB_SUP
            IF (KEEP(375) .EQ. 1) THEN 
              VAR=1.0D0
              A=dble(NELIM)
              B=dble(NELIM)*(dble(NELIM)+dble(2*ACC+1))
              C=-(VAR*TMP)
            ELSE    
              A=1.0D0
              B=dble(ACC+NELIM)
              C=-TMP
            ENDIF
            DELTA=((B*B)-(4*A*C))
            NB_ROWS(i)=int((-B+sqrt(DELTA))/(2*A))
            IF(NCB-ACC-X.LT.NB_ROWS(i))THEN
               NB_ROWS(i)=NCB-ACC-X
               ACC=NCB-X
               EXIT
            ENDIF
            ACC=ACC+NB_ROWS(i)
         ENDDO
         ACC=ACC+X
         IF(ACC.NE.NCB)THEN
            IF(PROCS(SLAVEF+1).EQ.NB_SUP)THEN
               TMP_NROW=(NCB-ACC)/PROCS(SLAVEF+1)+1
               DO i=1,PROCS(SLAVEF+1)
                  NB_ROWS(i)=NB_ROWS(i)+TMP_NROW
                  IF(ACC+TMP_NROW.GT.NCB)THEN
                     NB_ROWS(i)=NB_ROWS(i)-TMP_NROW+NCB-ACC
                     ACC=NCB
                     EXIT
                  ENDIF
                  ACC=ACC+TMP_NROW
               ENDDO
            ELSE
               NB_ROWS(PROCS(SLAVEF+1)-NB_SUP)=
     &              NB_ROWS(PROCS(SLAVEF+1)
     &              -NB_SUP)+NCB-ACC
            ENDIF
         ENDIF
      ENDIF
 777  CONTINUE
      NSLAVES=0
      ACC=1
      J=1
      K=1
      DO i=1,PROCS(SLAVEF+1)
         IF(NB_ROWS(i).NE.0)THEN
            SLAVES_LIST(J)=IDWLOAD(i)
            TAB_POS(J)=ACC
            ACC=ACC+NB_ROWS(i)
            NB_ROW_MAX=max(NB_ROW_MAX,NB_ROWS(i))
            IF(K50.EQ.0)THEN
               MAX_SURF=max(int(NB_ROWS(i),8)*int(NCB,8),int(0,8))
            ELSE
               MAX_SURF=max(int(NB_ROWS(i),8)*int(ACC,8),int(0,8))
            ENDIF
            NSLAVES=NSLAVES+1
            J=J+1
         ELSE
            SLAVES_LIST(PROCS(SLAVEF+1)-K+1)=IDWLOAD(i)
            K=K+1
         ENDIF
      ENDDO
      TAB_POS(SLAVEF+2) = NSLAVES
      TAB_POS(NSLAVES+1)= NCB+1
      NSLAVES_NODE=NSLAVES
      END SUBROUTINE DMUMPS_SET_PARTI_REGULAR
