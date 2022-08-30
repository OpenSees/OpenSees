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
      MODULE DMUMPS_FAC_ASM_MASTER_ELT_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC_ASM_NIV1_ELT( COMM_LOAD, ASS_IRECV,
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, INFO, ND, 
     &    FILS, FRERE, DAD, MAXFRW, root,
     &    OPASSW, OPELIW, PTRIST, PTLUST, PTRFAC, PTRAST, 
     &    STEP, PIMASTER, PAMASTER,PTRARW, 
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, LRLUSM,
     &    ICNTL, KEEP,KEEP8,DKEEP,INTARR,LINTARR,DBLARR,LDBLARR,
     &
     &    NSTK_S,PROCNODE_STEPS, SLAVEF, COMM,MYID,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &    )
!$    USE OMP_LIB
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_PTR,
     &                                    DMUMPS_DM_IS_DYNAMIC,
     &                                    DMUMPS_DM_FREE_BLOCK
      USE MUMPS_BUILD_SORT_INDEX_ELT_M
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_LR_CORE, ONLY : IS_FRONT_BLR_CANDIDATE,
     &                           DMUMPS_BLR_ASM_NIV1
      USE DMUMPS_LR_DATA_M, ONLY : DMUMPS_BLR_INIT_FRONT,
     &                             DMUMPS_BLR_SAVE_NFS4FATHER
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_ANA_LR,    ONLY : GET_CUT
      USE DMUMPS_LR_CORE,   ONLY : MAX_CLUSTER
      USE MUMPS_LR_COMMON,  ONLY : COMPUTE_BLR_VCS
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER IZERO 
      PARAMETER (IZERO=0)
      INTEGER N, NSTEPS
      INTEGER NELT
      INTEGER(8) LA, LRLU, LRLUS, LRLUSM, IPTRLU, POSFAC
      INTEGER KEEP(500), ICNTL(60)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER, INTENT(INOUT) :: INFO(2)
      INTEGER INODE,MAXFRW,
     &        IWPOSCB, COMP
      INTEGER, TARGET :: IWPOS, LIW
      INTEGER IDUMMY(1)
      INTEGER, PARAMETER :: LIDUMMY = 1
      INTEGER, TARGET :: IW(LIW)
      INTEGER(8), INTENT(IN) :: PTRARW(NELT+1), PTRAIW(NELT+1)
      INTEGER ITLOC(N+KEEP(253)),
     &        ND(KEEP(28)), PERM(N), 
     &        FILS(N), FRERE(KEEP(28)), DAD(KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28)),
     &              PAMASTER(KEEP(28))
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER ETATASS
      LOGICAL SON_LEVEL2
      DOUBLE PRECISION, TARGET :: A(LA)
      INTEGER, INTENT(IN) :: LRGROUPS(N)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      DOUBLE PRECISION DBLARR(LDBLARR)
      INTEGER INTARR(LINTARR)
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK_S(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      LOGICAL PACKED_CB, IS_CB_LR
      INTEGER, EXTERNAL :: MUMPS_TYPENODE
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INCLUDE 'mpif.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
!$    INTEGER :: NOMP
      INCLUDE 'mumps_headers.h'
      INTEGER LP, HS, HF
      LOGICAL LPOK
      INTEGER NBPANELS_L, NBPANELS_U
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFS4FATHER, ESTIM_NFS4FATHER_ATSON
      INTEGER IFATH
      INTEGER PARPIV_T1
      INTEGER(8) NFRONT8, LAELL8, KEEP8TMPCOPY
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,ISTCHK_CB_RIGHT,LSTK,LREQ
      INTEGER LREQ_OOC
      INTEGER :: SON_XXS, SON_XXLR
      INTEGER(8) LSTK8, SIZFR8
      LOGICAL :: IS_DYNAMIC_CB
      INTEGER(8) :: DYN_SIZE
      INTEGER SIZFI, NCB
      INTEGER NCOLS, NROWS, LDA_SON
      INTEGER NELIM, IORG, IBROT
      INTEGER :: J253
#if ! defined(ZERO_TRIANGLE)
      INTEGER(8) :: NUMROWS, JJ3
#endif
      INTEGER :: TOPDIAG
!$    INTEGER    :: CHUNK  
!$    INTEGER(8) :: CHUNK8 
      INTEGER(8) APOS, APOS2, LAPOS2
      INTEGER(8) POSELT, POSEL1, ICT12, ICT21
      INTEGER(8) IACHK
      INTEGER(8) JJ2
      INTEGER(8) :: JJ8, J18, J28
      INTEGER(8) :: AINPUT8, AII8
      INTEGER :: K1, K2, K3, KK, KK1
      INTEGER JPOS,ICT11, IJROW
      INTEGER Pos_First_NUMORG,NUMORG,IOLDPS,
     &        NUMELT, ELBEG
      INTEGER :: J
      INTEGER NSLAVES, NSLSON, NPIVS, NPIV_ANA, NPIV
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER ELTI
      INTEGER(8) :: SIZE_ELTI8
      INTEGER(8) :: II8
      INTEGER :: I
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER LRSTATUS  
      LOGICAL COMPRESS_PANEL, LR_ACTIVATED, COMPRESS_CB,
     &        OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER, POINTER :: SON_IWPOS, SON_LIW
      INTEGER, POINTER, DIMENSION(:) :: SON_IW
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SON_A
      INTEGER NCBSON
      LOGICAL SAME_PROC
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR
      INTEGER :: NB_BLR, NPARTSCB, NPARTSASS, MAXI_CLUSTER, 
     &           IBCKSZ2, MINSIZE
      INTRINSIC real
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      LOGICAL  MUMPS_INSSARBR, SSARBR
      EXTERNAL MUMPS_INSSARBR
      DOUBLE PRECISION FLOP1,FLOP1_EFF
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
!$    NOMP = OMP_GET_MAX_THREADS()
      LP      = ICNTL(1)
      LPOK    = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      NFS4FATHER = -1
      ETATASS    = 0  
      PACKED_CB  = .FALSE.
      IS_CB_LR   = .FALSE.
      IN         = INODE
      LEVEL = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
      IF (LEVEL.NE.1) THEN
       WRITE(*,*) 'INTERNAL ERROR 1 in DMUMPS_FAC_ASM_NIV1_ELT '
       CALL MUMPS_ABORT()
      END IF
      NSLAVES = 0
      HF = 6 + NSLAVES + KEEP(IXSZ)
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .ne. 0 ) THEN
        ELBEG  = FRT_PTR(INODE)
      ELSE
        ELBEG  = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NPIV_ANA=NUMORG
      NSTEPS = NSTEPS + 1
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .NE. 0) THEN
        DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         SON_IW => IW
         NASS = NASS + SON_IW(PIMASTER(STEP(ISON))+1+KEEP(IXSZ))
         ISON = FRERE(STEP(ISON))
        END DO
      ENDIF
      NFRONT = ND(STEP(INODE)) + NASS + KEEP(253)
      NASS1 = NASS + NUMORG
      CALL IS_FRONT_BLR_CANDIDATE(INODE, 1, ND(STEP(INODE)), 
     &       NUMORG, KEEP(486),
     &       KEEP(489), KEEP(490), KEEP(491), KEEP(492), 
     &       KEEP(20), KEEP(60), DAD(STEP(INODE)), KEEP(38),
     &       LRSTATUS, N, LRGROUPS)
      IF (DAD(STEP(INODE)).NE.0) THEN
        IF ( MUMPS_PROCNODE(PROCNODE_STEPS(STEP(DAD(STEP(INODE)))),
     &                                     KEEP(199) )
     &    .NE. MYID
     &   .AND.
     &   MUMPS_TYPENODE(PROCNODE_STEPS(STEP(DAD(STEP(INODE)))),
     &                                     KEEP(199))
     &    .EQ.1
     &    ) THEN
          IF (LRSTATUS.EQ.1 .OR. LRSTATUS.EQ.3) THEN
            LRSTATUS = LRSTATUS-1
          ENDIF
        ENDIF
      ENDIF
      COMPRESS_PANEL = (LRSTATUS.GE.2)
      COMPRESS_CB    = ((LRSTATUS.EQ.1).OR.
     &                  (LRSTATUS.EQ.3))
      LR_ACTIVATED   = (LRSTATUS.GT.0)
      IF (COMPRESS_CB.AND.(.NOT.COMPRESS_PANEL)) THEN
        COMPRESS_PANEL = .TRUE.
        LRSTATUS      = 3 
      ENDIF
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( .NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN 
        CALL DMUMPS_OOC_GET_PP_SIZES(KEEP(50), NFRONT, NFRONT, NASS1,
     &       NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC   
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL DMUMPS_COMPRE_NEW(N, KEEP(28),
     &        IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216), LRLUS,
     &        KEEP(IXSZ), COMP, DKEEP(97), MYID, SLAVEF,
     &        KEEP(199), PROCNODE_STEPS, DAD)
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress DMUMPS_FAC_ASM_NIV1_ELT'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          END IF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      END IF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .TRUE.
          CALL MUMPS_ELT_BUILD_SORT( 
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, IWPOS,
     &        IFSON, STEP, PIMASTER, PTRIST, PTRAIW, NELT, IW, LIW,
     &        INTARR, LINTARR, ITLOC, FILS, FRERE,
     &        KEEP,
     &        SON_LEVEL2, NIV1, INFO(1),
     &        DAD,PROCNODE_STEPS, SLAVEF,
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG,
     &        IDUMMY, LIDUMMY )
      IF (INFO(1).LT.0) GOTO 300
      IF (NFRONT_EFF.NE.NFRONT) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
           IF(MUMPS_IN_OR_ROOT_SSARBR(PROCNODE_STEPS(STEP(INODE)),
     &          KEEP(199)))THEN
              NPIV=NASS1-(NFRONT_EFF-ND(STEP(INODE)))
              CALL MUMPS_GET_FLOPS_COST(ND(STEP(INODE))+KEEP(253),
     &                                 NPIV,NPIV,
     &                                 KEEP(50),1,FLOP1)             
              NPIV=NPIV_ANA
              CALL MUMPS_GET_FLOPS_COST(ND(STEP(INODE))+KEEP(253),
     &                                 NPIV,NPIV,
     &                                 KEEP(50),1,FLOP1_EFF)
              CALL DMUMPS_LOAD_UPDATE(0,.FALSE.,FLOP1-FLOP1_EFF,
     &             KEEP,KEEP8)
           ENDIF
        IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF))
        NFRONT = NFRONT_EFF
        LREQ = HF + 2 * NFRONT + LREQ_OOC   
        ELSE
         IF (LPOK) THEN
          WRITE(LP,*) 
     &     ' ERROR 1 during ass_niv1_ELT', NFRONT, NFRONT_EFF
         ENDIF
         GOTO 270
        ENDIF
      ENDIF
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1.AND.
     &     OOCWRITE_COMPATIBLE_WITH_BLR) THEN
        CALL DMUMPS_OOC_PP_SET_PTR(KEEP(50),
     &       NBPANELS_L, NBPANELS_U, NASS1, 
     &       IOLDPS + HF + 2 * NFRONT, IW, LIW)
      ENDIF
      NCB   = NFRONT - NASS1
      MAXFRW = max0(MAXFRW, NFRONT)
      ICT11 = IOLDPS + HF - 1 + NFRONT
      CALL DMUMPS_SET_PARPIVT1 ( INODE, NFRONT, NASS1, KEEP, 
     &                           LR_ACTIVATED, PARPIV_T1) 
      NFRONT8=int(NFRONT,8)
      LAELL8 = NFRONT8 * NFRONT8
      IF(PARPIV_T1.NE.0) THEN
       LAELL8 = LAELL8+int(NASS1,8)
      ENDIF
      CALL DMUMPS_GET_SIZE_NEEDED 
     &                   (0, LAELL8, .FALSE.,
     &                    KEEP(1), KEEP8(1),
     &                    N,KEEP(28),IW,LIW,A,LA,
     &                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     &                    PTRIST,PTRAST,
     &                    STEP, PIMASTER,PAMASTER,KEEP(216),LRLUS,
     &                    KEEP(IXSZ), COMP, DKEEP(97), MYID,
     &                    SLAVEF, PROCNODE_STEPS, DAD, 
     &                    INFO(1), INFO(2))
      IF (INFO(1).LT.0) GOTO 490
      LRLU = LRLU - LAELL8
      LRLUS = LRLUS - LAELL8
      LRLUSM = min( LRLUS, LRLUSM )
      IF (KEEP(405).EQ.0) THEN
        KEEP8(69) = KEEP8(69) + LAELL8 
        KEEP8(68) = max(KEEP8(69), KEEP8(68))
      ELSE
!$OMP   ATOMIC CAPTURE
        KEEP8(69) = KEEP8(69) + LAELL8 
        KEEP8TMPCOPY = KEEP8(69)
!$OMP   END ATOMIC
!$OMP   ATOMIC UPDATE
        KEEP8(68) = max(KEEP8(68), KEEP8TMPCOPY)
!$OMP   END ATOMIC
      ENDIF
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL8
      SSARBR=MUMPS_INSSARBR(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
      CALL DMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &     LA-LRLUS, 
     &     0_8,
     &     LAELL8,  
     &     KEEP,KEEP8,
     &     LRLUS)
      IF (KEEP(405).EQ.0) KEEP(429)= KEEP(429)+1
#if defined(ZERO_TRIANGLE)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
!$        CHUNK8=int(KEEP(361),8)
!$OMP   PARALLEL DO PRIVATE(JJ8) SCHEDULE(STATIC, CHUNK8)
!$OMP&  IF (LAPOS2 - POSELT + 1_8 > int(KEEP(361),8) .AND. NOMP .GT. 1)
        DO JJ8 = POSELT, LAPOS2
           A( JJ8 ) = ZERO
        ENDDO
!$OMP   END PARALLEL DO
      ELSE
        TOPDIAG = max(KEEP(7), KEEP(8), KEEP(218))-1
        IF  (LR_ACTIVATED) THEN
         NULLIFY(BEGS_BLR)
         CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NFRONT-1), NASS1,
     &        NCB, LRGROUPS, NPARTSCB, 
     &        NPARTSASS, BEGS_BLR)
            NB_BLR = NPARTSASS + NPARTSCB
         CALL MAX_CLUSTER(BEGS_BLR,NB_BLR,MAXI_CLUSTER)
         DEALLOCATE(BEGS_BLR)
         CALL COMPUTE_BLR_VCS(KEEP(472), IBCKSZ2, KEEP(488), NASS1)
         MINSIZE = int(IBCKSZ2 / 2)
         TOPDIAG = max(2*MINSIZE + MAXI_CLUSTER-1,TOPDIAG)
        ENDIF
        NUMROWS = NFRONT8
!$        CHUNK = max(KEEP(360)/2,
!$   &          ( ((int(NUMROWS)+NOMP-1) / NOMP + 2) / 3) )
!$OMP     PARALLEL DO PRIVATE(APOS,JJ3) SCHEDULE(STATIC, CHUNK )  
!$OMP&    IF (int(NUMROWS - 1) .GT. KEEP(360) .AND. NOMP .GT. 1)
          DO JJ8 = 0_8, NUMROWS - 1_8
             APOS = POSELT + JJ8 * NFRONT8
             JJ3 = min( NFRONT8 - 1_8, JJ8 + TOPDIAG )
             A(APOS:APOS + JJ3) = ZERO
          ENDDO
!$OMP     END PARALLEL DO
      END IF
#endif
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      PTLUST(STEP(INODE)) = IOLDPS
      IW(IOLDPS+XXI) = LREQ  
      CALL MUMPS_STOREI8(LAELL8,IW(IOLDPS+XXR))
      CALL MUMPS_STOREI8(0_8,IW(IOLDPS+XXD))
      IW(IOLDPS+XXS) = -9999
      IW(IOLDPS+XXN) = -99999
      IW(IOLDPS+XXP) = -99999
      IW(IOLDPS+XXA) = -99999
      IW(IOLDPS+XXF) = -99999
      IW(IOLDPS+XXLR) = LRSTATUS
      IW(IOLDPS + KEEP(IXSZ))   = NFRONT
      IW(IOLDPS + KEEP(IXSZ)+ 1) = 0
      IW(IOLDPS + KEEP(IXSZ) + 2) = -NASS1
      IW(IOLDPS + KEEP(IXSZ) + 3) = -NASS1
      IW(IOLDPS + KEEP(IXSZ) + 4) = STEP(INODE)
      IW(IOLDPS + KEEP(IXSZ) + 5)   = NSLAVES
      IF (LR_ACTIVATED.AND.
     &       (KEEP(480).NE.0
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       .OR.COMPRESS_CB
     &       )) THEN
         CALL DMUMPS_BLR_INIT_FRONT (IW(IOLDPS+XXF), INFO,
     &                               MTK405=KEEP(405))
         IF (INFO(1).LT.0) GOTO 500
      ENDIF
       ESTIM_NFS4FATHER_ATSON = -9999
       IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
         IFATH = DAD( STEP( INODE) )
         IF (IFATH.NE.0) THEN
           IF (COMPRESS_CB.AND.
     &      MUMPS_TYPENODE(PROCNODE_STEPS(STEP(IFATH)),KEEP(199))
     &      .EQ. 2 ) THEN
            IOLDPS = PTLUST(STEP(INODE))
            CALL DMUMPS_COMPUTE_ESTIM_NFS4FATHER (
     &        N, INODE, IFATH, FILS, PERM, KEEP, 
     &        IOLDPS, HF, IW, LIW, NFRONT, NASS1,
     &        ESTIM_NFS4FATHER_ATSON
     &         )
            CALL DMUMPS_BLR_SAVE_NFS4FATHER ( IW(IOLDPS+XXF),
     &              ESTIM_NFS4FATHER_ATSON )
           ENDIF
         ENDIF
      ENDIF
      IF (NUMSTK.NE.0) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK    = PIMASTER(STEP(ISON))
          SON_IW    => IW
          SON_LIW   => LIW
          SON_IWPOS => IWPOS
          SON_A     => A
          LSTK      = SON_IW(ISTCHK + KEEP(IXSZ))
          LSTK8     = int(LSTK,8)
          NELIM     = SON_IW(ISTCHK + KEEP(IXSZ) + 1)
          NPIVS     = SON_IW(ISTCHK + KEEP(IXSZ) + 3)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON    = SON_IW(ISTCHK + KEEP(IXSZ) + 5)
          HS        = 6 + KEEP(IXSZ) + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LT.SON_IWPOS)
          IF ( SAME_PROC ) THEN
            ISTCHK_CB_RIGHT = PTRIST(STEP(ISON))
          ELSE
            ISTCHK_CB_RIGHT = ISTCHK
          ENDIF
          SON_XXS  = SON_IW(ISTCHK_CB_RIGHT+XXS)
          SON_XXLR = SON_IW(ISTCHK_CB_RIGHT+XXLR)
          PACKED_CB = ( SON_XXS .EQ. S_CB1COMP )
          IS_CB_LR = ( SON_XXLR.EQ.1 .OR.  SON_XXLR.EQ.3 )
     &        .AND. (KEEP(489).EQ.1.OR.KEEP(489).EQ.3)
          LEVEL1    = NSLSON.EQ.0
          IF (.NOT.SAME_PROC) THEN
           NROWS = SON_IW( ISTCHK + KEEP(IXSZ) + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          SIZFI   = HS + NROWS + NCOLS 
          K1 = ISTCHK + HS + NROWS + NPIVS
          IF ( .NOT. LEVEL1 .AND. NELIM.EQ.0 ) GOTO 205
          IF (LEVEL1 .AND. .NOT. IS_CB_LR) THEN
           K2 = K1 + LSTK - 1
           IF (PACKED_CB) THEN
             SIZFR8 = (LSTK8*(LSTK8+1_8)/2_8)
           ELSE
             SIZFR8 = LSTK8*LSTK8
           ENDIF
          ELSE 
           IF ( KEEP(50).eq.0 ) THEN
             SIZFR8 = int(NELIM,8) * LSTK8
           ELSE
             IF (PACKED_CB) THEN
               SIZFR8 = int(NELIM,8) * int(NELIM+1,8)/2_8
             ELSE
               SIZFR8 = int(NELIM,8) * int(NELIM,8)
             ENDIF
           END IF
           K2 = K1 + NELIM - 1
          ENDIF
            IF (LEVEL1 .AND. .NOT. IS_CB_LR) THEN
              IF (KEEP(50).EQ.0) THEN
                OPASSW = OPASSW + LSTK8*LSTK8
              ELSE
                OPASSW = OPASSW + LSTK8*(LSTK8+1)/2_8
              ENDIF
            ELSE
              IF (KEEP(50).EQ.0) THEN
                OPASSW = OPASSW + int(NELIM,8)*LSTK8
              ELSE
                OPASSW = OPASSW + int(NELIM,8)*int(NELIM,8)/2_8
              ENDIF
            ENDIF
          CALL MUMPS_GETI8(DYN_SIZE, SON_IW(ISTCHK_CB_RIGHT+XXD))
          IS_DYNAMIC_CB = DYN_SIZE .GT. 0_8
          IF ( IS_DYNAMIC_CB ) THEN
            CALL DMUMPS_DM_SET_PTR( PAMASTER(STEP(ISON)), DYN_SIZE,
     &           SON_A )
            IACHK = 1_8 
          ELSE
            IACHK = PAMASTER(STEP(ISON))
          ENDIF
          IF (IS_CB_LR .AND. LEVEL1) THEN
            POSEL1 = PTRAST(STEP(INODE))
            CALL DMUMPS_BLR_ASM_NIV1 (A, LA,
     &           POSEL1, NFRONT, NASS1, SON_IW(ISTCHK+XXF),
     &           SON_IW, SON_LIW,
     &           LSTK, NELIM, K1, K1+LSTK-1, KEEP(50),
     &           KEEP, KEEP8, OPASSW)
          ENDIF
          IF ( KEEP(50) .eq. 0 ) THEN
            POSEL1 = PTRAST(STEP(INODE)) - NFRONT8
            IF (K2.GE.K1) THEN
              DO 170 KK = K1, K2
                APOS = POSEL1 + int(SON_IW(KK),8) * NFRONT8
                DO 160 KK1 = 1, LSTK
                  JJ2 = APOS + int(SON_IW(K1 + KK1 - 1) - 1,8)
                  A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1 - 1,8))
  160           CONTINUE
                IACHK = IACHK + LSTK8
  170         CONTINUE
            END IF
          ELSE
            IF (LEVEL1 .AND. .NOT. IS_CB_LR) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
              IF (SIZFR8 .GT. 0) THEN
                CALL DMUMPS_LDLT_ASM_NIV12(A, LA, SON_A(IACHK),
     &           PTRAST(STEP( INODE )), NFRONT, NASS1,
     &           LDA_SON, SIZFR8,
     &           SON_IW( K1 ), K2 - K1 + 1, NELIM, ETATASS, 
     &           PACKED_CB
!$   &           , KEEP(360)
     &          )
              ENDIF
          ENDIF
  205     IF (LEVEL1) THEN 
           IF (SAME_PROC) ISTCHK = PTRIST(STEP(ISON))
           IF (SAME_PROC) THEN
             IF (KEEP(50).NE.0) THEN
              K2 = K1 + LSTK - 1
              DO KK = K1, K2
               SON_IW(KK) = SON_IW(KK - NROWS)
              ENDDO
             ELSE
              K2 = K1 + LSTK - 1
              K3 = K1 + NELIM
              DO KK = K3, K2
               SON_IW(KK) = SON_IW(KK - NROWS)
              ENDDO
              IF (NELIM .NE. 0) THEN
                K3 = K3 - 1
                DO KK = K1, K3
                 JPOS = SON_IW(KK) + ICT11
                 SON_IW(KK) = IW(JPOS)
                ENDDO
              ENDIF
             ENDIF
           ENDIF
             IF ( SAME_PROC ) THEN
               PTRIST(STEP(ISON)) = -99999999
             ELSE
               PIMASTER(STEP( ISON )) = -99999999
             ENDIF
               CALL DMUMPS_FREE_BLOCK_CB_STATIC(
     &         SSARBR, MYID, N, ISTCHK,
     &         IW, LIW, LRLU, LRLUS, IPTRLU,
     &         IWPOSCB, LA, KEEP,KEEP8,
     &          .FALSE.
     &         )
             IF (IS_DYNAMIC_CB) THEN
               CALL DMUMPS_DM_FREE_BLOCK( SON_A, SIZFR8,
     &                                    KEEP(405).EQ.1, KEEP8 )
             ENDIF
          ELSE
           PDEST = ISTCHK + 6 + KEEP(IXSZ)
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP, KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
              CALL DMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV, 
     &             BUFR, LBUFR, LBUFR_BYTES,
     &             INODE, ISON, NSLAVES, IDUMMY,
     &             NFRONT, NASS1, NFS4FATHER,
     &             TROW_SIZE, IW( INDX ),  
     &         PROCNODE_STEPS,
     &         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &         LRLUS, N, IW, LIW, A, LA,
     &         PTRIST, PTLUST, PTRFAC, PTRAST, STEP,
     &         PIMASTER, PAMASTER, NSTK_S, COMP,
     &         INFO(1), INFO(2), MYID, COMM, PERM, IPOOL, LPOOL,
     &         LEAF, NBFIN, ICNTL, KEEP, KEEP8, DKEEP, root,
     &         OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD,
     &         PTRARW, PTRAIW,
     &         INTARR, DBLARR, ND, FRERE,
     &         NELT+1, NELT, FRT_PTR, FRT_ELT,
     &
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, LRGROUPS
     &         )
              IF ( INFO(1) .LT. 0 ) GOTO 500
               EXIT
             ENDIF
           ENDDO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  DMUMPS_BUF_SEND_MAPLIG( INODE, NFRONT, 
     &       NASS1, NFS4FATHER,ISON, MYID,
     &       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     &       COMM, IERR, IW(PDEST), NSLSON, 
     &       SLAVEF, 
     &       KEEP,KEEP8, STEP, N, 
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS,
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA,
     &         PTRIST, PTLUST, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         INFO(1), INFO(2), COMM,
     &         PERM,
     &         IPOOL, LPOOL, LEAF,
     &         NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, DAD, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ICNTL, KEEP, KEEP8,DKEEP, ND, FRERE,
     &         NELT+1, NELT, FRT_PTR, FRT_ELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE., LRGROUPS )
               IF ( INFO(1) .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
           ENDIF
          ENDIF
      ISON = FRERE(STEP(ISON))
  220 CONTINUE
      END IF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J18= PTRAIW(ELTI)
        J28= PTRAIW(ELTI+1)-1
        AII8 = PTRARW(ELTI)
        SIZE_ELTI8 = J28 - J18 + 1_8
        DO II8=J18,J28
         I = INTARR(II8)
         IF (KEEP(50).EQ.0) THEN
          AINPUT8    = AII8 + II8 - J18
          ICT12 = POSELT + int(I-1,8) * NFRONT8
          DO JJ8=J18,J28
           APOS2 = ICT12 + int(INTARR(JJ8) - 1,8)
           A(APOS2) = A(APOS2) + DBLARR(AINPUT8)
           AINPUT8 = AINPUT8 + SIZE_ELTI8
          END DO
         ELSE
          ICT12 = POSELT + int(- NFRONT + I - 1,8)
          ICT21 = POSELT + int(I-1,8)*NFRONT8 - 1_8
          DO JJ8=II8,J28
           J =  INTARR(JJ8)
           IF (I.LT.J) THEN
              APOS2 = ICT12 + int(J,8)*NFRONT8
           ELSE
              APOS2 = ICT21 + int(J,8)
           ENDIF
           A(APOS2) = A(APOS2) + DBLARR(AII8)
           AII8 = AII8 + 1_8
          END DO
         END IF
        END DO
      END DO
      IF (KEEP(253).GT.0) THEN
       POSELT = PTRAST(STEP(INODE))
       IBROT = INODE
       IJROW = Pos_First_NUMORG  
       DO IORG = 1, NUMORG
        IF (KEEP(50).EQ.0) THEN
          DO J253=1, KEEP(253)
            APOS = POSELT+
     &           int(IJROW-1,8)               * NFRONT8 +
     &           int(NFRONT-KEEP(253)+J253-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (J253-1)*KEEP(254)+IBROT )
          ENDDO
        ELSE
          DO J253=1, KEEP(253)
            APOS = POSELT+
     &           int(NFRONT-KEEP(253)+J253-1,8)  * NFRONT8 +
     &           int(IJROW-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (J253-1)*KEEP(254)+IBROT )
          ENDDO
        ENDIF
       IBROT = FILS(IBROT)
       IJROW = IJROW+1
       ENDDO
      ENDIF
      IF (PARPIV_T1.NE.0.AND.(.NOT.SON_LEVEL2)) THEN
        IOLDPS = PTLUST(STEP(INODE))
        CALL DMUMPS_PARPIVT1_SET_NVSCHUR_and_MAX (
     &     N, INODE, IW, LIW, A, LA, KEEP, PERM,
     &     IOLDPS, POSELT, 
     &     NFRONT, NASS1, LR_ACTIVATED, PARPIV_T1)
      ENDIF
      GOTO 500
  270 CONTINUE
      INFO(1) = -8
      INFO(2) = LREQ
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE IN INTEGER ALLOCATION DURING DMUMPS_ASM_NIV1_ELT'
      ENDIF
      GOTO 490
  290 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE, SEND BUFFER TOO SMALL DURING DMUMPS_ASM_NIV1_ELT'
      ENDIF
      INFO(1) = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE, RECV BUFFER TOO SMALL DURING DMUMPS_ASM_NIV1_ELT'
      ENDIF
      INFO(1) = -20
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF (INFO(1).EQ.-13) THEN
       IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION DURING DMUMPS_ASM_NIV1_ELT'
       ENDIF
       INFO(2)  = NUMSTK 
      ENDIF
  490 CONTINUE
      IF ( KEEP(405) .EQ. 0 ) THEN
        CALL  DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      ENDIF
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_ASM_NIV1_ELT
      SUBROUTINE DMUMPS_FAC_ASM_NIV2_ELT( COMM_LOAD, ASS_IRECV,
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, INFO,
     &    ND, FILS, FRERE, DAD,
     &    CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW, root,
     &    OPASSW, OPELIW, PTRIST, PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_S,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS, 
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP, KEEP8,DKEEP,INTARR,LINTARR,DBLARR,LDBLARR,
     &    PROCNODE_STEPS, SLAVEF, COMM,MYID,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, LEAF, IPOOL, LPOOL,
     &    PERM, MEM_DISTRIB
     &               , LRGROUPS
     &    )
!$    USE OMP_LIB
      USE MUMPS_BUILD_SORT_INDEX_ELT_M
      USE DMUMPS_BUF
      USE DMUMPS_LOAD
      USE DMUMPS_LR_CORE, ONLY : IS_FRONT_BLR_CANDIDATE
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : DMUMPS_DM_SET_PTR,
     &                                    DMUMPS_DM_IS_DYNAMIC
      USE DMUMPS_ANA_LR,    ONLY : GET_CUT
      USE DMUMPS_LR_CORE,   ONLY : MAX_CLUSTER
      USE MUMPS_LR_COMMON,  ONLY : COMPUTE_BLR_VCS
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N,LIW,NSTEPS, NBFIN
      INTEGER NELT
      INTEGER KEEP(500), ICNTL(60)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION       DKEEP(230)
      INTEGER(8) :: LRLUS, LRLU, IPTRLU, POSFAC, LA
      INTEGER, INTENT(INOUT) :: INFO(2)
      INTEGER INODE, MAXFRW, LPOOL, LEAF,
     &        IWPOS, IWPOSCB, COMP, SLAVEF
      DOUBLE PRECISION, TARGET :: A(LA)
      INTEGER, intent(in) :: LRGROUPS(N)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER(8), INTENT(IN) :: PTRARW(NELT+1), PTRAIW(NELT+1)
      INTEGER IW(LIW), ITLOC(N+KEEP(253)),
     &        ND(KEEP(28)),
     &        FILS(N), FRERE(KEEP(28)), DAD (KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     &        STEP(N),
     & PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), PERM(N)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER   CAND(SLAVEF+1, max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER(8), INTENT(IN) :: LINTARR,LDBLARR
      DOUBLE PRECISION DBLARR(LDBLARR)
      INTEGER INTARR(LINTARR)
      INTEGER MYID, COMM
      INTEGER IFATH
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mpif.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
!$    INTEGER :: NOMP
      INTEGER LP, HS, HF, HF_OLD, NCBSON, NSLAVES_OLD
      LOGICAL LPOK
      INTEGER NCBSON_MAX
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER :: IBC_SOURCE
      DOUBLE PRECISION, DIMENSION(:), POINTER :: SON_A
      INTEGER NFS4FATHER, ESTIM_NFS4FATHER_ATSON
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,ISTCHK_CB_RIGHT,LSTK,LREQ
      INTEGER :: SON_XXS
      INTEGER(8) :: LAELL8
      INTEGER LREQ_OOC
      INTEGER NBPANELS_L, NBPANELS_U
      LOGICAL PACKED_CB, IS_CB_LR
      INTEGER(8) :: LCB
      LOGICAL :: IS_DYNAMIC_CB
      INTEGER(8) :: DYN_SIZE
      INTEGER NCB
      INTEGER MP
      INTEGER :: K1, K2, KK, KK1
      INTEGER :: J253
      INTEGER(8) :: AII8, AINPUT8, II8
      INTEGER(8) :: J18,J28,JJ8
      INTEGER(8) :: LAPOS2, JJ2, JJ3
      INTEGER(8) :: NFRONT8, POSELT, POSEL1, LDAFS8,
     &           IACHK, ICT12, ICT21
      INTEGER(8) APOS, APOS2
#if ! defined(ZERO_TRIANGLE)
      INTEGER :: TOPDIAG
#endif
!$    INTEGER    :: CHUNK  
!$    INTEGER(8) :: CHUNK8  
      INTEGER NELIM,NPIVS,NCOLS,NROWS,
     &        IORG
      INTEGER LDAFS, LDA_SON, IJROW, IBROT
      INTEGER Pos_First_NUMORG,NBCOL,NUMORG,IOLDPS
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, PDEST
      INTEGER PDEST1(1)
      INTEGER :: ISLAVE
      INTEGER ELTI
      INTEGER(8) :: SIZE_ELTI8
      INTEGER :: I, J
      INTEGER :: ELBEG, NUMELT
      LOGICAL SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER LRSTATUS  
      LOGICAL COMPRESS_PANEL, LR_ACTIVATED, COMPRESS_CB,
     &        OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER IZERO
      INTEGER IDUMMY(1)
      PARAMETER( IZERO = 0 )
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION RZERO
      PARAMETER( RZERO = 0.0D0 )
      PARAMETER( ZERO = 0.0D0 )
      logical :: force_cand
      INTEGER ETATASS
      INTEGER(8) :: APOSMAX
      DOUBLE PRECISION  MAXARR
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok
      INTEGER NCB_SPLIT, SIZE_LIST_SPLIT, NBSPLIT,
     &        NUMORG_SPLIT, TYPESPLIT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST, COPY_CAND
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SONROWS_PER_ROW
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR
      INTEGER :: NB_BLR, NPARTSCB, NPARTSASS, MAXI_CLUSTER, 
     &           IBCKSZ2, MINSIZE
      INTEGER :: NMB_OF_CAND, NMB_OF_CAND_ORIG
      LOGICAL :: IS_ofType5or6, SPLIT_MAP_RESTART
!$    NOMP = OMP_GET_MAX_THREADS()
      MP      = ICNTL(2)
      LP      = ICNTL(1)
      LPOK    = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      PACKED_CB = .FALSE.
      ETATASS = 0  
      IN = INODE
      NSTEPS    = NSTEPS + 1
      KEEP(429) = KEEP(429)+1
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .NE. 0 ) THEN
        ELBEG = FRT_PTR(INODE)
      ELSE
        ELBEG = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      ENDDO
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      NCBSON_MAX = 0
      DO WHILE (ISON .GT. 0)
        NUMSTK = NUMSTK + 1
        IF ( KEEP(48)==5 .AND. 
     &       MUMPS_TYPENODE(PROCNODE_STEPS(STEP(ISON)),
     &       KEEP(199)) .EQ. 1) THEN
          NCBSON_MAX =
     &      max(NCBSON_MAX,IW(PIMASTER(STEP(ISON))+KEEP(IXSZ)))
        ENDIF
        NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 + KEEP(IXSZ))
        ISON = FRERE(STEP(ISON))
      ENDDO
      NFRONT = ND(STEP(INODE)) + NASS + KEEP(253)
      NASS1 = NASS + NUMORG
      NCB   = NFRONT - NASS1
      CALL IS_FRONT_BLR_CANDIDATE(INODE, 2, NFRONT, NASS1, KEEP(486),
     &       KEEP(489), KEEP(490), KEEP(491), KEEP(492), 
     &       KEEP(20), KEEP(60), DAD(STEP(INODE)), KEEP(38),
     &       LRSTATUS, N, LRGROUPS)
      COMPRESS_PANEL = (LRSTATUS.GE.2)
      COMPRESS_CB    = ((LRSTATUS.EQ.1).OR.
     &                  (LRSTATUS.EQ.3))
      LR_ACTIVATED   = (LRSTATUS.GT.0)
      IF (COMPRESS_CB.AND.(.NOT.COMPRESS_PANEL)) THEN
        COMPRESS_PANEL = .TRUE.
        LRSTATUS      = 3 
      ENDIF
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( .NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      IF((KEEP(24).eq.0).or.(KEEP(24).eq.1)) then
         force_cand=.FALSE.
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      end if
      TYPESPLIT =  MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              KEEP(199))
      IS_ofType5or6 =    (TYPESPLIT.EQ.5 .OR. TYPESPLIT.EQ.6)
      ISTCHK            = PIMASTER(STEP(IFSON))
      PDEST             = ISTCHK + 6 + KEEP(IXSZ)
      NSLSON            = IW(ISTCHK + KEEP(IXSZ) + 5)
      SPLIT_MAP_RESTART = .FALSE.
      IF (force_cand) THEN
         INIV2                = ISTEP_TO_INIV2( STEP( INODE ))
         NMB_OF_CAND          = CAND( SLAVEF+1, INIV2 )
         NMB_OF_CAND_ORIG     = NMB_OF_CAND
         SIZE_TMP_SLAVES_LIST = NMB_OF_CAND
         IF  (IS_ofType5or6) THEN
           DO I=NMB_OF_CAND+1,SLAVEF
            IF ( CAND( I, INIV2 ).LT.0) EXIT
            NMB_OF_CAND = NMB_OF_CAND +1
           ENDDO
           SIZE_TMP_SLAVES_LIST = NSLSON-1
          WRITE(6,*) "NMB_OF_CAND, SIZE_TMP_SLAVES_LIST ", 
     & NMB_OF_CAND, SIZE_TMP_SLAVES_LIST
           IF (INODE.EQ.-999999) THEN
              SPLIT_MAP_RESTART = .TRUE.
           ENDIF
         ENDIF
         IF (IS_ofType5or6.AND.SPLIT_MAP_RESTART) THEN
           TYPESPLIT     = 4
           IS_ofType5or6 = .FALSE.
           SIZE_TMP_SLAVES_LIST = NMB_OF_CAND 
           CAND (SLAVEF+1, INIV2) = SIZE_TMP_SLAVES_LIST
         ENDIF
      ELSE
         INIV2 = 1
         SIZE_TMP_SLAVES_LIST = SLAVEF - 1 
         NMB_OF_CAND          =  SLAVEF - 1
         NMB_OF_CAND_ORIG     =  SLAVEF - 1
      ENDIF
      ALLOCATE(TMP_SLAVES_LIST(SIZE_TMP_SLAVES_LIST),stat=allocok)
      IF (allocok > 0 ) THEN
        GOTO 265
      ENDIF
       TYPESPLIT =  MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              KEEP(199))
       IF  ( (TYPESPLIT.EQ.4) 
     &               .OR.(TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) 
     &     )  THEN
        IF (TYPESPLIT.EQ.4) THEN
         ALLOCATE(COPY_CAND(SLAVEF+1),stat=allocok)
         IF (allocok > 0 ) THEN
            GOTO 245
         ENDIF
         CALL DMUMPS_SPLIT_PREP_PARTITION (
     &      INODE, STEP, N, SLAVEF, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      CAND(1,INIV2), ICNTL, COPY_CAND,
     &      NBSPLIT, NUMORG_SPLIT, TMP_SLAVES_LIST(1),
     &      SIZE_TMP_SLAVES_LIST 
     &                                    )
         NCB_SPLIT = NCB-NUMORG_SPLIT
         SIZE_LIST_SPLIT = SIZE_TMP_SLAVES_LIST - NBSPLIT
         CALL DMUMPS_LOAD_SET_PARTITION( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, COPY_CAND,
     &     MEM_DISTRIB(0), NCB_SPLIT, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST(NBSPLIT+1),
     &     SIZE_LIST_SPLIT,INODE 
     &       ) 
         DEALLOCATE (COPY_CAND)
         CALL DMUMPS_SPLIT_POST_PARTITION (
     &      INODE, STEP, N, SLAVEF, NBSPLIT, NCB,
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      ICNTL, 
     &      TAB_POS_IN_PERE(1,INIV2),
     &      NSLAVES
     &                                    )
         IF (SPLIT_MAP_RESTART) THEN
          IS_ofType5or6 = .TRUE.
          TYPESPLIT =  MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              KEEP(199))
          CAND( SLAVEF+1, INIV2 ) = NMB_OF_CAND_ORIG
         ENDIF
        ELSE
         ISTCHK    = PIMASTER(STEP(IFSON))
         PDEST     = ISTCHK + 6 + KEEP(IXSZ)
         NSLSON    = IW(ISTCHK + KEEP(IXSZ) + 5)
         IF (KEEP(376) .EQ. 1) THEN
           NFRONT = IW( PIMASTER(STEP(IFSON)) + KEEP(IXSZ))
         ENDIF
         CALL DMUMPS_SPLIT_PROPAGATE_PARTI (
     &      INODE, TYPESPLIT, IFSON, 
     &      CAND(1,INIV2), NMB_OF_CAND_ORIG,
     &      IW(PDEST), NSLSON,
     &      STEP, N, SLAVEF, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      ICNTL, ISTEP_TO_INIV2, INIV2,
     &      TAB_POS_IN_PERE, NSLAVES, 
     &      TMP_SLAVES_LIST,
     &      SIZE_TMP_SLAVES_LIST
     &                                    )
        ENDIF
       ELSE
        CALL DMUMPS_LOAD_SET_PARTITION( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, CAND(1,INIV2),
     &     MEM_DISTRIB(0), NCB, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST,
     &     SIZE_TMP_SLAVES_LIST,INODE
     &       ) 
       ENDIF
      HF   = NSLAVES + 6 + KEEP(IXSZ)
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
        CALL DMUMPS_OOC_GET_PP_SIZES(KEEP(50), NASS1, NFRONT, NASS1,
     &                               NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL DMUMPS_COMPRE_NEW(N, KEEP(28),
     &        IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER,
     &        KEEP(216),LRLUS,KEEP(IXSZ),
     &        COMP, DKEEP(97), MYID, SLAVEF,
     &        KEEP(199), PROCNODE_STEPS, DAD)
          IF (LRLU .NE. LRLUS) THEN
           IF (LPOK) THEN
            WRITE(LP, * ) 'PB compress DMUMPS_FAC_ASM_NIV2_ELT',
     &                    'LRLU,LRLUS=',LRLU,LRLUS
           ENDIF
           GOTO 270
          ENDIF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      ENDIF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .FALSE.
      ALLOCATE(SONROWS_PER_ROW(NFRONT-NASS1), stat=allocok)
      IF (allocok > 0) THEN
         GOTO 275
      ENDIF
        CALL MUMPS_ELT_BUILD_SORT(
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, IWPOS,
     &        IFSON, STEP, PIMASTER, PTRIST, PTRAIW, NELT, IW, LIW,
     &        INTARR, LINTARR, ITLOC, FILS, FRERE,
     &        KEEP, SON_LEVEL2, NIV1, INFO(1),
     &        DAD,PROCNODE_STEPS, SLAVEF,
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG,
     &        SONROWS_PER_ROW, NFRONT - NASS1)
      IF (INFO(1).LT.0) GOTO 250
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (
     &        (TYPESPLIT.EQ.5) .OR. (TYPESPLIT.EQ.6)) THEN
          WRITE(6,*) ' Internal error 1 in fac_ass due to splitting ',
     &     ' INODE, NFRONT, NFRONT_EFF =', INODE, NFRONT, NFRONT_EFF 
          WRITE(6,*) ' SPLITTING NOT YET READY FOR THAT'
          CALL MUMPS_ABORT()
        ENDIF
        IF (NFRONT.GT.NFRONT_EFF) THEN
            NCB    = NFRONT_EFF - NASS1
            NSLAVES_OLD = NSLAVES
            HF_OLD      = HF
            IF (TYPESPLIT.EQ.4) THEN
             ALLOCATE(COPY_CAND(SLAVEF+1),stat=allocok)
             IF (allocok > 0 ) THEN
                GOTO 245
             ENDIF
             CALL DMUMPS_SPLIT_PREP_PARTITION (
     &          INODE, STEP, N, SLAVEF, 
     &          PROCNODE_STEPS, KEEP, DAD, FILS,
     &          CAND(1,INIV2), ICNTL, COPY_CAND,
     &          NBSPLIT, NUMORG_SPLIT, TMP_SLAVES_LIST(1),
     &          SIZE_TMP_SLAVES_LIST 
     &                                    )
             NCB_SPLIT = NCB-NUMORG_SPLIT
             SIZE_LIST_SPLIT = SIZE_TMP_SLAVES_LIST - NBSPLIT
             CALL DMUMPS_LOAD_SET_PARTITION( NCBSON_MAX, 
     &         SLAVEF, KEEP,KEEP8,
     &         ICNTL, COPY_CAND,
     &         MEM_DISTRIB(0), NCB_SPLIT, NFRONT_EFF, NSLAVES,
     &         TAB_POS_IN_PERE(1,INIV2),
     &         TMP_SLAVES_LIST(NBSPLIT+1),
     &         SIZE_LIST_SPLIT,INODE 
     &       ) 
             DEALLOCATE (COPY_CAND)
             CALL DMUMPS_SPLIT_POST_PARTITION (
     &          INODE, STEP, N, SLAVEF, NBSPLIT, NCB,
     &          PROCNODE_STEPS, KEEP, DAD, FILS,
     &          ICNTL, 
     &          TAB_POS_IN_PERE(1,INIV2),
     &          NSLAVES
     &                                    )
            ELSE
             CALL DMUMPS_LOAD_SET_PARTITION( NCBSON_MAX,
     &       SLAVEF, KEEP, KEEP8, ICNTL,
     &       CAND(1,INIV2),
     &       MEM_DISTRIB(0), NCB, NFRONT_EFF, NSLAVES,
     &       TAB_POS_IN_PERE(1,INIV2),
     &       TMP_SLAVES_LIST, SIZE_TMP_SLAVES_LIST,INODE 
     &       ) 
            ENDIF
            HF = NSLAVES + 6 + KEEP(IXSZ)
            IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF)) -
     &                   (NSLAVES_OLD - NSLAVES)
            IF (NSLAVES_OLD .NE. NSLAVES) THEN
              IF (NSLAVES_OLD > NSLAVES) THEN
               DO KK=0,2*NFRONT_EFF-1
                 IW(IOLDPS+HF+KK)=IW(IOLDPS+HF_OLD+KK)
               ENDDO
              ELSE
               IF (IWPOS - 1 > IWPOSCB ) GOTO 270
               DO KK=2*NFRONT_EFF-1, 0, -1
                 IW(IOLDPS+HF+KK) = IW(IOLDPS+HF_OLD+KK)
               ENDDO
              END IF
            END IF
            NFRONT = NFRONT_EFF
            LREQ = HF + 2 * NFRONT + LREQ_OOC
        ELSE
          IF (LPOK) THEN
           WRITE(LP,*) ' INTERNAL ERROR 2 during ass_niv2'
          ENDIF
          GOTO 270
        ENDIF
      ENDIF
      NFRONT8=int(NFRONT,8)
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1.AND.
     &    OOCWRITE_COMPATIBLE_WITH_BLR) THEN
        CALL DMUMPS_OOC_PP_SET_PTR(KEEP(50),
     &       NBPANELS_L, NBPANELS_U, NASS1, 
     &       IOLDPS + HF + 2 * NFRONT, IW, LIW)
      ENDIF
      MAXFRW = max0(MAXFRW, NFRONT)
      PTLUST(STEP(INODE)) = IOLDPS
      IW(IOLDPS+KEEP(IXSZ))     = NFRONT
      IW(IOLDPS + 1+KEEP(IXSZ)) = 0
      IW(IOLDPS + 2+KEEP(IXSZ)) = -NASS1
      IW(IOLDPS + 3+KEEP(IXSZ)) = -NASS1
      IW(IOLDPS + 4+KEEP(IXSZ)) = STEP(INODE)
      IW(IOLDPS+5+KEEP(IXSZ)) = NSLAVES
      IW(IOLDPS+6+KEEP(IXSZ):IOLDPS+5+NSLAVES+KEEP(IXSZ))=
     &                     TMP_SLAVES_LIST(1:NSLAVES)
      ESTIM_NFS4FATHER_ATSON = -9999
      IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
        IFATH = DAD( STEP( INODE) )
        IF (IFATH.NE.0) THEN
         IF (COMPRESS_CB.AND.
     &    MUMPS_TYPENODE(PROCNODE_STEPS(STEP(IFATH)),KEEP(199))
     &    .EQ. 2 ) THEN
          IOLDPS = PTLUST(STEP(INODE))
          CALL DMUMPS_COMPUTE_ESTIM_NFS4FATHER (
     &      N, INODE, IFATH, FILS, PERM, KEEP, 
     &      IOLDPS, HF, IW, LIW, NFRONT, NASS1,
     &      ESTIM_NFS4FATHER_ATSON
     &         )
         ENDIF
        ENDIF
      ENDIF
        CALL DMUMPS_LOAD_MASTER_2_ALL(MYID, SLAVEF, COMM_LOAD,
     &     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &     NASS1, KEEP, KEEP8, IW(IOLDPS+6+KEEP(IXSZ)), NSLAVES,INODE)
      IF(KEEP(86).EQ.1)THEN
         IF(mod(KEEP(24),2).eq.0)THEN
            CALL DMUMPS_LOAD_SEND_MD_INFO(SLAVEF,
     &           CAND(SLAVEF+1,INIV2),
     &           CAND(1,INIV2),
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, 
     &           NSLAVES,INODE)
         ELSEIF((KEEP(24).EQ.0).OR.(KEEP(24).EQ.1))THEN
            CALL DMUMPS_LOAD_SEND_MD_INFO(SLAVEF,
     &           SLAVEF-1,
     &           TMP_SLAVES_LIST,
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, 
     &           NSLAVES,INODE)
         ENDIF
      ENDIF
      DEALLOCATE(TMP_SLAVES_LIST)
      IF (KEEP(50).EQ.0) THEN
        LAELL8 = int(NASS1,8) * NFRONT8
        LDAFS = NFRONT
        LDAFS8 = NFRONT8
      ELSE
        LAELL8 = int(NASS1,8)*int(NASS1,8)
        IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2)
     &     LAELL8 = LAELL8+int(NASS1,8)
        LDAFS = NASS1
        LDAFS8 = int(NASS1,8)
      ENDIF
      CALL DMUMPS_GET_SIZE_NEEDED 
     &                   (0, LAELL8, .FALSE.,
     &                    KEEP(1), KEEP8(1),
     &                    N,KEEP(28),IW,LIW,A,LA,
     &                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     &                    PTRIST,PTRAST,
     &                    STEP, PIMASTER,PAMASTER,KEEP(216),LRLUS,
     &                    KEEP(IXSZ), COMP, DKEEP(97), MYID,
     &                    SLAVEF, PROCNODE_STEPS, DAD, 
     &                    INFO(1), INFO(2))
      IF (INFO(1).LT.0) GOTO 490
      LRLU = LRLU - LAELL8
      LRLUS = LRLUS - LAELL8
      KEEP8(67) = min(LRLUS, KEEP8(67))
      KEEP8(69) = KEEP8(69) + LAELL8 
      KEEP8(68) = max(KEEP8(69), KEEP8(68))
      POSELT = POSFAC
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      POSFAC = POSFAC + LAELL8
      IW(IOLDPS+XXI)   = LREQ  
      CALL MUMPS_STOREI8(LAELL8,IW(IOLDPS+XXR))
      CALL MUMPS_STOREI8(0_8,IW(IOLDPS+XXD))
      IW(IOLDPS+XXS) =  -9999
      IW(IOLDPS+XXN) = -99999   
      IW(IOLDPS+XXP) = -99999   
      IW(IOLDPS+XXA) = -99999
      IW(IOLDPS+XXF) = -99999
      IW(IOLDPS+XXLR)= LRSTATUS
      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,LA-LRLUS,0_8,LAELL8,
     & KEEP,KEEP8,
     & LRLUS)
      POSEL1 = POSELT - LDAFS8
#if defined(ZERO_TRIANGLE)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
!$      CHUNK8 = int(KEEP(361),8)
!$OMP   PARALLEL DO PRIVATE(JJ8) SCHEDULE(STATIC, CHUNK8)
!$OMP&  IF (LAPOS2 - POSELT > int(KEEP(361),8) .AND. NOMP .GT. 1)
        DO JJ8 = POSELT, LAPOS2
          A(JJ8) = ZERO
        ENDDO
!$OMP   END PARALLEL DO
      ELSE
        TOPDIAG = max(KEEP(7), KEEP(8))-1
        IF  (LR_ACTIVATED) THEN
         NULLIFY(BEGS_BLR)
         CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NFRONT-1), NASS1,
     &        0, LRGROUPS, NPARTSCB, 
     &        NPARTSASS, BEGS_BLR)
         NB_BLR = NPARTSASS + NPARTSCB
         CALL MAX_CLUSTER(BEGS_BLR,NB_BLR,MAXI_CLUSTER)
         DEALLOCATE(BEGS_BLR)
         CALL COMPUTE_BLR_VCS(KEEP(472), IBCKSZ2, KEEP(488), NASS1)
         MINSIZE = int(IBCKSZ2 / 2)
         TOPDIAG = max(2*MINSIZE + MAXI_CLUSTER-1, TOPDIAG)
        ENDIF
!$      CHUNK = max(KEEP(360)/2,
!$   &       ( (( LDAFS + NOMP -1 ) / NOMP + 2) / 3) )
        APOS = POSELT
!$OMP   PARALLEL DO PRIVATE(APOS,JJ3) SCHEDULE(STATIC, CHUNK)
!$OMP&  IF (LDAFS - 1 .GT. KEEP(360) .AND. NOMP .GT. 1)
        DO JJ8 = 0_8, int(LDAFS-1,8)
          APOS = POSELT + JJ8 * int(LDAFS,8)
          JJ3 = min( int(LDAFS,8)  - 1_8, JJ8 + TOPDIAG )
          A(APOS:APOS+JJ3) = ZERO
        END DO
!$OMP   END PARALLEL DO
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
          APOSMAX = POSELT + int(NASS1,8)*int(NASS1,8)
          A(APOSMAX:APOSMAX+int(LDAFS-1,8))=ZERO
        ENDIF
      END IF
#endif
      IF ((NUMSTK.NE.0).AND.(NASS.NE.0)) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK = PIMASTER(STEP(ISON))
          NELIM = IW(ISTCHK + KEEP(IXSZ) + 1)
          IF (NELIM.EQ.0) GOTO 210
          LSTK    = IW(ISTCHK + KEEP(IXSZ))
          NPIVS   = IW(ISTCHK + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS=0
          NSLSON  = IW(ISTCHK + 5+KEEP(IXSZ))
          HS      = 6 + NSLSON + KEEP(IXSZ)
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LT.IWPOS)
          IF ( SAME_PROC ) THEN
            ISTCHK_CB_RIGHT=PTRIST(STEP(ISON))
          ELSE
            ISTCHK_CB_RIGHT=ISTCHK
          ENDIF
          SON_XXS = IW(ISTCHK_CB_RIGHT + XXS)
          PACKED_CB = ( SON_XXS .EQ. S_CB1COMP )
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW(ISTCHK + KEEP(IXSZ) + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          IF (KEEP(50).EQ.0) THEN
            LDA_SON = LSTK
            LCB = int(NELIM,8)*int(LSTK,8)
          ELSE 
            IF (NSLSON.EQ.0) THEN
             IF (SAME_PROC) THEN
               IS_CB_LR = IW(ISTCHK_CB_RIGHT+XXLR).EQ. 1 .OR.
     &                    IW(ISTCHK_CB_RIGHT+XXLR).EQ. 3
               IF (IS_CB_LR) THEN
                 LDA_SON = NELIM
               ELSE
                 LDA_SON = LSTK
               ENDIF
             ELSE 
               LDA_SON = LSTK
             ENDIF
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (PACKED_CB) THEN
              LCB = (int(NELIM,8)*int(NELIM+1,8))/2_8
            ELSE
              LCB = int(LDA_SON,8)*int(NELIM,8)
            ENDIF
          ENDIF
          IF (KEEP(50) .EQ. 0) THEN
            OPASSW = OPASSW + dble(LCB)
          ELSE
            OPASSW = OPASSW + int(NELIM,8)*int(NELIM+1,8)/2_8
          ENDIF
          IS_DYNAMIC_CB =
     &           DMUMPS_DM_IS_DYNAMIC(IW(ISTCHK_CB_RIGHT+XXD:
     &                                   ISTCHK_CB_RIGHT+XXD+1))
          IF ( IS_DYNAMIC_CB ) THEN
            CALL MUMPS_GETI8(DYN_SIZE, IW(ISTCHK_CB_RIGHT+XXD))
            CALL DMUMPS_DM_SET_PTR( PAMASTER(STEP(ISON)), DYN_SIZE,
     &                              SON_A )
            IACHK = 1_8 
          ELSE
            IACHK = PAMASTER(STEP(ISON))
            SON_A=>A
          ENDIF
          K1 = ISTCHK + HS + NROWS + NPIVS
          K2 = K1 + NELIM - 1
          IF (KEEP(50).eq.0) THEN
           IF (IS_ofType5or6) THEN
            APOS = POSELT  
            DO JJ8 = 1_8, int(NELIM,8)*int(LSTK,8)
             A(APOS+JJ8-1_8) = A(APOS+JJ8-1_8) + SON_A(IACHK+JJ8-1_8)
            ENDDO
           ELSE
            DO 170 KK = K1, K2
             APOS = POSEL1 + int(IW(KK),8) * LDAFS8
             DO 160 KK1 = 1, LSTK
              JJ2 = APOS + int(IW(K1 + KK1 - 1),8) - 1_8
              A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1 - 1,8))
  160        CONTINUE
             IACHK = IACHK + int(LSTK,8)
  170       CONTINUE
           ENDIF
          ELSE
            IF (LCB .GT. 0) THEN
              CALL DMUMPS_LDLT_ASM_NIV12(A, LA, SON_A(IACHK),
     &           POSELT, LDAFS, NASS1,
     &           LDA_SON, LCB,
     &           IW( K1 ), NELIM, NELIM, ETATASS,
     &           PACKED_CB
!$   &           , KEEP(360)
     &          )
            ENDIF
          ENDIF
  210     ISON = FRERE(STEP(ISON))
  220   CONTINUE
      ENDIF
      APOSMAX = POSELT + int(NASS1,8)*int(NASS1,8)
      IF (KEEP(219).NE.0) THEN
        IF (KEEP(50).EQ.2) THEN
          A( APOSMAX: APOSMAX+int(NASS1-1,8))=ZERO
        ENDIF
      ENDIF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J18= PTRAIW(ELTI)
        J28= PTRAIW(ELTI+1) - 1_8
        AII8 = PTRARW(ELTI)
        SIZE_ELTI8 = J28 - J18 + 1_8
        DO II8=J18,J28
         I = INTARR(II8)
         IF (KEEP(50).EQ.0) THEN
          IF (I.LE.NASS1) THEN
           AINPUT8   = AII8 + II8 - J18
           ICT12 = POSELT + int(I-1,8) * LDAFS8
           DO JJ8=J18,J28
            APOS2 = ICT12 + int(INTARR(JJ8) - 1,8)
            A(APOS2) = A(APOS2) + DBLARR(AINPUT8)
            AINPUT8 = AINPUT8 + SIZE_ELTI8
           END DO
          ENDIF
         ELSE
          ICT12 = POSELT - LDAFS8 + int(I,8) - 1_8
          ICT21 = POSELT + int(I-1,8)*LDAFS8 - 1_8
          IF ( I .GT. NASS1 ) THEN
           IF (KEEP(219).NE.0 .AND. KEEP(50).EQ.2) THEN
              AINPUT8=AII8
              DO JJ8=II8,J28
               J=INTARR(JJ8)
               IF (J.LE.NASS1) THEN
                A(APOSMAX+int(J-1,8))=
     &              max(dble(A(APOSMAX+int(J-1,8))),
     &                  abs(DBLARR(AINPUT8)))
               ENDIF
               AINPUT8=AINPUT8+1_8
              ENDDO
           ENDIF
           AII8 = AII8 + J28 - II8 + 1_8
           CYCLE
          ELSE
            IF (KEEP(219).NE.0) THEN
              MAXARR = RZERO
            ENDIF
            DO JJ8=II8,J28
              J =  INTARR(JJ8)
              IF ( J .LE. NASS1) THEN
                IF (I.LT.J) THEN
                  APOS2 = ICT12 + int(J,8)*LDAFS8
                ELSE
                  APOS2 = ICT21 + int(J,8)
                ENDIF
                A(APOS2) = A(APOS2) + DBLARR(AII8)
              ELSE IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                MAXARR = max(MAXARR,abs(DBLARR(AII8)))
              ENDIF
              AII8 = AII8 + 1_8
            END DO
            IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
                A(APOSMAX+int(I-1,8)) = 
     &             max( MAXARR, dble(A(APOSMAX+int(I-1,8))))
            ENDIF
          ENDIF 
         END IF 
        END DO
      END DO
      IF (KEEP(253).GT.0) THEN
       POSELT = PTRAST(STEP(INODE))
       IBROT = INODE
       IJROW = Pos_First_NUMORG  
       DO IORG = 1, NUMORG
        IF (KEEP(50).EQ.0) THEN
          DO J253 = 1, KEEP(253)
            APOS = POSELT +
     &             int(IJROW-1,8) * int(LDAFS,8) +
     &             int(LDAFS-KEEP(253)+J253-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (J253-1)*KEEP(254)+IBROT )
          ENDDO
        ENDIF
        IBROT = FILS(IBROT)
        IJROW = IJROW+1
       ENDDO
      ENDIF
      PTRCOL = IOLDPS + HF + NFRONT 
      PTRROW = IOLDPS + HF + NASS1 
      PDEST  = IOLDPS + 6 + KEEP(IXSZ)
      IBC_SOURCE = MYID
      DO ISLAVE = 1, NSLAVES
              CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB,
     &                NSLAVES, 
     &                NBLIG, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
        IERR = -1
        DO WHILE (IERR .EQ.-1)
         IF ( KEEP(50) .eq. 0 ) THEN
           NBCOL =  NFRONT
           CALL DMUMPS_BUF_SEND_DESC_BANDE( INODE,
     &      sum(SONROWS_PER_ROW(FIRST_INDEX:FIRST_INDEX+NBLIG-1)),
     &      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     &      IZERO, IDUMMY, 
     &      ESTIM_NFS4FATHER_ATSON,   
     &      IW(PDEST), IBC_SOURCE, NFRONT, COMM, KEEP, IERR
     &      , IW(IOLDPS+XXLR)
     &      )
         ELSE
           NBCOL = NASS1+SHIFT_INDEX+NBLIG
           CALL DMUMPS_BUF_SEND_DESC_BANDE( INODE,
     &      sum(SONROWS_PER_ROW(FIRST_INDEX:FIRST_INDEX+NBLIG-1)),
     &      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     &      NSLAVES-ISLAVE,
     &      IW( PTLUST(STEP(INODE))+6+KEEP(IXSZ)+ISLAVE),
     &      ESTIM_NFS4FATHER_ATSON,   
     &      IW(PDEST), IBC_SOURCE, NFRONT, COMM, KEEP, IERR
     &      , IW(IOLDPS+XXLR)
     &      )
         ENDIF
         IF (IERR.EQ.-1) THEN
          BLOCKING  = .FALSE.
          SET_IRECV = .TRUE.
          MESSAGE_RECEIVED = .FALSE.
          CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &     MPI_ANY_SOURCE, MPI_ANY_TAG,
     &     STATUS, BUFR, LBUFR,
     &     LBUFR_BYTES,
     &     PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &     LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &     PTLUST, PTRFAC,
     &     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1),
     &     INFO(2), COMM,
     &     PERM,
     &     IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &     root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &     FILS, DAD, PTRARW, PTRAIW,
     &     INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &     NELT+1, NELT, FRT_PTR, FRT_ELT,
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &       )
          IF ( INFO(1) .LT. 0 ) GOTO 500
          IF (MESSAGE_RECEIVED) THEN
           IOLDPS = PTLUST(STEP(INODE))
           PTRCOL = IOLDPS + HF + NFRONT
           PTRROW = IOLDPS + HF + NASS1 + SHIFT_INDEX
          ENDIF
         ENDIF
        ENDDO
        IF (IERR .EQ. -2) GOTO 300
        IF (IERR .EQ. -3) GOTO 305
        PTRROW = PTRROW + NBLIG
        PDEST  = PDEST + 1
      ENDDO
      DEALLOCATE(SONROWS_PER_ROW)
      IF (NUMSTK.EQ.0) GOTO 500
      ISON = IFSON
      DO IELL = 1, NUMSTK
        ISTCHK = PIMASTER(STEP(ISON))
        NELIM = IW(ISTCHK + 1 + KEEP(IXSZ))
        LSTK    = IW(ISTCHK + KEEP(IXSZ))
        NPIVS   = IW(ISTCHK + 3 + KEEP(IXSZ))
        IF ( NPIVS .LT. 0 ) NPIVS = 0
        NSLSON  = IW(ISTCHK + 5 + KEEP(IXSZ))
        HS      = 6 + NSLSON + KEEP(IXSZ)
        NCOLS     = NPIVS + LSTK
        SAME_PROC     = (ISTCHK.LT.IWPOS)
        IF (.NOT.SAME_PROC) THEN
         NROWS = IW(ISTCHK + 2 + KEEP(IXSZ) )
        ELSE
         NROWS = NCOLS
        ENDIF
        PDEST   = ISTCHK + 6 + KEEP(IXSZ)
        NCBSON  = LSTK - NELIM
        PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
           NFS4FATHER = NCBSON
           DO I=0,NCBSON-1
              IF(IW(PTRCOL+I) .GT. NASS1) THEN
                 NFS4FATHER = I
                 EXIT
              ENDIF
           ENDDO
           NFS4FATHER = NFS4FATHER + NELIM
        ELSE
          NFS4FATHER = 0
        ENDIF
        IF (NSLSON.EQ.0) THEN
          NSLSON = 1
          PDEST1(1)  = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(ISON)),
     &                                KEEP(199))
          IF (PDEST1(1).EQ.MYID) THEN
            CALL DMUMPS_MAPLIG_FILS_NIV1( COMM_LOAD, ASS_IRECV, 
     &      BUFR, LBUFR, LBUFR_BYTES,
     &      INODE, ISON, NSLAVES, 
     &      IW( PTLUST(STEP(INODE)) + 6 +KEEP(IXSZ)),
     &      NFRONT, NASS1, NFS4FATHER, NCBSON, IW( PTRCOL ),
     &      PROCNODE_STEPS,
     &      SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &      LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST, PTRFAC, PTRAST, STEP,
     &      PIMASTER, PAMASTER, NSTK_S, COMP,
     &      INFO(1), INFO(2), MYID, COMM, PERM,
     &      IPOOL, LPOOL, LEAF,
     &      NBFIN, ICNTL, KEEP, KEEP8, DKEEP, root,
     &      OPASSW, OPELIW,
     &      ITLOC, RHS_MUMPS, FILS, DAD,
     &      PTRARW, PTRAIW, INTARR, DBLARR, ND, FRERE, NELT+1, NELT, 
     &      FRT_PTR, FRT_ELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &      LRGROUPS )
           IF ( INFO(1) .LT. 0 ) GOTO 500
          ELSE
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            CALL  DMUMPS_BUF_SEND_MAPLIG( 
     &           INODE, NFRONT,NASS1,NFS4FATHER,
     &           ISON, MYID,
     &      NSLAVES, IW( PTLUST(STEP(INODE)) + 6 +KEEP(IXSZ)),
     &      IW(PTRCOL), NCBSON,
     &      COMM, IERR, PDEST1, NSLSON, SLAVEF, 
     &      KEEP,KEEP8, STEP, N, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &        MPI_ANY_SOURCE, MPI_ANY_TAG,
     &        STATUS, BUFR, LBUFR, LBUFR_BYTES,
     &        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &        LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &        PTLUST, PTRFAC,
     &        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1),
     &        INFO(2), COMM,
     &        PERM,
     &        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &        root,OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD,
     &        PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &        )
              IF ( INFO(1) .LT. 0 ) GOTO 500
            ENDIF
           ENDDO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
          ENDIF
        ELSE
          IF (PIMASTER(STEP(ISON)).GT.0) THEN
          IERR = -1
          DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  DMUMPS_BUF_SEND_MAPLIG( 
     &           INODE, NFRONT, NASS1, NFS4FATHER,
     &           ISON, MYID,
     &      NSLAVES, IW(PTLUST(STEP(INODE))+6+KEEP(IXSZ)),
     &      IW(PTRCOL), NCBSON,
     &      COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     &      KEEP,KEEP8, STEP, N, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL DMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &        MPI_ANY_SOURCE, MPI_ANY_TAG,
     &        STATUS, BUFR, LBUFR,
     &        LBUFR_BYTES,
     &        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &        LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &        PTLUST, PTRFAC,
     &        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1),
     &        INFO(2), COMM,
     &        PERM,
     &        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &        root,OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &        FILS, DAD, PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &        )
             IF ( INFO(1) .LT. 0 ) GOTO 500
            ENDIF
          ENDDO
          IF (IERR .EQ. -2) GOTO 290
          IF (IERR .EQ. -3) GOTO 295
          ENDIF
          DO ISLAVE = 0, NSLSON-1
            IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
               CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON,
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX        = PTRCOL + SHIFT_INDEX
              CALL DMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV, 
     &        BUFR, LBUFR, LBUFR_BYTES,
     &        INODE, ISON, NSLAVES, 
     &        IW( PTLUST(STEP(INODE))+6+KEEP(IXSZ)),
     &        NFRONT, NASS1,NFS4FATHER,
     &        TROW_SIZE, IW( INDX ),
     &        PROCNODE_STEPS,
     &        SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &        LRLUS, N, IW, LIW, A, LA,
     &        PTRIST, PTLUST, PTRFAC, PTRAST, STEP,
     &        PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1), INFO(2),
     &        MYID, COMM, PERM, IPOOL, LPOOL, LEAF,
     &        NBFIN, ICNTL, KEEP,KEEP8,DKEEP, root,
     &        OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD,
     &        PTRARW, PTRAIW, INTARR, DBLARR, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     & 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, LRGROUPS)
              IF ( INFO(1) .LT. 0 ) GOTO 500
              EXIT
            ENDIF
          ENDDO
        ENDIF
       ISON = FRERE(STEP(ISON))
      ENDDO
      GOTO 500
  250 CONTINUE
      IF (INFO(1).EQ.-13) THEN
       IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE IN INTEGER DYNAMIC ALLOCATION DURING
     & DMUMPS_FAC_ASM_NIV2_ELT'
       ENDIF
       INFO(2)   = NUMSTK + 1
      ENDIF
      GOTO 490
  245 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING COPY_CAND',
     &                 ' DURING DMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1)  = -13
      INFO(2)  = SLAVEF+1
      GOTO 490
  265 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING TMP_SLAVES_LIST',
     &                 ' DURING DMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1)   = -13
      INFO(2)   = SIZE_TMP_SLAVES_LIST
      GOTO 490
  270 CONTINUE
      INFO(1) = -8
      INFO(2) = LREQ
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE IN INTEGER ALLOCATION DURING DMUMPS_ASM_NIV2_ELT'
      ENDIF
      GOTO 490
  275 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING SONROWS_PER_ROW',
     &                 ' DURING DMUMPS_ASM_NIV2_ELT'
      ENDIF
      INFO(1)  = -13
      INFO(2)  = NFRONT-NASS1
      GOTO 490
  290 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (1) DURING DMUMPS_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -17
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (1) DURING DMUMPS_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -20
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (2)',
     &' DURING DMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -17
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  305 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (2)',
     &' DURING DMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -20
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  490 CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_ASM_NIV2_ELT
      END MODULE DMUMPS_FAC_ASM_MASTER_ELT_M
