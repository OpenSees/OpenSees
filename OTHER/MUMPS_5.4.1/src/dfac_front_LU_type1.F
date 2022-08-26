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
      MODULE DMUMPS_FAC1_LU_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC1_LU(
     &                           N, INODE, IW, LIW, A, 
     &                           LA,
     &                           IOLDPS, POSELT, IFLAG, IERROR, UU, 
     &                           NOFFW, NPVW, NBTINYW,
     &                           DET_EXPW, DET_MANTW, DET_SIGNW,
     &                           KEEP, KEEP8, STEP,
     &                           PROCNODE_STEPS, MYID, SLAVEF, SEUIL,
     &                           AVOID_DELAYED, ETATASS,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, 
     &     IWPOS 
     &               , LRGROUPS
     &               , PERM
     &     )
      USE DMUMPS_FAC_FRONT_AUX_M
      USE DMUMPS_OOC
      USE DMUMPS_FAC_LR
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_ANA_LR, ONLY : GET_CUT
      USE DMUMPS_LR_DATA_M
#if defined(BLR_MT)          
#endif
!$    USE OMP_LIB
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, IERROR
      INTEGER, INTENT(INOUT) :: NOFFW, NPVW, NBTINYW
      INTEGER, INTENT(INOUT) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, INTENT(INOUT) :: DET_MANTW
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER MYID, SLAVEF, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) ), STEP(N)
      DOUBLE PRECISION UU, SEUIL
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(230)
      INTEGER :: LRGROUPS(N), PERM(N)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      INTEGER Inextpiv
      INTEGER LAST_ROW, LAST_COL, FIRST_COL
      LOGICAL CALL_LTRSM, CALL_UTRSM
      DOUBLE PRECISION UUTEMP
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      INTEGER PIVOT_OPTION
      INTEGER LRTRSM_OPTION
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, LNextPiv2beWritten, 
     &        UNextPiv2beWritten, IFLAG_OOC,
     &        PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_LastPIVRPTRFilled_U
      INTEGER TYPEF_LOC
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER PARPIV_T1
      INTEGER CURRENT_BLR
      LOGICAL LR_ACTIVATED
      LOGICAL COMPRESS_CB, COMPRESS_PANEL
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR,
     &        OOC_EFFECTIVE_ON_FRONT, 
     &        OOC_EFF_AND_WRITE_BYPANEL
      INTEGER :: K473_LOC
      INTEGER FIRST_BLOCK, LAST_BLOCK
      INTEGER INFO_TMP(2), MAXI_RANK
      INTEGER HF, NPARTSASS, NPARTSCB, NB_BLR
      INTEGER MAXI_CLUSTER, LWORK, NELIM, NELIM_LOC
      INTEGER :: IROW_L, NVSCHUR
      INTEGER, POINTER, DIMENSION(:)          :: PTDummy
      INTEGER, POINTER, DIMENSION(:)          :: BEGS_BLR
      TYPE(LRB_TYPE), POINTER, DIMENSION(:,:) :: CB_LRB
      TYPE(LRB_TYPE), POINTER, DIMENSION(:)   :: ACC_LUA
      TYPE(LRB_TYPE), POINTER, DIMENSION(:)   :: BLR_U, BLR_L
      INTEGER, POINTER, DIMENSION(:)          :: BEGS_BLR_TMP
      TYPE(LRB_TYPE), POINTER, DIMENSION(:)   :: BLR_PANEL
      DOUBLE PRECISION, POINTER, DIMENSION(:)          :: DIAG
      INTEGER :: DIAGSIZ_STA, DIAGSIZ_DYN, DPOS, LorU, I, MEM, MEM_TOT
      INTEGER(8) :: POSELT_DIAG
      CHARACTER(len=1) :: DIR
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), TAU(:)
      INTEGER, ALLOCATABLE :: JPVT(:)
      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: BLOCK(:,:)
      INTEGER :: allocok,J
      INTEGER :: OMP_NUM
      INTEGER :: IP
      INTEGER(8) :: UPOS, LPOS
      INTEGER :: MY_NUM
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: NEXT_BLR_U, NEXT_BLR_L 
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_STATIC
      DOUBLE PRECISION :: ZERO
      PARAMETER (ZERO=0.0D0)
      INCLUDE 'mumps_headers.h'
      INTEGER(8):: KEEP8TMPCOPY, KEEP873COPY
      FIRST_BLOCK = -99999
      LAST_BLOCK  = -99999
      IP=0
      IF (KEEP(206).GE.1) THEN
        Inextpiv = 1   
      ELSE 
        Inextpiv = 0   
      ENDIF
      INOPV = 0
      SEUIL_LOC = SEUIL
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
        UUTEMP=UU
      ENDIF
      PIVOT_OPTION  = KEEP(468)
      LRTRSM_OPTION = KEEP(475)
      LAFAC  = -9999_8  
      XSIZE      = KEEP(IXSZ)
      NFRONT     = IW(IOLDPS+XSIZE)
      NASS       = iabs(IW(IOLDPS+2+XSIZE))
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED = .FALSE.        
      COMPRESS_PANEL = .FALSE.
      COMPRESS_CB = .FALSE.
      NULLIFY(PTDummy)
      NULLIFY(BEGS_BLR)
      NULLIFY(CB_LRB)
      NULLIFY(ACC_LUA)
      NULLIFY(BLR_U)
      NULLIFY(BLR_L)
      NULLIFY(BEGS_BLR_TMP)
      NULLIFY(BLR_PANEL)
      NULLIFY(DIAG)
      COMPRESS_PANEL = (IW(IOLDPS+XXLR).GE.2)
      COMPRESS_CB    = ((IW(IOLDPS+XXLR).EQ.1).OR.
     &                  (IW(IOLDPS+XXLR).EQ.3))
      LR_ACTIVATED   = (IW(IOLDPS+XXLR).GT.0)
      IF (COMPRESS_CB.AND.(.NOT.COMPRESS_PANEL)) THEN
        K473_LOC = 1
      ELSE
        K473_LOC = KEEP(473)
      ENDIF
      K473_LOC = KEEP(473)
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( .NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      OOC_EFFECTIVE_ON_FRONT= ((KEEP(201).EQ.1).AND. 
     &                         OOCWRITE_COMPATIBLE_WITH_BLR)
      CALL DMUMPS_SET_PARPIVT1 ( INODE, NFRONT, NASS, KEEP, 
     &                           LR_ACTIVATED, PARPIV_T1)
      IF (UUTEMP.EQ.ZERO) THEN 
        PIVOT_OPTION=0
      ELSE IF (PARPIV_T1.NE.0) THEN
         PIVOT_OPTION = min(PIVOT_OPTION,2)
      ENDIF
      IF (LR_ACTIVATED) THEN
        IF (LRTRSM_OPTION.EQ.3) THEN
          PIVOT_OPTION = MIN(PIVOT_OPTION,1)
        ELSEIF (LRTRSM_OPTION.EQ.2) THEN
          PIVOT_OPTION = MIN(PIVOT_OPTION, 2)
        ENDIF
      ENDIF
      IF (PIVOT_OPTION.LE.1) THEN
         PARPIV_T1 = 0
      ENDIF
      IF (NASS.LT.KEEP(4)) THEN
        NBKJIB_ORIG = NASS
      ELSE IF (NASS .GT. KEEP(3)) THEN
        NBKJIB_ORIG = min( KEEP(6), NASS )
      ELSE
        NBKJIB_ORIG = min( KEEP(5), NASS )
      ENDIF
      IF (.not.LR_ACTIVATED) THEN
          NBLR_ORIG     = KEEP(420)
      ELSE
          NBLR_ORIG  = -9999 
      ENDIF
      IF ((KEEP(114).EQ.1) .AND. 
     &    (KEEP(116).GT.0) .AND. ((NFRONT-NASS-KEEP(253)).GT.0) 
     &   ) THEN
         IROW_L = IOLDPS+6+XSIZE+NASS
         CALL DMUMPS_COMPUTE_SIZE_SCHUR_IN_FRONT ( 
     &     N, 
     &     NFRONT-NASS-KEEP(253), 
     &     KEEP(116), 
     &     IW(IROW_L), PERM, 
     &     NVSCHUR )
      ELSE
         NVSCHUR = 0
      ENDIF
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      CURRENT_BLR = 0
      CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
      LIWFAC    = IW(IOLDPS+XXI)
      IF ( OOC_EFFECTIVE_ON_FRONT ) THEN 
          LNextPiv2beWritten = 1 
          UNextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = LNextPiv2beWritten 
          PP_FIRST2SWAP_U = UNextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0
          MonBloc%LastPanelWritten_U = 0
          PP_LastPIVRPTRFilled_L = 0 
          PP_LastPIVRPTRFilled_U = 0 
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 1
          MonBloc%NROW     = NFRONT
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -88877    
          NULLIFY(MonBloc%INDICES)   
      ENDIF
      IF (LR_ACTIVATED) THEN
             IF (KEEP(405) .EQ. 1) THEN
!$OMP ATOMIC UPDATE
               CNT_NODES = CNT_NODES + 1 
!$OMP END ATOMIC
             ELSE
               CNT_NODES = CNT_NODES + 1 
             ENDIF
      ELSE IF (KEEP(486).NE.0) THEN
      ENDIF
      OOC_EFF_AND_WRITE_BYPANEL  = ( (PIVOT_OPTION.GE.3) .AND.
     &                                     OOC_EFFECTIVE_ON_FRONT )
      HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
      IF (LR_ACTIVATED) THEN
         CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NFRONT-1), NASS,
     &        NFRONT-NASS, LRGROUPS, NPARTSCB, 
     &        NPARTSASS, BEGS_BLR)
         CALL REGROUPING2(BEGS_BLR, NPARTSASS, NASS, NPARTSCB,
     &        NFRONT-NASS, KEEP(488), .FALSE., KEEP(472))     
         NB_BLR = NPARTSASS + NPARTSCB
         call MAX_CLUSTER(BEGS_BLR,NB_BLR,MAXI_CLUSTER)
         MAXI_RANK = KEEP(479)*MAXI_CLUSTER
         LWORK = MAXI_CLUSTER*MAXI_CLUSTER
         OMP_NUM = 1
#if defined(BLR_MT)
!$       OMP_NUM = OMP_GET_MAX_THREADS()
#endif
         ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),
     &             stat=allocok)
         IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
           GOTO 490
         ENDIF
         ALLOCATE(ACC_LUA(OMP_NUM),stat=allocok)
         IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = OMP_NUM
           GOTO 490
         ENDIF
         IF (KEEP(480).GE.3) THEN
           DO MY_NUM=1,OMP_NUM
             CALL ALLOC_LRB(ACC_LUA(MY_NUM), MAXI_RANK, 
     &                      MAXI_CLUSTER, MAXI_CLUSTER, .TRUE.,
     &                      IFLAG, IERROR, KEEP8)
             IF (IFLAG.LT.0) GOTO 490
             ACC_LUA(MY_NUM)%K = 0
           ENDDO
         ENDIF
      ENDIF
      IF (LR_ACTIVATED.AND.
     &       (KEEP(480).NE.0
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       .OR.COMPRESS_CB
     &       )) THEN
        INFO_TMP(1) = IFLAG
        INFO_TMP(2) = IERROR
        CALL DMUMPS_BLR_SAVE_INIT(IW(IOLDPS+XXF), 
     &              .FALSE., 
     &              .FALSE., 
     &              .FALSE., 
     &              NPARTSASS, 
     &              BEGS_BLR, PTDummy, 
     &              huge(NPARTSASS), 
     &              INFO_TMP)
        IFLAG  = INFO_TMP(1) 
        IERROR = INFO_TMP(2) 
        IF (IFLAG.LT.0) GOTO 500
      ENDIF
      IF (COMPRESS_CB.AND.NPARTSCB.GT.0) THEN
        allocate(CB_LRB(NPARTSCB,NPARTSCB),stat=allocok)
        IF (allocok > 0) THEN
          IFLAG  = -13
          IERROR = NPARTSCB*NPARTSCB
          GOTO 490
        ENDIF
        CALL DMUMPS_BLR_SAVE_CB_LRB(IW(IOLDPS+XXF),CB_LRB)
      ENDIF
      DO WHILE (IEND_BLR < NASS ) 
        CURRENT_BLR = CURRENT_BLR + 1
        IBEG_BLR = IW(IOLDPS+1+KEEP(IXSZ)) + 1 
        IF (.NOT. LR_ACTIVATED) THEN
          IEND_BLR = min(IEND_BLR + NBLR_ORIG, NASS)
        ELSE
          IEND_BLR = min(BEGS_BLR(CURRENT_BLR+1)-1, NASS)
          BEGS_BLR( CURRENT_BLR ) = IBEG_BLR
          IF ( IEND_BLR - IBEG_BLR + 1 .GT. MAXI_CLUSTER ) THEN
            MAXI_CLUSTER = IEND_BLR - IBEG_BLR + 1
            LWORK = MAXI_CLUSTER*MAXI_CLUSTER
            DEALLOCATE(BLOCK, WORK, RWORK, TAU, JPVT)
            ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
              GOTO 490
            ENDIF
            IF (KEEP(480).GE.3) THEN
              DO MY_NUM=1,OMP_NUM
                CALL DEALLOC_LRB(ACC_LUA(MY_NUM),KEEP8)
                CALL ALLOC_LRB(ACC_LUA(MY_NUM), MAXI_RANK,
     &                         MAXI_CLUSTER, MAXI_CLUSTER, .TRUE.,
     &                         IFLAG, IERROR, KEEP8)
                IF (IFLAG.LT.0) GOTO 500
                ACC_LUA(MY_NUM)%K = 0
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        IF (LR_ACTIVATED) THEN
          IF (KEEP(480).GE.5) THEN
            IF (CURRENT_BLR.EQ.1) THEN
              ALLOCATE(BLR_U(NB_BLR-CURRENT_BLR),stat=allocok)
              IF (allocok > 0) THEN
                 IFLAG  = -13
                 IERROR = NB_BLR-CURRENT_BLR
                 GOTO 490
              ENDIF  
              ALLOCATE(BLR_L(NB_BLR-CURRENT_BLR),stat=allocok)
              IF (allocok > 0) THEN
                 IFLAG  = -13
                 IERROR = NB_BLR-CURRENT_BLR
                 GOTO 490
              ENDIF 
              IF (NB_BLR.GT.CURRENT_BLR) THEN
                BLR_U(1:NB_BLR-CURRENT_BLR)%ISLR=.FALSE.
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              1, 
     &              CURRENT_BLR, BLR_U)
                BLR_L(1:NB_BLR-CURRENT_BLR)%ISLR=.FALSE.
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR, BLR_L)
              ENDIF
            ELSE
              IF (NB_BLR.GT.CURRENT_BLR) THEN
                CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &              IW(IOLDPS+XXF),
     &              1, 
     &              CURRENT_BLR, BLR_U)
                CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR, BLR_L)
              ENDIF
            ENDIF
            IF (CURRENT_BLR.LT.NPARTSASS) THEN
              ALLOCATE(NEXT_BLR_U(NB_BLR-CURRENT_BLR-1),stat=allocok)
              IF (allocok > 0) THEN
                 IFLAG  = -13
                 IERROR = NB_BLR-CURRENT_BLR-1
                 GOTO 490
              ENDIF 
              ALLOCATE(NEXT_BLR_L(NB_BLR-CURRENT_BLR-1),stat=allocok)
              IF (allocok > 0) THEN
                 IFLAG  = -13
                 IERROR = NB_BLR-CURRENT_BLR-1
                 GOTO 490
              ENDIF 
              IF (NB_BLR.GT.CURRENT_BLR+1) THEN
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              1, 
     &              CURRENT_BLR+1, NEXT_BLR_U)
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR+1, NEXT_BLR_L)
              ENDIF
            ENDIF
          ELSE
            ALLOCATE(BLR_U(NB_BLR-CURRENT_BLR),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = NB_BLR-CURRENT_BLR
              GOTO 490
            ENDIF 
            ALLOCATE(BLR_L(NB_BLR-CURRENT_BLR),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = NB_BLR-CURRENT_BLR
              GOTO 490
            ENDIF
          ENDIF
        ENDIF
        DO WHILE (IEND_BLOCK < IEND_BLR ) 
          IBEG_BLOCK = IW(IOLDPS+1+KEEP(IXSZ)) + 1
          IF (KEEP(405).EQ.0) THEN
            KEEP(425)=max(KEEP(425),IEND_BLOCK-IBEG_BLOCK)
          ELSE
!$OMP       ATOMIC UPDATE
            KEEP(425)=max(KEEP(425),IEND_BLOCK-IBEG_BLOCK)
!$OMP       END ATOMIC
          ENDIF
          IEND_BLOCK = min(IEND_BLOCK + NBKJIB_ORIG, IEND_BLR)
  50      CONTINUE  
            CALL DMUMPS_FAC_I(NFRONT,NASS,NFRONT,
     &      IBEG_BLOCK,IEND_BLOCK,N,INODE,
     &      IW,LIW,A,LA,INOPV,NOFFW,NBTINYW,
     &      DET_EXPW, DET_MANTW, DET_SIGNW,
     &      IFLAG,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     &      DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &      PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled_L,
     &      PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &      PP_LastPIVRPTRFilled_U,
     &      PIVOT_OPTION, LR_ACTIVATED, IEND_BLR,
     &      Inextpiv, OOC_EFFECTIVE_ON_FRONT, 
     &      NVSCHUR, PARPIV_T1
     &      )
            IF (IFLAG.LT.0) GOTO 500  
          IF (INOPV.EQ.1) THEN
            IF(STATICMODE) THEN
              INOPV = -1
              GOTO 50 
            ENDIF
          ELSE IF ( INOPV.LE.0 ) THEN 
            IF (PIVOT_OPTION.GE.3) THEN
              LAST_COL = NFRONT
            ELSEIF (PIVOT_OPTION.EQ.2) THEN
              LAST_COL = NASS
            ELSE
              LAST_COL = IEND_BLR
            ENDIF
            CALL DMUMPS_FAC_MQ(IBEG_BLOCK, IEND_BLOCK,
     &              NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &              LAST_COL, A, LA, POSELT, IFINB,
     &              LR_ACTIVATED
     &              )
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ENDIF
          ENDIF
          IF ( OOC_EFF_AND_WRITE_BYPANEL ) THEN
            MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
            STRAT          = STRAT_TRY_WRITE
            LAST_CALL      = .FALSE.
            CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_U,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
            IF (IFLAG_OOC < 0 ) THEN
              IFLAG=IFLAG_OOC
              GOTO 500
            ENDIF
          ENDIF
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
            IF (PIVOT_OPTION.GE.3) THEN
              LAST_COL = NFRONT
            ELSEIF (PIVOT_OPTION.EQ.2) THEN
              LAST_COL = NASS
            ELSE
              LAST_COL = IEND_BLR
            ENDIF
            CALL DMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK,
     &            NPIV, NFRONT, IEND_BLR, LAST_COL,
     &            A, LA, POSELT, 
     &            -66666, 
     &            .TRUE., .FALSE., .TRUE., 
     &            .FALSE., 
     &            LR_ACTIVATED      
     &            )
          ENDIF
        END DO 
        NPIV   =  IW(IOLDPS+1+XSIZE)
        IF (.NOT. LR_ACTIVATED
     &      .OR. (.NOT. COMPRESS_PANEL)
     &     ) THEN
          IF (PIVOT_OPTION.EQ.4) THEN
            LAST_ROW = NFRONT
          ELSE
            LAST_ROW = NASS
          ENDIF
          IF (PIVOT_OPTION.GE.3) THEN
            LAST_COL = NFRONT
          ELSE
            LAST_COL = NASS
          ENDIF
          IF (IEND_BLR.LT.LAST_ROW) THEN
            CALL DMUMPS_FAC_SQ(IBEG_BLR, IEND_BLR,
     &            NPIV, NFRONT, LAST_ROW, LAST_COL, 
     &            A, LA, POSELT, IEND_BLR, .TRUE., (PIVOT_OPTION.LT.2),
     &            .TRUE., .FALSE.,  
     &            LR_ACTIVATED)       
          ENDIF
        ELSE
          NELIM = IEND_BLR - NPIV
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            IF (KEEP(480).GE.2
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &         ) THEN
              DO J=1,NB_BLR-CURRENT_BLR
                 BLR_U(J)%M=0
                 BLR_U(J)%N=0
                 BLR_U(J)%K=0
                 BLR_U(J)%ISLR=.FALSE.
                 NULLIFY(BLR_U(J)%Q)
                 NULLIFY(BLR_U(J)%R)
              ENDDO
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              1, 
     &             CURRENT_BLR, BLR_U)
              DO J=1,NB_BLR-CURRENT_BLR
                 BLR_L(J)%M=0
                 BLR_L(J)%N=0
                 BLR_L(J)%K=0
                 BLR_L(J)%ISLR=.FALSE.
                 NULLIFY(BLR_L(J)%Q)
                 NULLIFY(BLR_L(J)%R)
              ENDDO
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR, BLR_L)
              NULLIFY(BLR_L)
              NULLIFY(BLR_U)
            IF (KEEP(480).GE.2 .AND. IEND_BLR.LT.NASS) THEN
              IF (LRTRSM_OPTION.EQ.3) THEN
                FIRST_BLOCK = 1
              ELSE
                FIRST_BLOCK = NPARTSASS-CURRENT_BLR
              ENDIF
#if defined(BLR_MT)
!$OMP PARALLEL
#endif
              CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 0, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          1, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK, 
     &          KEEP(474), 0, BLR_U, 
     &          KEEP8, 
     &          FIRST_BLOCK=FIRST_BLOCK)
              IF (IFLAG.LT.0) GOTO 900
              CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 1, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          1, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK,
     &          KEEP(474), 0, BLR_U, 
     &          KEEP8, 
     &          FIRST_BLOCK=FIRST_BLOCK)
 900            CONTINUE              
#if defined(BLR_MT)
!$OMP END PARALLEL
#endif
              IF (IFLAG.LT.0) GOTO 500
            ENDIF
            ENDIF
            IF (KEEP(486).EQ.3) THEN
              IF (KEEP(480).EQ.0) THEN
                DEALLOCATE(BLR_U,BLR_L)
                NULLIFY(BLR_L)
                NULLIFY(BLR_U)
              ENDIF
            ENDIF
            GOTO 100
          ENDIF
          IF (PIVOT_OPTION.GE.3) THEN
            FIRST_COL = NFRONT
          ELSEIF (PIVOT_OPTION.EQ.2) THEN
            FIRST_COL = NASS
          ELSE
            FIRST_COL = IEND_BLR
          ENDIF
          IF (LRTRSM_OPTION.EQ.3) THEN
            LAST_COL = IEND_BLR
          ELSEIF (LRTRSM_OPTION.EQ.2) THEN
            LAST_COL = NASS
          ELSE
            LAST_COL = NFRONT
          ENDIF
          CALL_LTRSM = (LRTRSM_OPTION.EQ.0)
          CALL_UTRSM = (LAST_COL-FIRST_COL.GT.0)
          IF ((IEND_BLR.LT.NFRONT) .AND. 
     &          (CALL_LTRSM.OR.CALL_UTRSM)) THEN
          CALL DMUMPS_FAC_SQ(IBEG_BLR, IEND_BLR,
     &            NPIV, NFRONT, NFRONT, 
     &            LAST_COL,
     &            A, LA, POSELT, 
     &            FIRST_COL, CALL_LTRSM,
     &            CALL_UTRSM, .FALSE.,
     &            .FALSE.,  
     &            LR_ACTIVATED)    
        ENDIF
#if defined(BLR_MT)          
#endif
#if defined(BLR_MT)
!$OMP PARALLEL PRIVATE(UPOS,LPOS) FIRSTPRIVATE(FIRST_BLOCK,LAST_BLOCK)
#endif
              CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, 
     &              IERROR, 
     &              NFRONT,
     &              BEGS_BLR, NB_BLR, DKEEP(8), KEEP(466), K473_LOC,
     &              BLR_U, CURRENT_BLR,
     &              'H', WORK, TAU, JPVT, LWORK, RWORK,
     &              BLOCK, MAXI_CLUSTER, NELIM,
     &              .FALSE., 0, 0,
     &              1, KEEP(483), KEEP8,
     &              K480=KEEP(480)
     &        )
#if defined(BLR_MT)
!$OMP BARRIER
#endif
              IF (IFLAG.LT.0) GOTO 400
            CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &          NFRONT,
     &          BEGS_BLR, NB_BLR, DKEEP(8), KEEP(466), K473_LOC, BLR_L, 
     &          CURRENT_BLR,
     &          'V', WORK, TAU, JPVT, LWORK, RWORK,
     &          BLOCK, MAXI_CLUSTER, NELIM,
     &          .FALSE., 0, 0,
     &          1, KEEP(483), KEEP8,
     &          K480=KEEP(480)
     &          )
#if defined(BLR_MT)
!$OMP BARRIER
!$OMP MASTER
#endif
          IF (KEEP(480).NE.0
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       ) THEN
            IF (KEEP(480).LT.5) THEN
              CALL DMUMPS_BLR_SAVE_PANEL_LORU(
     &            IW(IOLDPS+XXF),
     &            1, 
     &            CURRENT_BLR, BLR_U)
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &             IW(IOLDPS+XXF),
     &             0, 
     &             CURRENT_BLR, BLR_L)
            ENDIF
          ENDIF
#if defined(BLR_MT)          
!$OMP END MASTER
!$OMP BARRIER
#endif          
          IF (IFLAG.LT.0) GOTO 400
          IF (LRTRSM_OPTION.GT.0) THEN
            CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSELT, NFRONT,
     &                IBEG_BLR, 
     &                NB_BLR, BLR_L, CURRENT_BLR, CURRENT_BLR+1, 
     &                NB_BLR, 1, 0, 0, .FALSE.)     
            IF (PIVOT_OPTION.LT.3.AND.LRTRSM_OPTION.GE.2) THEN
              IF (PIVOT_OPTION.LE.1.AND.LRTRSM_OPTION.EQ.3) THEN
                FIRST_BLOCK = CURRENT_BLR+1
              ELSE
                FIRST_BLOCK = NPARTSASS+1
              ENDIF
              CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSELT, NFRONT,
     &                IBEG_BLR, NB_BLR, BLR_U, 
     &                CURRENT_BLR, FIRST_BLOCK, NB_BLR,
     &                1, 0, 1, .FALSE.)
#if defined(BLR_MT)          
!$OMP BARRIER
#endif          
              CALL DMUMPS_BLR_UPD_NELIM_VAR_U(
     &                A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &                BEGS_BLR, CURRENT_BLR, BLR_U, NB_BLR, 
     &                FIRST_BLOCK, IBEG_BLR, NPIV, NELIM)
            ENDIF
          ENDIF
#if defined(BLR_MT)          
!$OMP BARRIER
#endif          
          IF (IFLAG.LT.0) GOTO 400
          IF (KEEP(480).GE.2) THEN
            UPOS = POSELT+int(BEGS_BLR(CURRENT_BLR)-1,8)*int(NFRONT,8)
     &                   +int(BEGS_BLR(CURRENT_BLR+1)-NELIM-1,8)
            LPOS = POSELT+int(BEGS_BLR(CURRENT_BLR+1)-1,8)*int(NFRONT,8)
     &                   +int(BEGS_BLR(CURRENT_BLR+1)-NELIM-1,8)
            CALL DMUMPS_BLR_UPD_NELIM_VAR_L(A, LA, UPOS, A, LA, 
     &        LPOS, IFLAG, IERROR, NFRONT, NFRONT,
     &        BEGS_BLR, CURRENT_BLR, BLR_L, NB_BLR, 
     &        CURRENT_BLR+1, NELIM, 'N')
            IF (IFLAG.LT.0) GOTO 444
            IF (IEND_BLR.LT.NASS) THEN
              IF (LRTRSM_OPTION.EQ.3) THEN
                FIRST_BLOCK = 1
              ELSE
                FIRST_BLOCK = NPARTSASS-CURRENT_BLR
              ENDIF
              CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 0, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          1, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK,
     &          KEEP(474), 0, BLR_U, 
     &          KEEP8, 
     &          FIRST_BLOCK=FIRST_BLOCK)
              IF (IFLAG.LT.0) GOTO 442
              CALL DMUMPS_BLR_UPD_PANEL_LEFT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 1, 
     &          BEGS_BLR, BEGS_BLR, CURRENT_BLR, ACC_LUA,
     &          NB_BLR, NPARTSASS, NELIM,
     &          1, 0, 
     &          .FALSE., IFLAG, IERROR, 0,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), MAXI_CLUSTER, MAXI_RANK,
     &          KEEP(474), 0, BLR_U, 
     &          KEEP8, 
     &          FIRST_BLOCK=FIRST_BLOCK)
 442          CONTINUE
            ENDIF
 444        CONTINUE
          ELSE
            CALL DMUMPS_BLR_UPDATE_TRAILING(A, LA, POSELT, 
     &        IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, BEGS_BLR, CURRENT_BLR, BLR_L, NB_BLR, 
     &        BLR_U, NB_BLR, 
     &        NELIM,.FALSE., 0,
     &        1, 0, 
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
          ENDIF
#if defined(BLR_MT)
!$OMP BARRIER
#endif
          IF (IFLAG.LT.0) GOTO 400
          IF (KEEP(486).NE.2) THEN
            LAST_BLOCK = NB_BLR
          ELSEIF(UU.GT.0) THEN
            LAST_BLOCK = NPARTSASS
          ELSE
            LAST_BLOCK = CURRENT_BLR
          ENDIF
          IF (LRTRSM_OPTION.GT.0) THEN
            FIRST_BLOCK = CURRENT_BLR+1
            CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &            NFRONT, .TRUE.,
     &             BEGS_BLR(CURRENT_BLR),
     &             BEGS_BLR(CURRENT_BLR+1), 
     &             NB_BLR, BLR_L, CURRENT_BLR, 'V', 1,
     &             BEG_I_IN=FIRST_BLOCK, END_I_IN=LAST_BLOCK)
#if defined(BLR_MT)          
#endif          
          ENDIF
          IF (LRTRSM_OPTION.GE.2) THEN
            IF (LRTRSM_OPTION.EQ.2) THEN
              FIRST_BLOCK = NPARTSASS+1
            ELSE
              FIRST_BLOCK = CURRENT_BLR+1
            ENDIF
            CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT,
     &            NFRONT, .TRUE.,
     &             BEGS_BLR(CURRENT_BLR),
     &             BEGS_BLR(CURRENT_BLR+1), 
     &             NB_BLR, BLR_U, CURRENT_BLR, 'H', 1,
     &             BEG_I_IN=FIRST_BLOCK, END_I_IN=LAST_BLOCK)
          ENDIF
 400      CONTINUE
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 500
          CALL UPD_MRY_LU_LRGAIN(BLR_U,
     &               NB_BLR-CURRENT_BLR-NPARTSCB,
     &               NPARTSCB, 'H')
          CALL UPD_MRY_LU_LRGAIN(BLR_L,
     &               NB_BLR-CURRENT_BLR-NPARTSCB,
     &               NPARTSCB, 'V')
          IF (KEEP(486).EQ.3) THEN
            IF (KEEP(480).EQ.0) THEN
              CALL DEALLOC_BLR_PANEL (BLR_U, NB_BLR-CURRENT_BLR, KEEP8)
              CALL DEALLOC_BLR_PANEL (BLR_L, NB_BLR-CURRENT_BLR, KEEP8)
              DEALLOCATE(BLR_U,BLR_L)
            ENDIF
          ENDIF
          NULLIFY(BLR_L)
          NULLIFY(BLR_U)
        ENDIF
        IF ( OOC_EFF_AND_WRITE_BYPANEL ) THEN
             IF (PIVOT_OPTION.LT.4) THEN
               TYPEF_LOC = TYPEF_U
             ELSE
               TYPEF_LOC = TYPEF_BOTH_LU
             ENDIF
             MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
             STRAT          = STRAT_TRY_WRITE
             LAST_CALL      = .FALSE.
             CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_LOC,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
             IF (IFLAG_OOC < 0 ) THEN
                IFLAG=IFLAG_OOC
                GOTO 500
             ENDIF
        ENDIF
 100    CONTINUE
      END DO 
      IF (LR_ACTIVATED) THEN
          IBEG_BLR = IW(IOLDPS+1+KEEP(IXSZ)) + 1 
          BEGS_BLR( CURRENT_BLR + 1 ) = IBEG_BLR
          IF (
     &         (KEEP(486).EQ.2) 
     &       ) THEN
            CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW(IOLDPS+XXF),
     &                        BEGS_BLR_STATIC)
            IF (UU.GT.0) THEN
              allocate(BEGS_BLR_TMP(NB_BLR+1),stat=allocok)
              IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = NB_BLR+1
                GOTO 500
              ENDIF
              DO IP=1,NB_BLR+1
                 BEGS_BLR_TMP(IP) = BEGS_BLR_STATIC(IP)
              ENDDO
            ENDIF
          ENDIF
          MEM_TOT = 0
#if defined(BLR_MT)          
!$OMP PARALLEL 
!$OMP& PRIVATE(IP, LorU, DIR, NELIM_LOC)
#endif
          IF (
     &         (KEEP(486).EQ.2) 
     &       ) THEN
#if defined(BLR_MT)
!$OMP DO PRIVATE(DIAG, DIAGSIZ_STA, DIAGSIZ_DYN, DPOS, POSELT_DIAG, MEM,
!$OMP&           allocok)
!$OMP&   REDUCTION(+:MEM_TOT)
#endif
            DO IP=1,NPARTSASS
              IF (IFLAG.LT.0) CYCLE
              DIAGSIZ_DYN = BEGS_BLR(IP+1)-BEGS_BLR(IP)
              DIAGSIZ_STA = BEGS_BLR_STATIC(IP+1)-BEGS_BLR(IP)
              MEM = DIAGSIZ_DYN*(2*DIAGSIZ_STA-DIAGSIZ_DYN)
              MEM_TOT = MEM_TOT + MEM
              ALLOCATE(DIAG(MEM), stat=allocok)
              IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = MEM
                CYCLE
              ENDIF 
              DPOS = 1
              POSELT_DIAG = POSELT + int(BEGS_BLR(IP)-1,8)*int(NFRONT,8)
     &                             + int(BEGS_BLR(IP)-1,8)
              DO I=1,DIAGSIZ_STA
                IF (I.LE.DIAGSIZ_DYN) THEN
                  DIAG(DPOS:DPOS+DIAGSIZ_STA-1) =
     &                   A(POSELT_DIAG:POSELT_DIAG+int(DIAGSIZ_STA-1,8))
                  DPOS = DPOS + DIAGSIZ_STA
                ELSE
                  DIAG(DPOS:DPOS+DIAGSIZ_DYN-1) =
     &                   A(POSELT_DIAG:POSELT_DIAG+int(DIAGSIZ_DYN-1,8))
                  DPOS = DPOS + DIAGSIZ_DYN
                ENDIF
                POSELT_DIAG = POSELT_DIAG + int(NFRONT,8)
              ENDDO
              CALL DMUMPS_BLR_SAVE_DIAG_BLOCK(
     &              IW(IOLDPS+XXF),
     &              IP, DIAG)
            ENDDO
#if defined(BLR_MT)
!$OMP ENDDO
!$OMP SINGLE
#endif
            IF (KEEP(405) .EQ. 0) THEN
              KEEP8(69)    = KEEP8(69) + int(MEM_TOT,8)
              KEEP8TMPCOPY = KEEP8(69)
              KEEP8(68)    = max(KEEP8TMPCOPY, KEEP8(68))
              KEEP8(71)    = KEEP8(71) + int(MEM_TOT,8) 
              KEEP8TMPCOPY = KEEP8(71)
              KEEP8(70)    = max(KEEP8(71), KEEP8(70))
              KEEP8(73)    = KEEP8(73) + int(MEM_TOT,8)
              KEEP873COPY  = KEEP8(73)
              KEEP8(74)    = max(KEEP8(74), KEEP873COPY)
            ELSE
!$OMP         ATOMIC CAPTURE
              KEEP8(69)    = KEEP8(69) + int(MEM_TOT,8)
              KEEP8TMPCOPY = KEEP8(69)
!$OMP         END ATOMIC
!$OMP         ATOMIC UPDATE
              KEEP8(68)    = max(KEEP8TMPCOPY, KEEP8(68))
!$OMP         END ATOMIC
!$OMP         ATOMIC CAPTURE
              KEEP8(71)    = KEEP8(71) + int(MEM_TOT,8) 
              KEEP8TMPCOPY = KEEP8(71)
!$OMP         END ATOMIC
!$OMP         ATOMIC UPDATE
              KEEP8(70)    = max(KEEP8TMPCOPY, KEEP8(70))
!$OMP         END ATOMIC
!$OMP         ATOMIC CAPTURE
              KEEP8(73)    = KEEP8(73) + int(MEM_TOT,8)
              KEEP873COPY  = KEEP8(73)
!$OMP         END ATOMIC
!$OMP         ATOMIC UPDATE
              KEEP8(74)    = max(KEEP8(74), KEEP873COPY)
!$OMP         END ATOMIC
            ENDIF
            IF ( KEEP873COPY .GT. KEEP8(75) ) THEN
             IFLAG = -19
             CALL MUMPS_SET_IERROR(
     &             (KEEP873COPY-KEEP8(75)), IERROR)
            ENDIF
#if defined(BLR_MT)
!$OMP END SINGLE
#endif
            IF (IFLAG.LT.0) GOTO 447
            IF (UU.GT.0) THEN
              DO IP=1,NPARTSASS
                NELIM_LOC = BEGS_BLR_TMP(IP+1)-BEGS_BLR(IP+1)
                DO LorU=0,1
#if defined(BLR_MT)
!$OMP SINGLE
#endif
                  CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IW(IOLDPS+XXF), LorU, IP, BLR_PANEL)
                  CALL DEALLOC_BLR_PANEL(BLR_PANEL, NPARTSASS-IP, KEEP8)
#if defined(BLR_MT)
!$OMP END SINGLE
#endif
                  IF (LorU.EQ.0) THEN 
                    DIR = 'V'
                  ELSE
                    DIR = 'H'
                  ENDIF
                  CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG,
     &              IERROR, NFRONT, BEGS_BLR_TMP,
     &              NB_BLR, DKEEP(8), KEEP(466), K473_LOC,
     &              BLR_PANEL, IP,
     &              DIR, WORK, TAU, JPVT, LWORK, RWORK,
     &              BLOCK, MAXI_CLUSTER, NELIM_LOC,
     &              .FALSE., 0, 0,
     &              1, KEEP(483), KEEP8,
     &              END_I_IN=NPARTSASS, FRSWAP=.TRUE.
     &             )
#if defined(BLR_MT)
!$OMP BARRIER
#endif
                  IF (IFLAG.LT.0) GOTO 445
                ENDDO
#if defined(BLR_MT)
!$OMP BARRIER
!$OMP SINGLE
#endif
                BEGS_BLR_TMP(IP+1) = BEGS_BLR(IP+1)
#if defined(BLR_MT)
!$OMP END SINGLE
#endif
              ENDDO
#if defined(BLR_MT)
!$OMP BARRIER
#endif
 445          CONTINUE
            ENDIF 
 447        CONTINUE
          ENDIF 
          IF (IFLAG .LT. 0) GOTO 450
          IF (KEEP(480) .GE. 2) THEN
#if defined(BLR_MT)
!$OMP SINGLE
#endif
            CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW(IOLDPS+XXF),
     &                        BEGS_BLR_STATIC)
#if defined(BLR_MT)
!$OMP END SINGLE
#endif
            CALL DMUMPS_BLR_UPD_CB_LEFT(A, LA, POSELT, NFRONT,
     &          BEGS_BLR_STATIC, BEGS_BLR_STATIC, 
     &          NPARTSCB, NPARTSCB, NPARTSASS, NASS, 
     &          IW(IOLDPS+XXF),
     &          1, .FALSE., IFLAG, IERROR,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          ACC_LUA, KEEP(480),KEEP(479),KEEP(478),KEEP(476), 
     &          KEEP(484), MAXI_CLUSTER, MAXI_RANK,
     &          KEEP(474), 0, BLR_U, 
     &          .FALSE.,
     &          CB_LRB, KEEP8)
#if defined(BLR_MT)
!$OMP BARRIER
#endif
          ENDIF
          IF (IFLAG.LT.0) GOTO 450
#if defined(BLR_MT)
!$OMP MASTER
#endif
          IF (COMPRESS_CB
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       ) THEN
            CALL DMUMPS_BLR_SAVE_BEGS_BLR_DYN(IW(IOLDPS+XXF),
     &        BEGS_BLR)
          ENDIF
          IF (COMPRESS_CB) THEN
            IEND_BLR = BEGS_BLR(CURRENT_BLR+2)
            IF ( IEND_BLR - IBEG_BLR + 1 .GT. MAXI_CLUSTER ) THEN
              MAXI_CLUSTER = IEND_BLR - IBEG_BLR + 1
              LWORK = MAXI_CLUSTER*MAXI_CLUSTER
              DEALLOCATE(BLOCK, WORK, RWORK, TAU, JPVT)
              ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
              IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
              ENDIF
            ENDIF
          ENDIF
#if defined(BLR_MT)
!$OMP END MASTER
!$OMP BARRIER
#endif
          IF (IFLAG.LT.0) GOTO 450
          IF (COMPRESS_CB) THEN
            CALL DMUMPS_COMPRESS_CB(A, LA, POSELT, NFRONT,
     &      BEGS_BLR, BEGS_BLR, NPARTSCB, NPARTSCB, NPARTSASS, 
     &      NFRONT-NASS, NFRONT-NASS, INODE,
     &      IW(IOLDPS+XXF), 0, 1, IFLAG, IERROR,
     &      DKEEP(12), KEEP(466), KEEP(484), KEEP(489), CB_LRB,
     &      WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &      MAXI_CLUSTER, KEEP8, 
     &      -9999, -9999, -9999, KEEP(1), 
     &      NELIM=NELIM)
#if defined(BLR_MT)
!$OMP BARRIER
#endif
          ENDIF
 450      CONTINUE          
#if defined(BLR_MT)          
!$OMP END PARALLEL
#endif
          IF (
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND.UU.GT.0
     &       ) THEN
            deallocate(BEGS_BLR_TMP)
          ENDIF
          IF (IFLAG.LT.0) GOTO 500
         CALL UPD_MRY_LU_FR(NASS, NFRONT-NASS, 0, NASS-NPIV)
         CALL UPD_FLOP_FACTO_FR(NFRONT, NASS, NPIV, 0, 1)
       ENDIF
       IF ( (PIVOT_OPTION.LT.4) .AND. (.NOT.LR_ACTIVATED) ) THEN
         CALL  DMUMPS_FAC_FR_UPDATE_CBROWS( INODE,
     &     NFRONT, NASS, (PIVOT_OPTION.LT.3), A, LA, LAFAC, POSELT, 
     &     IW, LIW, IOLDPS, MonBloc, MYID, NOFFW,
     &     DET_EXPW, DET_MANTW, DET_SIGNW,
     &     LIWFAC, 
     &     PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &     LNextPiv2beWritten, UNextPiv2beWritten, 
     &     PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &     
     &     XSIZE, SEUIL, UU, DKEEP, KEEP8, KEEP, IFLAG,
     &     OOC_EFFECTIVE_ON_FRONT, NVSCHUR )
       ENDIF
       IF (KEEP(486).NE.0) THEN
         IF (.NOT.LR_ACTIVATED) THEN
           CALL UPD_FLOP_FRFRONTS(NFRONT, NPIV, NASS, 0, 1)
         ENDIF
       ENDIF
      IF ( OOC_EFFECTIVE_ON_FRONT ) THEN
          STRAT            = STRAT_WRITE_MAX   
          MonBloc%Last     = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_BOTH_LU,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
          IF (IFLAG_OOC < 0 ) THEN
            IFLAG=IFLAG_OOC
            GOTO 500
          ENDIF
          CALL DMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 600
 490  CONTINUE
 500  CONTINUE
 600  CONTINUE
      IF (LR_ACTIVATED) THEN
        IF (allocated(WORK))  deallocate(WORK)
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(TAU))   deallocate(TAU)
        IF (allocated(JPVT))  deallocate(JPVT)
        IF (allocated(BLOCK)) deallocate(BLOCK)
        IF (associated(ACC_LUA)) THEN 
         IF (KEEP(480).GE.3) THEN
          DO MY_NUM=1,OMP_NUM
            CALL DEALLOC_LRB(ACC_LUA(MY_NUM),KEEP8)
          ENDDO
         ENDIF
         DEALLOCATE(ACC_LUA)
         NULLIFY(ACC_LUA)
        ENDIF
        IF (associated(BEGS_BLR)) THEN
          DEALLOCATE(BEGS_BLR)
          NULLIFY(BEGS_BLR)
        ENDIF
      ENDIF
      IF (LR_ACTIVATED.AND.(KEEP(480).NE.0)) THEN
        IF (.NOT.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &     ) THEN
          CALL DMUMPS_BLR_FREE_ALL_PANELS(IW(IOLDPS+XXF), 2, 
     &                       KEEP8)
        ENDIF
      ENDIF
      IF (LR_ACTIVATED) THEN
        IF (.NOT.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &    .AND..NOT.COMPRESS_CB) THEN
          CALL DMUMPS_BLR_END_FRONT(IW(IOLDPS+XXF), IFLAG, KEEP8,
     &    MTK405=KEEP(405))
        ENDIF
      ENDIF
      NPVW = NPVW + IW(IOLDPS+1+XSIZE)
      END SUBROUTINE DMUMPS_FAC1_LU
      END MODULE DMUMPS_FAC1_LU_M
