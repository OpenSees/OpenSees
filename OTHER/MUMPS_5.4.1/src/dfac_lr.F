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
      MODULE DMUMPS_FAC_LR
      USE DMUMPS_LR_CORE
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING_LDLT(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, CURRENT_BLR, BLR_L,
     &        NELIM, IW2, BLOCK,
     &        MAXI_CLUSTER, NPIV, NIV, 
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR,
     &   NELIM, MAXI_CLUSTER, NPIV, NIV, TOL_OPT
      INTEGER, intent(inout)         :: IFLAG, IERROR
      DOUBLE PRECISION, intent(inout)    :: A(LA)
      TYPE(LRB_TYPE),intent(in) :: BLR_L(:)
      DOUBLE PRECISION, INTENT(INOUT) :: BLOCK(MAXI_CLUSTER,*)
      INTEGER, intent(in) :: IW2(*)
      INTEGER, DIMENSION(:) :: BEGS_BLR
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER :: I, NB_BLOCKS_PANEL, J, MID_RANK
      LOGICAL :: BUILDQ
      INTEGER :: OMP_NUM
      INTEGER :: IBIS
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELTT, POSELTD
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      POSELTD = POSELT + int(NFRONT,8) * int(BEGS_BLR(CURRENT_BLR)-1,8)
     &          + int(BEGS_BLR(CURRENT_BLR) - 1,8)
      OMP_NUM = 0 
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, OMP_NUM,
!$OMP&         MID_RANK, BUILDQ)
#endif
      DO IBIS = 1, (NB_BLOCKS_PANEL*(NB_BLOCKS_PANEL+1)/2) 
        IF (IFLAG.LT.0) CYCLE
        I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
        J = IBIS - I*(I-1)/2
#if defined(BLR_MT)         
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
#endif
            POSELTT = POSELT + int(NFRONT,8) *
     &                int(BEGS_BLR(CURRENT_BLR+I)-1,8)
     &           + int(BEGS_BLR(CURRENT_BLR+J) - 1, 8)
            CALL DMUMPS_LRGEMM4(MONE,
     &            BLR_L(J), BLR_L(I), ONE, A, LA, 
     &            POSELTT, NFRONT, 1, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &            MID_RANK, BUILDQ,
     &            .FALSE., MAXI_CLUSTER=MAXI_CLUSTER,
     &            DIAG=A(POSELTD), LD_DIAG=NFRONT, IW2=IW2,
     &            BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
            IF (IFLAG.LT.0) CYCLE
            CALL UPD_FLOP_UPDATE(BLR_L(J), BLR_L(I),
     &          MIDBLK_COMPRESS, MID_RANK, BUILDQ,
     &          (I.EQ.J), .FALSE.)
         ENDDO
#if defined(BLR_MT) 
!$OMP END DO
#endif
      END SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING_LDLT
      SUBROUTINE DMUMPS_BLR_SLV_UPD_TRAIL_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NCOL, NROW,
     &        A_BLOCFACTO, LA_BLOCFACTO, LD_BLOCFACTO, 
     &        BEGS_BLR_LM, NB_BLR_LM, BLR_LM, ISHIFT_LM,
     &        BEGS_BLR_LS, NB_BLR_LS, BLR_LS, ISHIFT_LS,
     &        CURRENT_BLR_LM, CURRENT_BLR_LS,
     &        IW2, BLOCK,
     &        MAXI_CLUSTER,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT
     &        )
!$    USE OMP_LIB      
      INTEGER(8), intent(in)  :: LA, LA_BLOCFACTO
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      DOUBLE PRECISION, intent(in)  :: A_BLOCFACTO(LA_BLOCFACTO)
      INTEGER(8), intent(in)  :: POSELT 
      INTEGER, intent(inout)    :: IFLAG, IERROR
      INTEGER, intent(in)     :: NCOL, NROW, IW2(*), TOL_OPT,
     &                           MAXI_CLUSTER, LD_BLOCFACTO
      INTEGER, intent(in)     :: NB_BLR_LM, NB_BLR_LS, 
     &                           ISHIFT_LM, ISHIFT_LS, 
     &                           CURRENT_BLR_LM, CURRENT_BLR_LS
      DOUBLE PRECISION, INTENT(INOUT) :: BLOCK(MAXI_CLUSTER,*)
      INTEGER, DIMENSION(:) :: BEGS_BLR_LM, BEGS_BLR_LS
      TYPE(LRB_TYPE),intent(in) :: BLR_LM(NB_BLR_LM-CURRENT_BLR_LM),
     &                             BLR_LS(NB_BLR_LS-CURRENT_BLR_LS)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER :: I, NB_BLOCKS_PANEL_LM, NB_BLOCKS_PANEL_LS, J, MID_RANK
      LOGICAL :: BUILDQ
      INTEGER :: OMP_NUM
      INTEGER :: IBIS
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELTT, POSELTD
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL_LM = NB_BLR_LM-CURRENT_BLR_LM
      NB_BLOCKS_PANEL_LS = NB_BLR_LS-CURRENT_BLR_LS
      POSELTD = 1_8 
      OMP_NUM = 0
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, OMP_NUM, MID_RANK, BUILDQ) 
#endif
      DO IBIS = 1, (NB_BLOCKS_PANEL_LS*NB_BLOCKS_PANEL_LM) 
        IF (IFLAG.LT.0) CYCLE
        I = (IBIS-1)/NB_BLOCKS_PANEL_LM+1
        J = IBIS - (I-1)*NB_BLOCKS_PANEL_LM
#if defined(BLR_MT)         
        OMP_NUM = 0 
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
#endif
            POSELTT = POSELT 
     &           + int(NCOL,8) * 
     &             int((BEGS_BLR_LS(CURRENT_BLR_LS+I)+ISHIFT_LS-1),8)
     &           + int((BEGS_BLR_LM(CURRENT_BLR_LM+J)+ISHIFT_LM-1),8)
            CALL DMUMPS_LRGEMM4(MONE,
     &            BLR_LM(J), BLR_LS(I), ONE, A, LA, 
     &            POSELTT, NCOL, 
     &            1, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &            MID_RANK, BUILDQ, 
     &            .FALSE., MAXI_CLUSTER=MAXI_CLUSTER,
     &            DIAG=A_BLOCFACTO, LD_DIAG=LD_BLOCFACTO, IW2=IW2,
     &            BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
            IF (IFLAG.LT.0) CYCLE
            CALL UPD_FLOP_UPDATE(BLR_LM(J), BLR_LS(I),
     &           MIDBLK_COMPRESS, MID_RANK, BUILDQ,
     &           .FALSE., .FALSE.)
         ENDDO
#if defined(BLR_MT) 
!$OMP END DO
         IF (IFLAG.LT.0) RETURN
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, POSELTT, MID_RANK, OMP_NUM, BUILDQ) 
#endif
         DO IBIS = 1, (NB_BLOCKS_PANEL_LS*(NB_BLOCKS_PANEL_LS+1)/2)
          IF (IFLAG.LT.0) CYCLE
          I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
          J = IBIS - I*(I-1)/2
#if defined(BLR_MT)         
          OMP_NUM = 0
!$        OMP_NUM = OMP_GET_THREAD_NUM() 
#endif
          POSELTT = POSELT 
     &        + int(NCOL,8) * 
     &          int((BEGS_BLR_LS(CURRENT_BLR_LS+I)+ISHIFT_LS-1),8)
     &        + int((NCOL-NROW+(BEGS_BLR_LS(CURRENT_BLR_LS+J)-1)),8)
          CALL DMUMPS_LRGEMM4(MONE,
     &            BLR_LS(J),BLR_LS(I), ONE, A, LA, 
     &            POSELTT, NCOL, 
     &            1, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &            MID_RANK, BUILDQ,
     &            .FALSE., MAXI_CLUSTER=MAXI_CLUSTER,
     &            DIAG=A_BLOCFACTO, LD_DIAG=LD_BLOCFACTO, IW2=IW2,
     &            BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
            IF (IFLAG.LT.0) CYCLE
            CALL UPD_FLOP_UPDATE(BLR_LS(J), BLR_LS(I),
     &            MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &            (I.EQ.J), .FALSE.)
      ENDDO
#if defined(BLR_MT) 
!$OMP END DO
#endif
      RETURN
      END SUBROUTINE DMUMPS_BLR_SLV_UPD_TRAIL_LDLT
      SUBROUTINE DMUMPS_BLR_UPD_NELIM_VAR_U(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, CURRENT_BLR, BLR_U, NB_BLR, 
     &        FIRST_BLOCK, IBEG_BLR, NPIV, NELIM)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR,
     &                                IBEG_BLR, NPIV, NELIM, FIRST_BLOCK
      INTEGER, intent(inout)         :: IFLAG, IERROR
      DOUBLE PRECISION, TARGET, intent(inout) :: A(LA)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_U(:)
      INTEGER, DIMENSION(:) :: BEGS_BLR
      TYPE(LRB_TYPE), POINTER :: LRB
      INTEGER :: IP
      INTEGER :: allocok
      INTEGER(8) :: LPOS, UPOS
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      IF (NELIM.NE.0) THEN 
        LPOS = POSELT + int(NFRONT,8)*int(NPIV,8) + int(IBEG_BLR-1,8)
#if defined(BLR_MT) 
!$OMP DO PRIVATE(LRB, UPOS)
#endif
        DO IP = FIRST_BLOCK, NB_BLR 
          IF (IFLAG.LT.0) CYCLE
          LRB => BLR_U(IP-CURRENT_BLR)
          UPOS = POSELT + int(NFRONT,8)*int(NPIV,8) 
     &                  + int(BEGS_BLR(IP)-1,8)
          IF (LRB%ISLR) THEN
             IF (LRB%K.GT.0) THEN
               allocate(TEMP_BLOCK( LRB%K, NELIM ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * LRB%K
                 GOTO 100
               ENDIF
               CALL dgemm('N', 'N', LRB%K, NELIM, LRB%N, ONE,
     &              LRB%R(1,1), LRB%K, A(LPOS), NFRONT,
     &              ZERO, TEMP_BLOCK, LRB%K) 
               CALL dgemm('N', 'N', LRB%M, NELIM, LRB%K, MONE,
     &              LRB%Q(1,1), LRB%M, TEMP_BLOCK, LRB%K,
     &              ONE, A(UPOS), NFRONT) 
               deallocate(TEMP_BLOCK)
             ENDIF
          ELSE
            CALL dgemm('N', 'N', LRB%M, NELIM, LRB%N, MONE,
     &              LRB%Q(1,1), LRB%M, A(LPOS), NFRONT,
     &              ONE, A(UPOS), NFRONT) 
          ENDIF
 100    CONTINUE
        ENDDO
#if defined(BLR_MT) 
!$OMP ENDDO
#endif
      ENDIF
      END SUBROUTINE DMUMPS_BLR_UPD_NELIM_VAR_U
      SUBROUTINE DMUMPS_BLR_UPD_NELIM_VAR_L(
     &        A_U, LA_U, UPOS, A_L, LA_L, LPOS, IFLAG, IERROR, LDU, LDL,
     &        BEGS_BLR_L, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        FIRST_BLOCK, NELIM, UTRANS)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA_U, LA_L
      INTEGER(8), intent(in)       :: UPOS, LPOS
      INTEGER, intent(in)          :: LDU, LDL, NB_BLR_L, CURRENT_BLR,
     &                                NELIM,  FIRST_BLOCK
      CHARACTER(len=1),INTENT(IN)  :: UTRANS
      INTEGER, intent(inout)         :: IFLAG, IERROR
      DOUBLE PRECISION, TARGET, intent(inout) :: A_L(LA_L), A_U(LA_U)
      TYPE(LRB_TYPE),intent(in) :: BLR_L(:)
      INTEGER                          :: BEGS_BLR_L(:)
      INTEGER :: I, NB_BLOCKS_PANEL_L, KL, ML, NL
      INTEGER :: allocok
      INTEGER(8) :: IPOS
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL_L = NB_BLR_L-CURRENT_BLR
      IF (NELIM.NE.0) THEN 
#if defined(BLR_MT) 
!$OMP DO PRIVATE(KL, ML, NL, IPOS)
#endif
        DO I = FIRST_BLOCK-CURRENT_BLR, NB_BLOCKS_PANEL_L
          IF (IFLAG.LT.0) CYCLE
          KL = BLR_L(I)%K 
          ML = BLR_L(I)%M 
          NL = BLR_L(I)%N 
          IPOS = LPOS + int(LDL,8) * 
     &        int(BEGS_BLR_L(CURRENT_BLR+I)-BEGS_BLR_L(CURRENT_BLR+1),8)
          IF (BLR_L(I)%ISLR) THEN
             IF (KL.GT.0) THEN
               allocate(TEMP_BLOCK( NELIM, KL ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * KL
                 write(*,*) 'Allocation problem in BLR routine 
     &         DMUMPS_BLR_UPD_NELIM_VAR_L: ',
     &         'not enough memory? memory requested = ', IERROR
                 GOTO 100
               ENDIF
               CALL dgemm(UTRANS , 'T' , NELIM, KL, NL , ONE ,
     &              A_U(UPOS) , LDU , BLR_L(I)%R(1,1) , KL ,
     &              ZERO , TEMP_BLOCK , NELIM) 
               CALL dgemm('N' , 'T' , NELIM , ML , KL , MONE ,
     &              TEMP_BLOCK , NELIM , BLR_L(I)%Q(1,1) , ML ,
     &              ONE , A_L(IPOS) , LDL) 
               deallocate(TEMP_BLOCK)
             ENDIF
          ELSE
            CALL dgemm(UTRANS , 'T' , NELIM, ML, NL , MONE ,
     &              A_U(UPOS) , LDU , BLR_L(I)%Q(1,1) , ML ,
     &              ONE , A_L(IPOS) , LDL) 
          ENDIF
 100      CONTINUE
        ENDDO
#if defined(BLR_MT) 
!$OMP ENDDO
#endif
      ENDIF
      END SUBROUTINE DMUMPS_BLR_UPD_NELIM_VAR_L
      SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR_L, BEGS_BLR_U, CURRENT_BLR, BLR_L, NB_BLR_L, 
     &        BLR_U,
     &        NB_BLR_U, NELIM, LBANDSLAVE, ISHIFT, NIV, SYM,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR_L, NB_BLR_U, 
     &                                CURRENT_BLR,
     &                                NELIM, NIV, SYM, TOL_OPT
      INTEGER, intent(inout)       :: IFLAG, IERROR
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER, intent(in)          :: ISHIFT
      DOUBLE PRECISION, intent(inout) :: A(LA)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_U(:)
      TYPE(LRB_TYPE),TARGET,intent(in) :: BLR_L(:)
      INTEGER :: BEGS_BLR_L(:), BEGS_BLR_U(:)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER :: I, NB_BLOCKS_PANEL_L, NB_BLOCKS_PANEL_U, 
     &           KL, ML, NL, J, IS, MID_RANK
      INTEGER :: allocok
      LOGICAL :: BUILDQ
      INTEGER :: OMP_NUM
      INTEGER :: IBIS
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELT_INCB, POSELT_TOP
      DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: TEMP_BLOCK
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL_L = NB_BLR_L-CURRENT_BLR
      NB_BLOCKS_PANEL_U = NB_BLR_U-CURRENT_BLR
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS = 0
      ENDIF
#if defined(BLR_MT)
!$OMP SINGLE
#endif
      IF (NELIM.NE.0) THEN 
         DO I = 1, NB_BLOCKS_PANEL_L
            KL = BLR_L(I)%K 
            ML = BLR_L(I)%M 
            NL = BLR_L(I)%N 
            IF (BLR_L(I)%ISLR) THEN
               IF (KL.GT.0) THEN
               allocate(TEMP_BLOCK( NELIM, KL ), stat=allocok )
               IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NELIM * KL
                 GOTO 100
               ENDIF
               POSELT_TOP  = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_U(CURRENT_BLR)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1,8)
               POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
               CALL dgemm('N' , 'T' , NELIM, KL, NL , ONE ,
     &                A(POSELT_TOP) , NFRONT , BLR_L(I)%R(1,1) , KL ,
     &                ZERO , TEMP_BLOCK , NELIM) 
               CALL dgemm('N' , 'T' , NELIM , ML , KL , MONE ,
     &                TEMP_BLOCK , NELIM , BLR_L(I)%Q(1,1) , ML ,
     &                ONE , A(POSELT_INCB) , NFRONT) 
               deallocate(TEMP_BLOCK)
               ENDIF
            ELSE
              POSELT_TOP  = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1)+IS-NELIM-1,8)
              POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+1) + IS - NELIM - 1, 8)
               CALL dgemm('N' , 'T' , NELIM, ML, NL , MONE ,
     &                A(POSELT_TOP) , NFRONT , BLR_L(I)%Q(1,1) , ML ,
     &                ONE , A(POSELT_INCB) , NFRONT) 
            ENDIF
         ENDDO
      ENDIF
 100  CONTINUE
#if defined(BLR_MT) 
!$OMP END SINGLE
#endif
      IF (IFLAG.LT.0) GOTO 200
      OMP_NUM = 0
#if defined(BLR_MT)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, POSELT_INCB, MID_RANK, BUILDQ)
#endif
      DO IBIS = 1, (NB_BLOCKS_PANEL_L*NB_BLOCKS_PANEL_U) 
        IF (IFLAG.LT.0) CYCLE
        I = (IBIS-1)/NB_BLOCKS_PANEL_U+1
        J = IBIS - (I-1)*NB_BLOCKS_PANEL_U
        POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR_L(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR_U(CURRENT_BLR+J) +IS - 1,8)
         CALL DMUMPS_LRGEMM4(MONE, BLR_U(J),
     &           BLR_L(I), ONE, A, LA, POSELT_INCB,
     &           NFRONT, 0, IFLAG, IERROR, 
     &           MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &           KPERCENT, MID_RANK, BUILDQ, .FALSE.)
         IF (IFLAG.LT.0) CYCLE
         CALL UPD_FLOP_UPDATE(BLR_U(J), BLR_L(I), 
     &           MIDBLK_COMPRESS, MID_RANK, BUILDQ,
     &           .FALSE., .FALSE.)
       ENDDO
#if defined(BLR_MT) 
!$OMP END DO
#endif
 200  CONTINUE
      END SUBROUTINE DMUMPS_BLR_UPDATE_TRAILING
      SUBROUTINE DMUMPS_BLR_UPD_PANEL_LEFT_LDLT(
     &        A, LA, POSELT, NFRONT, IWHANDLER,
     &        BEGS_BLR, CURRENT_BLR, NB_BLR, NPARTSASS,
     &        NELIM, IW2, BLOCK, ACC_LUA,
     &        MAXI_CLUSTER, MAXI_RANK, NIV, IFLAG, IERROR,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        K480, K479, K478, KPERCENT_LUA, KPERCENT,
     &        KEEP8,
     &        FIRST_BLOCK
     &        )
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, NPARTSASS,
     &                                CURRENT_BLR, IWHANDLER, TOL_OPT,
     &                                NELIM, NIV, K480, K479, K478,
     &                                MAXI_CLUSTER, MAXI_RANK,
     &                                KPERCENT_LUA, KPERCENT
      DOUBLE PRECISION, intent(inout)    :: A(LA)
      INTEGER, intent(in) :: IW2(*)
      DOUBLE PRECISION :: BLOCK(MAXI_CLUSTER,*)
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:)
      INTEGER(8) :: KEEP8(150)
      INTEGER, DIMENSION(:) :: BEGS_BLR
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER,intent(inout) :: IFLAG, IERROR
      INTEGER,OPTIONAL,intent(in) :: FIRST_BLOCK
      TYPE(LRB_TYPE), POINTER :: BLR_L(:), NEXT_BLR_L(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK, FRFR_UPDATES,
     &           I, II, J, JJ, NB_BLOCKS_PANEL, IND_U, IND_L, K_MAX,
     &           MAXRANK, NB_DEC, FR_RANK
      INTEGER :: MID_RANK, allocok
      INTEGER :: J_ORDER(CURRENT_BLR), J_RANK(CURRENT_BLR)
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      LOGICAL :: BUILDQ, COMPRESSED_FR
      INTEGER :: OFFSET_IW
      INTEGER :: OMP_NUM
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELT_INCB, POSELTD
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      ACC_LRB => ACC_LUA(1)
      IF (K480.GE.5) THEN
        IF (NB_BLOCKS_PANEL.GT.1) THEN
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         0, 
     &         CURRENT_BLR+1, NEXT_BLR_L)
        ENDIF
        IF (.not.(present(FIRST_BLOCK))) THEN
          write(*,*) "Internal error in 
     &      DMUMPS_BLR_UPD_PANEL_LEFT_LDLT: KEEP(480)=",K480,
     &      ">= 5, but FIRST_BLOCK argument is missing"
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      OMP_NUM = 0
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, JJ, POSELT_INCB, MID_RANK, BUILDQ, K_MAX,
!$OMP&         BLR_L, OMP_NUM, J_ORDER, J_RANK,
!$OMP&         IND_U, IND_L, ACC_LRB, POSELTD, NB_DEC,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK,
!$OMP&         FRFR_UPDATES, COMPRESSED_FR, FR_RANK, II, OFFSET_IW)
#endif
      DO I = 1, NB_BLOCKS_PANEL
#if defined(BLR_MT)         
        IF (IFLAG.LT.0) CYCLE
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
        POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR(CURRENT_BLR+1)-1,8)
        ACC_LRB%N = BEGS_BLR(CURRENT_BLR+I+1)-BEGS_BLR(CURRENT_BLR+I)
        ACC_LRB%M = BEGS_BLR(CURRENT_BLR+2)-BEGS_BLR(CURRENT_BLR+1)
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        COMPRESSED_FR = .FALSE.
        IF (K480.EQ.2) THEN
          DO J = 1, CURRENT_BLR
            J_ORDER(J) = J 
          ENDDO
        ELSE
          CALL DMUMPS_GET_LUA_ORDER(CURRENT_BLR, J_ORDER, J_RANK,
     &                              IWHANDLER, 
     &                              1, 0, I, 0, 
     &                              FRFR_UPDATES)
        ENDIF
        FR_RANK = 0
        IF ((K480.GE.5).AND.(I.NE.1)) THEN
          IF (I.GT.FIRST_BLOCK) THEN
            IF (FRFR_UPDATES.EQ.0) THEN
              CALL DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &               COMPRESSED_FR, 0, .FALSE.)
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              FR_RANK = ACC_LRB%K
            ENDIF
          ENDIF
        ENDIF
        NB_DEC = FRFR_UPDATES
        DO JJ = 1, CURRENT_BLR
          J = J_ORDER(JJ)
          K_MAX = J_RANK(JJ)
          POSELTD = POSELT + int(NFRONT,8) * int(BEGS_BLR(J)-1,8)
     &          + int(BEGS_BLR(J) - 1,8)
          OFFSET_IW = BEGS_BLR(J)
          IND_L = CURRENT_BLR+I-J
          IND_U = CURRENT_BLR+1-J
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IWHANDLER,
     &           0, 
     &           J, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (K480.GE.3) THEN
            IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
              NB_DEC = JJ-1
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, 0)
              COMPRESSED_FR = .FALSE.
              MAX_ACC_RANK = 0
            ENDIF
            OLD_ACC_RANK = ACC_LRB%K
          ENDIF
          CALL DMUMPS_LRGEMM4(MONE,
     &              BLR_L(IND_U), BLR_L(IND_L), ONE,
     &              A, LA, POSELT_INCB,
     &              NFRONT, 1, IFLAG, IERROR, 
     &              MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &              KPERCENT_RMB, MID_RANK, BUILDQ, 
     &              (K480.GE.3), LorU=0, 
     &              LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &              MAXI_CLUSTER=MAXI_CLUSTER,
     &              DIAG=A(POSELTD), LD_DIAG=NFRONT, 
     &              IW2=IW2(OFFSET_IW),
     &              BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
          IF (IFLAG.LT.0) GOTO 100
          CALL UPD_FLOP_UPDATE(BLR_L(IND_U), 
     &      BLR_L(IND_L), MIDBLK_COMPRESS, 
     &      MID_RANK, BUILDQ, (I.EQ.1), (K480.GE.3))
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            J_RANK(JJ) = MID_RANK
          ENDIF
          IF (K480.GE.3) THEN
            NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
            MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
            IF (K480.EQ.4) THEN
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
                IF (ACC_LRB%K.GT.0) THEN
                  CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &                 NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &                 TOL_OPT,
     &                 KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  MAX_ACC_RANK = ACC_LRB%K
                ENDIF
              ENDIF
            ENDIF
            IF ((K480.GE.5).AND.(I.NE.1)) THEN
              IF (I.GT.FIRST_BLOCK) THEN
                IF (JJ.EQ.FRFR_UPDATES) THEN
                  CALL DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &                NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                COMPRESSED_FR, 0, .FALSE.)
                  MAX_ACC_RANK = ACC_LRB%K
                  NEW_ACC_RANK = ACC_LRB%K
                  IF (COMPRESSED_FR) THEN
                    J_RANK(JJ) = ACC_LRB%K 
                    NB_DEC = FRFR_UPDATES-1 
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF (K480.GE.3) THEN
          IF ((K480.GE.5)) THEN
            IF (COMPRESSED_FR.OR.(K480.GE.6)) THEN  
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
                  IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                    CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
                  IF (FRFR_UPDATES.GT.0) THEN
                    allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR-NB_DEC
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'DMUMPS_BLR_UPD_PANEL_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,CURRENT_BLR-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
                  ELSE
                     allocate(POS_LIST(CURRENT_BLR+1),stat=allocok)
                     IF (allocok .GT. 0) THEN
                        IFLAG  = -13
                        IERROR = CURRENT_BLR+1
                        write(*,*) 'Allocation problem in BLR routine ',
     &                       'DMUMPS_BLR_UPD_PANEL_LEFT_LDLT: ',
     &                       'not enough memory? memory requested = ',
     &                       IERROR
                        GOTO 100
                     ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,CURRENT_BLR
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(CURRENT_BLR+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'DMUMPS_BLR_UPD_PANEL_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,CURRENT_BLR+1
                      RANK_LIST(II) = J_RANK(II-1)
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               CURRENT_BLR+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = floor(dble(ACC_LRB%M*ACC_LRB%N)/dble(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN  
              CALL ALLOC_LRB_FROM_ACC(ACC_LRB, NEXT_BLR_L(I-1), 
     &                       ACC_LRB%K, ACC_LRB%M, ACC_LRB%N, 0,
     &                       IFLAG, IERROR, KEEP8)
              IF (IFLAG.LT.0) CYCLE
              ACC_LRB%K = 0
            ELSE
              IF (I.NE.1) NEXT_BLR_L(I-1)%ISLR=.FALSE.
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, 0)
            ENDIF
          ELSE 
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
              IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &        THEN
               allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
               IF (allocok .GT. 0) THEN
                  IFLAG  = -13
                  IERROR = CURRENT_BLR-NB_DEC
                  GOTO 100
               ENDIF
              POS_LIST(1) = 1
              DO II = 1,CURRENT_BLR-NB_DEC-1
                POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
              ENDDO
              CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, 0)
          ENDIF                               
        ENDIF
 100    CONTINUE        
      ENDDO                                 
#if defined(BLR_MT) 
!$OMP END DO                                
#endif
      END SUBROUTINE DMUMPS_BLR_UPD_PANEL_LEFT_LDLT
      SUBROUTINE DMUMPS_BLR_UPD_PANEL_LEFT(
     &        A, LA, POSELT, NFRONT, IWHANDLER, LorU,
     &        BEGS_BLR, BEGS_BLR_U, CURRENT_BLR, ACC_LUA,
     &        NB_BLR, NPARTSASS, NELIM, NIV, SYM, 
     &        LBANDSLAVE, IFLAG, IERROR, ISHIFT,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        K480, K479, K478, KPERCENT_LUA, KPERCENT,
     &        MAXI_CLUSTER, MAXI_RANK, 
     &        K474, FSorCB, BLR_U_COL, KEEP8,
     &        FIRST_BLOCK, BEG_I_IN, END_I_IN)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, NPARTSASS,
     &                                CURRENT_BLR, IWHANDLER, LorU,
     &                                NELIM, NIV, SYM, K480, K479, K478,
     &                                MAXI_CLUSTER, MAXI_RANK,
     &                                KPERCENT_LUA, KPERCENT, ISHIFT,
     &                                K474, FSorCB
      LOGICAL, intent(in)          :: LBANDSLAVE
      DOUBLE PRECISION, TARGET, intent(inout) :: A(LA)
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:), BLR_U_COL(:)
      INTEGER(8) :: KEEP8(150)
      INTEGER, DIMENSION(:) :: BEGS_BLR, BEGS_BLR_U
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB, TOL_OPT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER,intent(inout) :: IFLAG, IERROR
      INTEGER,OPTIONAL,intent(in) :: FIRST_BLOCK
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN
      TYPE(LRB_TYPE), POINTER :: BLR_U(:), BLR_L(:), NEXT_BLR(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK, FRFR_UPDATES,
     &           NB_DEC, FR_RANK, MAXRANK, BEG_I, END_I
      INTEGER :: I,II,J,JJ, NB_BLOCKS_PANEL, IND_U, IND_L, K_MAX
      INTEGER :: MID_RANK, allocok
      INTEGER :: J_ORDER(CURRENT_BLR), J_RANK(CURRENT_BLR)
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      LOGICAL :: BUILDQ, COMPRESSED_FR
#if defined(BLR_MT)
      INTEGER :: OMP_NUM
      INTEGER :: CHUNK
#endif
      INTEGER(8) :: POSELT_INCB
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      IF (NIV.EQ.2.AND.LorU.EQ.0) THEN
        IF (LBANDSLAVE) THEN
          NB_BLOCKS_PANEL = NB_BLR
        ELSE
          NB_BLOCKS_PANEL = NPARTSASS-CURRENT_BLR
        ENDIF
      ELSE
        NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      ENDIF
      ACC_LRB => ACC_LUA(1)
      IF (K480.GE.5) THEN
        IF (NB_BLOCKS_PANEL.GT.1) THEN
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         LorU, 
     &         CURRENT_BLR+1, NEXT_BLR)
        ENDIF
        IF (.not.(present(FIRST_BLOCK))) THEN
          write(*,*) "Internal error in 
     &      DMUMPS_BLR_UPD_PANEL_LEFT: KEEP(480)=",K480,
     &      ">=5, but FIRST_BLOCK argument is missing"
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      IF (LorU.EQ.0) THEN 
          BEG_I = 1
      ELSE 
          BEG_I = 2
      ENDIF
      END_I = NB_BLOCKS_PANEL
      IF (K474.EQ.3) THEN
        IF(present(BEG_I_IN)) THEN
          BEG_I = BEG_I_IN - CURRENT_BLR
        ENDIF
        IF(present(END_I_IN)) THEN
          END_I = END_I_IN - CURRENT_BLR
        ENDIF
      ENDIF
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!$OMP& PRIVATE(I, J, JJ, POSELT_INCB, MID_RANK, BUILDQ,
!$OMP&         BLR_U, BLR_L, J_ORDER, J_RANK, K_MAX,
!$OMP&         IND_U, IND_L, OMP_NUM, ACC_LRB,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK,
!$OMP&         FRFR_UPDATES, FR_RANK, COMPRESSED_FR)
#endif
      DO I = BEG_I, END_I
        IF (IFLAG.LT.0) CYCLE
#if defined(BLR_MT)         
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
        IF (LorU.EQ.0) THEN 
          IF (LBANDSLAVE) THEN
            POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(I+1)-1),8)
     &           + int(BEGS_BLR_U(2)+ISHIFT-1,8)
            ACC_LRB%N = BEGS_BLR(I+2)-BEGS_BLR(I+1)
            ACC_LRB%M = BEGS_BLR_U(3)-BEGS_BLR_U(2)
            IF (K474.GE.2) THEN
              BLR_U => BLR_U_COL
            ENDIF
          ELSE
            POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(CURRENT_BLR+I)-1),8)
     &           + int(BEGS_BLR(CURRENT_BLR+1)-1,8)
            ACC_LRB%N = BEGS_BLR(CURRENT_BLR+I+1)
     &                 -BEGS_BLR(CURRENT_BLR+I)
            ACC_LRB%M = BEGS_BLR(CURRENT_BLR+2)-BEGS_BLR(CURRENT_BLR+1)
          ENDIF
        ELSE 
          POSELT_INCB = POSELT 
     &           + int(NFRONT,8) * int((BEGS_BLR(CURRENT_BLR+1)-1),8)
     &           + int(BEGS_BLR(CURRENT_BLR+I)-1,8)
          ACC_LRB%N = BEGS_BLR(CURRENT_BLR+2)-BEGS_BLR(CURRENT_BLR+1)
          ACC_LRB%M = BEGS_BLR(CURRENT_BLR+I+1)-BEGS_BLR(CURRENT_BLR+I)
        ENDIF
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        COMPRESSED_FR = .FALSE.
        IF (K480.EQ.2) THEN
          DO J = 1, CURRENT_BLR
            J_ORDER(J) = J 
          ENDDO
        ELSE
          CALL DMUMPS_GET_LUA_ORDER(CURRENT_BLR, J_ORDER, J_RANK, 
     &                              IWHANDLER, 
     &                              0, 0, I, LorU, 
     &                              FRFR_UPDATES, 
     &                              LBANDSLAVE, K474, BLR_U_COL)
        ENDIF
        FR_RANK = 0
        IF ((K480.GE.5).AND.(I.NE.1)) THEN
          IF (I.GT.FIRST_BLOCK) THEN
            IF (FRFR_UPDATES.EQ.0) THEN
              CALL DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &                NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                COMPRESSED_FR, LorU, .FALSE.)
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              FR_RANK = ACC_LRB%K
            ENDIF
          ENDIF
        ENDIF
        NB_DEC = FRFR_UPDATES
        DO JJ = 1, CURRENT_BLR
          J = J_ORDER(JJ)
          K_MAX = J_RANK(JJ)
          IF (LorU.EQ.0) THEN 
            IF (LBANDSLAVE) THEN
              IND_L = I
              IF (K474.LT.2) THEN
                IND_U = CURRENT_BLR+1-J
              ELSE
                IND_U = J
              ENDIF
            ELSE
              IND_L = CURRENT_BLR+I-J
              IND_U = CURRENT_BLR+1-J
            ENDIF
          ELSE 
            IND_L = CURRENT_BLR+1-J
            IND_U = CURRENT_BLR+I-J
          ENDIF
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IWHANDLER,
     &           0, 
     &           J, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (.NOT.LBANDSLAVE.OR.K474.LT.2) THEN
            CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &           IWHANDLER,
     &           1, 
     &           J, BLR_U)
          ENDIF
        IF (K480.GE.3) THEN
          IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
            NB_DEC = JJ-1
            CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB, MAXI_CLUSTER, 
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, LorU)
            COMPRESSED_FR = .FALSE.
            MAX_ACC_RANK = 0
          ENDIF
          OLD_ACC_RANK = ACC_LRB%K
        ENDIF
          CALL DMUMPS_LRGEMM4(MONE,
     &            BLR_U(IND_U), BLR_L(IND_L), ONE,
     &            A, LA, POSELT_INCB,
     &            NFRONT, 0, IFLAG, IERROR, 
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &            KPERCENT_RMB, MID_RANK, BUILDQ,
     &            (K480.GE.3), LorU=LorU, 
     &            LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &            MAXI_CLUSTER=MAXI_CLUSTER
     &            )
        IF (IFLAG.LT.0) GOTO 100
        CALL UPD_FLOP_UPDATE(BLR_U(IND_U), BLR_L(IND_L),
     &            MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &            .FALSE., (K480.GE.3))
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            J_RANK(JJ) = MID_RANK
          ENDIF
        IF (K480.GE.3) THEN
          NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
          MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
          IF (K480.EQ.4) THEN
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
              CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV,
     &            MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &            KPERCENT_LUA, NEW_ACC_RANK)
              MAX_ACC_RANK = ACC_LRB%K
            ENDIF
          ENDIF
        ENDIF
        IF ((K480.GE.5).AND.(I.NE.1)) THEN
          IF (I.GT.FIRST_BLOCK) THEN
            IF (JJ.EQ.FRFR_UPDATES) THEN
              CALL DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT,
     &               COMPRESSED_FR, LorU, .FALSE.)
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              IF (COMPRESSED_FR) THEN
                J_RANK(JJ) = ACC_LRB%K 
                NB_DEC = FRFR_UPDATES-1 
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        ENDDO
        IF (K480.GE.3) THEN
          IF ((K480.GE.5)) THEN
            IF (COMPRESSED_FR.OR.(K480.GE.6)) THEN  
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
                  IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                    CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
                  IF (FRFR_UPDATES.GT.0) THEN
                    allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR-NB_DEC
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'DMUMPS_BLR_UPD_PANEL_LEFT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,CURRENT_BLR-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
                  ELSE
                     allocate(POS_LIST(CURRENT_BLR+1),stat=allocok)
                     IF (allocok .GT. 0) THEN
                        IFLAG  = -13
                        IERROR = CURRENT_BLR+1
                        write(*,*) 'Allocation problem in BLR routine ',
     &                       'DMUMPS_BLR_UPD_PANEL_LEFT: ',
     &                       'not enough memory? memory requested = ',
     &                       IERROR
                        GOTO 100
                     ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,CURRENT_BLR
                      POS_LIST(II+1)=POS_LIST(II)+J_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(CURRENT_BLR+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = CURRENT_BLR+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'DMUMPS_BLR_UPD_PANEL_LEFT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,CURRENT_BLR+1
                      RANK_LIST(II) = J_RANK(II-1)
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               CURRENT_BLR+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = FLOOR(dble(ACC_LRB%M*ACC_LRB%N)/dble(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN  
              CALL ALLOC_LRB_FROM_ACC(ACC_LRB, NEXT_BLR(I-1), 
     &                    ACC_LRB%K, ACC_LRB%M, ACC_LRB%N, LorU,
     &                    IFLAG, IERROR, KEEP8)
              IF (IFLAG.LT.0) CYCLE
              ACC_LRB%K = 0
            ELSE
              IF (I.NE.1) NEXT_BLR(I-1)%ISLR=.FALSE.
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, LorU)
            ENDIF
          ELSE 
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
              IF (CURRENT_BLR-FRFR_UPDATES.GT.1) THEN
                CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &      THEN
              allocate(POS_LIST(CURRENT_BLR-NB_DEC),stat=allocok)
              IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = CURRENT_BLR-NB_DEC
                 GOTO 100
              ENDIF
              POS_LIST(1) = 1
              DO II = 1,CURRENT_BLR-NB_DEC-1
              POS_LIST(II+1)=POS_LIST(II)+J_RANK(NB_DEC+II) 
              ENDDO
              CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_INCB, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               J_RANK(NB_DEC+1:CURRENT_BLR), POS_LIST,
     &               CURRENT_BLR-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_INCB, NFRONT, NIV, LorU)
          ENDIF
        ENDIF
 100    CONTINUE        
      ENDDO
#if defined(BLR_MT) 
!$OMP END DO NOWAIT
#endif
      END SUBROUTINE DMUMPS_BLR_UPD_PANEL_LEFT
      SUBROUTINE DMUMPS_BLR_UPD_CB_LEFT_LDLT(A, LA, POSELT, NFRONT,
     &        BEGS_BLR, BEGS_BLR_DYN, NB_INCB, NB_INASM, NASS,      
     &        IWHANDLER, 
     &        IW2, BLOCK, ACC_LUA,
     &        MAXI_CLUSTER, MAXI_RANK, NIV, IFLAG, IERROR,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        K480, K479, K478, KPERCENT_LUA, KPERCENT, KEEP8)
!$    USE OMP_LIB 
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_INCB, NB_INASM
      INTEGER, INTENT(IN)          :: NIV, IWHANDLER, MAXI_CLUSTER,
     &                                MAXI_RANK, K480, K479, K478, NASS,
     &                                KPERCENT_LUA, KPERCENT
      INTEGER, intent(inout)         :: IFLAG, IERROR
      INTEGER(8) :: KEEP8(150)
      INTEGER, DIMENSION(:) :: BEGS_BLR
      INTEGER, DIMENSION(:) :: BEGS_BLR_DYN
      DOUBLE PRECISION, INTENT(INOUT) :: BLOCK(MAXI_CLUSTER,*)
      INTEGER, intent(in) :: IW2(*)
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB, TOL_OPT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER :: M, N, allocok
      INTEGER :: I, II, J, K, KK, IND_L, IND_U, K_MAX, IBIS,
     &           K_ORDER(NB_INASM), K_RANK(NB_INASM), NB_DEC
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      INTEGER(8) :: POSELT_BLOCK, POSELTD
      INTEGER :: NCB, MID_RANK, FRFR_UPDATES, MAXRANK, FR_RANK
      LOGICAL :: BUILDQ, COMPRESSED_FR
      TYPE(LRB_TYPE), POINTER :: BLR_L(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK
      INTEGER :: OFFSET_IW
      INTEGER :: OMP_NUM
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      NCB = NFRONT - NASS 
      ACC_LRB => ACC_LUA(1)
      OMP_NUM = 0
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, K, KK, POSELT_BLOCK, MID_RANK, BUILDQ,
!$OMP&         BLR_L, IND_U, IND_L, M, N, K_ORDER, K_RANK,
!$OMP&         K_MAX, OMP_NUM, ACC_LRB, POSELTD,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK, 
!$OMP&         FRFR_UPDATES, FR_RANK, NB_DEC, II)
#endif
      DO IBIS = 1,NB_INCB*(NB_INCB+1)/2
        IF (IFLAG.LT.0) CYCLE     
        I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
        J = IBIS - I*(I-1)/2
        I = I+NB_INASM
        J = J+NB_INASM
#if defined(BLR_MT)         
        OMP_NUM = 0
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        M = BEGS_BLR(I+1)-BEGS_BLR(I)
        N = BEGS_BLR(J+1)-BEGS_BLR(J)
        POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I)-1,8) + 
     &           int(BEGS_BLR(J)-1,8)
        ACC_LRB%M = N
        ACC_LRB%N = M
        IF (K480.EQ.2) THEN
          DO K = 1, NB_INASM
            K_ORDER(K) = K 
          ENDDO
        ELSE
          CALL DMUMPS_GET_LUA_ORDER(NB_INASM, K_ORDER, K_RANK, 
     &                             IWHANDLER, 
     &                             1, 1, I, J, 
     &                             FRFR_UPDATES)
        ENDIF
        FR_RANK = 0
        IF ((K480.GE.5).AND.(I.NE.J)) THEN
          IF (FRFR_UPDATES.EQ.0) THEN
            CALL DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                COMPRESSED_FR, 0, .TRUE.)
            FR_RANK = ACC_LRB%K
            MAX_ACC_RANK = ACC_LRB%K
            NEW_ACC_RANK = ACC_LRB%K
          ENDIF
        ENDIF
        NB_DEC = FRFR_UPDATES
        DO KK = 1, NB_INASM
          K = K_ORDER(KK)
          K_MAX = K_RANK(KK)
          POSELTD = POSELT + int(NFRONT,8) * int(BEGS_BLR_DYN(K)-1,8)
     &         + int(BEGS_BLR_DYN(K) - 1,8)
          OFFSET_IW = BEGS_BLR_DYN(K)
          IND_L = I-K
          IND_U = J-K
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IWHANDLER,
     &             0, 
     &             K, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (K480.GE.3) THEN 
            IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
              NB_DEC = KK-1
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB, 
     &              MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK, 
     &              NFRONT, NIV, 2)
              COMPRESSED_FR = .FALSE.
              MAX_ACC_RANK = 0
            ENDIF
            OLD_ACC_RANK = ACC_LRB%K
          ENDIF
          CALL DMUMPS_LRGEMM4(MONE,
     &              BLR_L(IND_U), BLR_L(IND_L), ONE,
     &              A, LA, POSELT_BLOCK,
     &              NFRONT, 1, IFLAG, IERROR, 
     &              MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &              KPERCENT_RMB, MID_RANK, BUILDQ,
     &              (K480.GE.3), LorU=2, 
     &              LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &              MAXI_CLUSTER=MAXI_CLUSTER,
     &              DIAG=A(POSELTD), LD_DIAG=NFRONT, 
     &              IW2=IW2(OFFSET_IW),
     &              BLOCK=BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1))
          IF (IFLAG.LT.0) GOTO 100
          CALL UPD_FLOP_UPDATE(BLR_L(IND_U), BLR_L(IND_L), 
     &           MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &           (I.EQ.J), (K480.GE.3))
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            K_RANK(KK) = MID_RANK
          ENDIF
          IF (K480.GE.3) THEN
            NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
            MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
            IF (K480.EQ.4) THEN
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
                IF (ACC_LRB%K.GT.0) THEN
                  CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &               TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  MAX_ACC_RANK = ACC_LRB%K
                ENDIF
              ENDIF
            ENDIF
            IF ((K480.GE.5).AND.(I.NE.J)) THEN
              IF (KK.EQ.FRFR_UPDATES) THEN
                CALL DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                   MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                   NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                   COMPRESSED_FR, 0, .TRUE.)
                IF (COMPRESSED_FR) THEN
                  K_RANK(KK) = ACC_LRB%K 
                  NB_DEC = FRFR_UPDATES-1 
                ENDIF
                MAX_ACC_RANK = ACC_LRB%K
                NEW_ACC_RANK = ACC_LRB%K
              ENDIF
            ENDIF
          ENDIF
        END DO
        IF (K480.GE.3) THEN
          IF ((K480.GE.5)) THEN
            IF (COMPRESSED_FR.OR.(K480.GE.6)) THEN  
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
                  IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                    CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &               TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
                  IF (FRFR_UPDATES.GT.0) THEN
                    allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM-NB_DEC
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'DMUMPS_BLR_UPD_CB_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,NB_INASM-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
                  ELSE
                    allocate(POS_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'DMUMPS_BLR_UPD_CB_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,NB_INASM
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       write(*,*) 'Allocation problem in BLR routine ',
     &                      'DMUMPS_BLR_UPD_CB_LEFT_LDLT: ',
     &                      'not enough memory? memory requested = ',
     &                      IERROR
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,NB_INASM+1
                      RANK_LIST(II) = K_RANK(II-1)
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               NB_INASM+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = FLOOR(dble(ACC_LRB%M*ACC_LRB%N)/dble(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN  
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2,
     &              COUNT_FLOPS=.FALSE.)
            ELSE
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
            ENDIF
          ELSE 
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
              IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &           THEN
              allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
              IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NB_INASM-NB_DEC
                 GOTO 100
              ENDIF
              POS_LIST(1) = 1
              DO II = 1,NB_INASM-NB_DEC-1
                POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
              ENDDO
              CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               KEEP8, NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, 
     &               TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
          ENDIF                               
        ENDIF
 100    CONTINUE     
      END DO
#if defined(BLR_MT) 
!$OMP END DO 
#endif
      END SUBROUTINE DMUMPS_BLR_UPD_CB_LEFT_LDLT
      SUBROUTINE DMUMPS_BLR_UPD_CB_LEFT(A, LA, POSELT, NFRONT,
     &        BEGS_BLR, BEGS_BLR_U, NB_ROWS, NB_INCB, NB_INASM, NASS,
     &        IWHANDLER, NIV, LBANDSLAVE, IFLAG, IERROR,
     &        MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB,
     &        ACC_LUA, K480, K479, K478, KPERCENT_LUA, KPERCENT,
     &        MAXI_CLUSTER, MAXI_RANK,
     &        K474, FSorCB, BLR_U_COL, COMPRESS_CB, CB_LRB, KEEP8)
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_ROWS, NB_INCB, NB_INASM
      INTEGER, INTENT(IN)          :: NIV, IWHANDLER, MAXI_CLUSTER,
     &                                MAXI_RANK, KPERCENT_LUA, KPERCENT
      INTEGER, INTENT(IN)          :: K480, K479, K478, NASS, K474,
     &                                FSorCB
      INTEGER, intent(inout)         :: IFLAG, IERROR
      INTEGER, DIMENSION(:) :: BEGS_BLR, BEGS_BLR_U
#if defined(MUMPS_F2003)
      TYPE(LRB_TYPE), POINTER, intent(inout) :: CB_LRB(:,:)
#else
      TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:)
#endif
      TYPE(LRB_TYPE), POINTER :: ACC_LUA(:), BLR_U_COL(:)
      INTEGER(8) :: KEEP8(150)
      INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB, TOL_OPT
      DOUBLE PRECISION,intent(in) :: TOLEPS
      LOGICAL, intent(in) :: LBANDSLAVE, COMPRESS_CB
      INTEGER :: M, N, allocok
      INTEGER :: I, II, J, K, KK, IND_L, IND_U, IBIS,
     &           K_ORDER(NB_INASM), K_RANK(NB_INASM)
      INTEGER, ALLOCATABLE :: POS_LIST(:), RANK_LIST(:)
      INTEGER(8) :: POSELT_BLOCK
      INTEGER :: MID_RANK, K_MAX, FRFR_UPDATES, NB_DEC
      INTEGER :: FRONT_CB_BLR_SAVINGS
      LOGICAL :: BUILDQ, COMPRESSED_FR
      TYPE(LRB_TYPE), POINTER :: BLR_U(:), BLR_L(:)
      TYPE(LRB_TYPE), POINTER :: ACC_LRB, LRB
      INTEGER :: OLD_ACC_RANK, MAX_ACC_RANK, NEW_ACC_RANK, MAXRANK,
     &           FR_RANK
#if defined(BLR_MT)
      INTEGER :: OMP_NUM
      INTEGER :: CHUNK
#endif
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      ACC_LRB => ACC_LUA(1)
      FRONT_CB_BLR_SAVINGS = 0
#if defined(BLR_MT)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, K, KK, POSELT_BLOCK, MID_RANK, BUILDQ,
!$OMP&         BLR_U, BLR_L, IND_U, IND_L, M, N, 
!$OMP&         ACC_LRB, OMP_NUM, K_MAX, K_ORDER, K_RANK,
!$OMP&         MAX_ACC_RANK, OLD_ACC_RANK, NEW_ACC_RANK,
!$OMP&         FRFR_UPDATES, LRB)
#endif
      DO IBIS = 1,NB_ROWS*NB_INCB
        IF (IFLAG.LT.0) CYCLE     
        I = (IBIS-1)/NB_INCB+1
        J = IBIS - (I-1)*NB_INCB
        IF (.NOT.LBANDSLAVE) THEN
          I = I+NB_INASM
        ENDIF
        J = J+NB_INASM
#if defined(BLR_MT)         
        OMP_NUM=0 
!$      OMP_NUM = OMP_GET_THREAD_NUM() 
        ACC_LRB => ACC_LUA(OMP_NUM+1)
#endif
        MAX_ACC_RANK = 0
        NEW_ACC_RANK = 0
        IF (LBANDSLAVE) THEN
          M = BEGS_BLR(I+2)-BEGS_BLR(I+1)
          IF (K474.EQ.1) THEN
            POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I+1)-1,8)
     &           +int(NASS,8) + int(BEGS_BLR_U(J-NB_INASM+1)-1,8)
            N = BEGS_BLR_U(J-NB_INASM+2)-BEGS_BLR_U(J-NB_INASM+1)
          ELSEIF (K474.GE.2) THEN
            BLR_U => BLR_U_COL
            POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I+1)-1,8)
     &          + int(NASS-1,8)
            N = BEGS_BLR_U(3)-BEGS_BLR_U(2)
          ELSE
            write(*,*) 'Internal error in DMUMPS_BLR_UPD_CB_LEFT',
     &        LBANDSLAVE,K474       
            CALL MUMPS_ABORT()
          ENDIF
        ELSE
          M = BEGS_BLR(I+1)-BEGS_BLR(I)
          POSELT_BLOCK = POSELT + int(NFRONT,8)*int(BEGS_BLR(I)-1,8) + 
     &           int(BEGS_BLR_U(J)-1,8)
          N = BEGS_BLR_U(J+1)-BEGS_BLR_U(J)
        ENDIF
        ACC_LRB%M = N
        ACC_LRB%N = M
        IF (K480.EQ.2) THEN
          DO K = 1, NB_INASM
            K_ORDER(K) = K 
          ENDDO
        ELSE
          CALL DMUMPS_GET_LUA_ORDER(NB_INASM, K_ORDER, K_RANK,
     &                             IWHANDLER, 
     &                             0, 1, I, J, 
     &                             FRFR_UPDATES,
     &                             LBANDSLAVE, K474, BLR_U_COL)
        ENDIF
        COMPRESSED_FR = .FALSE.
        FR_RANK = 0
        DO KK = 1, NB_INASM
          IF ((K480.GE.5.OR.COMPRESS_CB).AND.I.NE.J) THEN
            IF (KK-1.EQ.FRFR_UPDATES) THEN
              CALL DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                 NFRONT, NIV, TOLEPS, TOL_OPT, KPERCENT, 
     &                 COMPRESSED_FR, 0, .TRUE.)
              IF (COMPRESSED_FR) THEN
                K_RANK(KK) = ACC_LRB%K 
                NB_DEC = FRFR_UPDATES-1 
              ENDIF
              MAX_ACC_RANK = ACC_LRB%K
              NEW_ACC_RANK = ACC_LRB%K
              FR_RANK = ACC_LRB%K
            ENDIF
          ENDIF
          K = K_ORDER(KK)
          K_MAX = K_RANK(KK)
          IF (LBANDSLAVE) THEN
            IND_L = I
            IF (K474.LT.2) THEN
              IND_U = J-K
            ELSE
              IND_U = K
            ENDIF
          ELSE
            IND_L = I-K
            IND_U = J-K
          ENDIF
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IWHANDLER,
     &             0, 
     &             K, BLR_L)
          IF (BLR_L(IND_L)%M.EQ.0) THEN
            CYCLE
          ENDIF
          IF (.NOT.LBANDSLAVE.OR.K474.LT.2) THEN
            CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IWHANDLER,
     &             1, 
     &             K, BLR_U)
          ENDIF
          IF (K480.GE.3) THEN
            IF (ACC_LRB%K+K_MAX.GT.MAXI_RANK) THEN
              COMPRESSED_FR = .FALSE.
              NB_DEC = KK-1
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB, 
     &              MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK, 
     &              NFRONT, NIV, 2)
              MAX_ACC_RANK = 0
            ENDIF
            OLD_ACC_RANK = ACC_LRB%K
          ENDIF
            CALL DMUMPS_LRGEMM4(MONE,
     &              BLR_U(IND_U), BLR_L(IND_L), ONE,
     &              A, LA, POSELT_BLOCK,
     &              NFRONT, 0, IFLAG, IERROR,
     &              MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &              KPERCENT_RMB, MID_RANK, BUILDQ,
     &              (K480.GE.3), LorU=2,
     &              LRB3=ACC_LRB, MAXI_RANK=MAXI_RANK,
     &              MAXI_CLUSTER=MAXI_CLUSTER)
          IF (IFLAG.LT.0) GOTO 100
          CALL UPD_FLOP_UPDATE(BLR_U(IND_U), BLR_L(IND_L),
     &           MIDBLK_COMPRESS, MID_RANK, BUILDQ, 
     &           .FALSE., (K480.GE.3))
          IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
            K_RANK(KK) = MID_RANK
          ENDIF
          IF (K480.GE.3) THEN
            NEW_ACC_RANK = NEW_ACC_RANK + ACC_LRB%K - OLD_ACC_RANK
            MAX_ACC_RANK = MAX(MAX_ACC_RANK, ACC_LRB%K - OLD_ACC_RANK)
            IF (K480.EQ.4) THEN
              IF ((K478.GT.0).AND.((ACC_LRB%K-MAX_ACC_RANK).GE.K478)) 
     &          THEN
                CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &                 NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &                 TOL_OPT,
     &                 KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                MAX_ACC_RANK = ACC_LRB%K
              ENDIF
            ENDIF
          ENDIF
        END DO
        IF (K480.GE.3) THEN
          IF (K480.GE.5.OR.COMPRESS_CB) THEN
            IF (K480.GE.5.AND.(COMPRESSED_FR.OR.K480.GE.6)) THEN  
              IF (ACC_LRB%K.GT.0) THEN
                IF (K478.EQ.-1) THEN
                  IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                    CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
                  ENDIF
                ELSEIF (K478.LE.-2) THEN
                  IF (FRFR_UPDATES.GT.0) THEN
                    allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM-NB_DEC
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    DO II = 1,NB_INASM-NB_DEC-1
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
                  ELSE
                    allocate(POS_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       GOTO 100
                    ENDIF
                    POS_LIST(1) = 1
                    POS_LIST(2) = 1 + FR_RANK
                    DO II = 2,NB_INASM
                      POS_LIST(II+1)=POS_LIST(II)+K_RANK(II-1) 
                    ENDDO
                    allocate(RANK_LIST(NB_INASM+1),stat=allocok)
                    IF (allocok .GT. 0) THEN
                       IFLAG  = -13
                       IERROR = NB_INASM+1
                       GOTO 100
                    ENDIF
                    RANK_LIST(1) = FR_RANK
                    DO II = 2,NB_INASM+1
                      RANK_LIST(II) = K_RANK(II-1)
                    ENDDO
                    CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST, POS_LIST,
     &               NB_INASM+1, 0)
                    deallocate(RANK_LIST)
                  ENDIF
                  deallocate(POS_LIST)
                ENDIF
              ENDIF
            ENDIF
            MAXRANK = FLOOR(dble(ACC_LRB%M*ACC_LRB%N)/dble(ACC_LRB%M+
     &                                                     ACC_LRB%N))
            IF (COMPRESSED_FR.AND.(ACC_LRB%K.LE.MAXRANK)) THEN
              LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
              CALL ALLOC_LRB_FROM_ACC(ACC_LRB, LRB,
     &                       ACC_LRB%K, ACC_LRB%M, ACC_LRB%N, 0,
     &                       IFLAG, IERROR, KEEP8)
              FRONT_CB_BLR_SAVINGS = FRONT_CB_BLR_SAVINGS + 
     &                LRB%M*LRB%N - LRB%M*LRB%K - LRB%N*LRB%K
              ACC_LRB%K = 0
              IF (IFLAG.LT.0) GOTO 100
            ELSE
              CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &              MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
              LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
              CALL ALLOC_LRB(LRB, ACC_LRB%K, ACC_LRB%N, ACC_LRB%M, 
     &                 .FALSE., IFLAG, IERROR, KEEP8)
              IF (IFLAG.LT.0) GOTO 100
              DO II=1,ACC_LRB%N
                LRB%Q(II,1:ACC_LRB%M) =
     &          A( POSELT_BLOCK+int((II-1),8)*int(NFRONT,8) :
     &            POSELT_BLOCK+int((II-1),8)*int(NFRONT,8)
     &                        +int(ACC_LRB%M-1,8) )
              END DO  
            ENDIF
          ELSE 
            IF ((K480.EQ.4).AND.(K478.EQ.-1).AND.(ACC_LRB%K.GT.0)) THEN
              IF (NB_INASM-FRFR_UPDATES.GT.1) THEN               
                CALL DMUMPS_RECOMPRESS_ACC(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
              ENDIF
            ELSEIF ((K480.EQ.4).AND.(K478.LE.-2).AND.(ACC_LRB%K.GT.0)) 
     &        THEN
              allocate(POS_LIST(NB_INASM-NB_DEC),stat=allocok)
              IF (allocok .GT. 0) THEN
                 IFLAG  = -13
                 IERROR = NB_INASM-NB_DEC
                 GOTO 100
              ENDIF
              POS_LIST(1) = 1
              DO II = 1,NB_INASM-NB_DEC-1
                POS_LIST(II+1)=POS_LIST(II)+K_RANK(NB_DEC+II) 
              ENDDO
              CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELT_BLOCK,
     &               KEEP8,NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, 
     &               TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               K_RANK(NB_DEC+1:NB_INASM), POS_LIST,
     &               NB_INASM-NB_DEC, 0)
              deallocate(POS_LIST)
            ENDIF
            CALL DMUMPS_DECOMPRESS_ACC(ACC_LRB,MAXI_CLUSTER,
     &            MAXI_RANK, A, LA, POSELT_BLOCK, NFRONT, NIV, 2)
          ENDIF
        ENDIF                               
 100    CONTINUE       
      END DO
#if defined(BLR_MT)
!$OMP END DO 
#endif
      IF (COMPRESS_CB) THEN
#if defined(BLR_MT)
!$      OMP_NUM = OMP_GET_THREAD_NUM()
!$      IF (OMP_NUM.EQ.0) THEN
#endif
        CALL UPD_MRY_CB(NFRONT-NASS, NFRONT-NASS, 0, 1,
     &                        FRONT_CB_BLR_SAVINGS)
#if defined(BLR_MT)
!$      ELSE
!$        CALL UPD_MRY_CB(0, 0, 0, 1,
!$   &                        FRONT_CB_BLR_SAVINGS)
!$      ENDIF
#endif
      ENDIF
      END SUBROUTINE DMUMPS_BLR_UPD_CB_LEFT
      SUBROUTINE DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, LDA11,
     &        LDA21, COPY_DENSE_BLOCKS,
     &        BEGS_BLR_DIAG, BEGS_BLR_FIRST_OFFDIAG,
     &        NB_BLR, BLR_PANEL, CURRENT_BLR, DIR, DECOMP_TIMER,
     &        BEG_I_IN, END_I_IN, ONLY_NELIM_IN, CBASM_TOFIX_IN)
!$    USE OMP_LIB 
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      LOGICAL, intent(in)          :: COPY_DENSE_BLOCKS  
      INTEGER, intent(in)          :: NB_BLR, CURRENT_BLR
      INTEGER, intent(in)          :: BEGS_BLR_DIAG, 
     &                                BEGS_BLR_FIRST_OFFDIAG
      TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(:)
      CHARACTER(len=1) :: DIR
      INTEGER, intent(in) :: LDA11, LDA21
      INTEGER, intent(in) :: DECOMP_TIMER
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN, ONLY_NELIM_IN
      LOGICAL,OPTIONAL,intent(in) :: CBASM_TOFIX_IN
      INTEGER :: IP, M, N, BIP, BEG_I, END_I, ONLY_NELIM
      LOGICAL :: CBASM_TOFIX
#if defined(BLR_MT)
      INTEGER :: LAST_IP, CHUNK
#endif
      INTEGER :: K, I
      DOUBLE PRECISION :: PROMOTE_COST
      INTEGER(8) :: POSELT_BLOCK, LD_BLK_IN_FRONT
      DOUBLE PRECISION :: ONE, ALPHA, ZERO
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO = 0.0D0)
      IF(present(BEG_I_IN)) THEN
        BEG_I = BEG_I_IN
      ELSE
        BEG_I = CURRENT_BLR+1
      ENDIF
      IF(present(END_I_IN)) THEN
        END_I = END_I_IN
      ELSE
        END_I = NB_BLR
      ENDIF
      IF(present(ONLY_NELIM_IN)) THEN
        ONLY_NELIM = ONLY_NELIM_IN
      ELSE
        ONLY_NELIM = 0
      ENDIF
      IF (present(CBASM_TOFIX_IN)) THEN
        CBASM_TOFIX = CBASM_TOFIX_IN
      ELSE
        CBASM_TOFIX = .FALSE.
      ENDIF
      LD_BLK_IN_FRONT = int(LDA11,8)
      BIP             = BEGS_BLR_FIRST_OFFDIAG
#if !defined(BLR_MT)
      IF (BEG_I .NE. CURRENT_BLR+1) THEN
        DO I = 1, BEG_I - CURRENT_BLR - 1
          IF (CBASM_TOFIX) THEN
            BIP  = BIP +  BLR_PANEL(I)%N
          ELSE
            BIP  = BIP +  BLR_PANEL(I)%M
          ENDIF
        ENDDO
      ENDIF
#endif
#if defined(BLR_MT)
      LAST_IP = CURRENT_BLR+1
      CHUNK = 1
!$OMP DO PRIVATE(POSELT_BLOCK, M, N, K, I) SCHEDULE(DYNAMIC, CHUNK)
#endif
      DO IP = BEG_I, END_I
#if defined(BLR_MT)
        DO I = 1, IP - LAST_IP
          IF (CBASM_TOFIX) THEN
            BIP  = BIP +  BLR_PANEL(LAST_IP-CURRENT_BLR+I-1)%N
          ELSE
            BIP  = BIP +  BLR_PANEL(LAST_IP-CURRENT_BLR+I-1)%M
          ENDIF
        ENDDO
        LAST_IP = IP
#endif
        IF (DIR .eq. 'V') THEN
           IF (BIP .LE. LDA21) THEN
             IF (CBASM_TOFIX) THEN
               POSELT_BLOCK = POSELT 
     &              + int(LDA11,8)*int(BEGS_BLR_DIAG-1,8) + int(BIP-1,8)
             ELSE
               POSELT_BLOCK = POSELT + int(LDA11,8)*int(BIP-1,8) + 
     &              int(BEGS_BLR_DIAG - 1,8)
             ENDIF
           ELSE
             POSELT_BLOCK = POSELT + int(LDA11,8)*int(LDA21,8)+
     &              int(BEGS_BLR_DIAG - 1,8)
             POSELT_BLOCK = POSELT_BLOCK +
     &                      int(LDA21,8)*int(BIP-1-LDA21,8)
             LD_BLK_IN_FRONT=int(LDA21,8)
           ENDIF
        ELSE 
         POSELT_BLOCK = POSELT + int(LDA11,8)*int(BEGS_BLR_DIAG-1,8)
     &              + int(BIP-1,8)
        ENDIF
        M = BLR_PANEL(IP-CURRENT_BLR)%M
        N = BLR_PANEL(IP-CURRENT_BLR)%N
        IF(present(ONLY_NELIM_IN)) THEN
          ONLY_NELIM = ONLY_NELIM_IN
        ELSE
          ONLY_NELIM = N
        ENDIF
        K = BLR_PANEL(IP-CURRENT_BLR)%K
        IF (BLR_PANEL(IP-CURRENT_BLR)%ISLR) THEN
          IF (K.EQ.0) THEN
            IF (DIR .eq. 'V') THEN
              DO I = 1, M
                IF (BIP+I-1.GT.LDA21) THEN
                  LD_BLK_IN_FRONT = int(LDA21,8)
                ENDIF
                A(POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT :
     &                   POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT
     &                          + int(N-1,8)) = ZERO
              ENDDO
            ELSE 
              DO I = N-ONLY_NELIM+1, N
              A(POSELT_BLOCK+int(I-1,8)*int(LDA11,8):
     &          POSELT_BLOCK+int(I-1,8)*int(LDA11,8) + int(M-1,8)) 
     &              = ZERO
              ENDDO
            ENDIF
            GOTO 1800
          ENDIF
          IF (DIR .eq. 'V') THEN
            IF (DIR .eq.'V' .AND. BIP .LE. LDA21
     &                 .AND. BIP + M - 1 .GT. LDA21
     &                 .AND..NOT.CBASM_TOFIX) THEN
              CALL dgemm('T', 'T', N, LDA21-BIP+1, K, ONE ,
     &            BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &            BLR_PANEL(IP-CURRENT_BLR)%Q(1,1) , M, 
     &            ZERO, A(POSELT_BLOCK), int(LD_BLK_IN_FRONT))
              CALL dgemm('T', 'T', N, BIP+M-LDA21-1, K, ONE ,
     &            BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &            BLR_PANEL(IP-CURRENT_BLR)%Q(LDA21-BIP+2,1) , M, 
     &            ZERO, A(POSELT_BLOCK+int(LDA21-BIP,8)*int(LDA11,8)),
     &            LDA21)
            ELSE
              CALL dgemm('T', 'T', N, M, K, ONE ,
     &            BLR_PANEL(IP-CURRENT_BLR)%R(1,1) , K, 
     &            BLR_PANEL(IP-CURRENT_BLR)%Q(1,1) , M, 
     &            ZERO, A(POSELT_BLOCK), int(LD_BLK_IN_FRONT))
            ENDIF
          ELSE 
             CALL dgemm('N', 'N', M, ONLY_NELIM, K, ONE,
     &          BLR_PANEL(IP-CURRENT_BLR)%Q(1,1), M, 
     &          BLR_PANEL(IP-CURRENT_BLR)%R(1,N-ONLY_NELIM+1), K, ZERO,
     &          A(POSELT_BLOCK+int(N-ONLY_NELIM,8)*int(LDA11,8)), LDA11)
          ENDIF
          PROMOTE_COST = 2.0D0*M*K*ONLY_NELIM
          IF (CBASM_TOFIX) THEN
            CALL UPD_FLOP_DECOMPRESS(PROMOTE_COST, .TRUE.)
          ELSEIF(present(ONLY_NELIM_IN)) THEN
            CALL UPD_FLOP_DECOMPRESS(PROMOTE_COST, .FALSE.)
          ENDIF
        ELSE  IF (COPY_DENSE_BLOCKS) THEN
          IF (DIR .eq. 'V') THEN
            DO I = 1, M
              IF (BIP+I-1.GT.LDA21) THEN
                LD_BLK_IN_FRONT = int(LDA21,8)
              ENDIF
                A(POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT :
     &            POSELT_BLOCK+int(I-1,8)*LD_BLK_IN_FRONT
     &                                           + int(N-1,8)) 
     &          = BLR_PANEL(IP-CURRENT_BLR)%Q(I,1:N)
            ENDDO
          ELSE 
            DO I = N-ONLY_NELIM+1, N
              A(POSELT_BLOCK+int(I-1,8)*int(LDA11,8):
     &           POSELT_BLOCK+int(I-1,8)*int(LDA11,8) + int(M-1,8))
     &        = BLR_PANEL(IP-CURRENT_BLR)%Q(1:M,I)
            ENDDO
          ENDIF
        ENDIF
 1800   CONTINUE
#if !defined(BLR_MT)
        IF (CBASM_TOFIX) THEN
          BIP  = BIP +  BLR_PANEL(IP-CURRENT_BLR)%N
        ELSE
          BIP  = BIP +  BLR_PANEL(IP-CURRENT_BLR)%M
        ENDIF
#endif
      ENDDO
#if defined(BLR_MT)
!$OMP END DO
#endif
      END SUBROUTINE DMUMPS_DECOMPRESS_PANEL
      SUBROUTINE DMUMPS_COMPRESS_CB(A, LA, POSELT, LDA,
     &        BEGS_BLR, BEGS_BLR_U, NB_ROWS, NB_COLS, NB_INASM,
     &        NROWS, NCOLS, INODE,   
     &        IWHANDLER, SYM, NIV, IFLAG, IERROR,
     &        TOLEPS, TOL_OPT, KPERCENT, K489, CB_LRB,
     &        WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, KEEP8, 
     &        NFS4FATHER, NPIV, NVSCHUR_K253, KEEP,
     &        M_ARRAY,
     &        NELIM, 
     &        NBROWSinF
     &        )
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: LDA, NB_ROWS, NB_COLS, NB_INASM
      INTEGER, INTENT(IN)          :: NIV, IWHANDLER, MAXI_CLUSTER, 
     &                                KPERCENT, TOL_OPT, LWORK
      INTEGER, INTENT(IN)          :: K489, NROWS, NCOLS, INODE, SYM
      INTEGER, intent(inout)         :: IFLAG, IERROR
      TYPE(LRB_TYPE), TARGET, intent(inout) :: CB_LRB(:,:)
      INTEGER, DIMENSION(:) :: BEGS_BLR, BEGS_BLR_U
      DOUBLE PRECISION, TARGET, DIMENSION(:) :: RWORK
      DOUBLE PRECISION, TARGET, DIMENSION(:,:) :: BLOCK
      DOUBLE PRECISION, TARGET, DIMENSION(:) :: WORK, TAU
      INTEGER, TARGET, DIMENSION(:) :: JPVT
      INTEGER(8) :: KEEP8(150)
      DOUBLE PRECISION,intent(in) :: TOLEPS
      INTEGER, INTENT(in) :: NFS4FATHER, NPIV, NVSCHUR_K253, KEEP(500)
      DOUBLE PRECISION, OPTIONAL :: M_ARRAY(max(NFS4FATHER,1))
      INTEGER, intent(in), OPTIONAL :: NELIM
      INTEGER, intent(in), OPTIONAL :: NBROWSinF
      INTEGER :: M, N, INFO, FRONT_CB_BLR_SAVINGS
      INTEGER :: I, J, IBIS, IBIS_END, RANK, MAXRANK, II, JJ
      INTEGER(8) :: POSELT_BLOCK
      LOGICAL :: ISLR
      TYPE(LRB_TYPE), POINTER :: LRB
      INTEGER :: OMP_NUM
      INTEGER(8) :: POSA, ASIZE
      INTEGER    :: NROWS_CM
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: RWORK_THR
      DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: BLOCK_THR
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: WORK_THR, TAU_THR
      INTEGER, POINTER, DIMENSION(:) :: JPVT_THR
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
#if defined(BLR_MT)
!$OMP MASTER
#endif
      IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &     (NFS4FATHER.GT.0) ) THEN
       IF (NIV.EQ.1) THEN
        NROWS_CM  = NROWS  - (NFS4FATHER-NELIM)
       ELSE
        NROWS_CM  = NROWS  - NBROWSinF
       ENDIF
       IF (NROWS_CM-NVSCHUR_K253.GT.0)  THEN
         IF (NIV.EQ.1) THEN
          POSA     = POSELT
     &             + int(LDA,8)*int(NPIV+NFS4FATHER,8)
     &             + int(NPIV,8)
          ASIZE    = int(LDA,8)*int(LDA,8)
     &             - int(LDA,8)*int(NPIV+NFS4FATHER,8)
     &             - int(NPIV,8) 
         ELSE
          POSA     = POSELT 
     &             + int(LDA,8)*int(NBROWSinF,8)
     &             + int(NPIV,8) 
          ASIZE    = int(NROWS,8)*int(LDA,8)
     &             - int(LDA,8)*int(NBROWSinF,8)
     &             - int(NPIV,8)
         ENDIF
         CALL DMUMPS_COMPUTE_MAXPERCOL (
     &      A(POSA), ASIZE, LDA, 
     &      NROWS_CM-NVSCHUR_K253,
     &      M_ARRAY(1), NFS4FATHER, .FALSE., 
     &      -9999)
       ELSE 
          DO I=1, NFS4FATHER
            M_ARRAY(I) = ZERO
          ENDDO
       ENDIF
      ENDIF
#if defined(BLR_MT)
!$OMP END MASTER
!$OMP BARRIER
#endif
      FRONT_CB_BLR_SAVINGS = 0
      OMP_NUM = 0 
      IF (SYM.EQ.0.OR.NIV.EQ.2) THEN
        IBIS_END = NB_ROWS*NB_COLS
      ELSE
        IBIS_END = NB_ROWS*(NB_COLS+1)/2
      ENDIF
#if defined(BLR_MT)
      CHUNK = 1
!$OMP DO SCHEDULE(DYNAMIC,CHUNK) 
!$OMP& PRIVATE(I, J, POSELT_BLOCK, M, N, OMP_NUM, INFO, RANK,
!$OMP&         MAXRANK, ISLR, II, JJ, LRB)
#endif
      DO IBIS = 1,IBIS_END
        IF (IFLAG.LT.0) CYCLE     
#if defined(BLR_MT)         
        OMP_NUM = 0 
!$      OMP_NUM = OMP_GET_THREAD_NUM()
#endif
        BLOCK_THR => BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*MAXI_CLUSTER) 
        JPVT_THR  => JPVT(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        TAU_THR   => TAU(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        WORK_THR  => WORK(OMP_NUM*LWORK+1:
     &                 (OMP_NUM+1)*LWORK) 
        RWORK_THR => RWORK(OMP_NUM*2*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*2*MAXI_CLUSTER) 
        IF (SYM.EQ.0.OR.NIV.EQ.2) THEN
          I = (IBIS-1)/NB_COLS+1
          J = IBIS - (I-1)*NB_COLS
        ELSE
          I = CEILING((1.0D0+SQRT(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
          J = IBIS - I*(I-1)/2
        ENDIF
        IF (NIV.EQ.1) THEN
          I = I+NB_INASM
          J = J+NB_INASM
        ELSE
          J = J+NB_INASM
          IF (SYM.NE.0) THEN
            IF (BEGS_BLR_U(J).GE.BEGS_BLR(I+2)+NCOLS-NROWS-1+
     &          BEGS_BLR_U(NB_INASM+1)) THEN 
              CYCLE
            ENDIF
          ENDIF
        ENDIF
        IF (NIV.EQ.1) THEN
          M = BEGS_BLR(I+1)-BEGS_BLR(I)
          POSELT_BLOCK = POSELT + int(LDA,8)*int(BEGS_BLR(I)-1,8) + 
     &           int(BEGS_BLR_U(J)-1,8)
          IF (I .EQ. NB_INASM+1 .AND. present(NELIM)) THEN
            POSELT_BLOCK = POSELT_BLOCK + int(NELIM,8)*int(LDA,8)
            M = M - NELIM
          ENDIF
          N = BEGS_BLR_U(J+1)-BEGS_BLR_U(J)
        ELSE
          M = BEGS_BLR(I+2)-BEGS_BLR(I+1)
          POSELT_BLOCK = POSELT + int(LDA,8)*int(BEGS_BLR(I+1)-1,8)
     &           + int(BEGS_BLR_U(J)-1,8)
          IF (SYM.EQ.0) THEN
            N = BEGS_BLR_U(J+1)-BEGS_BLR_U(J)
          ELSE
            N = min(BEGS_BLR_U(J+1), BEGS_BLR(I+2) + NCOLS - NROWS -1
     &              + BEGS_BLR_U(NB_INASM+1)) - BEGS_BLR_U(J) 
          ENDIF
        ENDIF
        JPVT_THR(1:MAXI_CLUSTER) = 0
        IF (NIV.EQ.1) THEN
          LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
        ELSE
          LRB => CB_LRB(I,J-NB_INASM)
        ENDIF
        IF (K489.EQ.3) THEN
            MAXRANK = 1
            RANK = MAXRANK+1
            INFO = 0
            GOTO 3800
        ENDIF
        DO II=1,M
          BLOCK_THR(II,1:N)=
     &    A( POSELT_BLOCK+int(II-1,8)*int(LDA,8) :
     &    POSELT_BLOCK+int(II-1,8)*int(LDA,8)+int(N-1,8) )
        ENDDO  
        MAXRANK = floor(dble(M*N)/dble(M+N))
        MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
        CALL DMUMPS_TRUNCATED_RRQR( M, N,
     &       BLOCK_THR(1,1),
     &       MAXI_CLUSTER, JPVT_THR(1), 
     &       TAU_THR(1), 
     &       WORK_THR(1), N, 
     &       RWORK_THR(1), 
     &       TOLEPS, TOL_OPT, RANK, MAXRANK, INFO)
 3800 CONTINUE
        IF (INFO < 0) THEN
           WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                 " OF TRUNCATED_RRQR WHILE COMPRESSING A CB BLOCK"
           CALL MUMPS_ABORT()
        END IF
        ISLR = ((RANK.LE.MAXRANK).AND.(M.NE.0).AND.(N.NE.0))     
        CALL ALLOC_LRB(LRB, RANK, M, N, ISLR, IFLAG, IERROR, KEEP8)
        IF (IFLAG.LT.0) CYCLE
        IF (ISLR) THEN 
           IF (RANK .GT. 0) THEN 
               DO JJ=1,N
                  DO II=1,MIN(RANK,JJ)
                     LRB%R(II,JPVT_THR(JJ)) = BLOCK_THR(II,JJ)
                  ENDDO 
                  IF(JJ.LT.RANK) LRB%R(MIN(RANK,JJ)+1:RANK,JPVT_THR(JJ))
     &                 = ZERO
               ENDDO
               CALL dorgqr 
     &           (M, RANK, RANK,
     &           BLOCK_THR(1,1), 
     &           MAXI_CLUSTER, TAU_THR(1), 
     &           WORK_THR(1), LWORK, INFO )
               DO II=1,RANK
                 DO JJ= 1, M
                  LRB%Q(JJ,II) = BLOCK_THR(JJ,II)
                 ENDDO
               END DO
               IF (INFO < 0) THEN
                 WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                     " OF CUNGQR WHILE COMPRESSING A CB BLOCK"
                 CALL MUMPS_ABORT()
               END IF
               IF (K489.NE.3) THEN
                 CALL UPD_FLOP_COMPRESS(LRB, CB_COMPRESS=.TRUE.)
               ENDIF
          END IF
          FRONT_CB_BLR_SAVINGS = FRONT_CB_BLR_SAVINGS + 
     &                (M-RANK)*(N-RANK)-RANK*RANK
        ELSE 
           DO II=1,M
             LRB%Q(II,1:N) =
     &       A( POSELT_BLOCK+int((II-1),8)*int(LDA,8) :
     &         POSELT_BLOCK+int((II-1),8)*int(LDA,8)
     &                     +int(N-1,8) )
           END DO  
           IF (K489.NE.3) THEN
               CALL UPD_FLOP_COMPRESS(LRB, CB_COMPRESS=.TRUE.)
           ENDIF
           LRB%K = -1
        END IF
      END DO
#if defined(BLR_MT)
!$OMP END DO 
#endif
#if defined(BLR_MT)
!$      OMP_NUM = OMP_GET_THREAD_NUM()
!$      IF (OMP_NUM.EQ.0) THEN
#endif
        CALL UPD_MRY_CB(NROWS, NCOLS, SYM, NIV,
     &                        FRONT_CB_BLR_SAVINGS)
#if defined(BLR_MT)
!$      ELSE
!$        CALL UPD_MRY_CB(0, 0, SYM, NIV,
!$   &                        FRONT_CB_BLR_SAVINGS)
!$      ENDIF
#endif
      END SUBROUTINE DMUMPS_COMPRESS_CB
      SUBROUTINE DMUMPS_COMPRESS_PANEL(
     &        A, LA, POSELT, IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, TOLEPS, TOL_OPT, K473, BLR_PANEL, 
     &        CURRENT_BLR,
     &        DIR, WORK, TAU, JPVT, 
     &        LWORK, RWORK, BLOCK,
     &        MAXI_CLUSTER, NELIM, 
     &        LBANDSLAVE, NPIV, ISHIFT, NIV, KPERCENT, 
     &        KEEP8, K480,
     &        BEG_I_IN, END_I_IN, FRSWAP
     &        )
!$    USE OMP_LIB
      INTEGER(8), intent(in)       :: LA
      DOUBLE PRECISION, intent(inout)       :: A(LA)
      INTEGER(8), intent(in)       :: POSELT 
      INTEGER, intent(in)          :: NFRONT, NB_BLR, CURRENT_BLR, NIV
      INTEGER, intent(inout)          :: IFLAG, IERROR
      TYPE(LRB_TYPE), intent(inout) :: BLR_PANEL(:)
      DOUBLE PRECISION, TARGET, DIMENSION(:) :: RWORK
      DOUBLE PRECISION, TARGET, DIMENSION(:,:) :: BLOCK
      DOUBLE PRECISION, TARGET, DIMENSION(:) :: WORK, TAU
      INTEGER, TARGET, DIMENSION(:) :: JPVT
      INTEGER :: BEGS_BLR(:)
      INTEGER(8) :: KEEP8(150)
      INTEGER, OPTIONAL, intent(in) :: K480
      INTEGER,OPTIONAL,intent(in) :: BEG_I_IN, END_I_IN
      LOGICAL, OPTIONAL, intent(in) :: FRSWAP
      INTEGER, intent(in)          :: NPIV, ISHIFT, KPERCENT, K473,
     &                                TOL_OPT
      LOGICAL, intent(in)          :: LBANDSLAVE
      INTEGER                      :: MAXI_CLUSTER, LWORK, NELIM
      DOUBLE PRECISION,intent(in)              :: TOLEPS
      CHARACTER(len=1) :: DIR
      INTRINSIC maxval
      INTEGER :: IP, NB_BLOCKS_PANEL, M, N, RANK, MAXRANK
      INTEGER :: INFO, I, J, IS, BEG_I, END_I
      INTEGER(8) :: POSELT_BLOCK
      LOGICAL :: ISLR
      DOUBLE PRECISION :: ONE, ALPHA, ZERO
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO = 0.0D0)
      INTEGER :: OMP_NUM
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: RWORK_THR
      DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: BLOCK_THR
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: WORK_THR, TAU_THR
      INTEGER, POINTER, DIMENSION(:) :: JPVT_THR
#if defined(BLR_MT) 
      INTEGER :: CHUNK
#endif
      IF(present(BEG_I_IN)) THEN
        BEG_I = BEG_I_IN
      ELSE
        BEG_I = CURRENT_BLR+1
      ENDIF
      IF(present(END_I_IN)) THEN
        END_I = END_I_IN
      ELSE
        END_I = NB_BLR
      ENDIF
      IF (LBANDSLAVE) THEN
       IS = ISHIFT
      ELSE
       IS=0
      ENDIF
      IF (DIR .eq. 'V') THEN
         IF (LBANDSLAVE) THEN
          N = NPIV
         ELSE
          N = BEGS_BLR(CURRENT_BLR+1)-BEGS_BLR(CURRENT_BLR)-NELIM
         ENDIF
      ELSE IF (DIR .eq. 'H') THEN
        N = BEGS_BLR(CURRENT_BLR+1)-BEGS_BLR(CURRENT_BLR)-NELIM
      ELSE
         WRITE(*,*) " WRONG ARGUMENT IN DMUMPS_COMPRESS_PANEL "
         CALL MUMPS_ABORT()
      END IF
      NB_BLOCKS_PANEL = NB_BLR-CURRENT_BLR
      OMP_NUM = 0 
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO PRIVATE(INFO, POSELT_BLOCK, RANK, MAXRANK, I, J, OMP_NUM)
!$OMP&   SCHEDULE(DYNAMIC,CHUNK)
#endif
      DO IP = BEG_I, END_I
        IF (IFLAG.LT.0) CYCLE
#if defined(BLR_MT)         
        OMP_NUM = 0 
!$      OMP_NUM = OMP_GET_THREAD_NUM()
#endif
        BLOCK_THR => BLOCK(1:MAXI_CLUSTER,OMP_NUM*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*MAXI_CLUSTER) 
        JPVT_THR  => JPVT(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        TAU_THR   => TAU(OMP_NUM*MAXI_CLUSTER+1:
     &                 (OMP_NUM+1)*MAXI_CLUSTER) 
        WORK_THR  => WORK(OMP_NUM*LWORK+1:
     &                 (OMP_NUM+1)*LWORK) 
        RWORK_THR => RWORK(OMP_NUM*2*MAXI_CLUSTER+1:
     &                  (OMP_NUM+1)*2*MAXI_CLUSTER) 
        RANK = 0
        M = BEGS_BLR(IP+1)-BEGS_BLR(IP)
        IF (DIR .eq. 'V') THEN
          POSELT_BLOCK = POSELT + 
     &              int(NFRONT,8) * int(BEGS_BLR(IP)-1,8) + 
     &              int(BEGS_BLR(CURRENT_BLR) + IS - 1,8)
        ELSE 
          POSELT_BLOCK = POSELT +
     &              int(NFRONT,8)*int(BEGS_BLR(CURRENT_BLR)-1,8) + 
     &              int( BEGS_BLR(IP) - 1,8)
        ENDIF
        IF (present(K480)) then 
        IF (K480.GE.5) THEN
          IF (BLR_PANEL(IP-CURRENT_BLR)%ISLR) THEN
             IF (M.NE.BLR_PANEL(IP-CURRENT_BLR)%M) THEN
              write(*,*) 'Internal error in DMUMPS_COMPRESS_PANEL',
     &                    ' M size inconsistency',M,
     &                    BLR_PANEL(IP-CURRENT_BLR)%M
              CALL MUMPS_ABORT()
            ENDIF
            IF (N.NE.BLR_PANEL(IP-CURRENT_BLR)%N) THEN
              write(*,*) 'Internal error in DMUMPS_COMPRESS_PANEL',
     &                    ' N size inconsistency',N,
     &                    BLR_PANEL(IP-CURRENT_BLR)%N
              CALL MUMPS_ABORT()
            ENDIF
            MAXRANK = floor(dble(M*N)/dble(M+N))
            IF (BLR_PANEL(IP-CURRENT_BLR)%K.GT.MAXRANK) THEN
              write(*,*) 'Internal error in DMUMPS_COMPRESS_PANEL',
     &                    ' MAXRANK inconsistency',MAXRANK,
     &                    BLR_PANEL(IP-CURRENT_BLR)%K
              CALL MUMPS_ABORT()
            ENDIF
            GOTO 3000
          ENDIF
        ENDIF
        ENDIF
        JPVT_THR(1:MAXI_CLUSTER) = 0
        IF (K473.EQ.1) THEN
            MAXRANK = 1
            RANK = MAXRANK+1
            INFO = 0
            GOTO 3800
        ENDIF
        IF (DIR .eq. 'V') THEN
            DO I=1,M
                BLOCK_THR(I,1:N)=
     &          A( POSELT_BLOCK+int(I-1,8)*int(NFRONT,8) :
     &          POSELT_BLOCK+int(I-1,8)*int(NFRONT,8)+int(N-1,8) )
            END DO  
        ELSE 
            DO I=1,N
                BLOCK_THR(1:M,I)=
     &          A( POSELT_BLOCK+int(I-1,8)*int(NFRONT,8) :
     &          POSELT_BLOCK+int(I-1,8)*int(NFRONT,8)+int(M-1,8) )
            END DO  
        END IF
        MAXRANK = floor(dble(M*N)/dble(M+N))
        MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
        CALL DMUMPS_TRUNCATED_RRQR( M, N,
     &       BLOCK_THR(1,1),
     &       MAXI_CLUSTER, JPVT_THR(1), 
     &       TAU_THR(1), 
     &       WORK_THR(1), N, 
     &       RWORK_THR(1), 
     &       TOLEPS, TOL_OPT, RANK, MAXRANK, INFO)
 3800 CONTINUE
      IF (INFO < 0) THEN
           WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                 " OF TRUNCATED_RRQR WHILE COMPRESSING A BLOCK "
           CALL MUMPS_ABORT()
        END IF
        ISLR = ((RANK.LE.MAXRANK).AND.(M.NE.0).AND.(N.NE.0))     
        CALL ALLOC_LRB(BLR_PANEL(IP-CURRENT_BLR), RANK,
     &                 M, N, ISLR, IFLAG, IERROR, KEEP8)
        IF (IFLAG.LT.0) CYCLE
        IF ((M.EQ.0).OR.(N.EQ.0)) GOTO 3000
        IF (ISLR) THEN 
           IF (RANK .EQ. 0) THEN 
           ELSE 
               DO J=1,N
                 BLR_PANEL(IP-CURRENT_BLR)%R(1:MIN(RANK,J),
     &               JPVT_THR(J)) =
     &               BLOCK_THR(1:MIN(RANK,J),J)
                 IF(J.LT.RANK) BLR_PANEL(IP-CURRENT_BLR)%
     &               R(MIN(RANK,J)+1:RANK,JPVT_THR(J))= ZERO
               ENDDO
               CALL dorgqr 
     &           (M, RANK, RANK,
     &           BLOCK_THR(1,1), 
     &           MAXI_CLUSTER, TAU_THR(1), 
     &           WORK_THR(1), LWORK, INFO )
               DO I=1,RANK
                 BLR_PANEL(IP-CURRENT_BLR)%Q(1:M,I) = BLOCK_THR(1:M,I)
               END DO
               IF (INFO < 0) THEN
                 WRITE(*,*) " PROBLEM IN ARGUMENT NUMBER ",INFO,
     &                     " OF CUNGQR WHILE COMPRESSING A BLOCK "
                 CALL MUMPS_ABORT()
               END IF
               IF (present(FRSWAP)) THEN
                 CALL UPD_FLOP_COMPRESS(
     &               BLR_PANEL(IP-CURRENT_BLR), FRSWAP=FRSWAP)
               ELSE
                 CALL UPD_FLOP_COMPRESS(BLR_PANEL(IP-CURRENT_BLR))
               ENDIF
          END IF
        ELSE 
           IF (DIR .eq. 'V') THEN
               DO I=1,M
                   BLR_PANEL(IP-CURRENT_BLR)%Q(I,1:N) =
     &             A( POSELT_BLOCK+int((I-1),8)*int(NFRONT,8) :
     &               POSELT_BLOCK+int((I-1),8)*int(NFRONT,8)
     &                           +int(N-1,8) )
               END DO  
           ELSE 
               DO I=1,N
                   BLR_PANEL(IP-CURRENT_BLR)%Q(1:M,I) =
     &             A( POSELT_BLOCK+int((I-1),8)*int(NFRONT,8) :
     &               POSELT_BLOCK+int((I-1),8)*int(NFRONT,8)
     &                           +int(M-1,8) )
               END DO  
           END IF
           IF (K473.EQ.0) THEN
             IF (present(FRSWAP)) THEN
               CALL UPD_FLOP_COMPRESS(BLR_PANEL(IP-CURRENT_BLR), 
     &                                   FRSWAP=FRSWAP)
             ELSE
               CALL UPD_FLOP_COMPRESS(BLR_PANEL(IP-CURRENT_BLR))
             ENDIF
           ENDIF
           BLR_PANEL(IP-CURRENT_BLR)%K = -1
        END IF
 3000   CONTINUE
      END DO 
#if defined(BLR_MT) 
!$OMP END DO NOWAIT
#endif
      END SUBROUTINE DMUMPS_COMPRESS_PANEL
      SUBROUTINE DMUMPS_BLR_PANEL_LRTRSM(
     &                A,
     &                LA, POSELT, NFRONT,
     &                IBEG_BLOCK, NB_BLR,
     &                BLR_LorU,
     &                CURRENT_BLR, FIRST_BLOCK, LAST_BLOCK,
     &                NIV, SYM, LorU, LBANDSLAVE,
     &                IW, OFFSET_IW, NASS) 
!$    USE OMP_LIB
      INTEGER(8), intent(in)  :: LA
      INTEGER, intent(in)     :: NFRONT, NB_BLR, CURRENT_BLR,
     &                           NIV, SYM, LorU
      LOGICAL, intent(in)     :: LBANDSLAVE
      INTEGER(8), intent(in)  :: POSELT 
      INTEGER, intent(in)     :: IBEG_BLOCK, FIRST_BLOCK, LAST_BLOCK
      INTEGER, OPTIONAL, intent(in)     :: NASS
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      TYPE(LRB_TYPE), intent(inout)   :: BLR_LorU(:)
      INTEGER, OPTIONAL :: OFFSET_IW
      INTEGER, OPTIONAL :: IW(*)
      INTEGER(8) :: POSELT_LOCAL
      INTEGER    :: IP, LDA
#if defined(BLR_MT)
      INTEGER :: CHUNK
#endif
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      LDA = NFRONT
      IF (LorU.EQ.0.AND.SYM.NE.0.AND.NIV.EQ.2
     &             .AND.(.NOT.LBANDSLAVE)) THEN
        IF (present(NASS)) THEN
          LDA = NASS
       ELSE
          write(*,*) 'Internal error in DMUMPS_BLR_PANEL_LRTRSM'
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      IF (LBANDSLAVE) THEN
        POSELT_LOCAL = POSELT
      ELSE
        POSELT_LOCAL = POSELT + 
     &      int(IBEG_BLOCK-1,8)*int(LDA,8) + int(IBEG_BLOCK - 1,8)
      ENDIF
#if defined(BLR_MT) 
      CHUNK = 1
!$OMP DO 
!$OMP& SCHEDULE(DYNAMIC,CHUNK)
#endif  
      DO IP = FIRST_BLOCK, LAST_BLOCK
        CALL DMUMPS_LRTRSM(A, LA, POSELT_LOCAL, NFRONT, LDA,
     &             BLR_LorU(IP-CURRENT_BLR), NIV, SYM, LorU,
     &             IW, OFFSET_IW) 
      END DO
#if defined(BLR_MT) 
!$OMP END DO NOWAIT
#endif          
      END SUBROUTINE DMUMPS_BLR_PANEL_LRTRSM
      END MODULE DMUMPS_FAC_LR
