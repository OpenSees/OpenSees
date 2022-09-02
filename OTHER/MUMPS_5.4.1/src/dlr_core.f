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
C Note: the last routine of this file, xMUMPS_TRUNCATED_RRQR is derived from
C the LAPACK package, for which BSD 3-clause license applies 
C (see header of the routine).
      MODULE DMUMPS_LR_CORE
      USE MUMPS_LR_COMMON
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_LR_DATA_M
!$    USE OMP_LIB
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE INIT_LRB(LRB_OUT,K,M,N,ISLR)
C This routine simply initializes a LR block but does NOT allocate it
C (allocation occurs somewhere else)
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K,M,N
        LOGICAL,INTENT(IN) :: ISLR
        LRB_OUT%M = M
        LRB_OUT%N = N
        LRB_OUT%K = K
        LRB_OUT%ISLR = ISLR
        NULLIFY(LRB_OUT%Q)
        NULLIFY(LRB_OUT%R)
      END SUBROUTINE INIT_LRB
C
C
      SUBROUTINE IS_FRONT_BLR_CANDIDATE(INODE, NIV, NFRONT, NASS, 
     &                    BLRON, K489,
     &                    K490, K491, K492, K20, K60, IDAD, K38,
     &                    LRSTATUS, N, LRGROUPS)
        INTEGER,INTENT(IN) :: INODE, NFRONT, NASS, BLRON, K489, K490,
     &                        K491, K492, NIV, K20, K60, IDAD, K38
        INTEGER,INTENT(OUT):: LRSTATUS 
        INTEGER, INTENT(IN):: N
        INTEGER, INTENT(IN), OPTIONAL :: LRGROUPS(N)
C
C     Local variables
        LOGICAL :: COMPRESS_PANEL, COMPRESS_CB
        LRSTATUS = 0
        COMPRESS_PANEL = .FALSE.
        IF ((BLRON.NE.0).and.( 
     &        ((K492.LT.0).and.INODE.EQ.abs(K492))
     &        .or.
     &        ( (K492.GT.0).and.(K491.LE.NFRONT)
     &        .and.(K490.LE.NASS)))) THEN
          COMPRESS_PANEL = .TRUE.
C         Compression for NASS =1 is useless
          IF (NASS.LE.1) THEN
            COMPRESS_PANEL =.FALSE. 
          ENDIF
          IF (present(LRGROUPS)) THEN
           IF (LRGROUPS (INODE) .LT. 0) COMPRESS_PANEL = .FALSE.
          ENDIF
        ENDIF
        COMPRESS_CB = .FALSE.
        IF ((BLRON.NE.0).and.
     &        (K489.GT.0.AND.(K489.NE.2.OR.NIV.EQ.2))
     &        .and.( 
     &        ((K492.LT.0).and.INODE.EQ.abs(K492))
     &        .or.
     &      ((K492.GT.0).AND.(NFRONT-NASS.GT.K491)))) 
     &     THEN
          COMPRESS_CB = .TRUE.
        ENDIF
        IF (.NOT.COMPRESS_PANEL) COMPRESS_CB=.FALSE.
        IF (COMPRESS_PANEL.OR.COMPRESS_CB) THEN
          IF (COMPRESS_CB.AND.(.NOT.COMPRESS_PANEL)) THEN
            LRSTATUS = 1 
          ELSE IF (COMPRESS_PANEL.AND.(.NOT.COMPRESS_CB)) THEN
            LRSTATUS = 2 
          ELSE
            LRSTATUS = 3 
          ENDIF
        ELSE 
         LRSTATUS = 0
        ENDIF
C
C       Schur complement cannot be BLR for now
C
        IF ( INODE .EQ. K20 .AND. K60 .NE. 0 ) THEN
          LRSTATUS = 0
        ENDIF
C
C       Do not compress CB of children of root
C
        IF ( IDAD .EQ. K38 .AND. K38 .NE.0 ) THEN
         COMPRESS_CB = .FALSE.
         IF (LRSTATUS.GE.2) THEN
            LRSTATUS = 2
         ELSE
            LRSTATUS = 0
         ENDIF
        ENDIF
      RETURN
      END SUBROUTINE IS_FRONT_BLR_CANDIDATE
      SUBROUTINE ALLOC_LRB(LRB_OUT,K,M,N,ISLR,IFLAG,IERROR,KEEP8)
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K,M,N
        INTEGER,INTENT(INOUT) :: IFLAG, IERROR
        LOGICAL,INTENT(IN) :: ISLR
        INTEGER(8) :: KEEP8(150)
        INTEGER :: MEM, allocok
        DOUBLE PRECISION :: ZERO
        PARAMETER (ZERO = 0.0D0)
        INTEGER(8) :: KEEP8TMPCOPY, KEEP873COPY
        LRB_OUT%M = M
        LRB_OUT%N = N
        LRB_OUT%K = K
        LRB_OUT%ISLR = ISLR
        IF ((M.EQ.0).OR.(N.EQ.0)) THEN 
         nullify(LRB_OUT%Q)
         nullify(LRB_OUT%R)
         RETURN
        ENDIF
        IF (ISLR) THEN
          IF (K.EQ.0) THEN
            nullify(LRB_OUT%Q)
            nullify(LRB_OUT%R)
          ELSE
            allocate(LRB_OUT%Q(M,K),LRB_OUT%R(K,N),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = K*(M+N)
              RETURN
            ENDIF
          ENDIF
        ELSE
          nullify(LRB_OUT%R)
          allocate(LRB_OUT%Q(M,N),stat=allocok)
          IF (allocok > 0) THEN
            IFLAG  = -13
            IERROR = M*N
            RETURN
          ENDIF
        ENDIF
        IF (ISLR) THEN
          MEM = M*K + N*K
        ELSE
          MEM = M*N
        ENDIF
!$OMP   ATOMIC CAPTURE
        KEEP8(69)    = KEEP8(69) + int(MEM,8)
        KEEP8TMPCOPY = KEEP8(69)
!$OMP   END ATOMIC
!$OMP   ATOMIC UPDATE
        KEEP8(68)    = max(KEEP8TMPCOPY, KEEP8(68))
!$OMP   END ATOMIC
!$OMP   ATOMIC CAPTURE
        KEEP8(71)    = KEEP8(71) + int(MEM,8)
        KEEP8TMPCOPY = KEEP8(71)
!$OMP   END ATOMIC
!$OMP   ATOMIC UPDATE
        KEEP8(70)    = max(KEEP8TMPCOPY, KEEP8(70))
!$OMP   END ATOMIC
!$OMP   ATOMIC CAPTURE
        KEEP8(73)    = KEEP8(73) + int(MEM,8)
        KEEP873COPY  = KEEP8(73)
!$OMP   END ATOMIC
!$OMP   ATOMIC UPDATE
        KEEP8(74)    = max(KEEP8(74), KEEP873COPY)
!$OMP   END ATOMIC
        IF ( KEEP873COPY .GT. KEEP8(75) ) THEN
             IFLAG = -19
             CALL MUMPS_SET_IERROR(
     &             (KEEP873COPY-KEEP8(75)), IERROR)
        ENDIF
      RETURN
      END SUBROUTINE ALLOC_LRB
      SUBROUTINE ALLOC_LRB_FROM_ACC(ACC_LRB, LRB_OUT, K, M, N, LorU,
     &                          IFLAG, IERROR, KEEP8)
        TYPE(LRB_TYPE), INTENT(IN) :: ACC_LRB
        TYPE(LRB_TYPE), INTENT(OUT) :: LRB_OUT
        INTEGER,INTENT(IN) :: K, M, N, LorU
        INTEGER,INTENT(INOUT) :: IFLAG, IERROR
        INTEGER(8) :: KEEP8(150)
        INTEGER :: I
        IF (LorU.EQ.1) THEN
          CALL ALLOC_LRB(LRB_OUT,K,M,N,.TRUE.,IFLAG,IERROR,KEEP8)
          IF (IFLAG.LT.0) RETURN
          DO I=1,K
            LRB_OUT%Q(1:M,I) = ACC_LRB%Q(1:M,I)
            LRB_OUT%R(I,1:N) = -ACC_LRB%R(I,1:N)
          ENDDO
        ELSE
          CALL ALLOC_LRB(LRB_OUT,K,N,M,.TRUE.,IFLAG,IERROR,KEEP8)
          IF (IFLAG.LT.0) RETURN
          DO I=1,K
            LRB_OUT%Q(1:N,I) = ACC_LRB%R(I,1:N)
            LRB_OUT%R(I,1:M) = -ACC_LRB%Q(1:M,I)
          ENDDO
        ENDIF
      END SUBROUTINE ALLOC_LRB_FROM_ACC
      SUBROUTINE REGROUPING2(CUT, NPARTSASS, NASS,
     &                   NPARTSCB, NCB, IBCKSZ, ONLYCB, K472)
        INTEGER, INTENT(IN) :: IBCKSZ, NASS, NCB
        INTEGER, INTENT(INOUT) :: NPARTSCB, NPARTSASS
        INTEGER, POINTER, DIMENSION(:) :: CUT
        INTEGER, POINTER, DIMENSION(:) :: NEW_CUT
        INTEGER :: I, INEW, MINSIZE, NEW_NPARTSASS, allocok
        LOGICAL :: ONLYCB, TRACE
        INTEGER, INTENT(IN) :: K472
        INTEGER :: IBCKSZ2,IFLAG,IERROR
        ALLOCATE(NEW_CUT(max(NPARTSASS,1)+NPARTSCB+1),stat=allocok)
        IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = max(NPARTSASS,1)+NPARTSCB+1
           write(*,*) 'Allocation problem in BLR routine REGROUPING2:',
     &          ' not enough memory? memory requested = ' , IERROR
           RETURN
        ENDIF
        CALL COMPUTE_BLR_VCS(K472, IBCKSZ2, IBCKSZ, NASS)
        MINSIZE = int(IBCKSZ2 / 2)
        NEW_NPARTSASS = max(NPARTSASS,1)
        IF (.NOT. ONLYCB) THEN
           NEW_CUT(1) = 1
           INEW = 2
           I = 2
           DO WHILE (I .LE. NPARTSASS + 1)
              NEW_CUT(INEW) = CUT(I)
              TRACE = .FALSE.
              IF (NEW_CUT(INEW) - NEW_CUT(INEW-1) .GT. MINSIZE) THEN
                 INEW = INEW + 1
                 TRACE = .TRUE.
              ENDIF
              I = I + 1
           END DO
           IF (TRACE) THEN
              INEW = INEW - 1 
           ELSE
              IF (INEW .NE. 2) THEN
                 NEW_CUT(INEW-1) = NEW_CUT(INEW)
                 INEW = INEW - 1
              ENDIF
           ENDIF
           NEW_NPARTSASS = INEW - 1
        ENDIF
        IF (ONLYCB) THEN
           DO I=1,max(NPARTSASS,1)+1
              NEW_CUT(I) = CUT(I)
           ENDDO
        ENDIF
        IF (NCB .EQ. 0) GO TO 50
        INEW = NEW_NPARTSASS+2
        I = max(NPARTSASS,1) + 2
        DO WHILE (I .LE. max(NPARTSASS,1) + NPARTSCB + 1)
              NEW_CUT(INEW) = CUT(I)
              TRACE = .FALSE.
              IF (NEW_CUT(INEW) - NEW_CUT(INEW-1) .GT. MINSIZE) THEN
                 INEW = INEW + 1
                 TRACE = .TRUE.
              ENDIF
              I = I + 1
        END DO
        IF (TRACE) THEN
           INEW = INEW - 1 
        ELSE
           IF (INEW .NE.  NEW_NPARTSASS+2) THEN
           NEW_CUT(INEW-1) = NEW_CUT(INEW)
              INEW = INEW - 1
           ENDIF
        ENDIF
        NPARTSCB = INEW - 1 - NEW_NPARTSASS
 50     CONTINUE       
        NPARTSASS = NEW_NPARTSASS
        DEALLOCATE(CUT)
        ALLOCATE(CUT(NPARTSASS+NPARTSCB+1),stat=allocok)
        IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = NPARTSASS+NPARTSCB+1
           write(*,*) 'Allocation problem in BLR routine REGROUPING2:',
     &          ' not enough memory? memory requested = ' , IERROR
           RETURN
        ENDIF
        DO I=1,NPARTSASS+NPARTSCB+1
           CUT(I) = NEW_CUT(I)
        ENDDO
        DEALLOCATE(NEW_CUT)
      END SUBROUTINE REGROUPING2
      SUBROUTINE DMUMPS_LRTRSM(A, LA, POSELT_LOCAL, NFRONT, LDA, LRB, 
     &                NIV, SYM, LorU, IW, OFFSET_IW) 
C     -----------
C     Parameters
C     -----------
      INTEGER(8), intent(in)  :: LA
      INTEGER, intent(in)     :: NFRONT, NIV, SYM, LorU, LDA
      INTEGER(8), intent(in)  :: POSELT_LOCAL
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      TYPE(LRB_TYPE), intent(inout)   :: LRB
      INTEGER, OPTIONAL:: OFFSET_IW
      INTEGER, OPTIONAL :: IW(*)
C     -----------
C     Local variables
C     -----------
      INTEGER(8) :: DPOS, POSPV1, POSPV2, OFFDAG
      INTEGER    :: M, N, I, J
      DOUBLE PRECISION, POINTER  :: LR_BLOCK_PTR(:,:)
      DOUBLE PRECISION :: ONE, MONE, ZERO
      DOUBLE PRECISION :: MULT1, MULT2, A11, DETPIV, A22, A12
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      N  = LRB%N
      IF (LRB%ISLR) THEN
        M  = LRB%K
        LR_BLOCK_PTR => LRB%R
      ELSE
        M  = LRB%M
        LR_BLOCK_PTR => LRB%Q
      END IF
      IF (M.NE.0) THEN
C           Why is it Right, Lower, Tranpose? 
C           Because A is stored by rows 
C           but BLR_L is stored by columns          
        IF (SYM.EQ.0.AND.LorU.EQ.0) THEN
          CALL dtrsm('R', 'L', 'T', 'N', M, N, ONE,
     &                 A(POSELT_LOCAL), NFRONT,
     &                 LR_BLOCK_PTR(1,1), M)
        ELSE
          CALL dtrsm('R', 'U', 'N', 'U', M, N, ONE,
     &                 A(POSELT_LOCAL), LDA,
     &                 LR_BLOCK_PTR(1,1), M)
          IF (LorU.EQ.0) THEN
C           Now apply D scaling              
            IF (.NOT.present(OFFSET_IW)) THEN
              write(*,*) 'Internal error in ',
     &             'DMUMPS_LRTRSM' 
              CALL MUMPS_ABORT()
            ENDIF
            DPOS = POSELT_LOCAL
            I = 1
            DO
            IF(I .GT. N) EXIT
            IF(IW(OFFSET_IW+I-1) .GT. 0) THEN
C             1x1 pivot
              A11 = ONE/A(DPOS)
              CALL dscal(M, A11, LR_BLOCK_PTR(1,I), 1)
              DPOS = DPOS + int(LDA + 1,8)
              I = I+1
            ELSE
C             2x2 pivot
              POSPV1 = DPOS
              POSPV2 = DPOS+ int(LDA + 1,8)
              OFFDAG = POSPV1+1_8
              A11 = A(POSPV1)
              A22 = A(POSPV2)
              A12 = A(OFFDAG)
              DETPIV = A11*A22 - A12**2
              A22 = A11/DETPIV
              A11 = A(POSPV2)/DETPIV
              A12 = -A12/DETPIV
              DO J = 1,M
                MULT1 = A11*LR_BLOCK_PTR(J,I)+A12*LR_BLOCK_PTR(J,I+1)
                MULT2 = A12*LR_BLOCK_PTR(J,I)+A22*LR_BLOCK_PTR(J,I+1)
                LR_BLOCK_PTR(J,I)   = MULT1
                LR_BLOCK_PTR(J,I+1) = MULT2
              ENDDO
              DPOS = POSPV2 + int(LDA + 1,8)
              I = I+2
            ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      CALL UPD_FLOP_TRSM(LRB, LorU)
      END SUBROUTINE DMUMPS_LRTRSM
      SUBROUTINE DMUMPS_LRGEMM_SCALING(LRB, SCALED, A, LA, DIAG, 
     &          LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, MAXI_CLUSTER) 
C This routine does the scaling (for the symmetric case) before 
C computing the LR product (done in DMUMPS_LRGEMM4)        
        TYPE(LRB_TYPE),INTENT(IN) :: LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        DOUBLE PRECISION, intent(inout), DIMENSION(:,:)  :: SCALED
        INTEGER,INTENT(IN) :: LD_DIAG, NFRONT, IW2(*)
        INTEGER(8), INTENT(IN) :: POSELTT
        DOUBLE PRECISION, INTENT(IN), OPTIONAL :: DIAG(*)
        INTEGER, INTENT(IN) :: MAXI_CLUSTER
        DOUBLE PRECISION, intent(inout)  :: BLOCK(MAXI_CLUSTER)
        INTEGER :: J, NROWS
        DOUBLE PRECISION :: PIV1, PIV2, OFFDIAG
        IF (LRB%ISLR) THEN
            NROWS = LRB%K
        ELSE 
            NROWS = LRB%M
        ENDIF
        J = 1
        DO WHILE (J <= LRB%N)
            IF (IW2(J) > 0) THEN
                SCALED(1:NROWS,J) = DIAG(1+LD_DIAG*(J-1)+J-1) 
     &           * SCALED(1:NROWS,J)
                J = J+1
            ELSE !2x2 pivot
                PIV1    = DIAG(1+LD_DIAG*(J-1)+J-1)
                PIV2    = DIAG(1+LD_DIAG*J+J)
                OFFDIAG = DIAG(1+LD_DIAG*(J-1)+J)
                BLOCK(1:NROWS)    = SCALED(1:NROWS,J)
                SCALED(1:NROWS,J) = PIV1 * SCALED(1:NROWS,J)
     &            + OFFDIAG * SCALED(1:NROWS,J+1)
                SCALED(1:NROWS,J+1) = OFFDIAG * BLOCK(1:NROWS)
     &            + PIV2 * SCALED(1:NROWS,J+1)
                 J=J+2
            ENDIF
        END DO
      END SUBROUTINE DMUMPS_LRGEMM_SCALING
      SUBROUTINE DMUMPS_LRGEMM4(ALPHA,
     &           LRB1, LRB2, BETA,
     &           A, LA, POSELTT, NFRONT, SYM, 
     &           IFLAG, IERROR,
     &           MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT,
     &           RANK, BUILDQ, 
     &           LUA_ACTIVATED, 
C Start of OPTIONAL arguments        
     &           LorU,
     &           LRB3, MAXI_RANK, 
     &           MAXI_CLUSTER,
     &           DIAG, LD_DIAG, IW2, BLOCK
     &           )
C
CC
        TYPE(LRB_TYPE),INTENT(IN) :: LRB1,LRB2
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, SYM, TOL_OPT
        INTEGER,INTENT(INOUT) :: IFLAG, IERROR
        INTEGER(8), INTENT(IN) :: POSELTT
        DOUBLE PRECISION, INTENT(IN), OPTIONAL :: DIAG(*)
        INTEGER,INTENT(IN), OPTIONAL :: LD_DIAG, IW2(*)
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT
        DOUBLE PRECISION, intent(in) :: TOLEPS
        DOUBLE PRECISION :: ALPHA,BETA
        LOGICAL, INTENT(OUT) :: BUILDQ
        DOUBLE PRECISION, intent(inout), OPTIONAL  :: BLOCK(*)
        INTEGER, INTENT(IN), OPTIONAL :: LorU
        LOGICAL, INTENT(IN) :: LUA_ACTIVATED
        INTEGER, INTENT(IN), OPTIONAL :: MAXI_CLUSTER
        INTEGER, INTENT(IN), OPTIONAL :: MAXI_RANK
        TYPE(LRB_TYPE), INTENT(INOUT), OPTIONAL :: LRB3
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: XY_YZ
        DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:) :: XQ, R_Y
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: X, Y, Y1, Y2, Z
        CHARACTER(len=1) :: SIDE, TRANSY
        INTEGER :: K_XY, K_YZ, LDY, LDY1, LDY2, K_Y
        INTEGER :: LDXY_YZ, SAVE_K
        INTEGER :: I, J, RANK, MAXRANK, INFO, LWORK
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:), 
     &                          Y_RRQR(:,:)
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION, EXTERNAL ::dnrm2
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
        IF (LRB1%M.EQ.0) THEN
          RETURN
        ENDIF
        IF (LRB2%M.EQ.0) THEN
        ENDIF
        RANK = 0
        BUILDQ = .FALSE.
        IF (LRB1%ISLR.AND.LRB2%ISLR) THEN
            IF ((LRB1%K.EQ.0).OR.(LRB2%K.EQ.0)) THEN
              GOTO 1200
            ENDIF
                allocate(Y(LRB1%K,LRB2%K),stat=allocok)
                IF (allocok > 0) THEN
                    MREQ = LRB1%K*LRB2%K
                    GOTO 1570
                ENDIF
            X    => LRB1%Q
            K_Y  =  LRB1%N
            IF (SYM .EQ. 0) THEN
                Y1  => LRB1%R
            ELSE
                allocate(Y1(LRB1%K,LRB1%N),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%K*LRB1%N
                  GOTO 1570
                ENDIF
                DO J=1,LRB1%N
                    DO I=1,LRB1%K
                        Y1(I,J) = LRB1%R(I,J)
                    ENDDO
                ENDDO
                CALL DMUMPS_LRGEMM_SCALING(LRB1, Y1, A, LA, DIAG,
     &                 LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                 MAXI_CLUSTER) 
            ENDIF
            LDY1 =  LRB1%K
            Z    => LRB2%Q
            Y2   => LRB2%R
            LDY2 =  LRB2%K
             CALL dgemm('N', 'T', LRB1%K, LRB2%K, K_Y, ONE,
     &            Y1(1,1), LDY1, Y2(1,1), LDY2, ZERO, Y(1,1), LRB1%K )
            IF (MIDBLK_COMPRESS.GE.1) THEN 
                LWORK = LRB2%K*(LRB2%K+1)
                allocate(Y_RRQR(LRB1%K,LRB2%K),
     &               WORK_RRQR(LWORK), RWORK_RRQR(2*LRB2%K), 
     &               TAU_RRQR(MIN(LRB1%K,LRB2%K)),
     &               JPVT_RRQR(LRB2%K),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%K*LRB2%K + LWORK + 2*LRB2%K +
     &                   MIN(LRB1%K,LRB2%K) + LRB2%K
                  GOTO 1570
                ENDIF
                DO J=1,LRB2%K
                    DO I=1,LRB1%K
                        Y_RRQR(I,J) = Y(I,J)
                    ENDDO
                ENDDO
                MAXRANK = MIN(LRB1%K, LRB2%K)-1
                MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
                JPVT_RRQR = 0
                CALL DMUMPS_TRUNCATED_RRQR(LRB1%K, LRB2%K, Y_RRQR(1,1),
     &               LRB1%K, JPVT_RRQR, TAU_RRQR, WORK_RRQR,
     &               LRB2%K, RWORK_RRQR, TOLEPS, TOL_OPT, RANK, 
     &               MAXRANK, INFO)
                IF (RANK.GT.MAXRANK) THEN 
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR,
     &                     JPVT_RRQR)
                    BUILDQ = .FALSE.
                ELSE
                    BUILDQ = .TRUE.
                ENDIF
                IF (BUILDQ) THEN
                  IF (RANK.EQ.0) THEN
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR, 
     &                         JPVT_RRQR) 
                      deallocate(Y)
                    nullify(Y)
C                   GOTO 1580 not ok because BUILDQ .EQV. true
C                   would try to free XQ and R_Y that are not allocated
C                   in that case. So we free Y1 now if it was allocated. 
                    IF (SYM .NE. 0) deallocate(Y1)
                    GOTO 1200
                  ELSE
                    allocate(XQ(LRB1%M,RANK), R_Y(RANK,LRB2%K),
     &                       stat=allocok)
                    IF (allocok > 0) THEN
                      MREQ = LRB1%M*RANK + RANK*LRB2%K
                      GOTO 1570
                    ENDIF
                    DO J=1, LRB2%K
                       R_Y(1:MIN(RANK,J),JPVT_RRQR(J)) =
     &                   Y_RRQR(1:MIN(RANK,J),J)
                       IF(J.LT.RANK) R_Y(MIN(RANK,J)+1:
     &                   RANK,JPVT_RRQR(J))= ZERO
                    END DO
C                   LWORK=LRB2%K*(LRB2%K+1), with LRB2%K>RANK
C                   large enough for dorgqr
                    CALL dorgqr 
     &                  (LRB1%K, RANK, RANK, Y_RRQR(1,1),
     &                  LRB1%K, TAU_RRQR(1),  
     &                  WORK_RRQR(1), LWORK, INFO )
                    CALL dgemm('N', 'N', LRB1%M, RANK, LRB1%K, ONE,
     &                    X(1,1), LRB1%M, Y_RRQR(1,1), LRB1%K, ZERO, 
     &                    XQ(1,1), LRB1%M)
                    deallocate(Y_RRQR, WORK_RRQR, RWORK_RRQR, TAU_RRQR, 
     &                         JPVT_RRQR)   
                    nullify(X)
                    X      => XQ
                    K_XY   =  RANK
                      deallocate(Y)
                    nullify(Y)
                    Y      => R_Y
                    LDY    =  RANK
                    K_YZ   =  LRB2%K
                    TRANSY =  'N'
                    SIDE   = 'R'
                  ENDIF
                ENDIF
            ENDIF
            IF (.NOT.BUILDQ) THEN 
                  LDY    = LRB1%K
                K_XY   = LRB1%K
                K_YZ   = LRB2%K
                TRANSY = 'N'
                IF (LRB1%K .GE. LRB2%K) THEN
                    SIDE = 'L'
                ELSE 
                    SIDE = 'R'
                ENDIF
            ENDIF
        ENDIF
        IF (LRB1%ISLR.AND.(.NOT.LRB2%ISLR)) THEN 
            IF (LRB1%K.EQ.0) THEN
              GOTO 1200
            ENDIF
            SIDE   =  'R'
            K_XY   =  LRB1%K
            TRANSY =  'N' 
            Z      => LRB2%Q 
            X   => LRB1%Q
            LDY =  LRB1%K
            IF (SYM .EQ. 0) THEN
              Y   => LRB1%R
            ELSE
              allocate(Y(LRB1%K,LRB1%N),stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LRB1%K*LRB1%N
                GOTO 1570
              ENDIF
              DO J=1,LRB1%N
                  DO I=1,LRB1%K
                      Y(I,J) = LRB1%R(I,J)
                  ENDDO
              ENDDO
              CALL DMUMPS_LRGEMM_SCALING(LRB1, Y, A, LA, DIAG,
     &               LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &               MAXI_CLUSTER) 
            ENDIF
            K_YZ = LRB2%N
        ENDIF
        IF ((.NOT.LRB1%ISLR).AND.LRB2%ISLR) THEN 
            IF (LRB2%K.EQ.0) THEN
              GOTO 1200
            ENDIF
            SIDE   =  'L'
            K_YZ   =  LRB2%K
            X      => LRB1%Q 
            TRANSY =  'T' 
            K_XY = LRB1%N
            IF (SYM .EQ. 0) THEN
                Y  => LRB2%R
            ELSE
                allocate(Y(LRB2%K,LRB2%N),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB2%K*LRB2%N
                  GOTO 1570
                ENDIF
                DO J=1,LRB2%N
                    DO I=1,LRB2%K
                        Y(I,J) = LRB2%R(I,J)
                    ENDDO
                ENDDO
                CALL DMUMPS_LRGEMM_SCALING(LRB2, Y, A, LA, DIAG,
     &                 LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                 MAXI_CLUSTER) 
            ENDIF
            LDY =  LRB2%K
            Z   => LRB2%Q
        ENDIF
        IF ((.NOT.LRB1%ISLR).AND.(.NOT.LRB2%ISLR)) THEN 
            IF (SYM .EQ. 0) THEN
                X => LRB1%Q
            ELSE
                allocate(X(LRB1%M,LRB1%N),stat=allocok)
                IF (allocok > 0) THEN
                  MREQ = LRB1%M*LRB1%N
                  GOTO 1570
                ENDIF
                DO J=1,LRB1%N
                    DO I=1,LRB1%M
                        X(I,J) = LRB1%Q(I,J)
                    ENDDO
                ENDDO
                CALL DMUMPS_LRGEMM_SCALING(LRB1, X, A, LA, DIAG,
     &                   LD_DIAG, IW2, POSELTT, NFRONT, BLOCK, 
     &                   MAXI_CLUSTER) 
            ENDIF
            SIDE   =  'N'
            Z      => LRB2%Q
            K_XY = LRB1%N  
        ENDIF
        IF (LUA_ACTIVATED) THEN
          SAVE_K = LRB3%K
          IF (SIDE == 'L') THEN
            LRB3%K = LRB3%K+K_YZ
          ELSEIF (SIDE == 'R') THEN 
            LRB3%K = LRB3%K+K_XY
          ENDIF
        ENDIF
        IF (SIDE == 'L') THEN ! LEFT: XY_YZ = X*Y; A = XY_YZ*Z
            IF (.NOT.LUA_ACTIVATED) THEN
              allocate(XY_YZ(LRB1%M,K_YZ),stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LRB1%M*K_YZ
                GOTO 1570
              ENDIF
              LDXY_YZ = LRB1%M
            ELSE
               IF (SAVE_K+K_YZ.GT.MAXI_RANK) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 1a',
     &           'K_ACC+K_CUR>K_MAX:',SAVE_K,K_YZ,MAXI_RANK
                CALL MUMPS_ABORT()
              ENDIF
              IF (LRB3%M.NE.LRB1%M) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 1b',
     &           'LRB1%M =/= LRB3%M',LRB1%M,LRB3%M
                CALL MUMPS_ABORT()
              ENDIF
              XY_YZ => LRB3%Q(1:LRB1%M,SAVE_K+1:SAVE_K+K_YZ)
              LDXY_YZ = MAXI_CLUSTER
              DO I=1,K_YZ
                LRB3%R(SAVE_K+I,1:LRB2%M) = Z(1:LRB2%M,I)
              ENDDO
            ENDIF
            CALL dgemm('N', TRANSY, LRB1%M, K_YZ, K_XY, ONE,
     &             X(1,1), LRB1%M, Y(1,1), LDY, ZERO, XY_YZ(1,1), 
     &             LDXY_YZ)
            IF (.NOT.LUA_ACTIVATED) THEN
              CALL dgemm('N', 'T', LRB1%M, LRB2%M, K_YZ, ALPHA,
     &               XY_YZ(1,1), LRB1%M, Z(1,1), LRB2%M, BETA, 
     &               A(POSELTT), NFRONT)
              deallocate(XY_YZ)
            ENDIF
        ELSEIF (SIDE == 'R') THEN ! RIGHT: XY_YZ = Y*Z; A = X*XY_YZ
            IF (.NOT.LUA_ACTIVATED) THEN
              allocate(XY_YZ(K_XY,LRB2%M),stat=allocok)
              IF (allocok > 0) THEN
                MREQ = K_XY*LRB2%M
                GOTO 1570
              ENDIF
              LDXY_YZ = K_XY
            ELSE
               IF (SAVE_K+K_XY.GT.MAXI_RANK) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 2a',
     &           'K_ACC+K_CUR>K_MAX:',SAVE_K,K_XY,MAXI_RANK
                CALL MUMPS_ABORT()
              ENDIF
              IF (LRB3%N.NE.LRB2%M) THEN
                write(*,*) 'Internal error in DMUMPS_LRGEMM4 2b',
     &           'LRB2%M =/= LRB3%N',LRB2%M,LRB3%N
                CALL MUMPS_ABORT()
              ENDIF
              XY_YZ => LRB3%R(SAVE_K+1:SAVE_K+K_XY,1:LRB2%M)
              LDXY_YZ = MAXI_RANK
              DO I=1,K_XY
                LRB3%Q(1:LRB1%M,SAVE_K+I) = X(1:LRB1%M,I)
              ENDDO
            ENDIF
            CALL dgemm(TRANSY, 'T', K_XY, LRB2%M, K_YZ, ONE,
     &             Y(1,1), LDY, Z(1,1), LRB2%M, ZERO, XY_YZ(1,1), 
     &             LDXY_YZ)
            IF (.NOT.LUA_ACTIVATED) THEN
              CALL dgemm('N', 'N', LRB1%M, LRB2%M, K_XY, ALPHA,
     &               X(1,1), LRB1%M, XY_YZ(1,1), K_XY, BETA, A(POSELTT),
     &               NFRONT)
              deallocate(XY_YZ)
            ENDIF
        ELSE ! SIDE == 'N' : NONE; A = X*Z
          CALL dgemm('N', 'T', LRB1%M, LRB2%M, K_XY, ALPHA,
     &             X(1,1), LRB1%M, Z(1,1), LRB2%M, BETA, A(POSELTT),
     &             NFRONT)
        ENDIF
        GOTO 1580
 1570 CONTINUE        
C       Alloc NOT ok!!        
        IFLAG  = -13
        IERROR = MREQ
        RETURN
 1580 CONTINUE       
C       Alloc ok!!        
        IF ((.NOT.LRB1%ISLR).AND.(.NOT.LRB2%ISLR)) THEN 
            IF (SYM .NE. 0) deallocate(X)
        ELSEIF ((.NOT.LRB1%ISLR).AND.LRB2%ISLR) THEN 
            IF (SYM .NE. 0) deallocate(Y)
        ELSEIF (LRB1%ISLR.AND.(.NOT.LRB2%ISLR)) THEN 
            IF (SYM .NE. 0) deallocate(Y)
        ELSE
            IF (SYM .NE. 0) deallocate(Y1)
            IF ((MIDBLK_COMPRESS.GE.1).AND.BUILDQ) THEN
                deallocate(XQ)
                deallocate(R_Y)
            ELSE
                  deallocate(Y)
            ENDIF
        ENDIF
 1200 CONTINUE       
      END SUBROUTINE DMUMPS_LRGEMM4
      SUBROUTINE DMUMPS_DECOMPRESS_ACC(ACC_LRB, MAXI_CLUSTER, 
     &           MAXI_RANK, A, LA, POSELTT, NFRONT, NIV, LorU, 
     &           COUNT_FLOPS)
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK
        INTEGER(8), INTENT(IN) :: POSELTT
        LOGICAL, OPTIONAL :: COUNT_FLOPS
        LOGICAL  :: COUNT_FLOPS_LOC
        INTEGER  :: LorU
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
        IF (present(COUNT_FLOPS)) THEN
          COUNT_FLOPS_LOC=COUNT_FLOPS
        ELSE
          COUNT_FLOPS_LOC=.TRUE.
        ENDIF
        CALL dgemm('N', 'N', ACC_LRB%M, ACC_LRB%N, ACC_LRB%K,
     &         MONE, ACC_LRB%Q(1,1), MAXI_CLUSTER, ACC_LRB%R(1,1), 
     &         MAXI_RANK, ONE, A(POSELTT), NFRONT)
        ACC_LRB%K = 0
      END SUBROUTINE DMUMPS_DECOMPRESS_ACC
      SUBROUTINE DMUMPS_COMPRESS_FR_UPDATES(ACC_LRB, MAXI_CLUSTER, 
     &          MAXI_RANK, A, LA, POSELTT, NFRONT, NIV,
     &          TOLEPS, TOL_OPT, KPERCENT, BUILDQ, LorU, CB_COMPRESS)
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, LorU, TOL_OPT
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT
        INTEGER(8), INTENT(IN) :: POSELTT
        DOUBLE PRECISION, intent(in) :: TOLEPS
        LOGICAL, INTENT(OUT) :: BUILDQ
        LOGICAL, INTENT(IN) :: CB_COMPRESS
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:)
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        INTEGER :: INFO, RANK, MAXRANK, LWORK
        INTEGER :: I, J, M, N
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
              M = ACC_LRB%M
              N = ACC_LRB%N
              MAXRANK = floor(dble(M*N)/dble(M+N))
              MAXRANK = max (1, int((MAXRANK*KPERCENT/100)))
              LWORK = N*(N+1)
              allocate(WORK_RRQR(LWORK), RWORK_RRQR(2*N), 
     &             TAU_RRQR(N),
     &             JPVT_RRQR(N), stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LWORK +4 *N
                GOTO 100
              ENDIF
                DO I=1,N
                  ACC_LRB%Q(1:M,I)= 
     &            - A(POSELTT+int(I-1,8)*int(NFRONT,8) :
     &            POSELTT+int(I-1,8)*int(NFRONT,8) + int(M-1,8) )
                END DO  
              JPVT_RRQR = 0
              CALL DMUMPS_TRUNCATED_RRQR(M, N, ACC_LRB%Q(1,1),
     &               MAXI_CLUSTER, JPVT_RRQR(1), TAU_RRQR(1), 
     &               WORK_RRQR(1),
     &               N, RWORK_RRQR(1), TOLEPS, TOL_OPT, 
     &               RANK, MAXRANK, INFO)
              BUILDQ = (RANK.LE.MAXRANK)
              IF (BUILDQ) THEN
                DO J=1, N
                   ACC_LRB%R(1:MIN(RANK,J),JPVT_RRQR(J)) =
     &               ACC_LRB%Q(1:MIN(RANK,J),J)
                   IF(J.LT.RANK) ACC_LRB%R(MIN(RANK,J)+1:
     &               RANK,JPVT_RRQR(J))= ZERO
                END DO
                CALL dorgqr 
     &              (M, RANK, RANK, ACC_LRB%Q(1,1),
     &              MAXI_CLUSTER, TAU_RRQR(1),  
     &              WORK_RRQR(1), LWORK, INFO )
                DO I=1,N
                  A( POSELTT+int(I-1,8)*int(NFRONT,8) :
     &              POSELTT+int(I-1,8)*int(NFRONT,8)+int(M-1,8) ) = ZERO
                END DO  
                ACC_LRB%K = RANK
                CALL UPD_FLOP_COMPRESS(ACC_LRB, CB_COMPRESS=CB_COMPRESS)
              ELSE
                ACC_LRB%K = RANK
                ACC_LRB%ISLR = .FALSE.
                CALL UPD_FLOP_COMPRESS(ACC_LRB, CB_COMPRESS=CB_COMPRESS)
                ACC_LRB%ISLR = .TRUE.
                ACC_LRB%K = 0
              ENDIF
              deallocate(JPVT_RRQR, TAU_RRQR, WORK_RRQR, RWORK_RRQR)
              RETURN
 100        CONTINUE        
C     Alloc NOT ok!!
            write(*,*) 'Allocation problem in BLR routine 
     &        DMUMPS_COMPRESS_FR_UPDATES: ',
     &        'not enough memory? memory requested = ' , MREQ
            CALL MUMPS_ABORT()
            RETURN
      END SUBROUTINE DMUMPS_COMPRESS_FR_UPDATES
      SUBROUTINE DMUMPS_RECOMPRESS_ACC(ACC_LRB, MAXI_CLUSTER, 
     &          MAXI_RANK, A, LA, POSELTT, NFRONT, NIV,
     &          MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB, 
     &          KPERCENT_LUA, NEW_ACC_RANK)
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, TOL_OPT
        INTEGER :: IFLAG, IERROR
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT_LUA
        INTEGER,INTENT(INOUT) :: NEW_ACC_RANK
        INTEGER(8), INTENT(IN) :: POSELTT
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
        DOUBLE PRECISION, intent(in) :: TOLEPS
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), TARGET :: Q1, R1, 
     &                                                  Q2, R2
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        TYPE(LRB_TYPE)       :: LRB1, LRB2
        INTEGER :: INFO, RANK1, RANK2, RANK, MAXRANK, LWORK
        LOGICAL :: BUILDQ, BUILDQ1, BUILDQ2, SKIP1, SKIP2
        INTEGER :: I, J, M, N, K
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
              SKIP1 = .FALSE.
              SKIP2 = .FALSE.
              SKIP1 = .TRUE.
 1500 CONTINUE              
              M = ACC_LRB%M
              N = ACC_LRB%N
              K = ACC_LRB%K
              MAXRANK = K-1
              MAXRANK = max (1, int((MAXRANK*KPERCENT_LUA/100)))
              LWORK = K*(K+1)
              IF (.FALSE.) THEN
                CALL DMUMPS_RECOMPRESS_ACC_V2(ACC_LRB,
     &                 MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT,
     &                 NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS,
     &                 TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, 
     &                 NEW_ACC_RANK)
                K = ACC_LRB%K
                MAXRANK = K-1
                MAXRANK = max (1, int((MAXRANK*KPERCENT_LUA/100)))
                LWORK = K*(K+1)
                SKIP1 = .TRUE.
                SKIP2 = K.EQ.0
              ENDIF
              IF (SKIP1.AND.SKIP2) GOTO 1600
              allocate(Q1(M,K), Q2(N,K), 
     &             WORK_RRQR(LWORK), 
     &             RWORK_RRQR(2*K), 
     &             TAU_RRQR(K),
     &             JPVT_RRQR(K), stat=allocok)
              IF (allocok > 0) THEN
                MREQ = LWORK + M*N + N*K+ 4 * K
                GOTO 100
              ENDIF
              IF (SKIP1) THEN
                BUILDQ1 = .FALSE.
              ELSE
                DO J=1,K
                  DO I=1,M
                    Q1(I,J) = ACC_LRB%Q(I,J)
                  ENDDO
                ENDDO
                JPVT_RRQR = 0
                CALL DMUMPS_TRUNCATED_RRQR(M, K, Q1(1,1),
     &               M, JPVT_RRQR(1), TAU_RRQR(1), WORK_RRQR(1),
     &               K, RWORK_RRQR(1), TOLEPS, TOL_OPT, RANK1, 
     &               MAXRANK, INFO)
                BUILDQ1 = (RANK1.LE.MAXRANK)
              ENDIF
              IF (BUILDQ1) THEN
                  allocate(R1(RANK1,K), stat=allocok)
                  IF (allocok > 0) THEN
                     MREQ = RANK1*K
                     GOTO 100
                  ENDIF
                  DO J=1, K
                     R1(1:MIN(RANK1,J),JPVT_RRQR(J)) =
     &                 Q1(1:MIN(RANK1,J),J)
                     IF(J.LT.RANK1) R1(MIN(RANK1,J)+1:
     &                 RANK1,JPVT_RRQR(J))= ZERO
                  END DO
                  CALL dorgqr 
     &                (M, RANK1, RANK1, Q1(1,1),
     &                M, TAU_RRQR(1),  
     &                WORK_RRQR(1), LWORK, INFO )
              ENDIF
              IF (SKIP2) THEN
                BUILDQ2 = .FALSE.
              ELSE
                DO J=1,K
                  DO I=1,N
                    Q2(I,J) = ACC_LRB%R(J,I)
                  ENDDO
                ENDDO
                JPVT_RRQR = 0
                CALL DMUMPS_TRUNCATED_RRQR(N, K, Q2(1,1),
     &               N, JPVT_RRQR(1), TAU_RRQR(1), WORK_RRQR(1),
     &               K, RWORK_RRQR(1), TOLEPS, TOL_OPT, 
     &               RANK2, MAXRANK, INFO)
                BUILDQ2 = (RANK2.LE.MAXRANK)
              ENDIF
              IF (BUILDQ2) THEN
                  allocate(R2(RANK2,K), stat=allocok)
                  IF (allocok > 0) THEN
                     MREQ = RANK2*K
                     GOTO 100
                  ENDIF
                  DO J=1, K
                     R2(1:MIN(RANK2,J),JPVT_RRQR(J)) =
     &                 Q2(1:MIN(RANK2,J),J)
                     IF(J.LT.RANK2) R2(MIN(RANK2,J)+1:
     &                 RANK2,JPVT_RRQR(J))= ZERO
                  END DO
                  CALL dorgqr 
     &                (N, RANK2, RANK2, Q2(1,1),
     &                N, TAU_RRQR(1),  
     &                WORK_RRQR(1), LWORK, INFO )
              ENDIF
              CALL INIT_LRB(LRB1,RANK1,M,K,BUILDQ1)
              CALL INIT_LRB(LRB2,RANK2,N,K,BUILDQ2)
              IF (BUILDQ1.OR.BUILDQ2) THEN
                IF (BUILDQ1) THEN
                  LRB1%R => R1
                ELSE
                  DO J=1,K
                    DO I=1,M
                      Q1(I,J) = ACC_LRB%Q(I,J)
                    ENDDO
                  ENDDO
                ENDIF
                LRB1%Q => Q1
                IF (BUILDQ2) THEN
                  LRB2%R => R2
                ELSE
                  DO J=1,K
                    DO I=1,N
                      Q2(I,J) = ACC_LRB%R(J,I)
                    ENDDO
                  ENDDO
                ENDIF
                LRB2%Q => Q2
                ACC_LRB%K = 0
                CALL DMUMPS_LRGEMM4(MONE, LRB1, LRB2, ONE,
     &             A, LA, POSELTT, NFRONT, 0, IFLAG, IERROR,
     &             MIDBLK_COMPRESS-1, TOLEPS, TOL_OPT, 
     &             KPERCENT_RMB, 
     &             RANK, BUILDQ, .TRUE., LRB3=ACC_LRB,
     &             MAXI_RANK=MAXI_RANK, MAXI_CLUSTER=MAXI_CLUSTER)
                IF (IFLAG.LT.0) GOTO 100
                CALL UPD_FLOP_UPDATE(LRB1, LRB2, 
     &                  MIDBLK_COMPRESS-1, RANK, BUILDQ,
     &                  .TRUE., .FALSE., REC_ACC=.TRUE.) 
              ENDIF
              IF (.NOT. SKIP1)
     &          CALL UPD_FLOP_COMPRESS(LRB1, REC_ACC=.TRUE.)
              IF (.NOT. SKIP2)
     &          CALL UPD_FLOP_COMPRESS(LRB2, REC_ACC=.TRUE.)
              deallocate(Q1,Q2)
              IF (BUILDQ1) deallocate(R1)
              IF (BUILDQ2) deallocate(R2)
              deallocate(JPVT_RRQR, TAU_RRQR, WORK_RRQR, RWORK_RRQR)
              IF (SKIP1.AND.(RANK2.GT.0)) THEN
                SKIP1 = .FALSE.
                SKIP2 = .TRUE.
                GOTO 1500
              ENDIF
 1600         CONTINUE              
              NEW_ACC_RANK = 0
              RETURN
 100          CONTINUE
C     Alloc NOT ok!!
              write(*,*) 'Allocation problem in BLR routine 
     &          DMUMPS_RECOMPRESS_ACC: ',
     &          'not enough memory? memory requested = ' , MREQ        
              CALL MUMPS_ABORT()             
              RETURN
      END SUBROUTINE DMUMPS_RECOMPRESS_ACC
      RECURSIVE SUBROUTINE DMUMPS_RECOMPRESS_ACC_NARYTREE(
     &          ACC_LRB, MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT, KEEP8,
     &          NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT, 
     &          KPERCENT_RMB,
     &          KPERCENT_LUA, K478, RANK_LIST, POS_LIST, NB_NODES,
     &          LEVEL, ACC_TMP)
        TYPE(LRB_TYPE),TARGET,INTENT(INOUT) :: ACC_LRB
        TYPE(LRB_TYPE),TARGET,INTENT(INOUT),OPTIONAL :: ACC_TMP
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, TOL_OPT
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT_LUA
        INTEGER(8), INTENT(IN) :: POSELTT
        INTEGER(8) :: KEEP8(150)
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
        DOUBLE PRECISION, intent(in) :: TOLEPS
        INTEGER,INTENT(IN) :: K478, NB_NODES, LEVEL
        INTEGER,INTENT(INOUT) :: RANK_LIST(NB_NODES), POS_LIST(NB_NODES)
        TYPE(LRB_TYPE)       :: LRB, ACC_NEW
        TYPE(LRB_TYPE), POINTER :: LRB_PTR
        LOGICAL :: RESORT
        INTEGER :: I, J, M, N, L, NODE_RANK, NARY, IOFF, IMAX, CURPOS
        INTEGER :: NB_NODES_NEW, KTOT, NEW_ACC_RANK
        INTEGER, ALLOCATABLE :: RANK_LIST_NEW(:), POS_LIST_NEW(:)
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
        INTEGER :: allocok
        RESORT = .FALSE.
        M = ACC_LRB%M
        N = ACC_LRB%N
        NARY = -K478
        IOFF = 0
        NB_NODES_NEW = NB_NODES/NARY
        IF (NB_NODES_NEW*NARY.NE.NB_NODES) THEN
          NB_NODES_NEW = NB_NODES_NEW + 1
        ENDIF
        ALLOCATE(RANK_LIST_NEW(NB_NODES_NEW),POS_LIST_NEW(NB_NODES_NEW),
     &       stat=allocok)
        IF (allocok > 0) THEN
           write(*,*) 'Allocation error of RANK_LIST_NEW/POS_LIST_NEW ',
     &      'in DMUMPS_RECOMPRESS_ACC_NARYTREE'
           call MUMPS_ABORT()
        ENDIF
        DO J=1,NB_NODES_NEW
          NODE_RANK = RANK_LIST(IOFF+1)
          CURPOS = POS_LIST(IOFF+1)
          IMAX = MIN(NARY,NB_NODES-IOFF)
          IF (IMAX.GE.2) THEN
            DO I=2,IMAX
              IF (POS_LIST(IOFF+I).NE.CURPOS+NODE_RANK) THEN
                DO L=0,RANK_LIST(IOFF+I)-1
                  ACC_LRB%Q(1:M,CURPOS+NODE_RANK+L) = 
     &               ACC_LRB%Q(1:M,POS_LIST(IOFF+I)+L)
                  ACC_LRB%R(CURPOS+NODE_RANK+L,1:N) = 
     &               ACC_LRB%R(POS_LIST(IOFF+I)+L,1:N)
                ENDDO
                POS_LIST(IOFF+I) = CURPOS+NODE_RANK
              ENDIF
              NODE_RANK = NODE_RANK+RANK_LIST(IOFF+I)
            ENDDO
            CALL INIT_LRB(LRB,NODE_RANK,M,N,.TRUE.)
            IF (.NOT.RESORT.OR.LEVEL.EQ.0) THEN
              LRB%Q => ACC_LRB%Q(1:M,CURPOS:CURPOS+NODE_RANK)
              LRB%R => ACC_LRB%R(CURPOS:CURPOS+NODE_RANK,1:N)
            ELSE
              LRB%Q => ACC_TMP%Q(1:M,CURPOS:CURPOS+NODE_RANK)
              LRB%R => ACC_TMP%R(CURPOS:CURPOS+NODE_RANK,1:N)
            ENDIF
            NEW_ACC_RANK = NODE_RANK-RANK_LIST(IOFF+1)
            IF (NEW_ACC_RANK.GT.0) THEN
              CALL DMUMPS_RECOMPRESS_ACC(LRB,
     &             MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT,
     &             NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, 
     &             TOL_OPT, KPERCENT_RMB, KPERCENT_LUA, NEW_ACC_RANK)
            ENDIF
            RANK_LIST_NEW(J) = LRB%K
            POS_LIST_NEW(J) = CURPOS
          ELSE
            RANK_LIST_NEW(J) = NODE_RANK
            POS_LIST_NEW(J) = CURPOS
          ENDIF
          IOFF = IOFF+IMAX
        ENDDO
        IF (NB_NODES_NEW.GT.1) THEN
          IF (RESORT) THEN
            KTOT = SUM(RANK_LIST_NEW)
            CALL INIT_LRB(ACC_NEW,KTOT,M,N,.TRUE.)
            ALLOCATE(ACC_NEW%Q(MAXI_CLUSTER,MAXI_RANK),
     &           ACC_NEW%R(MAXI_RANK,MAXI_CLUSTER), stat=allocok)
            IF (allocok > 0) THEN
               write(*,*) 'Allocation error of ACC_NEW%Q/ACC_NEW%R ',
     &              'in DMUMPS_RECOMPRESS_ACC_NARYTREE'
               call MUMPS_ABORT()
            ENDIF
            CALL MUMPS_SORT_INT(NB_NODES_NEW, RANK_LIST_NEW,
     &                                         POS_LIST_NEW)
            CURPOS = 1
            IF (LEVEL.EQ.0) THEN
              LRB_PTR => ACC_LRB
            ELSE
              LRB_PTR => ACC_TMP
            ENDIF
            DO J=1,NB_NODES_NEW
              DO L=0,RANK_LIST_NEW(J)-1
                ACC_NEW%Q(1:M,CURPOS+L) = 
     &               LRB_PTR%Q(1:M,POS_LIST_NEW(J)+L)
                ACC_NEW%R(CURPOS+L,1:N) = 
     &               LRB_PTR%R(POS_LIST_NEW(J)+L,1:N)
              ENDDO
              POS_LIST_NEW(J) = CURPOS
              CURPOS = CURPOS + RANK_LIST_NEW(J)
            ENDDO
            IF (LEVEL.GT.0) THEN
              CALL DEALLOC_LRB(ACC_TMP, KEEP8)
            ENDIF
            CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST_NEW, POS_LIST_NEW, NB_NODES_NEW,
     &               LEVEL+1, ACC_NEW)
          ELSE
            CALL DMUMPS_RECOMPRESS_ACC_NARYTREE(ACC_LRB,
     &               MAXI_CLUSTER, MAXI_RANK, A, LA, POSELTT, KEEP8,
     &               NFRONT, NIV, MIDBLK_COMPRESS, TOLEPS, TOL_OPT,
     &               KPERCENT_RMB, KPERCENT_LUA, K478, 
     &               RANK_LIST_NEW, POS_LIST_NEW, NB_NODES_NEW, LEVEL+1)
          ENDIF
        ELSE
          IF (POS_LIST_NEW(1).NE.1) THEN
            write(*,*) 'Internal error in ',
     &      'DMUMPS_RECOMPRESS_ACC_NARYTREE', POS_LIST_NEW(1) 
          ENDIF
          ACC_LRB%K = RANK_LIST_NEW(1)
          IF (RESORT.AND.LEVEL.GT.0) THEN
             DO L=1,ACC_LRB%K
                DO I=1,M
                   ACC_LRB%Q(I,L) = ACC_TMP%Q(I,L)
                ENDDO
                DO I=1,N
                   ACC_LRB%R(L,I) = ACC_TMP%R(L,I)
                ENDDO
            ENDDO
            CALL DEALLOC_LRB(ACC_TMP, KEEP8)
          ENDIF
        ENDIF
        DEALLOCATE(RANK_LIST_NEW, POS_LIST_NEW)
      END SUBROUTINE DMUMPS_RECOMPRESS_ACC_NARYTREE
      SUBROUTINE DMUMPS_RECOMPRESS_ACC_V2(ACC_LRB, MAXI_CLUSTER,
     &          MAXI_RANK, A, LA, POSELTT, NFRONT, NIV,
     &          MIDBLK_COMPRESS, TOLEPS, TOL_OPT, KPERCENT_RMB, 
     &          KPERCENT_LUA, NEW_ACC_RANK)
        TYPE(LRB_TYPE),INTENT(INOUT) :: ACC_LRB
        INTEGER(8), intent(in)  :: LA
        DOUBLE PRECISION, intent(inout)  :: A(LA)
        INTEGER,INTENT(IN) :: NFRONT, NIV, TOL_OPT
        INTEGER,INTENT(IN) :: MAXI_CLUSTER, MAXI_RANK, KPERCENT_LUA
        INTEGER,INTENT(INOUT) :: NEW_ACC_RANK
        INTEGER(8), INTENT(IN) :: POSELTT
        INTEGER,intent(in) :: MIDBLK_COMPRESS, KPERCENT_RMB
        DOUBLE PRECISION, intent(in) :: TOLEPS
        DOUBLE PRECISION,    ALLOCATABLE :: RWORK_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE :: WORK_RRQR(:), TAU_RRQR(:)
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), TARGET ::
     &                                Q1, R1, Q2, PROJ
        INTEGER, ALLOCATABLE :: JPVT_RRQR(:)
        INTEGER :: INFO, RANK1, MAXRANK, LWORK
        LOGICAL :: BUILDQ1
        INTEGER :: I, J, M, N, K, K1
        INTEGER :: allocok, MREQ
        DOUBLE PRECISION :: ONE, MONE, ZERO
        PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
        PARAMETER (ZERO=0.0D0)
              M = ACC_LRB%M
              N = ACC_LRB%N
              K = NEW_ACC_RANK
              K1 = ACC_LRB%K - K
              MAXRANK = K-1
              MAXRANK = max (1, int((MAXRANK*KPERCENT_LUA/100)))
              LWORK = K*(K+1)
              allocate(Q1(M,K), PROJ(K1, K),
     &             WORK_RRQR(LWORK), RWORK_RRQR(2*K), 
     &             TAU_RRQR(K),
     &             JPVT_RRQR(K), stat=allocok)
              IF (allocok > 0) THEN
                MREQ = M * K + K1 * K + LWORK + 4 * K
                GOTO 100
              ENDIF
              DO J=1,K
                  DO I=1,M
                      Q1(I,J) = ACC_LRB%Q(I,J+K1)
                  ENDDO
              ENDDO
              CALL dgemm('T', 'N', K1, K, M, ONE, ACC_LRB%Q(1,1), 
     &                  MAXI_CLUSTER, Q1(1,1), M, ZERO, PROJ(1,1), K1)
              CALL dgemm('N', 'N', M, K, K1, MONE, ACC_LRB%Q(1,1),
     &                  MAXI_CLUSTER, PROJ(1,1), K1, ONE, Q1(1,1), M)
              JPVT_RRQR = 0
              CALL DMUMPS_TRUNCATED_RRQR(M, K, Q1(1,1),
     &               M, JPVT_RRQR(1), TAU_RRQR(1), WORK_RRQR(1),
     &               K, RWORK_RRQR(1), TOLEPS, TOL_OPT, 
     &               RANK1, MAXRANK, INFO)
              BUILDQ1 = (RANK1.LE.MAXRANK)
              IF (BUILDQ1) THEN
                allocate(Q2(N,K), stat=allocok)
                IF (allocok > 0) THEN
                   MREQ = N*K
                   GOTO 100
                ENDIF
                DO J=1,K
                    DO I=1,N
                        Q2(I,J) = ACC_LRB%R(J+K1,I)
                    ENDDO
                ENDDO
                CALL dgemm('N', 'T', K1, N, K, ONE, PROJ(1,1), K1,
     &                  Q2(1,1), N, ONE, ACC_LRB%R(1,1), MAXI_RANK)
                IF (RANK1.GT.0) THEN
                   allocate(R1(RANK1,K), stat=allocok)
                   IF (allocok > 0) THEN
                      MREQ = RANK1*K
                      GOTO 100
                   ENDIF
                  DO J=1, K
                     R1(1:MIN(RANK1,J),JPVT_RRQR(J)) =
     &                 Q1(1:MIN(RANK1,J),J)
                     IF(J.LT.RANK1) R1(MIN(RANK1,J)+1:
     &                 RANK1,JPVT_RRQR(J))= ZERO
                  END DO
                  CALL dorgqr 
     &                (M, RANK1, RANK1, Q1(1,1),
     &                M, TAU_RRQR(1),  
     &                WORK_RRQR(1), LWORK, INFO )
                  DO J=1,K
                    DO I=1,M
                        ACC_LRB%Q(I,J+K1) = Q1(I,J)
                    ENDDO
                  ENDDO
                  CALL dgemm('N', 'T', RANK1, N, K, ONE, R1(1,1), RANK1,
     &                   Q2(1,1), N, ZERO, ACC_LRB%R(K1+1,1), MAXI_RANK)
                  deallocate(R1)
                ENDIF
                deallocate(Q2)
                ACC_LRB%K = K1 + RANK1
              ENDIF
              deallocate(PROJ)
              deallocate(Q1, JPVT_RRQR, TAU_RRQR, WORK_RRQR, RWORK_RRQR)
              RETURN
 100          CONTINUE
C     Alloc NOT ok!!
              write(*,*) 'Allocation problem in BLR routine 
     &          DMUMPS_RECOMPRESS_ACC_V2: ',
     &          'not enough memory? memory requested = ' , MREQ        
              CALL MUMPS_ABORT()
              RETURN
      END SUBROUTINE DMUMPS_RECOMPRESS_ACC_V2
      SUBROUTINE MAX_CLUSTER(CUT,CUT_SIZE,MAXI_CLUSTER)
        INTEGER, intent(in) :: CUT_SIZE
        INTEGER, intent(out) :: MAXI_CLUSTER
        INTEGER, POINTER, DIMENSION(:) :: CUT
        INTEGER :: I
        MAXI_CLUSTER = 0
        DO I = 1, CUT_SIZE
          IF (CUT(I+1) - CUT(I) .GE. MAXI_CLUSTER) THEN
            MAXI_CLUSTER = CUT(I+1) - CUT(I)
          END IF
        END DO
      END SUBROUTINE MAX_CLUSTER
      SUBROUTINE DMUMPS_GET_LUA_ORDER(NB_BLOCKS, ORDER, RANK, IWHANDLER,
     &                       SYM, FS_OR_CB, I, J, FRFR_UPDATES, 
     &                       LBANDSLAVE_IN, K474, BLR_U_COL)
C     -----------
C     Parameters
C     -----------
        INTEGER, INTENT(IN) :: NB_BLOCKS, IWHANDLER, SYM, FS_OR_CB, I, J
        INTEGER, INTENT(OUT) :: ORDER(NB_BLOCKS), RANK(NB_BLOCKS), 
     &                          FRFR_UPDATES
        LOGICAL, OPTIONAL, INTENT(IN) :: LBANDSLAVE_IN
        INTEGER, OPTIONAL, INTENT(IN) :: K474
        TYPE(LRB_TYPE), POINTER, OPTIONAL :: BLR_U_COL(:)
C     -----------
C     Local variables
C     -----------
        INTEGER :: K, IND_L, IND_U
        LOGICAL :: LBANDSLAVE
        TYPE(LRB_TYPE), POINTER :: BLR_L(:), BLR_U(:)
        IF (PRESENT(LBANDSLAVE_IN)) THEN
          LBANDSLAVE = LBANDSLAVE_IN
        ELSE
          LBANDSLAVE = .FALSE.
        ENDIF
        IF ((SYM.NE.0).AND.(FS_OR_CB.EQ.0).AND.(J.NE.0)) THEN
          write(6,*) 'Internal error in DMUMPS_GET_LUA_ORDER',
     &     'SYM, FS_OR_CB, J = ',SYM,FS_OR_CB,J
          CALL MUMPS_ABORT()
        ENDIF
        FRFR_UPDATES = 0
        DO K = 1, NB_BLOCKS
          ORDER(K) = K
          IF (FS_OR_CB.EQ.0) THEN ! FS
            IF (J.EQ.0) THEN ! L panel
              IND_L = NB_BLOCKS+I-K
              IND_U = NB_BLOCKS+1-K
            ELSE ! U panel
              IND_L = NB_BLOCKS+1-K
              IND_U = NB_BLOCKS+I-K
            ENDIF
          ELSE ! CB
            IND_L = I-K
            IND_U = J-K
          ENDIF
          IF (LBANDSLAVE) THEN
            IND_L = I
            IF (K474.GE.2) THEN
              IND_U = K
            ENDIF
          ENDIF
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         0, ! L Panel
     &         K, BLR_L)
          IF (SYM.EQ.0) THEN
            IF (LBANDSLAVE.AND.K474.GE.2) THEN
              BLR_U => BLR_U_COL
            ELSE
              CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IWHANDLER,
     &         1, ! L Panel
     &         K, BLR_U)
            ENDIF
          ELSE
            BLR_U => BLR_L
          ENDIF
          IF (BLR_L(IND_L)%ISLR) THEN
            IF (BLR_U(IND_U)%ISLR) THEN
              RANK(K) = min(BLR_L(IND_L)%K, BLR_U(IND_U)%K)
            ELSE
              RANK(K) = BLR_L(IND_L)%K
            ENDIF
          ELSE
            IF (BLR_U(IND_U)%ISLR) THEN
              RANK(K) = BLR_U(IND_U)%K
            ELSE
              RANK(K) = -1
              FRFR_UPDATES = FRFR_UPDATES + 1
            ENDIF
          ENDIF
        ENDDO
        CALL MUMPS_SORT_INT(NB_BLOCKS, RANK, ORDER)
      END SUBROUTINE DMUMPS_GET_LUA_ORDER
      SUBROUTINE DMUMPS_BLR_ASM_NIV1 (A, LA, POSEL1, NFRONT, NASS1, 
     &      IWHANDLER, SON_IW, LIW, LSTK, NELIM, K1, K2, SYM,
     &      KEEP, KEEP8, OPASSW)
C
C Purpose
C =======
C     
C       Called by a level 1 master assembling the contribution
C       block of a level 1 son that has been BLR-compressed        
C
C
C Parameters
C ==========
C
        INTEGER(8) :: LA, POSEL1
        INTEGER :: LIW, NFRONT, NASS1, LSTK, NELIM, K1, K2, IWHANDLER
        DOUBLE PRECISION :: A(LA)
C       INTEGER :: SON_IW(LIW)
        INTEGER :: SON_IW(:) ! contiguity information lost but no copy
        INTEGER :: KEEP(500)
        INTEGER(8) :: KEEP8(150)
        INTEGER :: SYM
        DOUBLE PRECISION, INTENT(INOUT) :: OPASSW
C
C Local variables
C ===============
C
        DOUBLE PRECISION, ALLOCATABLE :: SON_A(:)
        INTEGER(8) :: APOS, SON_APOS, IACHK, JJ2, NFRONT8
        INTEGER :: KK, KK1, allocok, SON_LA
        TYPE(LRB_TYPE), POINTER :: CB_LRB(:,:), LRB
        INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_DYNAMIC
        INTEGER :: NB_INCB, NB_INASM, NB_BLR, I, J, M, N, II, NPIV,
     &        IBIS, IBIS_END, FIRST_ROW, LAST_ROW, FIRST_COL, LAST_COL,
     &        SON_LDA
        DOUBLE PRECISION :: PROMOTE_COST
        DOUBLE PRECISION :: ONE, ZERO
        PARAMETER (ONE = 1.0D0)
        PARAMETER (ZERO = 0.0D0)
        CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_DYN(IWHANDLER,
     &                       BEGS_BLR_DYNAMIC)
        CALL DMUMPS_BLR_RETRIEVE_CB_LRB(IWHANDLER, CB_LRB)
        NB_BLR   = size(BEGS_BLR_DYNAMIC)-1
        NB_INCB  = size(CB_LRB,1)
        NB_INASM = NB_BLR - NB_INCB
        NPIV = BEGS_BLR_DYNAMIC(NB_INASM+1)-1
        NFRONT8 = int(NFRONT,8)
        IF (SYM.EQ.0) THEN
          IBIS_END = NB_INCB*NB_INCB
        ELSE
          IBIS_END = NB_INCB*(NB_INCB+1)/2
        ENDIF
#if defined(BLR_MT)         
!$OMP PARALLEL         
!$OMP DO PRIVATE(IBIS, I, J, M, N, SON_LA, SON_LDA, FIRST_ROW,
!$OMP&         LAST_ROW, FIRST_COL, LAST_COL, LRB, SON_A, II, KK,
!$OMP&         APOS, IACHK, KK1, JJ2, PROMOTE_COST, allocok, SON_APOS)
#endif
        DO IBIS = 1,IBIS_END
C         Determining I,J from IBIS
          IF (SYM.EQ.0) THEN
            I = (IBIS-1)/NB_INCB+1
            J = IBIS - (I-1)*NB_INCB
          ELSE
            I = ceiling((1.0D0+sqrt(1.0D0+8.0D0*dble(IBIS)))/2.0D0)-1
            J = IBIS - I*(I-1)/2
          ENDIF
          I = I+NB_INASM
          J = J+NB_INASM
          IF (I.EQ.NB_INASM+1) THEN
C           first CB block, add NELIM because FIRST_ROW starts at NELIM+1
            FIRST_ROW = BEGS_BLR_DYNAMIC(I)-NPIV+NELIM
          ELSE
            FIRST_ROW = BEGS_BLR_DYNAMIC(I)-NPIV
          ENDIF
          LAST_ROW  = BEGS_BLR_DYNAMIC(I+1)-1-NPIV
          M=LAST_ROW-FIRST_ROW+1
          FIRST_COL = BEGS_BLR_DYNAMIC(J)-NPIV
          LAST_COL  = BEGS_BLR_DYNAMIC(J+1)-1-NPIV
          N = BEGS_BLR_DYNAMIC(J+1)-BEGS_BLR_DYNAMIC(J)
          SON_APOS  = 1_8
          SON_LA    = M*N
          SON_LDA   = N
          LRB => CB_LRB(I-NB_INASM,J-NB_INASM)
          IF (LRB%ISLR.AND.LRB%K.EQ.0) THEN
C           No need to perform extend-add              
            CALL DEALLOC_LRB(LRB, KEEP8)
            NULLIFY(LRB)
            CYCLE
          ENDIF
          allocate(SON_A(SON_LA),stat=allocok)
          IF (allocok.GT.0) THEN
            write(*,*) 'Not enough memory in DMUMPS_BLR_ASM_NIV1',
     &     ", Memory requested = ", SON_LA
            CALL MUMPS_ABORT()
          ENDIF
C         decompress block        
          IF (LRB%ISLR) THEN
            CALL dgemm('T', 'T', N, M, LRB%K, ONE, LRB%R(1,1), LRB%K,
     &                 LRB%Q(1,1), M, ZERO, SON_A(SON_APOS), SON_LDA)
            PROMOTE_COST = 2.0D0*M*N*LRB%K
            CALL UPD_FLOP_DECOMPRESS(PROMOTE_COST, .TRUE.)
          ELSE
            IF (I.EQ.J.AND.SYM.NE.0) THEN
C             Diag block and LDLT, copy only lower half              
             IF (J-NB_INASM.EQ.1.AND.NELIM.GT.0) THEN
C             The first diagonal block is rectangular !!
C             with NELIM more cols than rows
              DO II=1,M
                DO KK=1,II+NELIM
                  SON_A(SON_APOS+int(II-1,8)*int(SON_LDA,8) + 
     &                           int(KK-1,8))
     &                 = LRB%Q(II,KK)
                ENDDO
              ENDDO  
             ELSE
              DO II=1,M
                DO KK=1,II
                  SON_A(SON_APOS+int(II-1,8)*int(SON_LDA,8) + 
     &                           int(KK-1,8))
     &                 = LRB%Q(II,KK)
                ENDDO
              ENDDO  
             ENDIF
            ELSE
              DO II=1,M
                DO KK=1,N
                  SON_A(SON_APOS+int(II-1,8)*int(SON_LDA,8) + 
     &                           int(KK-1,8))
     &                 = LRB%Q(II,KK)
                ENDDO
              ENDDO  
            ENDIF
          ENDIF
C         Deallocate block
          CALL DEALLOC_LRB(LRB, KEEP8)
          NULLIFY(LRB)
C         extend add in father                           
          IF (SYM.NE.0.AND.J-NB_INASM.EQ.1.AND.NELIM.GT.0) THEN
C           Case of LDLT with NELIM: first-block column is treated
C           differently as the NELIM are assembled at the end of the
C           father
              DO KK = FIRST_ROW, LAST_ROW
                IACHK = 1_8 + int(KK-FIRST_ROW,8)*int(SON_LDA,8)
                IF (SON_IW(KK+K1-1).LE.NASS1) THEN
C                 Fully summed row of the father => permute destination in
C                 father, symmetric swap to be done 
C                 First NELIM columns
                  APOS = POSEL1 + int(SON_IW(KK+K1-1),8) - 1_8
                  DO KK1 = FIRST_COL, FIRST_COL+NELIM-1
                   JJ2 = APOS + int(SON_IW(K1+KK1-1)-1,8)*NFRONT8
                   A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                  ENDDO
C                 Remaining columns
                  APOS = POSEL1 + int(SON_IW(KK+K1-1)-1,8)*NFRONT8
C                 DO KK1 = FIRST_COL+NELIM, LAST_COL
C                 In case I=J and first block, one may have
C                 LAST_COL > KK, but only lower triangular part
C                 should be assembled. We use min(LAST_COL,KK)
C                 below index to cover this case.
                  DO KK1 = FIRST_COL+NELIM, min(LAST_COL,KK)
                   JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                   A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                  ENDDO
                ELSE
                  APOS = POSEL1 + int(SON_IW(KK+K1-1)-1,8)*NFRONT8
                  DO KK1 = FIRST_COL, LAST_COL
                   JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                   A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                  ENDDO
                ENDIF
              ENDDO
          ELSE
C           Case of LDLT without NELIM or LU: everything is simpler
            DO KK = FIRST_ROW, LAST_ROW
              APOS = POSEL1 + int(SON_IW(KK+K1-1)-1,8)*NFRONT8
              IACHK = 1_8 + int(KK-FIRST_ROW,8)*int(SON_LDA,8)
              IF (I.EQ.J.AND.SYM.NE.0) THEN
C               LDLT diag block: assemble only lower half
                DO KK1 = FIRST_COL, KK
                  JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                  A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                ENDDO
              ELSE
                DO KK1 = FIRST_COL, LAST_COL
                  JJ2 = APOS + int(SON_IW(K1+KK1-1),8) - 1_8
                  A(JJ2) = A(JJ2) + SON_A(IACHK + int(KK1-FIRST_COL,8))
                ENDDO
              ENDIF
            ENDDO
          ENDIF
C         Deallocate SON_A
          DEALLOCATE(SON_A)
        ENDDO
#if defined(BLR_MT)         
!$OMP END DO
!$OMP END PARALLEL         
#endif
        CALL DMUMPS_BLR_FREE_CB_LRB(IWHANDLER,
C          Only CB_LRB structure is left to deallocate
     &             .TRUE., 
     &             KEEP8)
        IF ((KEEP(486).EQ.3).OR.KEEP(486).EQ.0) THEN
C        Case of FR solve: the BLR structure could not be freed
C        in DMUMPS_END_FACTO_SLAVE and should be freed here
C        Not reachable in case of error: set INFO1 to 0
         CALL DMUMPS_BLR_END_FRONT(IWHANDLER, 0, KEEP8,
     &   MTK405=KEEP(405))
        ENDIF
      END SUBROUTINE DMUMPS_BLR_ASM_NIV1
      END MODULE DMUMPS_LR_CORE
C --------------------------------------------------------------------
      SUBROUTINE DMUMPS_TRUNCATED_RRQR( M, N, A, LDA, JPVT, TAU, WORK,
     &     LDW, RWORK, TOLEPS, TOL_OPT, RANK, MAXRANK, INFO)
C     This routine computes a Rank-Revealing QR factorization of a dense
C     matrix A. The factorization is truncated when the absolute value of
C     a diagonal coefficient of the R factor becomes smaller than a
C     prescribed threshold TOLEPS. The resulting partial Q and R factors
C     provide a rank-k approximation of the input matrix A with accuracy
C     TOLEPS.
C     
C     This routine is obtained by merging the LAPACK
C     (http://www.netlib.org/lapack/) CGEQP3 and CLAQPS routines and by
C     applying a minor modification to the outer factorization loop in
C     order to stop computations as soon as possible when the required
C     accuracy is reached.
C
C     Copyright (c) 1992-2017 The University of Tennessee and The 
C     University of Tennessee Research Foundation.  All rights reserved.
C     Copyright (c) 2000-2017 The University of California Berkeley. 
C     All rights reserved.
C     Copyright (c) 2006-2017 The University of Colorado Denver.  
C     All rights reserved.
C
C     Redistribution and use in source and binary forms, with or without
C     modification, are permitted provided that the following conditions
C     are met:
C
C      - Redistributions of source code must retain the above copyright
C        notice, this list of conditions and the following disclaimer.
C
C      - Redistributions in binary form must reproduce the above 
C        copyright notice, this list of conditions and the following 
C        disclaimer listed in this license in the documentation and/or 
C        other materials provided with the distribution.
C
C      - Neither the name of the copyright holders nor the names of its
C        contributors may be used to endorse or promote products derived from
C        this software without specific prior written permission.
C
C      The copyright holders provide no reassurances that the source code
C      provided does not infringe any patent, copyright, or any other
C      intellectual property rights of third parties.  The copyright holders
C      disclaim any liability to any recipient for claims brought against
C      recipient by any third party for infringement of that parties
C      intellectual property rights.
C
C      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C      "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C      LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C      A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C      OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C      SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C      LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C      DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C      THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
      IMPLICIT NONE
C
      INTEGER            ::  INFO, LDA, LDW, M, N, RANK, MAXRANK
C     TOL_OPT controls the tolerance option used      
C       >0 => use 2-norm (||.||_X = ||.||_2)
C       <0 => use Frobenius-norm (||.||_X = ||.||_F)
C     Furthermore, depending on abs(TOL_OPT):      
C       1 => absolute: ||B_{I(k+1:end),J(k+1:end)}||_X <= TOLEPS     
C       2 => relative to 2-norm of the compressed block: 
C        ||B_{I(k+1:end),J(k+1:end)}||_X <= TOLEPS*||B_{I,J}||_2
C       3 => relative to the max of the 2-norms of the row and column diagonal blocks 
C        ||B_{I(k+1:end),J{k+1:end}}||_X <= TOLEPS*max(||B_{I,I}||_2,||B_{J,J}||_2)
C       4 => relative to the sqrt of product of the 2-norms of the row and column diagonal blocks 
C        ||B_{I(k+1:end),J{k+1:end}}||_X <= TOLEPS*sqrt(||B_{I,I}||_2*||B_{J,J}||_2)
      INTEGER            ::  TOL_OPT
      DOUBLE PRECISION               ::  TOLEPS
      INTEGER            ::  JPVT(*)
      DOUBLE PRECISION               ::  RWORK(*)
      DOUBLE PRECISION            ::  A(LDA,*), TAU(*)
      DOUBLE PRECISION            ::  WORK(LDW,*)
      DOUBLE PRECISION               ::  TOLEPS_EFF, TRUNC_ERR
      INTEGER, PARAMETER ::  INB=1, INBMIN=2
      INTEGER            :: J, JB, MINMN, NB
      INTEGER            :: OFFSET, ITEMP
      INTEGER            :: LSTICC, PVT, K, RK
      DOUBLE PRECISION               :: TEMP, TEMP2, TOL3Z
      DOUBLE PRECISION            :: AKK
      DOUBLE PRECISION, PARAMETER    :: RZERO=0.0D+0, RONE=1.0D+0
      DOUBLE PRECISION :: ZERO
      DOUBLE PRECISION :: ONE
      PARAMETER          ( ONE = 1.0D+0 )
      PARAMETER          ( ZERO = 0.0D+0 ) 
      DOUBLE PRECISION               :: dlamch
      INTEGER            :: ilaenv, idamax
      EXTERNAL           :: idamax, dlamch
      EXTERNAL           dgeqrf, dormqr, xerbla
      EXTERNAL           ilaenv
      EXTERNAL           dgemm, dgemv, dlarfg, dswap
      DOUBLE PRECISION, EXTERNAL :: dnrm2
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDW.LT.N ) THEN
            INFO = -8
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         WRITE(*,999) -INFO
         RETURN
      END IF
      MINMN = MIN(M,N)
      IF( MINMN.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
      NB = ilaenv( INB, 'CGEQRF', ' ', M, N, -1, -1 )
      SELECT CASE(abs(TOL_OPT))
      CASE(1)
        TOLEPS_EFF = TOLEPS
      CASE(2)
C      TOLEPS_EFF will be computed at step K=1 below        
      CASE DEFAULT
        write(*,*) 'Internal error in DMUMPS_TRUNCATED_RRQR: TOL_OPT =',
     &        TOL_OPT
        CALL MUMPS_ABORT()
      END SELECT
      TOLEPS_EFF = TOLEPS
C
C     Avoid pointers (and TARGET attribute on RWORK/WORK)
C     because of implicit interface. An implicit interface
C     is needed to avoid intermediate array copies
C     VN1  => RWORK(1:N)
C     VN2  => RWORK(N+1:2*N)
C     AUXV => WORK(1:LDW,1:1)
C     F    => WORK(1:LDW,2:NB+1)
C     LDF  =  LDW
*     Initialize partial column norms. The first N elements of work
*     store the exact column norms.
      DO J = 1, N
C        VN1( J ) = dnrm2( M, A( 1, J ), 1 )
         RWORK( J ) = dnrm2( M, A( 1, J ), 1 )
C        VN2( J ) = VN1( J )
         RWORK( N + J ) = RWORK( J )
         JPVT(J) = J
      END DO
      IF (TOL_OPT.LT.0) THEN
C       Compute TRUNC_ERR for first step              
C       TRUNC_ERR = dnrm2( N, VN1( 1 ), 1 )
        TRUNC_ERR = dnrm2( N, RWORK( 1 ), 1 )
      ENDIF
      OFFSET = 0
      TOL3Z  = SQRT(dlamch('Epsilon'))
      DO 
         JB     = MIN(NB,MINMN-OFFSET)
         LSTICC = 0
         K      = 0
         DO 
            IF(K.EQ.JB) EXIT
            K   = K+1
            RK  = OFFSET+K
C           PVT = ( RK-1 ) + IDAMAX( N-RK+1, VN1( RK ), 1 )
            PVT = ( RK-1 ) + idamax( N-RK+1, RWORK( RK ), 1 )
            IF (RK.EQ.1) THEN 
C             IF (abs(TOL_OPT).EQ.2) TOLEPS_EFF = VN1(PVT)*TOLEPS
              IF (abs(TOL_OPT).EQ.2) TOLEPS_EFF = RWORK(PVT)*TOLEPS
            ENDIF
            IF (TOL_OPT.GT.0) THEN
C             TRUNC_ERR = VN1(PVT)
              TRUNC_ERR = RWORK(PVT)
C           ELSE
C             TRUNC_ERR has been already computed at previous step
            ENDIF
            IF(TRUNC_ERR.LT.TOLEPS_EFF) THEN
               RANK = RK-1
               RETURN
            END IF
            IF (RK.GT.MAXRANK) THEN
               RANK = RK
               INFO = RK
               RETURN
            END IF
            IF( PVT.NE.RK ) THEN
               CALL dswap( M, A( 1, PVT ), 1, A( 1, RK ), 1 )
c              CALL dswap( K-1, F( PVT-OFFSET, 1 ), LDF,
c    &              F( K, 1 ), LDF )
               CALL dswap( K-1, WORK( PVT-OFFSET, 2 ), LDW,
     &              WORK( K, 2 ), LDW )
               ITEMP     = JPVT(PVT)
               JPVT(PVT) = JPVT(RK)
               JPVT(RK)  = ITEMP
C              VN1(PVT)  = VN1(RK)
C              VN2(PVT)  = VN2(RK)
               RWORK(PVT)    = RWORK(RK)
               RWORK(N+PVT)  = RWORK(N+RK)
            END IF
*     Apply previous Householder reflectors to column K:
*     A(RK:M,RK) := A(RK:M,RK) - A(RK:M,OFFSET+1:RK-1)*F(K,1:K-1)**H.
            IF( K.GT.1 ) THEN
               CALL dgemv( 'No transpose', M-RK+1, K-1, -ONE,
C    &              A(RK,OFFSET+1), LDA, F(K,1), LDF,
     &              A(RK,OFFSET+1), LDA, WORK(K,2), LDW,
     &              ONE, A(RK,RK), 1 )
            END IF
*     Generate elementary reflector H(k).
            IF( RK.LT.M ) THEN
               CALL dlarfg( M-RK+1, A(RK,RK), A(RK+1,RK), 1, TAU(RK) )
            ELSE
               CALL dlarfg( 1, A(RK,RK), A(RK,RK), 1, TAU(RK) )
            END IF
            AKK      = A(RK,RK)
            A(RK,RK) = ONE
*     Compute Kth column of F:
*     F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).
            IF( RK.LT.N ) THEN
               CALL dgemv( 'Transpose', M-RK+1, N-RK, TAU(RK),
     &              A(RK,RK+1), LDA, A(RK,RK), 1, ZERO,
C    &              F( K+1, K ), 1 )
     &              WORK( K+1, K+1 ), 1 )
            END IF
*     Padding F(1:K,K) with zeros.
            DO J = 1, K
C              F( J, K ) = ZERO
               WORK( J, K+1 ) = ZERO
            END DO
*     Incremental updating of F:
*     F(1:N,K) := F(1:N-OFFSET,K) - 
*             tau(RK)*F(1:N,1:K-1)*A(RK:M,OFFSET+1:RK-1)**H*A(RK:M,RK).
            IF( K.GT.1 ) THEN
               CALL dgemv( 'Transpose', M-RK+1, K-1, -TAU(RK),
     &              A(RK,OFFSET+1), LDA, A(RK,RK), 1, ZERO,
     &              WORK(1,1), 1 )
C    &              AUXV(1,1), 1 )
               CALL dgemv( 'No transpose', N-OFFSET, K-1, ONE,
     &              WORK(1,2), LDW, WORK(1,1), 1, ONE, WORK(1,K+1), 1 )
C    &              F(1,1), LDF, AUXV(1,1), 1, ONE, F(1,K), 1 )
            END IF
*     Update the current row of A:
*     A(RK,RK+1:N) := A(RK,RK+1:N) - A(RK,OFFSET+1:RK)*F(K+1:N,1:K)**H.
            IF( RK.LT.N ) THEN
C              CALL dgemv( 'No Transpose', N-RK, K, -ONE, F( K+1, 1 ), 
               CALL dgemv( 'No Transpose', N-RK, K, -ONE, WORK( K+1,2 ),
     &              LDW,
     &              A( RK, OFFSET+1 ), LDA, ONE, A( RK, RK+1 ), LDA )
            END IF
*     Update partial column norms.
*     
            IF( RK.LT.MINMN ) THEN
               DO J = RK + 1, N
C                 IF( VN1( J ).NE.RZERO ) THEN
                  IF( RWORK( J ).NE.RZERO ) THEN
*     
*     NOTE: The following 4 lines follow from the analysis in
*     Lapack Working Note 176.
*
C                    TEMP = ABS( A( RK, J ) ) / VN1( J )
                     TEMP = ABS( A( RK, J ) ) / RWORK( J )
                     TEMP = MAX( RZERO, ( RONE+TEMP )*( RONE-TEMP ) )
C                    TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                     TEMP2 = TEMP*( RWORK( J ) / RWORK( N+J ) )**2
                     IF( TEMP2 .LE. TOL3Z ) THEN
C                       VN2( J ) = dble( LSTICC )
                        RWORK( N+J ) = dble( LSTICC )
                        LSTICC = J
                     ELSE
C                       VN1( J ) = VN1( J )*SQRT( TEMP )
                        RWORK( J ) = RWORK( J )*SQRT( TEMP )
                     END IF
                  END IF
               END DO
            END IF
            A( RK, RK ) = AKK
            IF (LSTICC.NE.0) EXIT
            IF (TOL_OPT.LT.0) THEN
C             Compute TRUNC_ERR for next step              
C             TRUNC_ERR = dnrm2( N-RK, VN1( RK+1 ), 1 )
              TRUNC_ERR = dnrm2( N-RK, RWORK( RK+1 ), 1 )
            ENDIF
         END DO
*     Apply the block reflector to the rest of the matrix:
*     A(RK+1:M,RK+1:N) := A(RK+1:M,RK+1:N) -
*     A(RK+1:M,OFFSET+1:RK)*F(K+1:N-OFFSET,1:K)**H.
         IF( RK.LT.MIN(N,M) ) THEN
            CALL dgemm( 'No transpose', 'Transpose', M-RK,
     &           N-RK, K, -ONE, A(RK+1,OFFSET+1), LDA,
C    &           F(K+1,1), LDF, ONE, A(RK+1,RK+1), LDA )
     &           WORK(K+1,2), LDW, ONE, A(RK+1,RK+1), LDA )
         END IF
*     Recomputation of difficult columns.
         DO WHILE( LSTICC.GT.0 ) 
C           ITEMP = NINT( VN2( LSTICC ) )
            ITEMP = NINT( RWORK( N + LSTICC ) )
C           VN1( LSTICC ) = dnrm2( M-RK, A( RK+1, LSTICC ), 1 )
            RWORK( LSTICC ) = dnrm2( M-RK, A( RK+1, LSTICC ), 1 )
*     
*     NOTE: The computation of RWORK( LSTICC ) relies on the fact that 
*     SNRM2 does not fail on vectors with norm below the value of
*     SQRT(DLAMCH('S')) 
*     
C           VN2( LSTICC ) = VN1( LSTICC )
            RWORK( N + LSTICC ) = RWORK( LSTICC )
            LSTICC = ITEMP
         END DO
         IF(RK.GE.MINMN) EXIT
         OFFSET = RK
         IF (TOL_OPT.LT.0) THEN
C          Compute TRUNC_ERR for next step              
C          TRUNC_ERR = dnrm2( N-RK, VN1( RK+1 ), 1 )
           TRUNC_ERR = dnrm2( N-RK, RWORK( RK+1 ), 1 )
         ENDIF
      END DO
      RANK = RK
      RETURN
 999  FORMAT ('On entry to DMUMPS_TRUNCATED_RRQR, parameter number',
     &            I2,' had an illegal value')
      END SUBROUTINE DMUMPS_TRUNCATED_RRQR
