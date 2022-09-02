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
      SUBROUTINE DMUMPS_FREETOPSO( N, KEEP28, IWCB, LIWW,
     &       W, LWC,
     &       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER(8), INTENT(IN) :: LWC
      INTEGER(8), INTENT(INOUT) :: POSWCB
      INTEGER N,LIWW,IWPOSCB, KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28)
      INTEGER(8) :: PTRACB(KEEP28)
      DOUBLE PRECISION W(LWC)
      INTEGER SIZFI, SIZFR
      IF ( IWPOSCB .eq. LIWW ) RETURN
      DO WHILE ( IWCB( IWPOSCB + 2 ) .eq. 0 )
        SIZFR = IWCB( IWPOSCB + 1 )
        SIZFI =  2  
        IWPOSCB = IWPOSCB + SIZFI
        POSWCB  = POSWCB  + SIZFR
        IF ( IWPOSCB .eq. LIWW ) RETURN
      END DO
      RETURN
      END SUBROUTINE DMUMPS_FREETOPSO
      SUBROUTINE DMUMPS_COMPSO(N,KEEP28,IWCB,LIWW,W,LWC,
     &       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER(8), INTENT(IN)    :: LWC
      INTEGER(8), INTENT(INOUT) :: POSWCB
      INTEGER N,LIWW,IWPOSCB,KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28)
      INTEGER(8) :: PTRACB(KEEP28)
      DOUBLE PRECISION W(LWC)
      INTEGER IPTIW,SIZFI,LONGI
      INTEGER(8) :: IPTA, LONGR, SIZFR, I8
      INTEGER    :: I
      IPTIW = IWPOSCB
      IPTA  = POSWCB
      LONGI = 0
      LONGR = 0_8
      IF ( IPTIW .EQ. LIWW ) RETURN
10    CONTINUE
       IF (IWCB(IPTIW+2).EQ.0) THEN
        SIZFR  = int(IWCB(IPTIW+1),8)
        SIZFI =  2  
        IF (LONGI.NE.0) THEN
          DO 20 I=0,LONGI-1
            IWCB(IPTIW + SIZFI - I) = IWCB (IPTIW - I)
 20       CONTINUE 
          DO 30 I8=0,LONGR-1
            W(IPTA + SIZFR - I8)   = W(IPTA - I8)
 30       CONTINUE
        ENDIF
        DO 40 I=1,KEEP28
          IF ((PTRICB(I).LE.(IPTIW+1)).AND.
     &        (PTRICB(I).GT.IWPOSCB) ) THEN
            PTRICB(I) = PTRICB(I) + SIZFI
            PTRACB(I) = PTRACB(I) + SIZFR
          ENDIF 
40      CONTINUE 
        IWPOSCB = IWPOSCB + SIZFI
        IPTIW   = IPTIW + SIZFI
        POSWCB = POSWCB + SIZFR
        IPTA   = IPTA + SIZFR     
       ELSE
        SIZFR  = int(IWCB(IPTIW+1),8)
        SIZFI  = 2
        IPTIW = IPTIW + SIZFI
        LONGI = LONGI + SIZFI
        IPTA  = IPTA + SIZFR
        LONGR = LONGR + SIZFR
       ENDIF
       IF (IPTIW.NE.LIWW) GOTO 10
       RETURN
       END SUBROUTINE DMUMPS_COMPSO
      SUBROUTINE DMUMPS_SOL_X(A, NZ8, N, IRN, ICN, Z, KEEP,KEEP8)
      INTEGER N, I, J, KEEP(500)
      INTEGER(8), INTENT(IN) :: NZ8
      INTEGER(8) KEEP8(150)
      INTEGER IRN(NZ8), ICN(NZ8)
      DOUBLE PRECISION A(NZ8)
      DOUBLE PRECISION Z(N)
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      INTEGER(8) :: K
      INTRINSIC     abs
      DO 10 I = 1, N
        Z(I) = ZERO
   10 CONTINUE
      IF (KEEP(264).EQ.0) THEN
       IF (KEEP(50) .EQ.0) THEN
         DO K = 1_8, NZ8
          I = IRN(K)
          J = ICN(K)
          IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
          IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
          Z(I) = Z(I) + abs(A(K))
         ENDDO
        ELSE
         DO K = 1_8, NZ8
          I = IRN(K)
          J = ICN(K)
          IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
          IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
          Z(I) = Z(I) + abs(A(K))
          IF (J.NE.I) THEN 
            Z(J) = Z(J) + abs(A(K))
          ENDIF
         ENDDO
        ENDIF
      ELSE
       IF (KEEP(50) .EQ.0) THEN
         DO K = 1_8, NZ8
          I = IRN(K)
          J = ICN(K)
          Z(I) = Z(I) + abs(A(K))
         ENDDO
        ELSE
         DO K = 1_8, NZ8
          I = IRN(K)
          J = ICN(K)
          Z(I) = Z(I) + abs(A(K))
          IF (J.NE.I) THEN 
            Z(J) = Z(J) + abs(A(K))
          ENDIF
         ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SOL_X
      SUBROUTINE DMUMPS_SCAL_X(A, NZ8, N, IRN, ICN, Z,
     &            KEEP, KEEP8, COLSCA)
      INTEGER,    INTENT(IN)  :: N, KEEP(500)
      INTEGER(8), INTENT(IN)  :: NZ8
      INTEGER(8), INTENT(IN)  :: KEEP8(150)
      INTEGER,    INTENT(IN)  :: IRN(NZ8), ICN(NZ8)
      DOUBLE PRECISION,    INTENT(IN)  :: A(NZ8)
      DOUBLE PRECISION,       INTENT(IN)  :: COLSCA(N)
      DOUBLE PRECISION,       INTENT(OUT) :: Z(N)
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      INTEGER         :: I, J
      INTEGER(8)      :: K
      DO 10 I = 1, N
        Z(I) = ZERO
   10 CONTINUE
      IF (KEEP(50) .EQ.0) THEN
       DO K = 1_8, NZ8
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + abs(A(K)*COLSCA(J))
       ENDDO
      ELSE
       DO K = 1, NZ8
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + abs(A(K)*COLSCA(J))
        IF (J.NE.I) THEN
          Z(J) = Z(J) + abs(A(K)*COLSCA(I))
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SCAL_X
      SUBROUTINE DMUMPS_SOL_Y(A, NZ8, N, IRN, ICN, RHS, X, R, W,
     &           KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER,    INTENT(IN)   :: N, KEEP(500)
      INTEGER(8), INTENT(IN)   :: NZ8
      INTEGER(8), INTENT(IN)   :: KEEP8(150)
      INTEGER,    INTENT(IN)   :: IRN(NZ8), ICN(NZ8)
      DOUBLE PRECISION,    INTENT(IN)   :: A(NZ8), RHS(N), X(N)
      DOUBLE PRECISION,       INTENT(OUT)  :: W(N)
      DOUBLE PRECISION,    INTENT(OUT)  :: R(N)
      INTEGER I, J
      INTEGER(8) :: K8
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      DOUBLE PRECISION D
      DO I = 1, N
        R(I) = RHS(I)
        W(I) = ZERO
      ENDDO
      IF (KEEP(264).EQ.0) THEN
       IF (KEEP(50) .EQ.0) THEN
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            IF ((I .GT. N) .OR. (J .GT. N) .OR. (I .LT. 1) .OR. 
     &       (J .LT. 1)) CYCLE
            D = A(K8) * X(J)
            R(I) = R(I) - D
            W(I) = W(I) + abs(D)
          ENDDO
       ELSE
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            IF ((I .GT. N) .OR. (J .GT. N) .OR. (I .LT. 1) .OR. 
     &       (J .LT. 1)) CYCLE
            D = A(K8) * X(J)
            R(I) = R(I) - D
            W(I) = W(I) + abs(D)
            IF (I.NE.J) THEN
              D = A(K8) * X(I)
              R(J) = R(J) - D
              W(J) = W(J) + abs(D)
            ENDIF
          ENDDO
       ENDIF
      ELSE
       IF (KEEP(50) .EQ.0) THEN
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            D = A(K8) * X(J)
            R(I) = R(I) - D
            W(I) = W(I) + abs(D)
          ENDDO
       ELSE
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            D = A(K8) * X(J)
            R(I) = R(I) - D
            W(I) = W(I) + abs(D)
            IF (I.NE.J) THEN
              D = A(K8) * X(I)
              R(J) = R(J) - D
              W(J) = W(J) + abs(D)
            ENDIF
          ENDDO
       ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SOL_Y
      SUBROUTINE DMUMPS_SOL_MULR(N, R, W)
      INTEGER, intent(in)  :: N
      DOUBLE PRECISION,    intent(in)  :: W(N)
      DOUBLE PRECISION, intent(inout) :: R(N)
      INTEGER I
      DO 10 I = 1, N
        R(I) = R(I) * W(I)
   10 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_SOL_MULR
      SUBROUTINE DMUMPS_SOL_B(N, KASE, X, EST, W, IW, GRAIN)
      INTEGER, intent(in)    :: N
      INTEGER, intent(inout) :: KASE
      INTEGER IW(N)
      DOUBLE PRECISION W(N), X(N)
      DOUBLE PRECISION, intent(inout)    :: EST
      INTEGER, intent(in)    :: GRAIN
      INTRINSIC abs, nint, real, sign
      INTEGER DMUMPS_IXAMAX
      EXTERNAL DMUMPS_IXAMAX
      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      INTEGER I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION ALTSGN
      DOUBLE PRECISION TEMP
      SAVE ITER, J, JLAST, JUMP
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO = 0.0D0 )
      PARAMETER( ONE = 1.0D0 )
      DOUBLE PRECISION, PARAMETER :: RZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER :: RONE = 1.0D0
      IF (KASE .EQ. 0) THEN
        DO 10 I = 1, N
          X(I) = ONE / dble(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        RETURN
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 20
      CASE(2)
        GOTO 40
      CASE(3)
        GOTO 70
      CASE(4)
        GOTO 120
      CASE(5)
        GOTO 160
      CASE DEFAULT
      END SELECT
   20 CONTINUE
      IF (N .EQ. 1) THEN
        W(1) = X(1)
        EST = abs(W(1))
        GOTO 190
      ENDIF
      DO 30 I = 1, N
        X(I)  = sign( RONE,dble(X(I)) )
        IW(I) = nint(dble(X(I)))
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
   40 CONTINUE
      J = DMUMPS_IXAMAX(N, X, 1, GRAIN)
      ITER = 2
   50 CONTINUE
      DO 60 I = 1, N
        X(I) = ZERO
   60 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
   70 CONTINUE
      DO 80 I = 1, N
        W(I) = X(I)
   80 CONTINUE
      DO 90 I = 1, N
        IF (nint(sign(RONE, dble(X(I)))) .NE. IW(I)) GOTO 100
   90 CONTINUE
      GOTO 130
  100 CONTINUE
      DO 110 I = 1, N
        X(I) = sign(RONE, dble(X(I)))
        IW(I) = nint(dble(X(I)))
  110 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
  120 CONTINUE
      JLAST = J
      J = DMUMPS_IXAMAX(N, X, 1, GRAIN)
      IF ((abs(X(JLAST)) .NE. abs(X(J))) .AND. (ITER .LT. ITMAX)) THEN
        ITER = ITER + 1
        GOTO 50
      ENDIF
  130 CONTINUE
      EST = RZERO
      DO 140 I = 1, N
        EST = EST + abs(W(I))
  140 CONTINUE
      ALTSGN = RONE
      DO 150 I = 1, N
        X(I) = ALTSGN * (RONE + dble(I - 1) / dble(N - 1))
        ALTSGN = -ALTSGN
  150 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
  160 CONTINUE
      TEMP = RZERO
      DO 170 I = 1, N
        TEMP = TEMP + abs(X(I))
  170 CONTINUE
      TEMP = 2.0D0 * TEMP / dble(3 * N)
      IF (TEMP .GT. EST) THEN
        DO 180 I = 1, N
          W(I) = X(I)
  180   CONTINUE
        EST = TEMP
      ENDIF
  190 KASE = 0
      RETURN
      END SUBROUTINE DMUMPS_SOL_B
      SUBROUTINE DMUMPS_QD2( MTYPE, N, NZ8, ASPK, IRN, ICN,
     &    LHS, WRHS, W, RHS, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER MTYPE, N
      INTEGER(8), INTENT(IN) :: NZ8
      INTEGER, INTENT(IN) :: IRN( NZ8 ), ICN( NZ8 )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(IN) :: ASPK( NZ8 )
      DOUBLE PRECISION, INTENT(IN) :: LHS( N ), WRHS( N )
      DOUBLE PRECISION, INTENT(OUT):: RHS( N )
      DOUBLE PRECISION,    INTENT(OUT):: W( N )
      INTEGER I, J
      INTEGER(8) :: K8
      DOUBLE PRECISION, PARAMETER :: DZERO = 0.0D0
      DO I = 1, N
        W(I) = DZERO
        RHS(I) = WRHS(I)
      ENDDO
      IF ( KEEP(50) .EQ. 0 ) THEN
       IF (MTYPE .EQ. 1) THEN
        IF (KEEP(264).EQ.0) THEN
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. 
     &        (J .GT. N)) CYCLE
            RHS(I) = RHS(I) - ASPK(K8) * LHS(J)
            W(I) = W(I) + abs(ASPK(K8))
          ENDDO
        ELSE
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            RHS(I) = RHS(I) - ASPK(K8) * LHS(J)
            W(I) = W(I) + abs(ASPK(K8))
          ENDDO
        ENDIF
       ELSE
        IF (KEEP(264).EQ.0) THEN
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. 
     &        (J .GT. N)) CYCLE
            RHS(J) = RHS(J) - ASPK(K8) * LHS(I)
            W(J) = W(J) + abs(ASPK(K8))
          ENDDO
        ELSE
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            RHS(J) = RHS(J) - ASPK(K8) * LHS(I)
            W(J) = W(J) + abs(ASPK(K8))
          ENDDO
        ENDIF
       ENDIF
      ELSE
        IF (KEEP(264).EQ.0) THEN
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. 
     &        (J .GT. N)) CYCLE
            RHS(I) = RHS(I) - ASPK(K8) * LHS(J)
            W(I) = W(I) + abs(ASPK(K8))
            IF (J.NE.I) THEN
                RHS(J) = RHS(J) - ASPK(K8) * LHS(I)
                W(J) = W(J) + abs(ASPK(K8))
            ENDIF
          ENDDO
        ELSE
          DO K8 = 1_8, NZ8
            I = IRN(K8)
            J = ICN(K8)
            RHS(I) = RHS(I) - ASPK(K8) * LHS(J)
            W(I) = W(I) + abs(ASPK(K8))
            IF (J.NE.I) THEN
                RHS(J) = RHS(J) - ASPK(K8) * LHS(I)
                W(J) = W(J) + abs(ASPK(K8))
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_QD2
      SUBROUTINE DMUMPS_ELTQD2( MTYPE, N,
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT8, A_ELT,
     &    LHS, WRHS, W, RHS, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR
      INTEGER(8), INTENT(IN) :: NA_ELT8
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION A_ELT(NA_ELT8)
      DOUBLE PRECISION LHS( N ), WRHS( N ), RHS( N )
      DOUBLE PRECISION W(N)
      CALL DMUMPS_MV_ELT(N, NELT, ELTPTR, ELTVAR, A_ELT,
     &                         LHS, RHS, KEEP(50), MTYPE )
      RHS = WRHS - RHS
      CALL DMUMPS_SOL_X_ELT( MTYPE, N, 
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT8, A_ELT,
     &    W, KEEP,KEEP8 )
      RETURN
      END SUBROUTINE DMUMPS_ELTQD2
      SUBROUTINE DMUMPS_SOL_X_ELT( MTYPE, N, 
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT8, A_ELT,
     &    W, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR
      INTEGER(8), INTENT(IN) :: NA_ELT8
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION A_ELT(NA_ELT8)
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION W(N)
      INTEGER I, J, IEL, SIZEI, IELPTR
      INTEGER(8) :: K8
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO = 0.0D0)
      W = DZERO
      K8 = 1_8
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) = 
     &           W( ELTVAR( IELPTR + I) )
     &           + abs(A_ELT( K8 ))
               K8 = K8 + 1_8
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
               TEMP = TEMP + abs( A_ELT(K8))
               K8 = K8 + 1_8
              END DO
              W(ELTVAR( IELPTR + J )) = 
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) = 
     &        W(ELTVAR( IELPTR + J )) + abs(A_ELT( K8 ))
          K8 = K8 + 1_8
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) = 
     &           W(ELTVAR( IELPTR + J )) + abs(A_ELT( K8 ))
              W(ELTVAR( IELPTR + I ) ) = 
     &           W(ELTVAR( IELPTR + I )) + abs(A_ELT( K8 ))
              K8 = K8 + 1_8
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_SOL_X_ELT
      SUBROUTINE DMUMPS_SOL_SCALX_ELT(MTYPE, N, 
     &    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT8, A_ELT,
     &    W, KEEP,KEEP8, COLSCA )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR
      INTEGER(8), INTENT(IN) :: NA_ELT8
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION COLSCA(N)
      DOUBLE PRECISION A_ELT(NA_ELT8)
      DOUBLE PRECISION W(N)
      DOUBLE PRECISION TEMP, TEMP2
      INTEGER I, J, IEL, SIZEI, IELPTR
      INTEGER(8) :: K8
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO = 0.0D0)
      W = DZERO
      K8 = 1_8
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              TEMP2 = abs(COLSCA(ELTVAR( IELPTR + J) ))
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) =
     &           W( ELTVAR( IELPTR + I) )
     &           + abs(A_ELT( K8 )) * TEMP2
               K8 = K8 + 1_8
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              TEMP2= abs(COLSCA(ELTVAR( IELPTR + J) ))
              DO I = 1, SIZEI
               TEMP = TEMP + abs(A_ELT( K8 )) * TEMP2
               K8 = K8 + 1_8
              END DO
              W(ELTVAR( IELPTR + J )) =
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) =
     &        W(ELTVAR( IELPTR + J )) + 
     &        abs( A_ELT( K8 )*COLSCA(ELTVAR( IELPTR + J)) )
          K8 = K8 + 1_8
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) =
     &           W(ELTVAR( IELPTR + J )) + 
     &           abs(A_ELT( K8 )*COLSCA(ELTVAR( IELPTR + J)))
              W(ELTVAR( IELPTR + I ) ) =
     &           W(ELTVAR( IELPTR + I )) + 
     &           abs(A_ELT( K8 )*COLSCA(ELTVAR( IELPTR + I)))
              K8 = K8 + 1_8
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_SOL_SCALX_ELT
      SUBROUTINE DMUMPS_ELTYD( MTYPE, N, NELT, ELTPTR, 
     &                     LELTVAR, ELTVAR, NA_ELT8, A_ELT,
     &                     SAVERHS, X, Y, W, K50 )
      IMPLICIT NONE
      INTEGER N, NELT, K50, MTYPE, LELTVAR
      INTEGER(8) :: NA_ELT8
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( LELTVAR )
      DOUBLE PRECISION A_ELT( NA_ELT8 ), X( N ), Y( N ), 
     &                 SAVERHS(N)
      DOUBLE PRECISION W(N)
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION TEMP2
      PARAMETER( ZERO = 0.0D0 )
      Y = SAVERHS
      W = ZERO
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( K50 .eq. 0 ) THEN
          IF ( MTYPE .eq. 1 ) THEN
            DO J = 1, SIZEI
              TEMP = X( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                Y( ELTVAR( IELPTR + I ) ) =
     &          Y( ELTVAR( IELPTR + I ) ) -
     &             A_ELT( K ) * TEMP
                W( ELTVAR( IELPTR + I ) ) =
     &          W( ELTVAR( IELPTR + I ) ) +
     &             abs( A_ELT( K ) * TEMP )
                K = K + 1
              END DO
            END DO
          ELSE
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              TEMP2 = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP - 
     &          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                TEMP2 = TEMP2 +  abs(
     &          A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
              W( ELTVAR( IELPTR + J ) ) = TEMP2
            END DO
          END IF
        ELSE
          DO J = 1, SIZEI
            Y( ELTVAR( IELPTR + J ) ) =
     &      Y( ELTVAR( IELPTR + J ) ) -
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            W( ELTVAR( IELPTR + J ) ) =
     &      W( ELTVAR( IELPTR + J ) ) + abs(
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
            K = K + 1
            DO I = J+1, SIZEI
              Y( ELTVAR( IELPTR + I ) ) =
     &        Y( ELTVAR( IELPTR + I ) ) -
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     &        Y( ELTVAR( IELPTR + J ) ) -
     &           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              W( ELTVAR( IELPTR + I ) ) =
     &        W( ELTVAR( IELPTR + I ) ) + abs(
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
              W( ELTVAR( IELPTR + J ) ) =
     &        W( ELTVAR( IELPTR + J ) ) + abs(
     &           A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_ELTYD
      SUBROUTINE DMUMPS_SOLVE_GET_OOC_NODE(
     &     INODE,PTRFAC,KEEP,A,LA,STEP,
     &     KEEP8,N,MUST_BE_PERMUTED,IERR)
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER INODE,KEEP(500),N
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER IERR
      DOUBLE PRECISION A(LA)      
      INTEGER RETURN_VALUE
      LOGICAL MUST_BE_PERMUTED
      RETURN_VALUE=DMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,PTRFAC,
     &     KEEP(28),A,LA,IERR)
      IF(RETURN_VALUE.EQ.OOC_NODE_NOT_IN_MEM)THEN
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
         CALL DMUMPS_SOLVE_ALLOC_FACTOR_SPACE(INODE,PTRFAC,
     &        KEEP,KEEP8,A,IERR)
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
         CALL DMUMPS_READ_OOC(
     &        A(PTRFAC(STEP(INODE))),
     &        INODE,IERR
     &        )
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ELSE
         IF(IERR.LT.0)THEN
            RETURN
         ENDIF
      ENDIF
      IF(RETURN_VALUE.NE.OOC_NODE_PERMUTED)THEN
         MUST_BE_PERMUTED=.TRUE.
         CALL DMUMPS_SOLVE_MODIFY_STATE_NODE(INODE)
      ELSE
         MUST_BE_PERMUTED=.FALSE.
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_GET_OOC_NODE
      SUBROUTINE DMUMPS_BUILD_MAPPING_INFO(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCAL_LIST
      INTEGER :: I,IERR,TMP,NSTEPS,N_LOCAL_LIST
      INTEGER :: MASTER,TAG_SIZE,TAG_LIST
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      LOGICAL :: I_AM_SLAVE
      PARAMETER(MASTER=0, TAG_SIZE=85,TAG_LIST=86)
      I_AM_SLAVE = (id%MYID .NE. MASTER
     &     .OR. ((id%MYID.EQ.MASTER).AND.(id%KEEP(46).EQ.1)))
      NSTEPS = id%KEEP(28)
      ALLOCATE(LOCAL_LIST(NSTEPS),STAT=IERR)
      IF(IERR.GT.0) THEN
         WRITE(*,*)'Problem in solve: error allocating LOCAL_LIST'
         CALL MUMPS_ABORT()
      END IF
      N_LOCAL_LIST = 0
      IF(I_AM_SLAVE) THEN
         DO I=1,NSTEPS
            IF(id%PTLUST_S(I).NE.0) THEN
               N_LOCAL_LIST = N_LOCAL_LIST + 1
               LOCAL_LIST(N_LOCAL_LIST) = I
            END IF
         END DO
         IF(id%MYID.NE.MASTER) THEN 
            CALL MPI_SEND(N_LOCAL_LIST, 1,
     &           MPI_INTEGER, MASTER, TAG_SIZE, id%COMM,IERR)
            CALL MPI_SEND(LOCAL_LIST, N_LOCAL_LIST,
     &           MPI_INTEGER, MASTER, TAG_LIST, id%COMM,IERR)
            DEALLOCATE(LOCAL_LIST)
            ALLOCATE(id%IPTR_WORKING(1),
     &           id%WORKING(1),
     &           STAT=IERR)
            IF(IERR.GT.0) THEN
               WRITE(*,*)'Problem in solve: error allocating ',
     &              'IPTR_WORKING and WORKING'
               CALL MUMPS_ABORT()
            END IF
         END IF
      END IF
      IF(id%MYID.EQ.MASTER) THEN
         ALLOCATE(id%IPTR_WORKING(id%NPROCS+1), STAT=IERR)
         IF(IERR.GT.0) THEN
            WRITE(*,*)'Problem in solve: error allocating IPTR_WORKING'
            CALL MUMPS_ABORT()
         END IF
         id%IPTR_WORKING = 0
         id%IPTR_WORKING(1) = 1
         id%IPTR_WORKING(MASTER+2) = N_LOCAL_LIST
         DO I=1, id%NPROCS-1
            CALL MPI_RECV(TMP, 1, MPI_INTEGER, MPI_ANY_SOURCE,
     &           TAG_SIZE, id%COMM, STATUS, IERR)
            id%IPTR_WORKING(STATUS(MPI_SOURCE)+2) = TMP
         END DO
         DO I=2, id%NPROCS+1
            id%IPTR_WORKING(I) = id%IPTR_WORKING(I)
     &           + id%IPTR_WORKING(I-1)
         END DO
         ALLOCATE(id%WORKING(id%IPTR_WORKING(id%NPROCS+1)-1),STAT=IERR)
         IF(IERR.GT.0) THEN
            WRITE(*,*)'Problem in solve: error allocating LOCAL_LIST'
            CALL MUMPS_ABORT()
         END IF
         TMP = MASTER + 1
         IF (I_AM_SLAVE) THEN
            id%WORKING(id%IPTR_WORKING(TMP):id%IPTR_WORKING(TMP+1)-1)
     &           = LOCAL_LIST(1:id%IPTR_WORKING(TMP+1)
     &           -id%IPTR_WORKING(TMP))
         ENDIF
         DO I=1,id%NPROCS-1
            CALL MPI_RECV(LOCAL_LIST, NSTEPS, MPI_INTEGER,
     &           MPI_ANY_SOURCE, TAG_LIST, id%COMM, STATUS, IERR)
            TMP = STATUS(MPI_SOURCE)+1
            id%WORKING(id%IPTR_WORKING(TMP):id%IPTR_WORKING(TMP+1)-1)
     &           = LOCAL_LIST(1:id%IPTR_WORKING(TMP+1)-
     &           id%IPTR_WORKING(TMP))
         END DO
         DEALLOCATE(LOCAL_LIST)
      END IF
      END SUBROUTINE DMUMPS_BUILD_MAPPING_INFO
      SUBROUTINE DMUMPS_SOL_OMEGA(N, RHS,
     &    X, Y, R_W, C_W, IW, IFLAG,
     &    OMEGA, NOITER, TESTConv, 
     &    LP, ARRET, GRAIN )
      IMPLICIT NONE
      INTEGER N,  IFLAG
      INTEGER IW(N,2)
      DOUBLE PRECISION RHS(N)
      DOUBLE PRECISION X(N), Y(N)
      DOUBLE PRECISION R_W(N,2)
      DOUBLE PRECISION C_W(N)
      INTEGER LP, NOITER
      LOGICAL TESTConv
      DOUBLE PRECISION OMEGA(2)
      DOUBLE PRECISION ARRET
      INTEGER, intent(in) :: GRAIN
      DOUBLE PRECISION, PARAMETER :: CGCE=0.2D0
      DOUBLE PRECISION, PARAMETER :: CTAU=1.0D3
      INTEGER I, IMAX
      DOUBLE PRECISION OM1, OM2, DXMAX
      DOUBLE PRECISION TAU, DD
      DOUBLE PRECISION OLDOMG(2)
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D0
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D0
      INTEGER DMUMPS_IXAMAX
      INTRINSIC  abs, max
      SAVE  OM1, OLDOMG
      IMAX = DMUMPS_IXAMAX(N, X, 1, GRAIN)
      DXMAX = abs(X(IMAX))
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO I = 1, N
        TAU = (R_W(I, 2) * DXMAX + abs(RHS(I))) * dble(N) * CTAU
        DD = R_W(I, 1) + abs(RHS(I))
        IF (DD .GT. TAU * epsilon(CTAU)) THEN
          OMEGA(1) = max(OMEGA(1), abs(Y(I)) / DD)
          IW(I, 1) = 1
        ELSE
          IF (TAU .GT. ZERO) THEN
            OMEGA(2) = max(OMEGA(2),
     &                     abs(Y(I)) / (DD + R_W(I, 2) * DXMAX))
          ENDIF
          IW(I, 1) = 2
        ENDIF
      ENDDO
      IF (TESTConv) THEN
        OM2 = OMEGA(1) + OMEGA(2)
        IF (OM2 .LT. ARRET ) THEN
           IFLAG = 1
           GOTO 70
        ENDIF
        IF (NOITER .GE. 1) THEN
           IF (OM2 .GT. OM1 * CGCE) THEN
             IF (OM2 .GT. OM1) THEN
               OMEGA(1) = OLDOMG(1)
               OMEGA(2) = OLDOMG(2)
               DO I = 1, N
                 X(I) = C_W(I)
               ENDDO
               IFLAG = 2
               GOTO 70
             ENDIF
             IFLAG = 3
             GOTO 70
           ENDIF
        ENDIF
        DO I = 1, N
             C_W(I) = X(I)
        ENDDO
        OLDOMG(1) = OMEGA(1)
        OLDOMG(2) = OMEGA(2)
        OM1 = OM2
      ENDIF
      IFLAG = 0
      RETURN
   70 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_SOL_OMEGA
      SUBROUTINE DMUMPS_SOL_LCOND(N, RHS,
     &    X, Y, D, R_W, C_W, IW, KASE,
     &    OMEGA, ERX, COND, 
     &    LP, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER N, KASE, KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IW(N,2)
      DOUBLE PRECISION RHS(N)
      DOUBLE PRECISION X(N), Y(N)
      DOUBLE PRECISION D(N)
      DOUBLE PRECISION R_W(N,2)
      DOUBLE PRECISION C_W(N)
      INTEGER LP
      DOUBLE PRECISION COND(2),OMEGA(2)
      LOGICAL LCOND1, LCOND2
      INTEGER JUMP, I, IMAX
      DOUBLE PRECISION ERX, DXMAX
      DOUBLE PRECISION DXIMAX
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER :: ONE  = 1.0D0
      INTEGER DMUMPS_IXAMAX
      INTRINSIC     abs, max
      SAVE LCOND1, LCOND2, JUMP,  DXIMAX, DXMAX
      IF (KASE .EQ. 0) THEN
        LCOND1 = .FALSE.
        LCOND2 = .FALSE.
        COND(1) = ONE
        COND(2) = ONE
        ERX = ZERO
        JUMP = 1
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 30
      CASE(2)
        GOTO 10
      CASE(3)
        GOTO 110
      CASE(4)
        GOTO 150
      CASE(5)
        GOTO 35
      CASE DEFAULT
      END SELECT
   10 CONTINUE
   30 CONTINUE
   35 CONTINUE
      IMAX = DMUMPS_IXAMAX(N, X, 1, KEEP(361))
      DXMAX = abs(X(IMAX))
      DO I = 1, N
        IF (IW(I, 1) .EQ. 1) THEN
          R_W(I, 1) = R_W(I, 1) + abs(RHS(I))
          R_W(I, 2) = ZERO
          LCOND1 = .TRUE.
        ELSE
          R_W(I, 2) = R_W(I, 2) * DXMAX + R_W(I, 1)
          R_W(I, 1) = ZERO
          LCOND2 = .TRUE.
        ENDIF
      ENDDO
      DO I = 1, N
        C_W(I) = X(I) * D(I)
      ENDDO
      IMAX = DMUMPS_IXAMAX(N, C_W(1), 1, KEEP(361))
      DXIMAX = abs(C_W(IMAX))
      IF (.NOT.LCOND1) GOTO 130
  100 CONTINUE
      CALL DMUMPS_SOL_B(N, KASE, Y, COND(1), C_W, IW(1, 2), KEEP(361))
      IF (KASE .EQ. 0) GOTO 120
      IF (KASE .EQ. 1) CALL DMUMPS_SOL_MULR(N, Y, D)
      IF (KASE .EQ. 2) CALL DMUMPS_SOL_MULR(N, Y, R_W)
      JUMP = 3
      RETURN
  110 CONTINUE
      IF (KASE .EQ. 1) CALL DMUMPS_SOL_MULR(N, Y, R_W)
      IF (KASE .EQ. 2) CALL DMUMPS_SOL_MULR(N, Y, D)
      GOTO 100
  120 CONTINUE
      IF (DXIMAX .GT. ZERO) COND(1) = COND(1) / DXIMAX
      ERX = OMEGA(1) * COND(1)
  130 CONTINUE
      IF (.NOT.LCOND2) GOTO 170
      KASE = 0
  140 CONTINUE
      CALL DMUMPS_SOL_B(N, KASE, Y, COND(2), C_W, IW(1, 2), KEEP(361))
      IF (KASE .EQ. 0) GOTO 160
      IF (KASE .EQ. 1) CALL DMUMPS_SOL_MULR(N, Y, D)
      IF (KASE .EQ. 2) CALL DMUMPS_SOL_MULR(N, Y, R_W(1, 2))
      JUMP = 4
      RETURN
  150 CONTINUE
      IF (KASE .EQ. 1) CALL DMUMPS_SOL_MULR(N, Y, R_W(1, 2))
      IF (KASE .EQ. 2) CALL DMUMPS_SOL_MULR(N, Y, D)
      GOTO 140
  160 IF (DXIMAX .GT. ZERO) THEN
        COND(2) = COND(2) / DXIMAX
      ENDIF
      ERX = ERX + OMEGA(2) * COND(2)
  170 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_SOL_LCOND
      SUBROUTINE DMUMPS_SOL_CPY_FS2RHSCOMP( JBDEB, JBFIN, NBROWS,
     &   KEEP, RHSCOMP, NRHS, LRHSCOMP, FIRST_ROW_RHSCOMP, W, LD_W,
     &   FIRST_ROW_W )
         INTEGER :: JBDEB, JBFIN, NBROWS
         INTEGER :: NRHS, LRHSCOMP
         INTEGER :: FIRST_ROW_RHSCOMP
         INTEGER, INTENT(IN) :: KEEP(500)
         DOUBLE PRECISION, INTENT(INOUT) :: RHSCOMP(LRHSCOMP,NRHS)
         INTEGER :: LD_W, FIRST_ROW_W
         DOUBLE PRECISION :: W(LD_W*(JBFIN-JBDEB+1))
         INTEGER :: JJ, K, ISHIFT
!$OMP    PARALLEL DO PRIVATE(ISHIFT, JJ), IF
!$OMP&   (JBFIN-JBDEB+1 > 2*KEEP(362) .AND.
!$OMP&   NBROWS * (JBFIN-JBDEB+1) > 2*KEEP(363))
         DO K = JBDEB, JBFIN
           ISHIFT = FIRST_ROW_W + LD_W * (K-JBDEB)
           DO JJ = 0, NBROWS-1
              RHSCOMP(FIRST_ROW_RHSCOMP+JJ,K) = W(ISHIFT+JJ)
           END DO
         END DO
!$OMP    END PARALLEL DO
      RETURN
      END SUBROUTINE DMUMPS_SOL_CPY_FS2RHSCOMP
      SUBROUTINE DMUMPS_SOL_BWD_GTHR( JBDEB, JBFIN, J1, J2,
     &   RHSCOMP, NRHS, LRHSCOMP, W, LD_W, FIRST_ROW_W,
     &   IW, LIW, KEEP, N, POSINRHSCOMP_BWD )
      INTEGER, INTENT(IN) :: JBDEB, JBFIN, J1, J2
      INTEGER, INTENT(IN) :: NRHS, LRHSCOMP
      INTEGER, INTENT(IN) :: FIRST_ROW_W, LD_W, LIW
      INTEGER, INTENT(IN) :: IW(LIW)
      INTEGER, INTENT(IN) :: KEEP(500)
      DOUBLE PRECISION, INTENT(INOUT) :: RHSCOMP(LRHSCOMP,NRHS)
      DOUBLE PRECISION :: W(LD_W*(JBFIN-JBDEB+1))
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: POSINRHSCOMP_BWD(N)
      INTEGER :: ISHIFT, JJ, K, IPOSINRHSCOMP
!$OMP PARALLEL DO PRIVATE(JJ,ISHIFT,IPOSINRHSCOMP), IF
!$OMP& ((JBFIN-JBDEB+1 > 2*KEEP(362) .AND.
!$OMP& (JBFIN-JBDEB+1)*(J2-KEEP(253)-J1+1)>2*KEEP(363)))
             DO K=JBDEB, JBFIN
               ISHIFT = FIRST_ROW_W+(K-JBDEB)*LD_W
               DO JJ = J1, J2-KEEP(253)   
                 IPOSINRHSCOMP =  abs(POSINRHSCOMP_BWD(IW(JJ)))
                 W(ISHIFT+JJ-J1)= RHSCOMP(IPOSINRHSCOMP,K)
               ENDDO
             ENDDO
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE DMUMPS_SOL_BWD_GTHR
      SUBROUTINE DMUMPS_SOL_Q(MTYPE, IFLAG, N,
     &    LHS, WRHS, W, RES, GIVNORM, ANORM, XNORM, SCLNRM,
     &    MPRINT, ICNTL, KEEP,KEEP8)
      INTEGER MTYPE,N,IFLAG,ICNTL(60), KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION RES(N),LHS(N)
      DOUBLE PRECISION WRHS(N)
      DOUBLE PRECISION W(N)
      DOUBLE PRECISION RESMAX,RESL2,XNORM, SCLNRM
      DOUBLE PRECISION ANORM,DZERO
      LOGICAL GIVNORM,PROK
      INTEGER MPRINT, MP
      INTEGER K
      INTRINSIC abs, max, sqrt
      MP = ICNTL(2)
      PROK = (MPRINT .GT. 0)
      DZERO = 0.0D0
      IF (.NOT.GIVNORM) ANORM = DZERO
      RESMAX = DZERO
      RESL2  = DZERO
      DO 40 K = 1, N
        RESMAX = max(RESMAX, abs(RES(K)))
        RESL2 = RESL2 + abs(RES(K)) * abs(RES(K))
        IF (.NOT.GIVNORM) ANORM = max(ANORM, W(K))
   40 CONTINUE
      XNORM = DZERO
      DO 50 K = 1, N
        XNORM = max(XNORM, abs(LHS(K)))
   50 CONTINUE
      IF ( XNORM .EQ. DZERO .OR. (exponent(XNORM) .LT.
     &      minexponent(XNORM) + KEEP(122) ) 
     &     .OR.
     &        ( exponent(ANORM)+exponent(XNORM) .LT.
     &           minexponent(XNORM) + KEEP(122) )
     &     .OR.
     &       ( exponent(ANORM) + exponent(XNORM) -exponent(RESMAX) 
     &       .LT. minexponent(XNORM) + KEEP(122) )
     &      ) THEN
            IF (mod(IFLAG/2,2) .EQ. 0) THEN
              IFLAG = IFLAG + 2
            ENDIF
            IF ((MP .GT. 0) .AND. (ICNTL(4) .GE. 2)) WRITE( MP, * )
     &    ' max-NORM of computed solut. is zero or close to zero. '
      ENDIF
      IF (RESMAX .EQ. DZERO) THEN
        SCLNRM = DZERO
      ELSE
        SCLNRM = RESMAX / (ANORM * XNORM)
      ENDIF
      RESL2 = sqrt(RESL2)
      IF (PROK) WRITE( MPRINT, 90 ) RESMAX, RESL2, ANORM, XNORM, 
     &      SCLNRM
   90  FORMAT (/' RESIDUAL IS ............ (MAX-NORM)        =',1PD9.2/
     &       '                       .. (2-NORM)          =',1PD9.2/
     &       ' RINFOG(4):NORM OF input  Matrix  (MAX-NORM)=',1PD9.2/
     &       ' RINFOG(5):NORM OF Computed SOLUT (MAX-NORM)=',1PD9.2/
     &       ' RINFOG(6):SCALED RESIDUAL ...... (MAX-NORM)=',1PD9.2)
      RETURN
      END SUBROUTINE DMUMPS_SOL_Q
      SUBROUTINE DMUMPS_SOLVE_FWD_TRSOLVE (A, LA, APOS, NPIV, LDADIAG, 
     &           NRHS_B, WCB, LWCB, LDA_WCB, PPIV_COURANT, MTYPE, KEEP)
       INTEGER, INTENT(IN) :: MTYPE, LDADIAG, NPIV, KEEP(500)
       INTEGER, INTENT(IN) :: NRHS_B, LDA_WCB
       INTEGER(8), INTENT(IN) ::  LA, APOS, LWCB, PPIV_COURANT
       DOUBLE PRECISION, INTENT(IN) :: A(LA)
       DOUBLE PRECISION, INTENT(INOUT) :: WCB(LWCB)
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0D0)
          IF (KEEP(50).NE.0 .OR. MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dtrsv( 'U', 'T', 'U', NPIV, A(APOS), LDADIAG,
     &                   WCB(PPIV_COURANT), 1 )
               ELSE
#endif
                  CALL dtrsm( 'L','U','T','U', NPIV, NRHS_B, ONE,
     &                   A(APOS), LDADIAG, WCB(PPIV_COURANT),
     &                   LDA_WCB )
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
          ELSE
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dtrsv( 'L', 'N', 'N', NPIV, A(APOS), LDADIAG,
     &                   WCB(PPIV_COURANT), 1 )
               ELSE
#endif
                  CALL dtrsm( 'L','L','N','N', NPIV, NRHS_B, ONE,
     &                   A(APOS), LDADIAG, WCB(PPIV_COURANT),
     &                   LDA_WCB )
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
          ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_FWD_TRSOLVE
      SUBROUTINE DMUMPS_SOLVE_BWD_TRSOLVE (A, LA, APOS, NPIV, LDADIAG, 
     &           NRHS_B, WCB, LWCB, LDA_WCB, PPIV_COURANT, MTYPE, KEEP)
       INTEGER, INTENT(IN) :: MTYPE, LDADIAG, NPIV, KEEP(500)
       INTEGER, INTENT(IN) :: NRHS_B, LDA_WCB
       INTEGER(8), INTENT(IN) ::  LA, APOS, LWCB, PPIV_COURANT
       DOUBLE PRECISION, INTENT(IN) :: A(LA)
       DOUBLE PRECISION, INTENT(INOUT) :: WCB(LWCB)
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0D0)
          IF (MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dtrsv( 'L', 'T', 'N', NPIV, A(APOS), LDADIAG,
     &                   WCB(PPIV_COURANT), 1 )
               ELSE
#endif
                  CALL dtrsm( 'L','L','T','N', NPIV, NRHS_B, ONE,
     &                   A(APOS), LDADIAG, WCB(PPIV_COURANT),
     &                   LDA_WCB )
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
          ELSE
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dtrsv( 'U', 'N', 'U', NPIV, A(APOS), LDADIAG,
     &                   WCB(PPIV_COURANT), 1 )
               ELSE
#endif
                  CALL dtrsm( 'L','U','N','U', NPIV, NRHS_B, ONE,
     &                   A(APOS), LDADIAG, WCB(PPIV_COURANT),
     &                   LDA_WCB )
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
          ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_BWD_TRSOLVE
      SUBROUTINE DMUMPS_SOLVE_GEMM_UPDATE 
     &           (A, LA, APOS1, NX, LDA, NY,
     &           NRHS_B, WCB, LWCB, PTRX, LDX,
     &           PTRY, LDY,
     &           MTYPE, KEEP, COEF_Y )
       INTEGER, INTENT(IN) :: MTYPE, NY, NX, KEEP(500)
       INTEGER, INTENT(IN) :: NRHS_B, LDY, LDA, LDX
       INTEGER(8), INTENT(IN) ::  LA, APOS1, LWCB, PTRX, 
     &                            PTRY
       DOUBLE PRECISION, INTENT(IN) :: A(LA)
       DOUBLE PRECISION, INTENT(INOUT) :: WCB(LWCB)
       DOUBLE PRECISION, INTENT(IN) :: COEF_Y
      DOUBLE PRECISION ALPHA, ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, ALPHA=-1.0D0)
         IF ( NX .NE. 0 .AND. NY.NE.0 ) THEN
            IF ( MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dgemv('T', NX, NY, ALPHA, A(APOS1),
     &            LDA,  WCB(PTRX), 1, COEF_Y,
     &            WCB(PTRY), 1)
               ELSE
#endif
                  CALL dgemm('T', 'N', NY, NRHS_B, NX, ALPHA,
     &            A(APOS1), LDA, WCB(PTRX), LDX, COEF_Y,
     &            WCB(PTRY), LDY)
#if defined(MUMPS_USE_BLAS2)
               END IF
#endif
            ELSE                
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL dgemv('N',NY, NX, ALPHA, A(APOS1),
     &                 LDA, WCB(PTRX), 1,
     &                 COEF_Y, WCB(PTRY), 1 )
               ELSE
#endif
                  CALL dgemm('N', 'N', NY, NRHS_B, NX, ALPHA,
     &                 A(APOS1), LDA, WCB(PTRX), LDX,
     &                 COEF_Y, WCB(PTRY), LDY)
#if defined(MUMPS_USE_BLAS2)
               END IF
#endif
            END IF
         END IF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_GEMM_UPDATE
      SUBROUTINE DMUMPS_SOLVE_LD_AND_RELOAD (
     &  INODE, N, NPIV, LIELL, NELIM, NSLAVES,
     &  PPIV_COURANT, 
     &  IW, IPOS, LIW, 
     &  A, LA, APOS,
     &  WCB, LWCB, LD_WCBPIV, 
     &  RHSCOMP, LRHSCOMP, NRHS, 
     &  POSINRHSCOMP_FWD, JBDEB, JBFIN, 
     &  MTYPE, KEEP, OOCWRITE_COMPATIBLE_WITH_BLR
     &  )
      USE DMUMPS_OOC 
      INTEGER, INTENT(IN) :: MTYPE, INODE, N, NPIV, LIELL,
     &                       NELIM, NSLAVES
      INTEGER, INTENT(IN) :: LRHSCOMP, NRHS, LIW, JBDEB, JBFIN
      INTEGER, INTENT(IN) :: IW(LIW), IPOS, POSINRHSCOMP_FWD(N)
      INTEGER(8), INTENT(IN) :: LWCB, APOS, LA, PPIV_COURANT
      INTEGER, INTENT(IN) :: LD_WCBPIV
      INTEGER, INTENT(IN) :: KEEP(500)
      DOUBLE PRECISION, INTENT(IN) ::  WCB( LWCB ), A( LA )
      DOUBLE PRECISION, INTENT(INOUT) :: RHSCOMP(LRHSCOMP, NRHS)
      LOGICAL, INTENT(IN) :: OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER    :: TempNROW, J1, J3, PANEL_SIZE, TYPEF
      INTEGER    :: IPOSINRHSCOMP, JJ, K, NBK, LDAJ, 
     &              LDAJ_ini, NBK_ini, LDAJ_FIRST_PANEL, NRHS_B
      INTEGER(8) :: IFR8 , APOS1, APOS2, APOSOFF, IFR_ini8, 
     &              POSWCB1, POSWCB2
      DOUBLE PRECISION    :: VALPIV, A11, A22, A12, DETPIV
!$    LOGICAL :: OMP_FLAG
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0D0)
      NRHS_B = JBFIN-JBDEB+1
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
         J1 = IPOS + 1
         J3 = IPOS + NPIV
      ELSE
         J1 = IPOS + LIELL + 1
         J3 = IPOS + LIELL + NPIV
      END IF
      IPOSINRHSCOMP =  POSINRHSCOMP_FWD(IW(J1)) 
      IF ( KEEP(50) .eq. 0 ) THEN
!$         OMP_FLAG=(NRHS_B.GE.KEEP(362).AND.NRHS_B*NPIV.GE.KEEP(363))
!$OMP PARALLEL DO PRIVATE(IFR8) IF (OMP_FLAG)
           DO K=JBDEB,JBFIN
             IFR8 =  PPIV_COURANT + (K-JBDEB)*LD_WCBPIV
             RHSCOMP(IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1, K) =
     &            WCB(IFR8:IFR8+int(NPIV-1,8))
           ENDDO
!$OMP END PARALLEL DO
      ELSE
         IFR8 = PPIV_COURANT - 1_8
         IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN 
          IF (MTYPE.EQ.1) THEN
            IF ((MTYPE.EQ.1).AND.NSLAVES.NE.0) THEN
              TempNROW= NPIV+NELIM
              LDAJ_FIRST_PANEL=TempNROW
            ELSE
              TempNROW= LIELL
              LDAJ_FIRST_PANEL=TempNROW
            ENDIF
            TYPEF=TYPEF_L
          ELSE 
            TempNROW= NPIV
            LDAJ_FIRST_PANEL=LIELL  
            TYPEF= TYPEF_U
          ENDIF
          PANEL_SIZE = DMUMPS_OOC_PANEL_SIZE( LDAJ_FIRST_PANEL )
          LDAJ = TempNROW  
         ELSE                
            LDAJ = NPIV 
         ENDIF
         APOS1 = APOS
         JJ    = J1
         IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN 
            NBK   = 0           
         ENDIF
         IFR_ini8 = PPIV_COURANT - 1_8
         LDAJ_ini = LDAJ   
         IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) 
     &                                      NBK_ini = NBK  
!$           OMP_FLAG = ( JBFIN-JBDEB+1.GE.KEEP(362) .AND.
!$   &                  ((J3-J1+1)*(JBFIN-JBDEB+1) .GE. KEEP(363)))
!$OMP  PARALLEL DO PRIVATE(JJ,IFR8,NBK,APOS1,APOS2,APOSOFF,VALPIV,
!$OMP&      POSWCB1, POSWCB2,A11,A22,A12,DETPIV,LDAJ)  IF(OMP_FLAG)
         DO K = JBDEB, JBFIN
            IFR8  = IFR_ini8 + int(K-JBDEB,8)*int(LD_WCBPIV,8)
            NBK   = NBK_ini
            APOS1 = APOS
            LDAJ  = LDAJ_ini
            JJ    = J1
          DO 
            IF (JJ .GT. J3) EXIT
            IFR8  = IFR8 + 1_8
            IF (IW(JJ+LIELL) .GT. 0) THEN
               VALPIV  = ONE/A( APOS1 )
                 RHSCOMP(IPOSINRHSCOMP+JJ-J1 , K ) = 
     &                 WCB( IFR8 ) * VALPIV
              IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR)
     &        THEN
                NBK = NBK+1
                IF (NBK.EQ.PANEL_SIZE) THEN
                  NBK = 0
                  LDAJ = LDAJ - PANEL_SIZE
                ENDIF
              ENDIF
              APOS1 = APOS1 + int(LDAJ + 1,8)
              JJ = JJ+1
            ELSE
              IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR)
     &        THEN
                NBK = NBK+1
              ENDIF
              APOS2 = APOS1+int(LDAJ+1,8)
              IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR)
     &        THEN
                APOSOFF = APOS1+int(LDAJ,8)
              ELSE
                APOSOFF=APOS1+1_8
              ENDIF
              A11 = A(APOS1)
              A22 = A(APOS2)
              A12 = A(APOSOFF)
              DETPIV = A11*A22 - A12**2
              A22 = A11/DETPIV
              A11 = A(APOS2)/DETPIV
              A12 = -A12/DETPIV
              POSWCB1 = IFR8
              POSWCB2 = POSWCB1+1_8
              RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) =
     &               WCB(POSWCB1)*A11
     &               + WCB(POSWCB2)*A12
              RHSCOMP(IPOSINRHSCOMP+JJ-J1+1,K) = 
     &                 WCB(POSWCB1)*A12
     &                 + WCB(POSWCB2)*A22
              IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR)
     &        THEN
                 NBK = NBK+1
                 IF (NBK.GE.PANEL_SIZE) THEN
                    LDAJ = LDAJ - NBK
                    NBK = 0
                 ENDIF
              ENDIF
              APOS1 = APOS2 + int(LDAJ + 1,8)
              JJ = JJ+2
              IFR8 = IFR8+1_8
            ENDIF  
           ENDDO   
         ENDDO     
!$OMP END PARALLEL DO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_SOLVE_LD_AND_RELOAD
      SUBROUTINE DMUMPS_SET_SCALING_LOC( scaling_data, N, ILOC, LILOC,
     &  COMM, MYID, I_AM_SLAVE, MASTER, NB_BYTES, NB_BYTES_MAX,
     &  K16_8, LP, LPOK, ICNTL, INFO )
      IMPLICIT NONE
      type scaling_data_t
        SEQUENCE
        DOUBLE PRECISION, dimension(:), pointer :: SCALING
        DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type (scaling_data_t), INTENT(INOUT) :: scaling_data
      INTEGER, INTENT(IN) :: N, LILOC, COMM, MYID, MASTER, LP
      INTEGER, INTENT(IN) :: ILOC(LILOC)
      INTEGER(8), INTENT(INOUT) :: NB_BYTES, NB_BYTES_MAX
      INTEGER(8), INTENT(IN) :: K16_8 
      LOGICAL, INTENT(IN) :: I_AM_SLAVE, LPOK
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER, INTENT(IN) :: ICNTL(60)
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: SCALING
      INTEGER :: I, IERR_MPI, allocok
      INCLUDE 'mpif.h'
      NULLIFY(scaling_data%SCALING_LOC)
      IF (I_AM_SLAVE) THEN
        ALLOCATE(scaling_data%SCALING_LOC(max(1,LILOC)),
     &  stat=allocok)
        IF (allocok > 0) THEN
          INFO(1)=-13
          INFO(2)=max(1,LILOC)
          GOTO 35
        ENDIF
        NB_BYTES     = NB_BYTES + int(max(1,LILOC),8)*K16_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
      ENDIF
      IF (MYID .NE. MASTER) THEN
        ALLOCATE(SCALING(N), stat=allocok)
        IF (allocok > 0) THEN
          IF (LPOK) THEN
            WRITE(LP,*) 'Error allocating temporary scaling array'
          ENDIF
          INFO(1)=-13
          INFO(2)=N
          GOTO 35
        ENDIF
        NB_BYTES     = NB_BYTES + int(N,8)*K16_8
        NB_BYTES_MAX = max(NB_BYTES_MAX,NB_BYTES)
      ELSE
        SCALING => scaling_data%SCALING
      ENDIF
 35   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                         COMM, MYID )
      IF (INFO(1) .LT. 0) GOTO 90
      CALL MPI_BCAST( SCALING(1), N, MPI_DOUBLE_PRECISION,
     &                MASTER, COMM, IERR_MPI)
      IF ( I_AM_SLAVE ) THEN
        DO I = 1, LILOC
          IF (ILOC(I) .GE. 1 .AND. ILOC(I) .LE. N) THEN
            scaling_data%SCALING_LOC(I) = SCALING(ILOC(I))
          ENDIF
        ENDDO
      ENDIF
 90   CONTINUE
      IF (MYID.NE. MASTER) THEN
        IF (associated(SCALING)) THEN
          DEALLOCATE(SCALING)
          NB_BYTES     = NB_BYTES - int(N,8)*K16_8
        ENDIF
      ENDIF
      NULLIFY(SCALING)
      IF (INFO(1) .LT. 0) THEN
        IF (associated(scaling_data%SCALING_LOC)) THEN
          DEALLOCATE(scaling_data%SCALING_LOC)
          NULLIFY(scaling_data%SCALING_LOC)
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SET_SCALING_LOC
