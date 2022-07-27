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
      SUBROUTINE DMUMPS_MV_ELT( N, NELT, ELTPTR, ELTVAR, A_ELT,
     &                          X, Y, K50, MTYPE )
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C  To perform the matrix vector product
C      A_ELT X = Y    if MTYPE = 1
C      A_ELT^T X = Y  if MTYPE = 0
C
C  If K50 is different from 0, then the elements are
C  supposed to be in symmetric packed storage; the
C  lower part is stored by columns.
C  Otherwise, the element is square, stored by columns.
C
C  Note
C  ====
C
C  A_ELT is processed entry by entry and this code is not
C  optimized. In particular, one could gather/scatter
C  X / Y for each element to improve performance.
C
C  Arguments
C  =========
C
      INTEGER N, NELT, K50, MTYPE
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( * )
      DOUBLE PRECISION A_ELT( * ), X( N ), Y( N )
C
C  Local variables
C  ===============
C
      INTEGER IEL, I , J, SIZEI, IELPTR
      INTEGER(8) :: K8
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
C
C
C     Executable statements
C     =====================
C
      Y = ZERO
      K8 = 1_8
C     --------------------
C     Process the elements
C     --------------------
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( K50 .eq. 0 ) THEN
C         -------------------
C         Unsymmetric element
C         stored by columns
C         -------------------
          IF ( MTYPE .eq. 1 ) THEN
C           -----------------
C           Compute A_ELT x X
C           -----------------
            DO J = 1, SIZEI
              TEMP = X( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                Y( ELTVAR( IELPTR + I ) ) =
     &          Y( ELTVAR( IELPTR + I ) ) +
     &             A_ELT( K8 ) * TEMP
                K8 = K8 + 1
              END DO
            END DO
          ELSE
C           -------------------
C           Compute A_ELT^T x X
C           -------------------
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP + 
     &          A_ELT( K8 ) * X( ELTVAR( IELPTR + I ) )
                K8 = K8 + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
            END DO
          END IF
        ELSE
C         -----------------
C         Symmetric element
C         L stored by cols
C         -----------------
          DO J = 1, SIZEI
C           Diagonal counted once
            Y( ELTVAR( IELPTR + J ) ) =
     &      Y( ELTVAR( IELPTR + J ) ) +
     &           A_ELT( K8 ) * X( ELTVAR( IELPTR + J ) )
            K8 = K8 + 1
            DO I = J+1, SIZEI
C             Off diagonal + transpose
              Y( ELTVAR( IELPTR + I ) ) =
     &        Y( ELTVAR( IELPTR + I ) ) +
     &           A_ELT( K8 ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     &        Y( ELTVAR( IELPTR + J ) ) +
     &           A_ELT( K8 ) * X( ELTVAR( IELPTR + I ) )
              K8 = K8 + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_MV_ELT
      SUBROUTINE DMUMPS_LOC_MV8
     &( N, NZ_loc8, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     &  LDLT, MTYPE)
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Perform a distributed matrix vector product.
C        Y_loc <- A X   if MTYPE = 1
C        Y_loc <- A^T X if MTYPE = 0
C
C     Notes:
C     =====
C
C     1) assembly of all Y_loc still has to be done on exit.
C     2) X should be available on all processors.
C
C     Arguments:
C     =========
C
      INTEGER N
      INTEGER(8) :: NZ_loc8
      INTEGER IRN_loc( NZ_loc8 ), JCN_loc( NZ_loc8 )
      DOUBLE PRECISION A_loc( NZ_loc8 ), X( N ), Y_loc( N )
      INTEGER LDLT, MTYPE
C
C     Locals variables:
C     ================
C
      INTEGER I, J
      INTEGER(8) :: K8
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      Y_loc = ZERO
      IF ( LDLT .eq. 0 ) THEN
C       Unsymmetric
        IF ( MTYPE .eq. 1 ) THEN
C         No transpose
          DO K8 = 1_8, NZ_loc8
            I = IRN_loc(K8)
            J = JCN_loc(K8)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &          (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K8) * X(J)
        ENDDO
        ELSE
C         Transpose
          DO K8 = 1_8, NZ_loc8
            I = IRN_loc(K8)
            J = JCN_loc(K8)
            IF ((I .LE. 0) .OR. (I .GT. N)
     &        .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(J) = Y_loc(J) + A_loc(K8) * X(I)
        ENDDO
        END IF
      ELSE
C       Lower (or upper) part of symmetric
C       matrix was provided (LDLT facto)
        DO K8 = 1_8, NZ_loc8
          I = IRN_loc(K8)
          J = JCN_loc(K8)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &        (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K8) * X(J)
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + A_loc(K8) * X(I)
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_LOC_MV8
      SUBROUTINE DMUMPS_MV8( N, NZ8, IRN, ICN, ASPK, X, Y,
     &                      LDLT, MTYPE, MAXTRANS, PERM,
     &                      IFLAG, IERROR )
C
C     Purpose:
C     =======
C
C     Perform matrix-vector product
C        Y <- A X if MTYPE = 1
C        Y <- A^T X if MTYPE = 0
C
C
C     Note:
C     ====
C
C     MAXTRANS should be set to 1 if a column permutation
C     was applied on A and we still want the matrix vector
C     product wrt the original matrix.
C
C     Arguments:
C     =========
C
      INTEGER N, LDLT, MTYPE, MAXTRANS
      INTEGER(8) :: NZ8
      INTEGER IRN( NZ8 ), ICN( NZ8 ) 
      INTEGER PERM( N )
      DOUBLE PRECISION ASPK( NZ8 ), X( N ), Y( N )
      INTEGER, intent(inout) :: IFLAG, IERROR
C
C     Local variables
C     ===============
C
      INTEGER I, J
      INTEGER(8) :: K8
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PX
      DOUBLE PRECISION ZERO
      INTEGER :: allocok
      PARAMETER( ZERO = 0.0D0 )
      Y = ZERO
      ALLOCATE(PX(N), stat=allocok)
      IF (allocok < 0) THEN
        IFLAG  = -13
        IERROR = N
        RETURN
      ENDIF
C
C     --------------------------------------
C     Permute X if A has been permuted
C     with some max-trans column permutation
C     --------------------------------------
      IF ( MAXTRANS .eq. 1 .and. MTYPE .eq. 1) THEN
        DO I = 1, N
          PX(I) = X( PERM( I ) )
        END DO
      ELSE
        PX = X
      END IF
      IF ( LDLT .eq. 0 ) THEN
C
C     Complete unsymmetric matrix was provided (LU facto)
       IF (MTYPE .EQ. 1) THEN
        DO K8 = 1_8, NZ8
          I = IRN(K8)
          J = ICN(K8)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(I) = Y(I) + ASPK(K8) * PX(J)
        ENDDO
       ELSE
        DO K8 = 1_8, NZ8
          I = IRN(K8)
          J = ICN(K8)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(J) = Y(J) + ASPK(K8) * PX(I)
        ENDDO
       ENDIF
C
      ELSE
C
C       Lower (or upper) part of symmetric
C       matrix was provided (LDLT facto)
        DO K8 = 1_8, NZ8
          I = IRN(K8)
          J = ICN(K8)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(I) = Y(I) + ASPK(K8) * PX(J)
          IF (J.NE.I) THEN
            Y(J) = Y(J) + ASPK(K8) * PX(I)
          ENDIF
        ENDDO
      END IF
      IF ( MAXTRANS .EQ. 1 .AND. MTYPE .eq. 0 ) THEN
      PX = Y
      DO I = 1, N
        Y( PERM( I ) ) = PX( I )
      END DO
      END IF
      DEALLOCATE(PX)
      RETURN
      END SUBROUTINE DMUMPS_MV8
C
C
      SUBROUTINE DMUMPS_LOC_OMEGA1
     &( N, NZ_loc8, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     &  LDLT, MTYPE)
      IMPLICIT NONE
C
C     Purpose:
C     =======
C     Compute
C        * If MTYPE = 1
C            Y_loc(i) = Sum | Aij | | Xj |
C                        j
C        * If MTYPE = 0
C            Y_loc(j) = Sum | Aij | | Xi |
C
C
C     Notes:
C     =====
C
C     1) assembly of all Y_loc still has to be done.
C     2) X should be available on all processors.
C
C     Arguments:
C     =========
C
      INTEGER N
      INTEGER(8) :: NZ_loc8
      INTEGER IRN_loc( NZ_loc8 ), JCN_loc( NZ_loc8 )
      DOUBLE PRECISION A_loc( NZ_loc8 ), X( N )
      DOUBLE PRECISION Y_loc( N )
      INTEGER LDLT, MTYPE
C
C     Local variables:
C     ===============
C
      INTEGER I, J
      INTEGER(8) :: K8
      DOUBLE PRECISION, PARAMETER :: RZERO=0.0D0
C
      Y_loc = RZERO
      IF ( LDLT .eq. 0 ) THEN
C       Unsymmetric
        IF ( MTYPE .eq. 1 ) THEN
C         No transpose
          DO K8 = 1_8, NZ_loc8
            I = IRN_loc(K8)
            J = JCN_loc(K8)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &          (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
            Y_loc(I) = Y_loc(I) + abs( A_loc(K8) * X(J) )
          ENDDO
        ELSE
C         Transpose
          DO K8 = 1_8, NZ_loc8
            I = IRN_loc(K8)
            J = JCN_loc(K8)
            IF ((I .LE. 0) .OR. (I .GT. N)
     &        .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(J) = Y_loc(J) + abs( A_loc(K8) * X(I) )
          ENDDO
        END IF
      ELSE
C       Lower (or upper) part of symmetric
C       matrix was provided (LDLT facto)
        DO K8 = 1_8, NZ_loc8
          I = IRN_loc(K8)
          J = JCN_loc(K8)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &        (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + abs( A_loc(K8) * X(J) )
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + abs( A_loc(K8) * X(I) )
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_LOC_OMEGA1
