      SUBROUTINE DLAPTM( N, NRHS, ALPHA, D, E, X, LDX, BETA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            LDB, LDX, N, NRHS
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DLAPTM multiplies an N by NRHS matrix X by a symmetric tridiagonal
*  matrix A and stores the result in a matrix B.  The operation has the
*  form
*
*     B := alpha * A * X + beta * B
*
*  where alpha may be either 1. or -1. and beta may be 0., 1., or -1.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices X and B.
*
*  ALPHA   (input) DOUBLE PRECISION
*          The scalar alpha.  ALPHA must be 1. or -1.; otherwise,
*          it is assumed to be 0.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix A.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal or superdiagonal elements of A.
*
*  X       (input) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          The N by NRHS matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(N,1).
*
*  BETA    (input) DOUBLE PRECISION
*          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,
*          it is assumed to be 1.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N by NRHS matrix B.
*          On exit, B is overwritten by the matrix expression
*          B := alpha * A * X + beta * B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(N,1).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Multiply B by BETA if BETA.NE.1.
*
      IF( BETA.EQ.ZERO ) THEN
         DO 20 J = 1, NRHS
            DO 10 I = 1, N
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE IF( BETA.EQ.-ONE ) THEN
         DO 40 J = 1, NRHS
            DO 30 I = 1, N
               B( I, J ) = -B( I, J )
   30       CONTINUE
   40    CONTINUE
      END IF
*
      IF( ALPHA.EQ.ONE ) THEN
*
*        Compute B := B + A*X
*
         DO 60 J = 1, NRHS
            IF( N.EQ.1 ) THEN
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) +
     $                     E( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) + E( N-1 )*X( N-1, J ) +
     $                     D( N )*X( N, J )
               DO 50 I = 2, N - 1
                  B( I, J ) = B( I, J ) + E( I-1 )*X( I-1, J ) +
     $                        D( I )*X( I, J ) + E( I )*X( I+1, J )
   50          CONTINUE
            END IF
   60    CONTINUE
      ELSE IF( ALPHA.EQ.-ONE ) THEN
*
*        Compute B := B - A*X
*
         DO 80 J = 1, NRHS
            IF( N.EQ.1 ) THEN
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) -
     $                     E( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) - E( N-1 )*X( N-1, J ) -
     $                     D( N )*X( N, J )
               DO 70 I = 2, N - 1
                  B( I, J ) = B( I, J ) - E( I-1 )*X( I-1, J ) -
     $                        D( I )*X( I, J ) - E( I )*X( I+1, J )
   70          CONTINUE
            END IF
   80    CONTINUE
      END IF
      RETURN
*
*     End of DLAPTM
*
      END
