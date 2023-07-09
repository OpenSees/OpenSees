      SUBROUTINE XLAENV( ISPEC, NVALUE )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC, NVALUE
*     ..
*
*  Purpose
*  =======
*
*  XLAENV sets certain machine- and problem-dependent quantities
*  which will later be retrieved by ILAENV.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be set in the COMMON array IPARMS.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form)
*          = 7: the number of processors
*          = 8: another crossover point, for the multishift QR and QZ
*               methods for nonsymmetric eigenvalue problems.
*
*  NVALUE  (input) INTEGER
*          The value of the parameter specified by ISPEC.
*
*  =====================================================================
*
*     .. Arrays in Common ..
      INTEGER            IPARMS( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / CLAENV / IPARMS
*     ..
*     .. Save statement ..
      SAVE               / CLAENV /
*     ..
*     .. Executable Statements ..
*
      IF( ISPEC.GE.1 .AND. ISPEC.LE.8 ) THEN
         IPARMS( ISPEC ) = NVALUE
      END IF
*
      RETURN
*
*     End of XLAENV
*
      END
