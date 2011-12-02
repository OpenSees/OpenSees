/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Sparse BLAS-3, using some dense BLAS-3 operations.
 */

#include "pdsp_defs.h"

int
sp_dgemm(char *trans, int m, int n, int k, double alpha, SuperMatrix *A,
	 double *b, int ldb, double beta, double *c, int ldc)
{
/*  Purpose   
    =======   

    sp_dgemm()  performs one of the matrix-vector operations   
       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   
    where alpha and beta are scalars, x and y are vectors and A is a
    sparse m by n matrix.   

    Parameters   
    ==========   

    TRANS  - (input) char*
             On entry, TRANS specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   
                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   
                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

    M      - int   
             On entry,  M  specifies  the number of rows of the matrix 
	     op( A ) and of the matrix C.  M must be at least zero. 
	     Unchanged on exit.   

    N      - int
             On entry,  N specifies the number of columns of the matrix 
	     op( B ) and the number of columns of the matrix C. N must be 
	     at least zero.
	     Unchanged on exit.   

    K      - int
             On entry, K specifies the number of columns of the matrix 
	     op( A ) and the number of rows of the matrix op( B ). K must 
	     be at least  zero.   
             Unchanged on exit.
	     
    ALPHA  - (input) double
             On entry, ALPHA specifies the scalar alpha.   

    A      - (input) SuperMatrix*
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   

    B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ),
             Before entry thhe leading  k by n part of the array B must
	     contain the matrix B.
             Unchanged on exit.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
             in the calling (sub) program. LDB must be at least max( 1, n ).  
             Unchanged on exit.   

    BETA   - (input) double
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   

    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry, the leading m by n part of the array C must 
             contain the matrix C,  except when beta is zero, in which 
             case C need not be set on entry.   
             On exit, the array C is overwritten by the m by n matrix 
	     ( alpha*op( A )*B + beta*C ).   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
             in the calling (sub)program. LDC must be at least max(1,m).   
             Unchanged on exit.   

    ==== Sparse Level 3 Blas routine.   
*/
    int    incx = 1, incy = 1;
    int    j;

    for (j = 0; j < n; ++j) {
	sp_dgemv(trans, alpha, A, &b[ldb*j], incx, beta, &c[ldc*j], incy);
    }
    
    return 0;
}
