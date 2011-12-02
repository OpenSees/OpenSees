/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
#include <stdlib.h>
#include "pdsp_defs.h"
#include "util.h"

/*
 * Generate a banded square matrix A, with dimension n and semi-bandwidth b.
 */
void
dband(int n, int b, int nonz, double **nzval, int **rowind, int **colptr)
{
    int iseed[] = {1992,1993,1994,1995};    
    register int i, j, ub, lb, ilow, ihigh, lasta = 0;
    double *a;
    int    *asub, *xa;
    double *val;
    int    *row;
    extern double dlaran_();
    
    printf("A banded matrix.");
    dallocateA(n, nonz, nzval, rowind, colptr); /* Allocate storage */
    a    = *nzval;
    asub = *rowind;
    xa   = *colptr;
    ub = lb = b;
    
    for (i = 0; i < 4; ++i) iseed[i] = abs( iseed[i] ) % 4096;
    if ( iseed[3] % 2 != 1 ) ++iseed[3];

    for (j = 0; j < n; ++j) {
	xa[j] = lasta;
	val = &a[lasta];
	row = &asub[lasta];
	ilow = MAX(0, j - ub);
	ihigh = MIN(n-1, j + lb);
	for (i = ilow; i <= ihigh; ++i) {
	    val[i-ilow] = dlaran_(iseed);
	    row[i-ilow] = i;
	}
	lasta += ihigh - ilow + 1;
    } /* for j ... */
    xa[n] = lasta;
}

/*
 * Generate a block diagonal matrix A.
 */
void
dblockdiag(int nb, /* number of blocks */
	   int bs, /* block size */
	   int nonz, double **nzval, int **rowind, int **colptr)
{
    int iseed[] = {1992,1993,1994,1995};    
    register int i, j, b, n, lasta = 0, cstart, rstart;
    double *a;
    int    *asub, *xa;
    double *val;
    int    *row;
    extern double dlaran_();
    
    n = bs * nb;
    printf("A block diagonal matrix: nb %d, bs %d, n %d\n", nb, bs, n);
    dallocateA(n, nonz, nzval, rowind, colptr); /* Allocate storage */
    a    = *nzval;
    asub = *rowind;
    xa   = *colptr;
    
    for (i = 0; i < 4; ++i) iseed[i] = abs( iseed[i] ) % 4096;
    if ( iseed[3] % 2 != 1 ) ++iseed[3];

    for (b = 0; b < nb; ++b) {
	cstart = b * bs; /* start of the col # of the current block */
	rstart = b * bs; /* start of the row # of the current block */
	for (j = cstart; j < cstart + bs; ++j) {
	    xa[j] = lasta;
	    val = &a[lasta];
	    row = &asub[lasta];
	    for (i = 0; i < bs; ++i) {
		val[i] = dlaran_(iseed);
		row[i] = i + rstart;
	    }
	    lasta += bs;
	} /* for j ... */
    } /* for b ... */
    
    xa[n] = lasta;
}

double dlaran_(int *iseed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   

    Purpose   
    =======   

    DLARAN returns a random real number from a uniform (0,1)   
    distribution.   

    Arguments   
    =========   

    ISEED   (input/output) INT array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    Further Details   
    ===============   

    This routine uses a multiplicative congruential method with modulus   
    2**48 and multiplier 33952834046453 (see G.S.Fishman,   
    'Multiplicative congruential random number generators with modulus   
    2**b: an exhaustive analysis for b = 32 and a partial analysis for   
    b = 48', Math. Comp. 189, pp 331-344, 1990).   

    48-bit integers are stored in 4 integer array elements with 12 bits   
    per element. Hence the routine is portable across machines with   
    integers of 32 bits or more.   

    ===================================================================== 
*/
    
    /* Local variables */
    int it1, it2, it3, it4;

    --iseed;

    /* multiply the seed by the multiplier modulo 2**48 */
    it4 = iseed[4] * 2549;
    it3 = it4 / 4096;
    it4 -= it3 << 12;
    it3 = it3 + iseed[3] * 2549 + iseed[4] * 2508;
    it2 = it3 / 4096;
    it3 -= it2 << 12;
    it2 = it2 + iseed[2] * 2549 + iseed[3] * 2508 + iseed[4] * 322;
    it1 = it2 / 4096;
    it2 -= it1 << 12;
    it1 = it1 + iseed[1] * 2549 + iseed[2] * 2508 + iseed[3] * 322 + iseed[4] 
	    * 494;
    it1 %= 4096;

   /* return updated seed */

    iseed[1] = it1;
    iseed[2] = it2;
    iseed[3] = it3;
    iseed[4] = it4;

   /* convert 48-bit integer to a real number in the interval (0,1) */

    return ((double) it1 +
	    ((double) it2 + ((double) it3 + (double) it4 * 2.44140625e-4) *
	     2.44140625e-4) * 2.44140625e-4) * 2.44140625e-4;

} /* dlaran_ */

