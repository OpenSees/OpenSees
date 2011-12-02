#include <math.h>
#include <stdlib.h>
#include "pdsp_defs.h"

int
pdgstrf_pivotL(
	       const int  pnum,  /* process number */
	       const int  jcol,  /* current column */
	       const double u,   /* diagonal pivoting threshold */
	       yes_no_t *usepr,  /* re-use the pivot sequence given by
				    perm_r[]/inv_perm_r[] */
	       int   *perm_r,    /* modified - row pivotings */
	       int   *inv_perm_r,/* modified - inverse of perm_r */
	       int   *inv_perm_c,/* in - used to find diagonal of Pc*A*Pc' */
	       int   *pivrow,    /* the pivot row for this column */
	       GlobalLU_t *Glu,  /* modified - global LU data structures */
	       Gstat_t *Gstat    /* modified */
       )
{
/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Purpose
 * =======
 *   Performs the numerical pivoting on the current column of L,
 *   and the CDIV operation.
 *
 * Pivot policy
 * ============
 *   (1) Compute thresh = u * max_(i>=j) abs(A_ij);
 *   (2) IF user specifies pivot row k and abs(A_kj) >= thresh THEN
 *           pivot row = k;
 *       ELSE IF abs(A_jj) >= thresh THEN
 *           pivot row = j;
 *       ELSE
 *           pivot row = m;
 * 
 *   Note: If you absolutely want to use a given pivot order, then set u=0.0.
 *
 * Return value
 * ============
 *   0      success;
 *   i > 0  U(i,i) is exactly zero.
 *
 */
    register int fsupc; /* first column in the supernode */
    register int nsupc; /* no of columns in the supernode */
    register int nsupr; /* no of rows in the supernode */
    register int lptr;  /* the starting subscript of the supernode */
    register int pivptr, old_pivptr, diag, diagind;
    register int isub, icol, k, itemp;
    register double pivmax, temp, thresh;
    double       *lu_sup_ptr; 
    double       *lu_col_ptr;
    int          *lsub_ptr;
    int          *lsub;
    double       *lusup;
    int          *xlusup;

    /* Initialize pointers */
    lsub       = Glu->lsub;
    lusup      = Glu->lusup;
    xlusup     = Glu->xlusup;
    fsupc      = Glu->xsup[Glu->supno[jcol]];
    nsupc      = jcol - fsupc;	        /* excluding jcol; nsupc >= 0 */
    lptr       = Glu->xlsub[fsupc];
    nsupr      = Glu->xlsub_end[fsupc] - lptr;
    lu_sup_ptr = &lusup[xlusup[fsupc]];	/* start of the current supernode */
    lu_col_ptr = &lusup[xlusup[jcol]];	/* start of jcol in the supernode */
    lsub_ptr   = &lsub[lptr];	/* start of row indices of the supernode */

#ifdef CHK_PIVOT
    printf("Before cdiv: col %d\n", jcol);
    for (k = nsupc; k < nsupr; k++) 
	printf("  lu[%d] %f\n", lsub_ptr[k], lu_col_ptr[k]);
#endif
    
    /* Determine the largest abs numerical value for partial pivoting;
       Also search for user-specified pivot, and diagonal element. */
    if ( *usepr == YES ) *pivrow = inv_perm_r[jcol];
    diagind = inv_perm_c[jcol];
    pivmax = 0.0;
    pivptr = nsupc;
    diag = EMPTY;
    old_pivptr = nsupc;
    for (isub = nsupc; isub < nsupr; ++isub) {
	temp = fabs (lu_col_ptr[isub]);
	if ( temp > pivmax ) {
	    pivmax = temp;
	    pivptr = isub;
	}
	if ( *usepr == YES && lsub_ptr[isub] == *pivrow ) old_pivptr = isub;
	if ( lsub_ptr[isub] == diagind ) diag = isub;
    }

    /* Test for singularity */
    if ( pivmax == 0.0 ) {
	*pivrow = lsub_ptr[pivptr];
	perm_r[*pivrow] = jcol;
	inv_perm_r[jcol] = *pivrow;
	*usepr = NO;
	return (jcol+1);
    }

    thresh = u * pivmax;
    
    /* Choose appropriate pivotal element by our policy. */
    if ( *usepr == YES ) {
	temp = lu_col_ptr[old_pivptr];
	if ( temp != 0.0 && fabs(temp) >= thresh )
	    pivptr = old_pivptr;
	else
	    *usepr = NO;
    }
    if ( *usepr == NO ) {
	/* Can we use diagonal as pivot? */
	if ( diag >= 0 ) { /* diagonal exists */
	    temp = lu_col_ptr[diag];
	    if ( temp != 0.0 && fabs(temp) >= thresh ) pivptr = diag;
	}
	*pivrow = lsub_ptr[pivptr];
    }
    
    /* Record pivot row */
    perm_r[*pivrow] = jcol;
    inv_perm_r[jcol] = *pivrow;
    
    /* Interchange row subscripts */
    if ( pivptr != nsupc ) {
	itemp = lsub_ptr[pivptr];
	lsub_ptr[pivptr] = lsub_ptr[nsupc];
	lsub_ptr[nsupc] = itemp;

	/* Interchange numerical values as well, for the whole supernode,
	 * such that L is indexed the same way as A.
 	 */
	k = 0;
	for (icol = 0; icol <= nsupc; ++icol, k += nsupr) {
	    itemp = pivptr + k;
	    temp = lu_sup_ptr[itemp];
	    lu_sup_ptr[itemp] = lu_sup_ptr[nsupc + k];
	    lu_sup_ptr[nsupc + k] = temp;
	}
    } /* if */

    
    /* CDIV operation */
    Gstat->procstat[pnum].fcops += nsupr - nsupc;
/*    ops[FACT] += nsupr - nsupc;*/
    temp = 1.0 / lu_col_ptr[nsupc];
    for (k = nsupc+1; k < nsupr; k++) 
	lu_col_ptr[k] *= temp;

#ifdef CHK_PIVOT
    printf("After cdiv: col %d\n", jcol);
    for (k = nsupc; k < nsupr; k++) 
	printf("  lu[%d] %f\n", lsub_ptr[k], lu_col_ptr[k]);
#endif

    return 0;
}

