#include "pdsp_defs.h"

int
pdgstrf_snode_dfs(
		  const int  pnum,      /* process number */
		  const int  jcol,	  /* in - start of the supernode */
		  const int  kcol, 	  /* in - end of the supernode */
		  const int  *asub,     /* in */
		  const int  *xa_begin, /* in */
		  const int  *xa_end,   /* in */
		  int        *xprune,   /* out */
		  int        *marker,   /* modified */
		  int        *col_lsub, /* values are irrelevant on entry 
					   and on return */
		  pxgstrf_shared_t *pxgstrf_shared /* modified */
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
 *    pdgstrf_snode_dfs() determines the union of the row structures of 
 *    those columns within the relaxed snode.
 *    Note: The relaxed snodes are leaves of the supernodal etree, 
 *    therefore, the portion outside the rectangular supernode must be zero.
 *
 * Return value
 * ============
 *     0   success;
 *    >0   number of bytes allocated when run out of memory.
 *
 */
    GlobalLU_t *Glu = pxgstrf_shared->Glu;
    register int i, k, ifrom, nextl, nsuper;
    int          ito;
    int          krow, kmark, mem_error;
    int          *supno, *lsub, *xlsub, *xlsub_end;
    
    supno                 = Glu->supno;
    xlsub                 = Glu->xlsub;
    xlsub_end             = Glu->xlsub_end;
    nsuper = NewNsuper(pnum, &pxgstrf_shared->lu_locks[NSUPER_LOCK],
		       &Glu->nsuper);
    Glu->xsup[nsuper]     = jcol;
    Glu->xsup_end[nsuper] = kcol + 1;
    
    nextl = 0;
    for (i = jcol; i <= kcol; i++) {
	/* for each nonzero in A[*,i] */
	for (k = xa_begin[i]; k < xa_end[i]; k++) {	
	    krow = asub[k];
	    kmark = marker[krow];
	    if ( kmark != kcol ) { /* First time visit krow */
		marker[krow] = kcol;
		col_lsub[nextl++] = krow;
	    }
    	}
	supno[i] = nsuper;
    }

    if ( (mem_error = Glu_alloc(pnum, jcol, 2*nextl, LSUB, &ito, 
				pxgstrf_shared)) )
	return mem_error;
    
    xlsub[jcol] = ito;
    lsub        = Glu->lsub;
    for (ifrom = 0; ifrom < nextl; ++ifrom)
	lsub[ito++] = col_lsub[ifrom];
    xlsub_end[jcol] = ito;

    return 0;
}
