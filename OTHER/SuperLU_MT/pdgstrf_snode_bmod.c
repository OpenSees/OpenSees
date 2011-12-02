#include "pdsp_defs.h"

int
pdgstrf_snode_bmod(
		   const int  pnum,   /* process number */
		   const int  jcol,   /* in - current column in the s-node */
		   const int  jsupno, /* in */
		   const int  fsupc,  /* in - first column in the s-node */
		   double     *dense, /* in */
		   double     *tempv, /* working array */
		   GlobalLU_t *Glu,   /* modified */
		   Gstat_t *Gstat     /* modified */
		   )
{
/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Performs numeric block updates within the relaxed supernode. 
 */
#if ( MACH==CRAY_PVP )
    _fcd ftcs1, ftcs2, ftcs3;
#endif
#ifdef USE_VENDOR_BLAS    
    int            incx = 1, incy = 1;
    double         alpha = -1.0, beta = 1.0;
#endif
    
    int            luptr, nsupc, nsupr, nrow;
    int            isub, irow, i, iptr; 
    register int   ufirst, nextlu;
    double         *lusup;
    int            *lsub, *xlsub, *xlsub_end, *xlusup, *xlusup_end;
    register float flopcnt;

    lsub       = Glu->lsub;
    xlsub      = Glu->xlsub;
    xlsub_end  = Glu->xlsub_end;
    lusup      = Glu->lusup;
    xlusup     = Glu->xlusup;
    xlusup_end = Glu->xlusup_end;

    nextlu = xlusup[jcol];
    
    /*
     *	Process the supernodal portion of L\U[*,j]
     */
    for (isub = xlsub[fsupc]; isub < xlsub_end[fsupc]; isub++) {
  	irow = lsub[isub];
	lusup[nextlu] = dense[irow];
	dense[irow] = 0;
	++nextlu;
    }

    xlusup_end[jcol] = nextlu;
    
    if ( fsupc < jcol ) {

	luptr = xlusup[fsupc];
	nsupr = xlsub_end[fsupc] - xlsub[fsupc];
	nsupc = jcol - fsupc;	/* Excluding jcol */
	ufirst = xlusup[jcol];	/* Points to the beginning of column
				   jcol in supernode L\U(jsupno). */
	nrow = nsupr - nsupc;
	
	flopcnt = nsupc * (nsupc - 1) + 2 * nrow * nsupc;
	Gstat->procstat[pnum].fcops += flopcnt;

/*	ops[TRSV] += nsupc * (nsupc - 1);
	ops[GEMV] += 2 * nrow * nsupc;    */

#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
	ftcs1 = _cptofcd("L", strlen("L"));
	ftcs2 = _cptofcd("N", strlen("N"));
	ftcs3 = _cptofcd("U", strlen("U"));
	STRSV( ftcs1, ftcs2, ftcs3, &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	SGEMV( ftcs2, &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
	      &lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#else
	dtrsv_( "L", "N", "U", &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	dgemv_( "N", &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
		&lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#endif
#else
	dlsolve ( nsupr, nsupc, &lusup[luptr], &lusup[ufirst] );
	dmatvec ( nsupr, nrow, nsupc, &lusup[luptr+nsupc], 
		 &lusup[ufirst], &tempv[0] );

        /* Scatter tempv[*] into lusup[*] */
	iptr = ufirst + nsupc;
	for (i = 0; i < nrow; i++) {
	    lusup[iptr++] -= tempv[i];
	    tempv[i] = 0.0;
	}
#endif

    }

    return 0;
}
