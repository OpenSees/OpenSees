#include <stdio.h>
#include <stdlib.h>
#include "pdsp_defs.h"
#include "util.h"


void dlsolve(int, int, double *, double *);
void dmatvec(int, int, int, double *, double *, double *);
void dmatvec2(int, int, int, double*, double*, double*, double*, double*);

void
pdgstrf_bmod1D_mv2(
		   const int pnum, /* process number */
		   const int n,    /* number of rows in the matrix */
		   const int w,    /* current panel width */
		   const int jcol, /* leading column of the current panel */
		   const int fsupc,/* leading column of the updating s-node */ 
		   const int krep, /* last column of the updating s-node */ 
		   const int nsupc,/* number of columns in the updating s-node */ 
		   int nsupr, /* number of rows in the updating supernode */  
		   int nrow,  /* number of rows below the diagonal block of
				 the updating supernode */ 
		   int *repfnz,    /* in */
		   int *panel_lsub,/* modified */
		   int *w_lsub_end,/* modified */
		   int *spa_marker,/* modified; size n-by-w */
		   double *dense,  /* modified */
		   double *tempv,  /* working array - zeros on entry/exit */
		   GlobalLU_t *Glu,/* modified */
		   Gstat_t *Gstat  /* modified */
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
 *
 *    Performs numeric block updates (sup-panel) in topological order.
 *    It features: col-col, 2cols-col, 3cols-col, and sup-col updates.
 *    Results are returned in SPA dense[*,w].
 *
 */
    double  zero = 0.0;

#if ( MACH==CRAY_PVP )
    _fcd ftcs1 = _cptofcd("L", strlen("L")),
         ftcs2 = _cptofcd("N", strlen("N")),
         ftcs3 = _cptofcd("U", strlen("U"));
#endif
#ifdef USE_VENDOR_BLAS
    int          incx = 1, incy = 1;
    double       alpha = 1.0, beta = zero;
#endif

    double       ukj, ukj1, ukj2;
    int          luptr, luptr1, luptr2;
    int          segsze;
    register int lptr; /* start of row subscripts of the updating supernode */
    register int i, j, kfnz, krep_ind, isub, irow, no_zeros, twocols;
    register int jj;	       /* index through each column in the panel */
    int      kfnz2[2], jj2[2]; /* detect two identical columns */
    int  *repfnz_col, *repfnz_col1; /* repfnz[] for a column in the panel */
    double *dense_col, *dense_col1;  /* dense[] for a column in the panel */
    double *tri[2], *matvec[2];
    int  *col_marker, *col_marker1; /* each column of the spa_marker[*,w] */
    int  *col_lsub, *col_lsub1;   /* each column of the panel_lsub[*,w] */
    int          *lsub, *xlsub_end;
    double       *lusup;
    int          *xlusup;
    register float flopcnt;
    
#ifdef TIMING
    double *utime = Gstat->utime;
    double f_time;
#endif    
    
    lsub      = Glu->lsub;
    xlsub_end = Glu->xlsub_end;
    lusup     = Glu->lusup;
    xlusup    = Glu->xlusup;
    lptr      = Glu->xlsub[fsupc];
    krep_ind  = lptr + nsupc - 1;
    twocols = 0;
    tri[0] = tempv;
    tri[1] = tempv + n;

#ifdef DEBUG
if (jcol == BADPAN && krep == BADREP) {
    printf("(%d) dbmod1D[1] jcol %d, fsupc %d, krep %d, nsupc %d, nsupr %d, nrow %d\n",
	   pnum, jcol, fsupc, krep, nsupc, nsupr, nrow);
    PrintInt10("lsub[xlsub[2774]", nsupr, &lsub[lptr]);
}    
#endif
    
    /* -----------------------------------------------
     * Sequence through each column in the panel ...
     * ----------------------------------------------- */
    repfnz_col= repfnz;
    dense_col = dense;
    col_marker= spa_marker;
    col_lsub  = panel_lsub;
    
    for (jj = jcol; jj < jcol + w; ++jj, col_marker += n, col_lsub += n,
	 repfnz_col += n, dense_col += n) {
	
	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue;	/* skip any zero segment */

	segsze = krep - kfnz + 1;
	luptr = xlusup[fsupc];

	flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
	Gstat->procstat[pnum].fcops += flopcnt;

	/* Case 1: Update U-segment of size 1 -- col-col update */
	if ( segsze == 1 ) {
#ifdef TIMING
	    f_time = SuperLU_timer_();
#endif	    
	    ukj = dense_col[lsub[krep_ind]];
	    luptr += nsupr*(nsupc-1) + nsupc;
#ifdef DEBUG
if (krep == BADCOL && jj == -1) {
    printf("(%d) dbmod1D[segsze=1]: k %d, j %d, ukj %.10e\n",
	   pnum, lsub[krep_ind], jj, ukj);
    PrintInt10("segsze=1", nsupr, &lsub[lptr]);
}
#endif	    
	    for (i = lptr + nsupc; i < xlsub_end[fsupc]; i++) {
		irow = lsub[i];
		dense_col[irow] -= ukj * lusup[luptr];
		++luptr;
#ifdef SCATTER_FOUND		
		if ( col_marker[irow] != jj ) {
		    col_marker[irow] = jj;
		    col_lsub[w_lsub_end[jj-jcol]++] = irow;
		}
#endif		
	    }
#ifdef TIMING
	    utime[FLOAT] += SuperLU_timer_() - f_time;
#endif
	    
	} else if ( segsze <= 3 ) {
	    
#ifdef TIMING
	    f_time = SuperLU_timer_();
#endif	    
	    ukj = dense_col[lsub[krep_ind]];
	    luptr += nsupr*(nsupc-1) + nsupc-1;
	    ukj1 = dense_col[lsub[krep_ind - 1]];
	    luptr1 = luptr - nsupr;
	    if ( segsze == 2 ) {
		ukj -= ukj1 * lusup[luptr1];
		dense_col[lsub[krep_ind]] = ukj;
/*#pragma ivdep*/
		for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    irow = lsub[i];
		    ++luptr;  ++luptr1;
		    dense_col[irow] -= (ukj * lusup[luptr]
					+ ukj1 * lusup[luptr1]);
#ifdef SCATTER_FOUND		
		    if ( col_marker[irow] != jj ) {
			col_marker[irow] = jj;
			col_lsub[w_lsub_end[jj-jcol]++] = irow;
		    }
#endif		
		}
	    } else {
		ukj2 = dense_col[lsub[krep_ind - 2]];
		luptr2 = luptr1 - nsupr;
		ukj1 -= ukj2 * lusup[luptr2-1];
		ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
		dense_col[lsub[krep_ind]] = ukj;
		dense_col[lsub[krep_ind-1]] = ukj1;
		for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    irow = lsub[i];
		    ++luptr; ++luptr1; ++luptr2;
		    dense_col[irow] -= (ukj * lusup[luptr]
				+ ukj1*lusup[luptr1] + ukj2*lusup[luptr2]);
#ifdef SCATTER_FOUND		
		    if ( col_marker[irow] != jj ) {
			col_marker[irow] = jj;
			col_lsub[w_lsub_end[jj-jcol]++] = irow;
		    }
#endif		
		}
	    }
#ifdef TIMING
	    utime[FLOAT] += SuperLU_timer_() - f_time;
#endif
	} else { /* segsze >= 4 */
	    if ( twocols == 1 ) {
		jj2[1] = jj; /* got two columns */
		twocols = 0;
		
		for (j = 0; j < 2; ++j) { /* Do two tri-solves */
		    i = n * (jj2[j] - jcol);
		    repfnz_col1 = &repfnz[i];
		    dense_col1  = &dense[i];
		    kfnz2[j] = repfnz_col1[krep];		    
		    no_zeros = kfnz2[j] - fsupc;
		    segsze = krep - kfnz2[j] + 1;
		    matvec[j] = tri[j] + segsze;

		    /* Gather U[*,j] segment from dense[*] to tri[*]. */
		    isub = lptr + no_zeros;
		    for (i = 0; i < segsze; ++i) {
			irow = lsub[isub];
			tri[j][i] = dense_col1[irow]; /* Gather */
			++isub;
		    }

#ifdef TIMING
		    f_time = SuperLU_timer_();
#endif
		    /* start effective triangle */
		    luptr = xlusup[fsupc] + nsupr * no_zeros + no_zeros;
		    
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
		    STRSV( ftcs1, ftcs2, ftcs3, &segsze, &lusup[luptr], 
			   &nsupr, tri[j], &incx );
#else
		    dtrsv_( "L", "N", "U", &segsze, &lusup[luptr], 
			   &nsupr, tri[j], &incx );
#endif
#else
		    dlsolve ( nsupr, segsze, &lusup[luptr], tri[j] );
		    
#endif
#ifdef TIMING
		    utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    
		} /* end for j ... two tri-solves */

#ifdef TIMING
		f_time = SuperLU_timer_();
#endif
		
		if ( kfnz2[0] < kfnz2[1] ) { /* First column is bigger */
		    no_zeros = kfnz2[0] - fsupc;
		    segsze = kfnz2[1] - kfnz2[0];
		    luptr = xlusup[fsupc] + nsupr * no_zeros + nsupc;
#ifdef USE_VENDOR_BLAS		    
#if ( MACH==CRAY_PVP )
		    SGEMV( ftcs2, &nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr, tri[0], &incx, &beta, matvec[0], &incy );
#else
		    dgemv_( "N", &nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr, tri[0], &incx, &beta, matvec[0], &incy );
#endif
#else
		    dmatvec (nsupr, nrow, segsze, &lusup[luptr],
			     tri[0], matvec[0]);
#endif
		    
		} else if ( kfnz2[0] > kfnz2[1] ) {
		    no_zeros = kfnz2[1] - fsupc;
		    segsze = kfnz2[0] - kfnz2[1];
		    luptr = xlusup[fsupc] + nsupr * no_zeros + nsupc;
#ifdef USE_VENDOR_BLAS		    
#if ( MACH==CRAY_PVP )
		    SGEMV( ftcs2, &nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr, tri[1], &incx, &beta, matvec[1], &incy );
#else
		    dgemv_( "N", &nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr, tri[1], &incx, &beta, matvec[1], &incy );
#endif
#else
		    dmatvec (nsupr, nrow, segsze, &lusup[luptr],
			     tri[1], matvec[1]);
#endif
		}
		
		/* Do matrix-vector multiply with two destinations */
		kfnz = MAX( kfnz2[0], kfnz2[1] );
		no_zeros = kfnz - fsupc;
		segsze = krep - kfnz + 1;
		luptr = xlusup[fsupc] + nsupr * no_zeros + nsupc;
#if ( MACH==DEC )
		dgemv2_(&nsupr, &nrow, &segsze, &lusup[luptr],
			&tri[0][kfnz-kfnz2[0]], &tri[1][kfnz-kfnz2[1]],
			matvec[0], matvec[1]);
		/*#elif ( MACH==CRAY_PVP )
		DGEMV2(&nsupr, &nrow, &segsze, &lusup[luptr],
		       &tri[0][kfnz-kfnz2[0]], &tri[1][kfnz-kfnz2[1]],
		       matvec[0], matvec[1]);*/
#else
		dmatvec2(nsupr, nrow, segsze, &lusup[luptr],
			 &tri[0][kfnz-kfnz2[0]], &tri[1][kfnz-kfnz2[1]],
			 matvec[0], matvec[1]);
#endif

#ifdef TIMING
		utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    

		for (j = 0; j < 2; ++j) {
		    i = n * (jj2[j] - jcol);
		    dense_col1  = &dense[i];
		    col_marker1 = &spa_marker[i];
		    col_lsub1   = &panel_lsub[i];
		    no_zeros = kfnz2[j] - fsupc;
		    segsze = krep - kfnz2[j] + 1;
		    
		    /* Scatter tri[*] into SPA dense[*]. */
		    isub = lptr + no_zeros;
		    for (i = 0; i < segsze; i++) {
			irow = lsub[isub];
			dense_col1[irow] = tri[j][i]; /* Scatter */
			tri[j][i] = zero;
			++isub;
#ifdef DEBUG
	if (jj == -1 && krep == 3423)
	    printf("(%d) dbmod1D[scatter] jj %d, dense_col[%d] %e\n",
		   pnum, jj, irow, dense_col[irow]);
#endif
		    }
		    
		    /* Scatter matvec[*] into SPA dense[*]. */
/*#pragma ivdep*/
		    for (i = 0; i < nrow; i++) {
			irow = lsub[isub];
			dense_col1[irow] -= matvec[j][i]; /* Scatter-add */
#ifdef SCATTER_FOUND		
			if ( col_marker1[irow] != jj2[j] ) {
			    col_marker1[irow] = jj2[j];
			    col_lsub1[w_lsub_end[jj2[j]-jcol]++] = irow;
			}
#endif		
			matvec[j][i] = zero;
			++isub;
		    }
		    
		} /* end for two destination update */
		
	    } else { /* wait for a second column */
		jj2[0] = jj;
		twocols = 1;
	    }
	} /* else segsze >= 4 */
	
    } /* for jj ... */

    
    if ( twocols == 1 ) { /* one more column left */
	i = n * (jj2[0] - jcol);
	repfnz_col1 = &repfnz[i];
	dense_col1  = &dense[i];
	col_marker1 = &spa_marker[i];
	col_lsub1   = &panel_lsub[i];
	kfnz = repfnz_col1[krep];		    
	no_zeros = kfnz - fsupc;
	segsze = krep - kfnz + 1;

	/* Gather U[*,j] segment from dense[*] to tri[*]. */
	isub = lptr + no_zeros;
	for (i = 0; i < segsze; ++i) {
	    irow = lsub[isub];
	    tri[0][i] = dense_col1[irow]; /* Gather */
	    ++isub;
	}

#ifdef TIMING
	f_time = SuperLU_timer_();
#endif
	/* start effective triangle */
	luptr = xlusup[fsupc] + nsupr * no_zeros + no_zeros;
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
	STRSV( ftcs1, ftcs2, ftcs3, &segsze, &lusup[luptr], 
	       &nsupr, tri[0], &incx );
#else
	dtrsv_( "L", "N", "U", &segsze, &lusup[luptr], 
	       &nsupr, tri[0], &incx );
#endif
#else
	dlsolve ( nsupr, segsze, &lusup[luptr], tri[0] );
#endif
	
	luptr += segsze;	/* Dense matrix-vector */
	matvec[0] = tri[0] + segsze;
		
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
	SGEMV( ftcs2, &nrow, &segsze, &alpha, &lusup[luptr], 
	       &nsupr, tri[0], &incx, &beta, matvec[0], &incy );
#else
	dgemv_( "N", &nrow, &segsze, &alpha, &lusup[luptr], 
	       &nsupr, tri[0], &incx, &beta, matvec[0], &incy );
#endif
#else
	dmatvec (nsupr, nrow, segsze, &lusup[luptr], tri[0], matvec[0]);
#endif
#ifdef TIMING
	utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    

	/* Scatter tri[*] into SPA dense[*]. */
	isub = lptr + no_zeros;
	for (i = 0; i < segsze; i++) {
	    irow = lsub[isub];
	    dense_col1[irow] = tri[0][i]; /* Scatter */
	    tri[0][i] = zero;
	    ++isub;
#ifdef DEBUG
	if (jj == -1 && krep == 3423)
	    printf("(%d) dbmod1D[scatter] jj %d, dense_col[%d] %e\n",
		   pnum, jj, irow, dense_col[irow]);
#endif
	}
		    
	/* Scatter matvec[*] into SPA dense[*]. */
	for (i = 0; i < nrow; i++) {
	    irow = lsub[isub];
	    dense_col1[irow] -= matvec[0][i]; /* Scatter-add */
#ifdef SCATTER_FOUND		
	    if ( col_marker1[irow] != jj2[0] ) {
		col_marker1[irow] = jj2[0];
		col_lsub1[w_lsub_end[jj2[0]-jcol]++] = irow;
	    }
#endif		
	    matvec[0][i] = zero;
	    ++isub;
	}
    
    } /* if twocols == 1 */

}


