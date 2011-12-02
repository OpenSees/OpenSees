#include "pdsp_defs.h"
#include "util.h"

int
pdgstrf_column_dfs(
		   const int  pnum,    /* process number */
		   const int  m,       /* number of rows in the matrix */
		   const int  jcol,    /* current column in the panel */
		   const int  fstcol,  /* first column in the panel */
		   int *perm_r,   /* row pivotings that are done so far */
		   int *ispruned, /* in */
		   int *col_lsub, /* the RHS vector to start the dfs */
		   int lsub_end,  /* size of col_lsub[] */
		   int *super_bnd,/* supernode partition by upper bound */
		   int *nseg,     /* modified - with new segments appended */
		   int *segrep,   /* modified - with new segments appended */
		   int *repfnz,   /* modified */
		   int *xprune,   /* modified */
		   int *marker2,  /* modified */
		   int *parent,   /* working array */
		   int *xplore,   /* working array */
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
 *   pdgstrf_column_dfs() performs a symbolic factorization on column jcol, 
 *   and detects whether column jcol belongs in the same supernode as jcol-1.
 *
 * Local parameters
 * ================
 *   A supernode representative is the last column of a supernode.
 *   The nonzeros in U[*,j] are segments that end at supernodal
 *   representatives. The routine returns a list of such supernodal 
 *   representatives in topological order of the dfs that generates them.
 *   The location of the first nonzero in each such supernodal segment
 *   (supernodal entry location) is also returned.
 *
 *   nseg: no of segments in current U[*,j]
 *   samesuper: samesuper=NO if column j does not belong in the same
 *	        supernode as j-1. Otherwise, samesuper=YES.
 *
 *   marker2: A-row --> A-row/col (0/1)
 *   repfnz: SuperA-col --> PA-row
 *   parent: SuperA-col --> SuperA-col
 *   xplore: SuperA-col --> index to L-structure
 *
 * Return value
 * ============
 *     0  success;
 *   > 0  number of bytes allocated when run out of space.
 *
 */
    GlobalLU_t *Glu = pxgstrf_shared->Glu; /* modified */
    register int jcolm1, jcolm1size, nextl, ifrom;
    register int k, krep, krow, kperm, samesuper, nsuper;
    register int no_lsub;
    int	    fsupc;		/* first column in a supernode */
    int     myfnz;		/* first nonz column in a U-segment */
    int	    chperm, chmark, chrep, kchild;
    int     xdfs, maxdfs, kpar;
    int     ito;	        /* Used to compress row subscripts */
    int     mem_error;
    int     *xsup, *xsup_end, *supno, *lsub, *xlsub, *xlsub_end;
    static  int  first = 1, maxsuper;

    if ( first ) {
	maxsuper = sp_ienv(3);
	first = 0;
    }

    /* Initialize pointers */
    xsup      = Glu->xsup;
    xsup_end  = Glu->xsup_end;
    supno     = Glu->supno;
    lsub      = Glu->lsub;
    xlsub     = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    jcolm1    = jcol - 1;
    nextl     = lsub_end;
    no_lsub   = 0;
    samesuper = YES;

    /* Test whether the row structure of column jcol is contained
       in that of column jcol-1. */
    for (k = 0; k < lsub_end; ++k) {
	krow = col_lsub[k];
	if ( perm_r[krow] == EMPTY ) { /* krow is in L */
	    ++no_lsub;
	    if (marker2[krow] != jcolm1) 
	        samesuper = NO; /* row subset test */
	    marker2[krow] = jcol;
	}
    }

#if ( DEBUGlevel>=2 )
  if (jcol == BADCOL)
    printf("(%d) pdgstrf_column_dfs[1] %d, fstcol %d, lsub_end %d, no_lsub %d, samesuper? %d\n",
	   pnum, jcol, fstcol, lsub_end, no_lsub, samesuper);
#endif
    
    /*
     * For each nonzero in A[fstcol:n,jcol] perform DFS ...
     */
    for (k = 0; k < lsub_end; ++k) {
	krow = col_lsub[k];
	
	/* if krow was visited before, go to the next nonzero */
	if ( marker2[krow] == jcol ) continue;
	marker2[krow] = jcol;
	kperm = perm_r[krow];
#if ( DEBUGlevel>=3 )
  if (jcol == BADCOL)
    printf("(%d) pdgstrf_column_dfs[inner]: perm_r[krow=%d] %d\n", pnum, krow, kperm);
#endif
	
	/* Ignore the nonzeros in U corresponding to the busy columns
	   during the panel DFS. */
	/*if ( lbusy[kperm] != fstcol ) {  xiaoye? */
	if ( kperm >= fstcol ) {
	    /*
	     * krow is in U: if its supernode representative krep
	     * has been explored, update repfnz[*].
	     */
	    krep = SUPER_REP(supno[kperm]);
	    myfnz = repfnz[krep];
	    
#if ( DEBUGlevel>=3 )
  if (jcol == BADCOL)
    printf("(%d) pdgstrf_column_dfs[inner-U]: krep %d, myfnz %d, kperm %d\n",
	   pnum, krep, myfnz, kperm);
#endif
	    if ( myfnz != EMPTY ) {	/* Visited before */
		if ( myfnz > kperm ) repfnz[krep] = kperm;
		/* continue; */
	    } else {
		/* Otherwise, perform dfs starting at krep */
		parent[krep] = EMPTY;
		repfnz[krep] = kperm;
		if ( ispruned[krep] ) {
		    if ( SINGLETON( supno[krep] ) )
			xdfs = xlsub_end[krep];
		    else xdfs = xlsub[krep];
		    maxdfs = xprune[krep];
#ifdef PROFILE
		    procstat[pnum].pruned++;
#endif		    
		} else {
		    fsupc = SUPER_FSUPC( supno[krep] );
		    xdfs = xlsub[fsupc] + krep-fsupc+1;
		    maxdfs = xlsub_end[fsupc];
#ifdef PROFILE
		    procstat[pnum].unpruned++;
#endif		    
		}
		
		do {
		    /* 
		     * For each unmarked kchild of krep ...
		     */
		    while ( xdfs < maxdfs ) {
			
			kchild = lsub[xdfs];
			xdfs++;
			chmark = marker2[kchild];
			
			if ( chmark != jcol ) { /* Not reached yet */
			    marker2[kchild] = jcol;
			    chperm = perm_r[kchild];
			    
			    if ( chperm == EMPTY ) {
				/* kchild is in L: place it in L[*,k]. */
				++no_lsub;
				col_lsub[nextl++] = kchild;
				if (chmark != jcolm1) samesuper = NO;
			    } else {
				/* kchild is in U: chrep = its supernode
				 * representative. If its rep has 
				 * been explored, update its repfnz[*].
				 */
				chrep = SUPER_REP( supno[chperm] );
				myfnz = repfnz[chrep];
				if ( myfnz != EMPTY ) { /* Visited before */
				    if ( myfnz > chperm )
					repfnz[chrep] = chperm;
				} else {
				    /* Continue dfs at super-rep of kchild */
				    xplore[krep] = xdfs;	
				    xplore[m + krep] = maxdfs;	
				    parent[chrep] = krep;
				    krep = chrep; /* Go deeper down G(L^t) */
				    repfnz[krep] = chperm;
				    if ( ispruned[krep] ) {
					if ( SINGLETON( supno[krep] ) )
					    xdfs = xlsub_end[krep];
					else xdfs = xlsub[krep];
					maxdfs = xprune[krep];
#ifdef PROFILE
					procstat[pnum].pruned++;
#endif		    
				    } else {
					fsupc = SUPER_FSUPC( supno[krep] );
					xdfs = xlsub[fsupc] + krep-fsupc+1;
					maxdfs = xlsub_end[fsupc];
#ifdef PROFILE
					procstat[pnum].unpruned++;
#endif		    
				    }
				}
			    } /* else */
			} /* if */
		    } /* while */
		    
		    /* krow has no more unexplored nbrs:
		     *    place supernode-rep krep in postorder DFS,
		     *    backtrack dfs to its parent.
		     */
		    segrep[*nseg] = krep;
		    ++(*nseg);
#if ( DEBUGlevel>=3 )
  if (jcol == BADCOL)
    printf("(%d) pdgstrf_column_dfs[inner-dfs] new nseg %d, repfnz[krep=%d] %d\n",
	   pnum, *nseg, krep, repfnz[krep]);
#endif
		    kpar = parent[krep]; /* Pop from stack, mimic recursion */
		    if ( kpar == EMPTY ) break; /* dfs done */
		    krep = kpar;
		    xdfs = xplore[krep];
		    maxdfs = xplore[m + krep];
		} while ( kpar != EMPTY ); /* Do ... until empty stack */
		
	    } /* else myfnz ... */
	} /* if kperm >= fstcol ... */
    } /* for each nonzero ... */
	
#if ( DEBUGlevel>=3 )
  if (jcol == BADCOL)
    printf("(%d) pdgstrf_column_dfs[2]: nextl %d, samesuper? %d\n",
	   pnum, nextl, samesuper);
#endif

    /* assert(no_lsub == nextl - no_usub);*/

    /* ---------------------------------------------------------
       Check to see if j belongs in the same supernode as j-1.
       --------------------------------------------------------- */
    
    /* Does it matter if jcol == 0? - xiaoye */
    if ( samesuper == YES ) {
	nsuper = supno[jcolm1];
	jcolm1size = xlsub_end[jcolm1] - xlsub[jcolm1];
#if ( DEBUGlevel>=3 )
  if (jcol == BADCOL)
    printf("(%d) pdgstrf_column_dfs[YES] jcol-1 %d, jcolm1size %d, supno[%d] %d\n",
	   pnum, jcolm1, jcolm1size, jcolm1, nsuper);
#endif	
	if ( no_lsub != jcolm1size-1 )
	    samesuper = NO;        /* Enforce T2 supernode */
	else {
	    /* Make sure the number of columns in a supernode does not
	       exceed threshold. */
	    fsupc = xsup[nsuper];
	    if ( jcol - fsupc >= maxsuper )
		samesuper = NO;
	    else {
		/* start of a supernode in H (coarser partition) */
		if ( super_bnd[jcol] != 0 ) samesuper = NO;
	    }
	}
    }
    
    /* If jcol starts a new supernode, allocate storage for 
     * the subscript set of both first and last column of
     * a previous supernode. (first for num values, last for pruning)
     */
    if ( samesuper == NO ) { /* starts a new supernode */
	nsuper = NewNsuper(pnum, &pxgstrf_shared->lu_locks[NSUPER_LOCK], 
			   &Glu->nsuper);
	xsup[nsuper] = jcol;
	
	/* Copy column jcol; also reserve space to store pruned graph */
	if ((mem_error = Glu_alloc(pnum, jcol, 2*no_lsub, LSUB, &ito, 
				  pxgstrf_shared)))
	    return mem_error;
	xlsub[jcol] = ito;
	lsub = Glu->lsub;
	for (ifrom = 0; ifrom < nextl; ++ifrom) {
	    krow = col_lsub[ifrom];
	    if ( perm_r[krow] == EMPTY ) /* Filter U-subscript */
		lsub[ito++] = krow;
	}
	k = ito;
	xlsub_end[jcol] = k;
	
	/* Make a copy in case it is a singleton supernode */
	for (ifrom = xlsub[jcol]; ifrom < ito; ++ifrom)
	    lsub[k++] = lsub[ifrom];
	
    } else { /* Supernode of size > 1: overwrite column jcol-1 */
	k = xlsub_end[fsupc];
	xlsub[jcol] = k;
	xprune[fsupc] = k;
	for (ifrom = 0; ifrom < nextl; ++ifrom) {
	    krow = col_lsub[ifrom];
	    if ( perm_r[krow] == EMPTY ) /* Filter U-subscript */
		lsub[k++] = krow;
	}
	xlsub_end[jcol] = k;
    }

#if ( DEBUGlevel>=3 )
  if (jcol == BADCOL) {
    printf("(%d) pdgstrf_column_dfs[3]: %d in prev s-node %d? %d\n",
	   pnum, jcol, fsupc, samesuper);
    PrintInt10("lsub", xlsub_end[jcol]-xlsub[jcol], &lsub[xlsub[jcol]]);
  }
#endif
    
    /* Tidy up the pointers before exit */
    xprune[jcol] = k;     /* upper bound for pruning */
    supno[jcol] = nsuper;
    xsup_end[nsuper] = jcol + 1;
    
    return 0;
}
