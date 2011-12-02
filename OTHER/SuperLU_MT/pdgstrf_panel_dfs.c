#include "pdsp_defs.h"

void
pdgstrf_panel_dfs(
		  const int  pnum,  /* process number */
		  const int  m,     /* number of rows in the matrix */
		  const int  w,     /* current panel width */
		  const int  jcol,  /* leading column of the current panel */
		  SuperMatrix *A,   /* original matrix */
		  int *perm_r, /* row pivotings that are done so far */
		  int *xprune, /* in */
		  int *ispruned,   /* in */
		  int *lbusy,      /* in; size n */
		  int *nseg,	   /* out */
		  int *panel_lsub, /* out */
		  int *w_lsub_end, /* out; values irrelevant on entry */
		  int *segrep,     /* out */
		  int *repfnz,     /* out */
		  int *marker,     /* modified */
		  int *spa_marker, /* modified; size n-by-w */
		  int        *parent,     /* working array */
		  int *xplore,     /* working array */
		  double *dense,      /* out; size n-by-w */
		  GlobalLU_t *Glu         /* modified */
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
 *   Performs a symbolic factorization on a panel of columns [jcol, jcol+w).
 *   It skips all those busy descendants that are worked on by other
 *   processors along the e-tree path.
 *
 * Notes
 * =====
 *
 * (1) panel_lsub[0:w*n-1]: temporary for the nonzero row indices below 
 *     the panel diagonal, which will be used later in the inner LU
 *     factorization. For the busy columns, some of the nonzeros in U
 *     may be mistakenly placed in this list, because "perm_r" is
 *     still "empty". Later, during dcolumn_dfs in the inner factorization,
 *     we must filter those nonzeros belonging in U.
 *
 * (2) A supernode representative is the last column of a supernode.
 *     The nonzeros in U[*,j] are segments that end at supernodal
 *     representatives.
 *
 * (3) The routine returns one list of the supernodal representatives
 *     in topological order of the DFS that generates them. This list is
 *     a superset of the topological order of each individual column within
 *     the panel. The location of the first nonzero in each supernodal
 *     segment (supernodal entry location) is also returned. Each column
 *     has a separate list for this purpose.
 *
 * (4) Two marker arrays are used to facilitate dfs:
 *     marker[i] == jj, if i was visited during dfs of current column jj;
 *     marker1[i] == jcol, if i was visited by earlier columns in this panel;
 *
 * (5) The dfs stack is the combination of xplore[2*m] and parent[m]:
 *     xplore[k]     - pointer to k's adjancency list where search begins
 *     xplore[m + k] - pointer to k's adjancency list where search ends
 *
 * (6) Array mappings
 *     marker: A-row --> A-row/col (0/1)
 *     repfnz: SuperA-col --> PA-row
 *     parent: SuperA-col --> SuperA-col
 *     xplore: SuperA-col --> index to L-structure
 *
 */
    NCPformat *Astore;
    double    *a;
    int       *asub;
    int       *xa_begin, *xa_end;
    register int krep, chperm, chmark, chrep, kchild, myfnz;
    register int k, krow, kmark, kperm, fsupc;
    register int xdfs, maxdfs, kpar, jj, nextp;
    register int nextl_col;/* next open position in panel_lsub[*,jj] */
    int       *marker1;	   /* marker1[jj] == jcol if vertex jj was visited 
			      by a previous column within this panel.   */
    int       *repfnz_col; /* start of each column in the panel */
    double    *dense_col;  /* start of each column in the panel */
    int       *xsup, *xsup_end, *supno, *lsub, *xlsub, *xlsub_end;

    int       *col_marker; /* marker array of each column in the panel */

    /* Initialize pointers */
    xsup       = Glu->xsup;
    xsup_end   = Glu->xsup_end;
    supno      = Glu->supno;
    lsub       = Glu->lsub;
    xlsub      = Glu->xlsub;
    xlsub_end  = Glu->xlsub_end;
    Astore     = A->Store;
    a          = Astore->nzval;
    asub       = Astore->rowind;
    xa_begin   = Astore->colbeg;
    xa_end     = Astore->colend;
    marker1    = marker + m;
    repfnz_col = repfnz;
    dense_col  = dense;
    nextp      = 0;
    *nseg      = 0;

#if ( DEBUGlevel>=2 )
if (jcol == BADPAN)    
    printf("(%d) pdgstrf_panel_dfs[begin] jcol %d, w %d\n", pnum, jcol, w);
#endif
    
    /*
     * For each column in the panel ...
     */
    for (jj = jcol; jj < jcol + w; ++jj, nextp += m) {
	nextl_col = nextp;
	col_marker = &spa_marker[nextp];

	/*
	 * For each nonz in A[*,jj] perform dfs ...
	 */
	for (k = xa_begin[jj]; k < xa_end[jj]; ++k) {
	    krow = asub[k];
	    dense_col[krow] = a[k];
	    kmark = col_marker[krow];
	    
	    /* if krow was visited before, go to the next nonzero */
	    if ( kmark == jj ) continue;

	    /*
	     * For each unmarked nbr krow of jj ...
	     */
	    col_marker[krow] = jj;
	    kperm = perm_r[krow];
	    
	    if ( kperm == EMPTY ) {
		/* krow is in L: place it in structure of L[*,jj].
		 * NOTE: some entries in U may get here, because "perm_r"
		 *       is not yet available from a preceeding busy column.
		 */
		panel_lsub[nextl_col++] = krow; /* krow is indexed into A */
	    } else {
		/* 
		 * krow is in U (0 <= kperm < jcol): if its supernode
		 * representative krep has been explored, update repfnz[*].
		 */
		if ( lbusy[kperm] == jcol ) { /* kperm is busy */
#if ( DEBUGlevel>=3 )
  if (jj == BADCOL)		    
    printf("(%d) pdgstrf_panel_dfs(%d) skip busy krow %d, kperm %d\n",
	   pnum, jj, krow, kperm);
#endif		    
		    continue;
		}

		/* Here, krep cannot possibly be "busy" */
		krep = SUPER_REP( supno[kperm] );
		myfnz = repfnz_col[krep];

#ifdef CHK_DFS
if (jj == BADCOL)		
    printf("(%d) pdgstrf_panel_dfs[1] %d, krep %d, fsupc %d, Pr[krow %d] %d, myfnz %d\n",
	   pnum, jj, krep, SUPER_FSUPC(supno[krep]), krow, kperm, myfnz);
#endif
		if ( myfnz != EMPTY ) {	/* Representative visited before */
		    if ( myfnz > kperm ) repfnz_col[krep] = kperm;
		    /* continue; */
		} else {
		    /* Otherwise, performs dfs starting from krep */
		    parent[krep] = EMPTY;
		    repfnz_col[krep] = kperm;
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
#ifdef CHK_DFS
if (jj == BADCOL)		    
{
    register int i;
    printf("(%d) pdgstrf_panel_dfs[2] %d, ispruned[%d] %d, xdfs %d, maxdfs %d\n",
	   pnum, jj, krep, ispruned[krep], xdfs, maxdfs);
    /*for (i = xdfs; i < maxdfs; i++) printf("(%d) lsub-%d", pnum, lsub[i]);*/
    printf("\n");
}
#endif
		    do {
			while ( xdfs < maxdfs ) {
			    /* for each unmarked kchild of krep ... */
			    kchild = lsub[xdfs];
			    xdfs++;
			    chmark = col_marker[kchild];
			    
			    if ( chmark != jj ) { /* Not reached yet */
				col_marker[kchild] = jj;
				chperm = perm_r[kchild];
				
				if ( chperm == EMPTY ) {
				    /* kchild is in L: place it in L[*,j]. */
				    panel_lsub[nextl_col++] = kchild;
				} else {
				    /* kchild is in U (0 <= chperm < jcol): 
				     * chrep = its supernode-rep. If its rep
				     * has been explored, update its repfnz[*].
				     */

				    if ( lbusy[chperm] == jcol ) {
#ifdef DEBUG
if (jj == BADCOL)					
    printf("(%d) pdgstrf_panel_dfs(%d) skip busy kchild %d, chperm %d\n",
	   pnum, jj, kchild, chperm);
#endif		    
	                                     continue;
                                    }
				    
				    chrep = SUPER_REP( supno[chperm] );
				    myfnz = repfnz_col[chrep];
#ifdef DEBUG
if (jj == BADCOL)				    
    printf("(%d) pdgstrf_panel_dfs[3] %d, krep %d, Pr[kchild %d] %d, chrep %d, fsupc %d, myfnz %d\n",
	   pnum, jj, krep, kchild, chperm, chrep,
	   SUPER_FSUPC(supno[chrep]), myfnz);
#endif
				    if ( myfnz != EMPTY ) {/* Visited before */
					if ( myfnz > chperm )
					    repfnz_col[chrep] = chperm;
				    } else {
					/* Cont. dfs at snode-rep of kchild */
					xplore[krep] = xdfs;
					xplore[m + krep] = maxdfs;
					parent[chrep] = krep;
					krep = chrep; /* Go deeper down G(L) */
					repfnz_col[krep] = chperm;
					if ( ispruned[krep] ) {
					    if ( SINGLETON( supno[krep] ) )
						xdfs = xlsub_end[krep];
					    else xdfs = xlsub[krep];
					    maxdfs = xprune[krep];
#ifdef PROFILE
					    procstat[pnum].pruned++;
#endif		    
					} else {
					    fsupc = SUPER_FSUPC(supno[krep]);
					    xdfs = xlsub[fsupc] + krep-fsupc+1;
					    maxdfs = xlsub_end[fsupc];
#ifdef PROFILE
					    procstat[pnum].unpruned++;
#endif		    
					}
#ifdef CHK_DFS
if (jj == BADCOL)
    printf("(%d) pdgstrf_panel_dfs[4] %d, ispruned[%d] %d, xdfs %d, maxdfs %d\n",
	   pnum, jj, krep, ispruned[krep], xdfs, maxdfs);
#endif
					
				    } /* else */
				} /* else */
			      
			    } /* if... */
			    
			} /* while xdfs < maxdfs */
			
			/* krow has no more unexplored nbrs:
			 *    Place snode-rep krep in postorder DFS, if this 
			 *    segment is seen for the first time. (Note that
			 *    "repfnz[krep]" may change later.)
			 *    Backtrack dfs to its parent.
			 */
			if ( marker1[krep] != jcol ) {
			    segrep[*nseg] = krep;
			    ++(*nseg);
			    marker1[krep] = jcol;
#ifdef CHK_DFS
if (jj == BADCOL)			    
    printf("(%d) pdgstrf_panel_dfs(%d) repfnz[%d] %d added to top.list by jj %d\n",
	   pnum, jj, krep, repfnz_col[krep], jj);
#endif			    
			}
			
			kpar = parent[krep]; /* Pop stack, mimic recursion */
			if ( kpar == EMPTY ) break; /* dfs done */
			krep = kpar;
			xdfs = xplore[krep];
			maxdfs = xplore[m + krep];
			
#ifdef CHK_DFS
if (jj == BADCOL)			
{
    register int i;
    printf("(%d) pdgstrf_panel_dfs[5] pop stack: %d, krep %d, xdfs %d, maxdfs %d\n",
	   pnum, jj, krep, xdfs, maxdfs);
    /* for (i = xdfs; i < maxdfs; i++) printf("(%d) lsub-%d", pnum, lsub[i]);*/
    printf("\n");
}
#endif

		    } while ( kpar != EMPTY ); /* until empty stack */
		    
		} /* else: myfnz == EMPTY */
		
	    } /* else: kperm != EMPTY */
	    
	} /* for each nonzero in A[*,jj] */

#if ( DEBUGlevel>=3 )
if (jj == BADCOL) {
#define REPCOL 0    
    krep = REPCOL;
    printf("(%d) pdgstrf_panel_dfs(end) w_lsub_end[jj=%d] %d, repfnz_col[%d] %d\n",
	   pnum, jj, nextl_col - nextp, krep, repfnz_col[krep]);
    PrintInt10("lsub", nextl_col - nextp, &panel_lsub[nextp]);
}
#endif
	
	w_lsub_end[jj-jcol] = nextl_col - nextp;
	repfnz_col += m;
        dense_col += m;
	
    } /* for jj ... */

}
