/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include "pdsp_defs.h"

#define XPAND_HINT(memtype, new_next, jcol, param) {\
fprintf(stderr, "Storage for %12s exceeded; Current column %d; Need at least %d;\n",\
	memtype, jcol, new_next); \
fprintf(stderr, "You may set it by the %d-th parameter in routine sp_ienv().\n", param); \
ABORT("Memory allocation failed"); \
}

/*
 * Set up pointers for integer working arrays.
 */
void
pxgstrf_SetIWork(int n, int panel_size, int *iworkptr, int **segrep,
		 int **parent, int **xplore, int **repfnz, int **panel_lsub,
		 int **marker, int **lbusy)
{
    *segrep = iworkptr;                                      /* n  */
    *parent = iworkptr + n;                                  /* n  */
    *xplore = iworkptr + 2*n;                                /* 2*n */
    *repfnz = iworkptr + 4*n;                                /* w*n */
    *panel_lsub = iworkptr + 4*n + panel_size*n;             /* w*n */
    *marker = iworkptr + 4*n + 2*panel_size*n;               /* 3*n */
    *lbusy  = iworkptr + (4+NO_MARKER)*n + 2*panel_size*n;   /* n   */
    ifill (*repfnz, n * panel_size, EMPTY);
}


void
copy_mem_int(int howmany, void *old, void *new)
{
    register int i;
    int *iold = old;
    int *inew = new;
    for (i = 0; i < howmany; i++) inew[i] = iold[i];
}


void
user_bcopy(char *src, char *dest, int bytes)
{
    char *s_ptr, *d_ptr;

    s_ptr = src + bytes - 1;
    d_ptr = dest + bytes - 1;
    for (; d_ptr >= dest; --s_ptr, --d_ptr ) *d_ptr = *s_ptr;
}



int *intMalloc(int n)
{
    int *buf;
    buf = (int *) SUPERLU_MALLOC(n * sizeof(int));
    if ( !buf ) {
	fprintf(stderr, "SUPERLU_MALLOC failed for buf in intMalloc()\n");
	exit (1);
    }
    return (buf);
}

int *intCalloc(int n)
{
    int *buf;
    register int i;
    buf = (int *) SUPERLU_MALLOC(n * sizeof(int));
    if ( !buf ) {
	fprintf(stderr, "SUPERLU_MALLOC failed for buf in intCalloc()\n");
	exit (1);
    }
    for (i = 0; i < n; ++i) buf[i] = 0;
    return (buf);
}

/*
 * Allocate n elements storage from a global array.
 * It uses lock for mutually exclusive access to the next position, so that
 * more than one processors can call aalloc on the same array correctly.
 * Return value: 0 - success
 *              >0 - number of bytes allocated when run out of space
 */
int
Glu_alloc(
	  const int pnum,     /* process number */
	  const int jcol,
	  const int num,      /* number of elements requested */
	  const MemType mem_type,
          int   *prev_next,   /* return "next" value before allocation */
	  pxgstrf_shared_t *pxgstrf_shared
	  )
{
    GlobalLU_t *Glu = pxgstrf_shared->Glu;
    register int fsupc, nextl, nextu, new_next;
#ifdef PROFILE
    double   t;
#endif
    
    switch ( mem_type ) {
	
      case LUSUP: 
	/* Storage for the H-supernode is already set aside, so we do
	   not use lock here. */
	if ( Glu->map_in_sup[jcol] < 0 )
	    fsupc = jcol + Glu->map_in_sup[jcol];
	else fsupc = jcol;
	*prev_next = Glu->map_in_sup[fsupc];
	Glu->map_in_sup[fsupc] += num;

#if 0
	{
	    register int i, j;
	    i = fsupc + part_super_h[fsupc];
	    if ( Glu->dynamic_snode_bound == YES )
	      new_next = Glu->nextlu;
	    else new_next = Glu->map_in_sup[i];
	    if (new_next < 0) { /* relaxed s-node */
		for (i = fsupc+1; Glu->map_in_sup[i] < 0; ++i) ;
		new_next = Glu->map_in_sup[i];
	    }
	    if ( Glu->map_in_sup[fsupc]>Glu->nzlumax 
		|| Glu->map_in_sup[fsupc]>new_next ) {
		printf("(%d) jcol %d, map_[%d]=%d, map_[%d]=new_next %d\n",
		       pnum, jcol, fsupc, Glu->map_in_sup[fsupc],
		       i, new_next);
		printf("(%d) snode type %d,size %d, |H-snode| %d\n",
		       pnum, Glu->pan_status[fsupc].type, 
		       Glu->pan_status[fsupc].size, part_super_h[fsupc]);
		for (j = fsupc; j < i; j += part_super_h[j])
		    printf("(%d) H snode %d, size %d\n",
			   pnum, j, part_super_h[j]);
		ABORT("LUSUP exceeded.");  /* xiaoye */
	    }
	}	    
#endif	
	break;

	
      case UCOL: case USUB:

#ifdef PROFILE
	t = SuperLU_timer_();
#endif
	
#if ( MACH==SUN )
	mutex_lock( &pxgstrf_shared->lu_locks[ULOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
	pthread_mutex_lock( &pxgstrf_shared->lu_locks[ULOCK] );
#elif ( MACH==SGI || MACH==ORIGIN )
#pragma critical lock( pxgstrf_shared->lu_locks[ULOCK] )
#elif ( MACH==CRAY_PVP )
#pragma _CRI guard ULOCK
#endif	
	{
	    nextu = Glu->nextu;
	    new_next = nextu + num;
	    if ( new_next > Glu->nzumax ) {
	        XPAND_HINT("U columns", new_next, jcol, 7);
	    }
	    *prev_next = nextu;
	    Glu->nextu = new_next;

	} /* end of critical region */
	
#if ( MACH==SUN )
	mutex_unlock( &pxgstrf_shared->lu_locks[ULOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
	pthread_mutex_unlock( &pxgstrf_shared->lu_locks[ULOCK] );
#elif ( MACH==CRAY_PVP )
#pragma _CRI endguard ULOCK
#endif

#ifdef PROFILE
	Gstat->procstat[pnum].cs_time += SuperLU_timer_() - t;
#endif
	
	break;

	
	case LSUB:

#ifdef PROFILE
	t = SuperLU_timer_();
#endif
	
#if ( MACH==SUN )
	mutex_lock( &pxgstrf_shared->lu_locks[LLOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
	pthread_mutex_lock( &pxgstrf_shared->lu_locks[LLOCK] );
#elif ( MACH==SGI || MACH==ORIGIN )
#pragma critical lock( pxgstrf_shared->lu_locks[LLOCK] )
#elif ( MACH==CRAY_PVP )
#pragma _CRI guard LLOCK
#endif	
	{
	  nextl = Glu->nextl;
	  new_next = nextl + num;
	  if ( new_next > Glu->nzlmax ) {
	      XPAND_HINT("L subscripts", new_next, jcol, 8);
	  }
	  *prev_next = nextl;
	  Glu->nextl = new_next;
	  
	} /* end of #pragama critical lock() */
	
#if ( MACH==SUN )
	mutex_unlock( &pxgstrf_shared->lu_locks[LLOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
	pthread_mutex_unlock( &pxgstrf_shared->lu_locks[LLOCK]);
#elif ( MACH==CRAY_PVP )
#pragma _CRI endguard LLOCK
#endif	

#ifdef PROFILE
	Gstat->procstat[pnum].cs_time += SuperLU_timer_() - t;
#endif
	
	  break;

    }
    
    return 0;
}

/*
 * Set up memory image in lusup[*], using the supernode boundaries in 
 * the Householder matrix.
 * 
 * In both static and dynamic scheme, the relaxed supernodes (leaves) 
 * are stored in the beginning of lusup[*]. In the static scheme, the
 * memory is also set aside for the internal supernodes using upper
 * bound information from H. In the dynamic scheme, however, the memory
 * for the internal supernodes is not allocated by this routine.
 *
 * Return value
 *   o Static scheme: number of nonzeros of all the supernodes in H.
 *   o Dynamic scheme: number of nonzeros of the relaxed supernodes. 
 */
int
PresetMap(
	  const int n,
	  SuperMatrix *A, /* original matrix permuted by columns */
	  pxgstrf_relax_t *pxgstrf_relax, /* relaxed supernodes */
	  pdgstrf_options_t *pdgstrf_options, /* input */
	  GlobalLU_t *Glu /* modified */
	  )
{
    register int i, j, k, w, rs, rs_lastcol, krow, kmark, maxsup, nextpos;
    register int rs_nrow; /* number of nonzero rows in a relaxed supernode */
    int          *marker, *asub, *xa_begin, *xa_end;
    NCPformat    *Astore;
    int *map_in_sup; /* memory mapping function; values irrelevant on entry. */
    int *colcnt;     /* column count of Lc or H */
    int *super_bnd;  /* supernodes partition in H */
    char *snode_env, *getenv();

    snode_env = getenv("SuperLU_DYNAMIC_SNODE_STORE");
    if ( snode_env != NULL ) {
	Glu->dynamic_snode_bound = YES;
#if ( PRNTlevel>=1 )
	printf(".. Use dynamic alg. to allocate storage for L supernodes.\n");
#endif
    } else  Glu->dynamic_snode_bound = NO;

    Astore   = A->Store;
    asub     = Astore->rowind;
    xa_begin = Astore->colbeg;
    xa_end   = Astore->colend;
    rs       = 1;
    marker   = intMalloc(n);
    ifill(marker, n, EMPTY);
    map_in_sup = Glu->map_in_sup = intCalloc(n+1);
    colcnt = pdgstrf_options->colcnt_h;
    super_bnd = pdgstrf_options->part_super_h;
    nextpos = 0;

    /* Split large supernode into smaller pieces */
    maxsup = sp_ienv(3);
    for (j = 0; j < n; ) {
	w = super_bnd[j];
	k = j + w;
	if ( w > maxsup ) {
	    w = w % maxsup;
	    if ( w == 0 ) w = maxsup;
	    while ( j < k ) {
		super_bnd[j] = w;
		j += w;
		w = maxsup;
	    }
	}
	j = k;
    }
    
    for (j = 0; j < n; j += w) {
        if ( Glu->dynamic_snode_bound == NO ) map_in_sup[j] = nextpos;

	if ( pxgstrf_relax[rs].fcol == j ) {
	    /* Column j starts a relaxed supernode. */
	    map_in_sup[j] = nextpos;
	    rs_nrow = 0;
	    w = pxgstrf_relax[rs++].size;
	    rs_lastcol = j + w;
	    for (i = j; i < rs_lastcol; ++i) {
		/* for each nonzero in A[*,i] */
		for (k = xa_begin[i]; k < xa_end[i]; k++) {	
		    krow = asub[k];
		    kmark = marker[krow];
		    if ( kmark != j ) { /* first time visit krow */
			marker[krow] = j;
			++rs_nrow;
		    }
		}
	    }
	    nextpos += w * rs_nrow;
	    
	    /* Find the next H-supernode, with leading column i, which is
	       outside the relaxed supernode, rs. */
	    for (i = j; i < rs_lastcol; k = i, i += super_bnd[i]);
	    if ( i > rs_lastcol ) {
		/* The w columns [rs_lastcol, i) may join in the
		   preceeding relaxed supernode; make sure we leave
		   enough room for the combined supernode. */
		w = i - rs_lastcol;
		nextpos += w * MAX(rs_nrow, colcnt[k]);
	    }
	    w = i - j;
	} else { /* Column j starts a supernode in H */
	    w = super_bnd[j];
	    if ( Glu->dynamic_snode_bound == NO ) nextpos += w * colcnt[j];
	}

	/* Set up the offset (negative) to the leading column j of a
	   supernode in H */ 
	for (i = 1; i < w; ++i) map_in_sup[j + i] = -i;
	
    } /* for j ... */

    if ( Glu->dynamic_snode_bound == YES ) Glu->nextlu = nextpos;
    else map_in_sup[n] = nextpos;

#if ( PRNTlevel>=1 )
    printf("** PresetMap() allocates %d reals to lusup[*]....\n", nextpos);
#endif

    free (marker);
    return nextpos;
}

/*
 * Dynamically set up storage image in lusup[*], using the supernode
 * boundaries in H.
 */
int
DynamicSetMap(
	      const int pnum,      /* process number */
	      const int jcol,      /* leading column of the s-node in H */
	      const int num,       /* number of elements requested */
	      pxgstrf_shared_t *pxgstrf_shared
	      )
{
    GlobalLU_t *Glu = pxgstrf_shared->Glu;
    register int nextlu, new_next;
    int *map_in_sup = Glu->map_in_sup; /* modified; memory mapping function */
    
#ifdef PROFILE
    double t = SuperLU_timer_();
#endif

#if ( MACH==SUN )
    mutex_lock( &pxgstrf_shared->lu_locks[LULOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_lock( &pxgstrf_shared->lu_locks[LULOCK] );
#elif ( MACH==SGI || MACH==ORIGIN )
#pragma critical lock ( pxgstrf_shared->lu_locks[LULOCK] )
#elif ( MACH==CRAY_PVP )
#pragma _CRI guard LULOCK
#endif
    {
	nextlu = Glu->nextlu;
	map_in_sup[jcol] = nextlu;
	new_next = nextlu + num;
	if ( new_next > Glu->nzlumax ) {
	    XPAND_HINT("L supernodes", new_next, jcol, 6);
	}
	Glu->nextlu = new_next;
    } /* end of critical region */

#if ( MACH==SUN )
    mutex_unlock( &pxgstrf_shared->lu_locks[LULOCK] );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_unlock( &pxgstrf_shared->lu_locks[LULOCK] );
#elif ( MACH==CRAY_PVP )
#pragma _CRI endguard LULOCK
#endif

#ifdef PROFILE
    Gstat->procstat[pnum].cs_time += SuperLU_timer_() - t;
#endif

    return 0;
}

