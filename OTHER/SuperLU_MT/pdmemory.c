/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
#include <stdio.h>
#include <malloc.h>
#include "pdsp_defs.h"

/* ------------------
   Constants & Macros
   ------------------ */
#define EXPAND      1.5
#define NO_MEMTYPE  4      /* 0: lusup;
			      1: ucol;
			      2: lsub;
			      3: usub */
#define GluIntArray(n)   (9 * (n) + 5)

/* -------------------
   Internal prototypes
   ------------------- */
void    *pxgstrf_expand (int *, MemType,int, int, GlobalLU_t *);
void    copy_mem_double (int, void *, void *);
void    pdgstrf_StackCompress(GlobalLU_t *);
void    pxgstrf_SetupSpace (void *, int);
void    *duser_malloc   (int, int);
void    duser_free      (int, int);

/* ----------------------------------------------
   External prototypes (in memory.c - prec-indep)
   ---------------------------------------------- */
extern void    copy_mem_int    (int, void *, void *);
extern void    user_bcopy      (char *, char *, int);

typedef struct {
    int  size;
    int  used;
    int  top1;  /* grow upward, relative to &array[0] */
    int  top2;  /* grow downward */
    void *array;
} LU_stack_t;

typedef enum {HEAD, TAIL}   stack_end_t;
typedef enum {SYSTEM, USER} LU_space_t;

ExpHeader *expanders; /* Array of pointers to 4 types of memory */
static LU_stack_t stack;
static int        no_expand;
static int        ndim;
static LU_space_t whichspace; /* 0 - system malloc'd; 1 - user provided */

/* Macros to manipulate stack */
#define StackFull(x)         ( x + stack.used >= stack.size )
#define NotDoubleAlign(addr) ( (long int)addr & 7 )
#define DoubleAlign(addr)    ( ((long int)addr + 7) & ~7L )

#define Reduce(alpha)        ((alpha + 1) / 2)     /* i.e. (alpha-1)/2 + 1 */

/* temporary space used by BLAS calls */
/*#define NUM_TEMPV(n,w,t,b)  (MAX(n, (t + b)*w))*/
#define NUM_TEMPV(n,w,t,b)  (MAX( 2*n, (t + b)*w ))


/*
 * Setup the memory model to be used for factorization.
 *    lwork = 0: use system malloc;
 *    lwork > 0: use user-supplied work[] space.
 */
void pxgstrf_SetupSpace(void *work, int lwork)
{
    if ( lwork == 0 ) {
	whichspace = SYSTEM; /* malloc/free */
    } else if ( lwork > 0 ) {
	whichspace = USER;   /* user provided space */
	stack.size = lwork;
	stack.used = 0;
	stack.top1 = 0;
	stack.top2 = lwork;
	stack.array = (void *) work;
    }
}


void *duser_malloc(int bytes, int which_end)
{
    void *buf;
    
    if ( StackFull(bytes) ) return (NULL);

    if ( which_end == HEAD ) {
	buf = (char*) stack.array + stack.top1;
	stack.top1 += bytes;
    } else {
	stack.top2 -= bytes;
	buf = (char*) stack.array + stack.top2;
    }
    
    stack.used += bytes;
    return buf;
}


void duser_free(int bytes, int which_end)
{
    if ( which_end == HEAD ) {
	stack.top1 -= bytes;
    } else {
	stack.top2 += bytes;
    }
    stack.used -= bytes;
}


/* Returns the working storage used during factorization */
int superlu_TempSpace(n, w, p)
{
    register float tmp, ptmp;
    register int iword = sizeof(int), dword = sizeof(double);
    int    maxsuper = sp_ienv(3),
           rowblk   = sp_ienv(4);

    /* globally shared */
    tmp = 14 * n * iword;

    /* local to each processor */
    ptmp = (2 * w + 5 + NO_MARKER) * n * iword;
    ptmp += (n * w + NUM_TEMPV(n,w,maxsuper,rowblk)) * dword;
#if ( PRNTlevel>=1 )
    printf("Per-processor work[] %.0f MB\n", ptmp/1024/1024);
#endif
    ptmp *= p;

    return (tmp + ptmp);
}


/*
 * superlu_memusage consists of the following fields:
 *    o for_lu (float)
 *      The amount of space used in bytes for L\U data structures.
 *    o total_needed (float)
 *      The amount of space needed in bytes to perform factorization.
 *    o expansions (int)
 *      The number of memory expansions during the LU factorization.
 */
int superlu_QuerySpace(int P, SuperMatrix *L, SuperMatrix *U, int panel_size,
		       superlu_memusage_t *superlu_memusage)
{
    SCPformat *Lstore;
    NCPformat *Ustore;
    register int n, iword, dword, lwork;

    Lstore = L->Store;
    Ustore = U->Store;
    n = L->ncol;
    iword = sizeof(int);
    dword = sizeof(double);

    /* L supernodes of type SCP */
    superlu_memusage->for_lu = (7*n + 3) * iword 
                             + Lstore->nzval_colend[n-1] * dword
	                     + Lstore->rowind_colend[n-1] * iword;

    /* U columns of type NCP */
    superlu_memusage->for_lu += (2*n + 1) * iword
	+ Ustore->colend[n-1] * (dword + iword);

    /* Working storage to support factorization */
    lwork = superlu_TempSpace(n, panel_size, P);
    superlu_memusage->total_needed = superlu_memusage->for_lu + lwork;

    superlu_memusage->expansions = --no_expand;

    return 0;
}


int pdgstrf_memory_use(const int nzlmax, const int nzumax, const int nzlumax)
{
    register int iword, dword;

    iword   = sizeof(int);
    dword   = sizeof(double);
    
    return (10 * ndim * iword +
	    nzlmax * iword + nzumax * (iword + dword) + nzlumax * dword);

}


/*
 * Allocate storage for the data structures common to all factor routines.
 * For those unpredictable size, make a guess as FILL * nnz(A).
 * Return value:
 *     If lwork = -1, return the estimated amount of space required;
 *     otherwise, return the amount of space actually allocated when
 *     memory allocation failure occurred.
 */
int
pdgstrf_MemInit(int n, int annz, pdgstrf_options_t *pdgstrf_options,
		SuperMatrix *L, SuperMatrix *U, GlobalLU_t *Glu)
{
    register int nprocs = pdgstrf_options->nprocs;
    yes_no_t refact = pdgstrf_options->refact;
    register int panel_size = pdgstrf_options->panel_size;
    register int lwork = pdgstrf_options->lwork;
    void     *work = pdgstrf_options->work;
    int      iword, dword, retries = 0;
    SCPformat *Lstore;
    NCPformat *Ustore;
    int      *xsup, *xsup_end, *supno;
    int      *lsub, *xlsub, *xlsub_end;
    double   *lusup;
    int      *xlusup, *xlusup_end;
    double   *ucol;
    int      *usub, *xusub, *xusub_end;
    int      nzlmax, nzumax, nzlumax;
    int      FILL_LUSUP = sp_ienv(6); /* Guess the fill-in growth for LUSUP */
    int      FILL_UCOL = sp_ienv(7); /* Guess the fill-in growth for UCOL */
    int      FILL_LSUB = sp_ienv(8); /* Guess the fill-in growth for LSUB */
    
    no_expand = 0;
    ndim      = n;
    iword     = sizeof(int);
    dword     = sizeof(double);

    expanders = (ExpHeader *) SUPERLU_MALLOC(NO_MEMTYPE * sizeof(ExpHeader));

    if ( refact == NO ) {

	/* Guess amount of storage needed by L\U factors. */
        if ( FILL_UCOL < 0 ) nzumax = -FILL_UCOL * annz;
	else nzumax = FILL_UCOL;
	if ( FILL_LSUB < 0 ) nzlmax = -FILL_LSUB * annz;
	else nzlmax = FILL_LSUB;

	if ( Glu->dynamic_snode_bound == YES ) {
	    if ( FILL_LUSUP < 0 ) nzlumax = -FILL_LUSUP * annz;
	    else nzlumax = FILL_LUSUP; /* estimate an upper bound */
	} else {
	    nzlumax = Glu->nzlumax; /* preset as static upper bound */
	}

	if ( lwork == -1 ) {
	    return (GluIntArray(n) * iword + 
		    superlu_TempSpace(n, panel_size, nprocs)
		    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword);
        } else {
	    pxgstrf_SetupSpace(work, lwork);
	}
	
	/* Integer pointers for L\U factors */
	if ( whichspace == SYSTEM ) {
	    xsup       = intMalloc(n+1);
	    xsup_end   = intMalloc(n);
	    supno      = intMalloc(n+1);
	    xlsub      = intMalloc(n+1);
	    xlsub_end  = intMalloc(n);
	    xlusup     = intMalloc(n+1);
	    xlusup_end = intMalloc(n);
	    xusub      = intMalloc(n+1);
	    xusub_end  = intMalloc(n);
	} else {
	    xsup       = (int *)duser_malloc((n+1) * iword, HEAD);
	    xsup_end   = (int *)duser_malloc((n) * iword, HEAD);
	    supno      = (int *)duser_malloc((n+1) * iword, HEAD);
	    xlsub      = (int *)duser_malloc((n+1) * iword, HEAD);
	    xlsub_end  = (int *)duser_malloc((n) * iword, HEAD);
	    xlusup     = (int *)duser_malloc((n+1) * iword, HEAD);
	    xlusup_end = (int *)duser_malloc((n) * iword, HEAD);
	    xusub      = (int *)duser_malloc((n+1) * iword, HEAD);
	    xusub_end  = (int *)duser_malloc((n) * iword, HEAD);
	}

	lusup = (double *) pxgstrf_expand( &nzlumax, LUSUP, 0, 0, Glu );
	ucol  = (double *) pxgstrf_expand( &nzumax, UCOL, 0, 0, Glu );
	lsub  = (int *)    pxgstrf_expand( &nzlmax, LSUB, 0, 0, Glu );
	usub  = (int *)    pxgstrf_expand( &nzumax, USUB, 0, 1, Glu );

	while ( !ucol || !lsub || !usub ) {
	    /*ABORT("Not enough core in LUMemInit()");*/
#if (PRNTlevel==1)
	    printf(".. pdgstrf_MemInit(): #retries %d\n", ++retries);
#endif
	    if ( whichspace == SYSTEM ) {
		SUPERLU_FREE(ucol);
		SUPERLU_FREE(lsub);
		SUPERLU_FREE(usub);
	    } else {
		duser_free(nzumax*dword+(nzlmax+nzumax)*iword, HEAD);
	    }
	    nzumax /= 2;    /* reduce request */
	    nzlmax /= 2;
	    if ( nzumax < annz/2 ) {
		printf("Not enough memory to perform factorization.\n");
		return (pdgstrf_memory_use(nzlmax, nzumax, nzlumax) + n);
	    }
	    ucol  = (double *) pxgstrf_expand( &nzumax, UCOL, 0, 0, Glu );
	    lsub  = (int *)    pxgstrf_expand( &nzlmax, LSUB, 0, 0, Glu );
	    usub  = (int *)    pxgstrf_expand( &nzumax, USUB, 0, 1, Glu );
	}
	
	if ( !lusup )  {
	    printf("Not enough memory to perform factorization.\n");
	    return (pdgstrf_memory_use(nzlmax, nzumax, nzlumax) + n);
	}
	
    } else { /* refact == YES */
	Lstore   = L->Store;
	Ustore   = U->Store;
	xsup     = Lstore->sup_to_colbeg;
	xsup_end = Lstore->sup_to_colend;
	supno    = Lstore->col_to_sup;
	xlsub    = Lstore->rowind_colbeg;
	xlsub_end= Lstore->rowind_colend;
	xlusup   = Lstore->nzval_colbeg;
	xlusup_end= Lstore->nzval_colend;
	xusub    = Ustore->colbeg;
	xusub_end= Ustore->colend;
	nzlmax   = Glu->nzlmax;    /* max from previous factorization */
	nzumax   = Glu->nzumax;
	nzlumax  = Glu->nzlumax;
	
	if ( lwork == -1 ) {
	    return (GluIntArray(n) * iword + superlu_TempSpace(n, panel_size, nprocs)
		    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword);
        } else if ( lwork == 0 ) {
	    whichspace = SYSTEM;
	} else {
	    whichspace = USER;
	    stack.size = lwork;
	    stack.top2 = lwork;
	}
	
	lsub  = expanders[LSUB].mem  = Lstore->rowind;
	lusup = expanders[LUSUP].mem = Lstore->nzval;
	usub  = expanders[USUB].mem  = Ustore->rowind;
	ucol  = expanders[UCOL].mem  = Ustore->nzval;;

	expanders[LSUB].size         = nzlmax;
	expanders[LUSUP].size        = nzlumax;
	expanders[USUB].size         = nzumax;
	expanders[UCOL].size         = nzumax;	
    }

    Glu->xsup       = xsup;
    Glu->xsup_end   = xsup_end;
    Glu->supno      = supno;
    Glu->lsub       = lsub;
    Glu->xlsub      = xlsub;
    Glu->xlsub_end  = xlsub_end;
    Glu->lusup      = lusup;
    Glu->xlusup     = xlusup;
    Glu->xlusup_end = xlusup_end;
    Glu->ucol       = ucol;
    Glu->usub       = usub;
    Glu->xusub      = xusub;
    Glu->xusub_end  = xusub_end;
    Glu->nzlmax     = nzlmax;
    Glu->nzumax     = nzumax;
    Glu->nzlumax    = nzlumax;
    ++no_expand;

#if ( PRNTlevel>=1 )
    printf(".. pdgstrf_MemInit() refact %d, space? %d, nzlumax %d, nzumax %d, nzlmax %d\n",
	   refact, whichspace, nzlumax, nzumax, nzlmax);
    fflush(stdout);
#endif

    return 0;
    
} /* pdgstrf_MemInit */

/* 
 * Allocate known working storage. Returns 0 if success, otherwise
 * returns the number of bytes allocated so far when failure occurred.
 */
int
pdgstrf_WorkInit(int n, int panel_size, int **iworkptr, double **dworkptr)
{
    int    isize, dsize, extra;
    double *old_ptr;
    int    maxsuper = sp_ienv(3),
           rowblk   = sp_ienv(4);

    isize = (2*panel_size + 5 + NO_MARKER) * n * sizeof(int);
    dsize = (n * panel_size +
	     NUM_TEMPV(n,panel_size,maxsuper,rowblk)) * sizeof(double);
    
    if ( whichspace == SYSTEM ) 
	*iworkptr = (int *) intCalloc(isize/sizeof(int));
    else
	*iworkptr = (int *) duser_malloc(isize, TAIL);
    if ( ! *iworkptr ) {
	fprintf(stderr, "pdgstrf_WorkInit: malloc fails for local iworkptr[]\n");
	return (isize + n);
    }

    if ( whichspace == SYSTEM )
	*dworkptr = (double *) SUPERLU_MALLOC(dsize);
    else {
	*dworkptr = (double *) duser_malloc(dsize, TAIL);
	if ( NotDoubleAlign(*dworkptr) ) {
	    old_ptr = *dworkptr;
	    *dworkptr = (double*) DoubleAlign(*dworkptr);
	    *dworkptr = (double*) ((double*)*dworkptr - 1);
	    extra = (char*)old_ptr - (char*)*dworkptr;
#ifdef CHK_EXPAND	    
	    printf("pdgstrf_WorkInit: not aligned, extra %d\n", extra);
#endif	    
	    stack.top2 -= extra;
	    stack.used += extra;
	}
    }
    if ( ! *dworkptr ) {
	fprintf(stderr, "malloc fails for local dworkptr[].");
	return (isize + dsize + n);
    }
	
    return 0;
}


/*
 * Set up pointers for real working arrays.
 */
void
pdgstrf_SetRWork(int n, int panel_size, double *dworkptr,
		 double **dense, double **tempv)
{
    int maxsuper = sp_ienv(3);
    int rowblk   = sp_ienv(4);
    *dense = dworkptr;
    *tempv = *dense + panel_size*n;
    dfill (*dense, n * panel_size, 0.0);
    dfill (*tempv, NUM_TEMPV(n,panel_size,maxsuper,rowblk), 0.0);     
}
	
/*
 * Free the working storage used by factor routines.
 */
void pdgstrf_WorkFree(int *iwork, double *dwork, GlobalLU_t *Glu)
{
    if ( whichspace == SYSTEM ) {
	SUPERLU_FREE (iwork);
	SUPERLU_FREE (dwork);
    } else {
	stack.used -= (stack.size - stack.top2);
	stack.top2 = stack.size;
/*	pdgstrf_StackCompress(Glu);  */
    }
}

/* 
 * Expand the data structures for L and U during the factorization.
 * Return value:   0 - successful return
 *               > 0 - number of bytes allocated when run out of space
 */
int
pdgstrf_MemXpand(
		 int jcol,
		 int next, /* number of elements currently in the factors */
		 MemType mem_type,/* which type of memory to expand  */
		 int *maxlen, /* modified - max. length of a data structure */
		 GlobalLU_t *Glu /* modified - global LU data structures */
		 )
{
    void   *new_mem;
    
#ifdef CHK_EXPAND    
    printf("pdgstrf_MemXpand(): jcol %d, next %d, maxlen %d, MemType %d\n",
	   jcol, next, *maxlen, mem_type);
#endif    

    if (mem_type == USUB) 
    	new_mem = pxgstrf_expand(maxlen, mem_type, next, 1, Glu);
    else
	new_mem = pxgstrf_expand(maxlen, mem_type, next, 0, Glu);
    
    if ( !new_mem ) {
	int    nzlmax  = Glu->nzlmax;
	int    nzumax  = Glu->nzumax;
	int    nzlumax = Glu->nzlumax;
    	fprintf(stderr, "Can't expand MemType %d: jcol %d\n", mem_type, jcol);
    	return (pdgstrf_memory_use(nzlmax, nzumax, nzlumax) + ndim);
    }

    switch ( mem_type ) {
      case LUSUP:
	Glu->lusup   = (double *) new_mem;
	Glu->nzlumax = *maxlen;
	break;
      case UCOL:
	Glu->ucol   = (double *) new_mem;
	Glu->nzumax = *maxlen;
	break;
      case LSUB:
	Glu->lsub   = (int *) new_mem;
	Glu->nzlmax = *maxlen;
	break;
      case USUB:
	Glu->usub   = (int *) new_mem;
	Glu->nzumax = *maxlen;
	break;
    }
    
    return 0;
    
}


void
copy_mem_double(int howmany, void *old, void *new)
{
    register int i;
    double *dold = old;
    double *dnew = new;
    for (i = 0; i < howmany; i++) dnew[i] = dold[i];
}

/*
 * Expand the existing storage to accommodate more fill-ins.
 */
void
*pxgstrf_expand(
		int *prev_len,   /* length used from previous call */
		MemType type,    /* which part of the memory to expand */
		int len_to_copy, /* size of memory to be copied to new store */
		int keep_prev,   /* = 1: use prev_len;
				    = 0: compute new_len to expand */
		GlobalLU_t *Glu  /* modified - global LU data structures */
		)
{
    double   alpha = EXPAND;
    void     *new_mem, *old_mem;
    int      new_len, tries, lword, extra, bytes_to_copy;

    if ( no_expand == 0 || keep_prev ) /* First time allocate requested */
        new_len = *prev_len;
    else {
	new_len = alpha * *prev_len;
    }
    
    if ( type == LSUB || type == USUB ) lword = sizeof(int);
    else lword = sizeof(double);

    if ( whichspace == SYSTEM ) {
	new_mem = (void *) SUPERLU_MALLOC(new_len * lword);
    
	if ( no_expand != 0 ) {
	    tries = 0;
	    if ( keep_prev ) {
		if ( !new_mem ) return (NULL);
	    } else {
		while ( !new_mem ) {
		    if ( ++tries > 10 ) return (NULL);
		    alpha = Reduce(alpha);
		    new_len = alpha * *prev_len;
		    new_mem = (void *) SUPERLU_MALLOC(new_len * lword); 
		}
	    }
	    if ( type == LSUB || type == USUB ) {
		copy_mem_int(len_to_copy, expanders[type].mem, new_mem);
	    } else {
		copy_mem_double(len_to_copy, expanders[type].mem, new_mem);
	    }
	    SUPERLU_FREE (expanders[type].mem);
	}
	expanders[type].mem = (void *) new_mem;
#if ( MACH==SGI || MACH==ORIGIN )
/*	bzero(new_mem, new_len*lword);*/
#endif	
	
    } else { /* whichspace == USER */
	if ( no_expand == 0 ) {
	    new_mem = duser_malloc(new_len * lword, HEAD);
	    if ( NotDoubleAlign(new_mem) &&
		(type == LUSUP || type == UCOL) ) {
		old_mem = new_mem;
		new_mem = (void *)DoubleAlign(new_mem);
		extra = (char*)new_mem - (char*)old_mem;
#ifdef CHK_EXPAND		
		printf("expand(): not aligned, extra %d\n", extra);
#endif		
		stack.top1 += extra;
		stack.used += extra;
	    }
	    expanders[type].mem = (void *) new_mem;
	}
	else {
	    tries = 0;
	    extra = (new_len - *prev_len) * lword;
	    if ( keep_prev ) {
		if ( StackFull(extra) ) return (NULL);
	    } else {
		while ( StackFull(extra) ) {
		    if ( ++tries > 10 ) return (NULL);
		    alpha = Reduce(alpha);
		    new_len = alpha * *prev_len;
		    extra = (new_len - *prev_len) * lword;	    
		}
	    }

	    if ( type != USUB ) {
		new_mem = (void*)((char*)expanders[type + 1].mem + extra);
		bytes_to_copy = (char*)stack.array + stack.top1
		    - (char*)expanders[type + 1].mem;
		user_bcopy(expanders[type+1].mem, new_mem, bytes_to_copy);

		if ( type < USUB ) {
		    Glu->usub = expanders[USUB].mem =
			(void*)((char*)expanders[USUB].mem + extra);
		}
		if ( type < LSUB ) {
		    Glu->lsub = expanders[LSUB].mem =
			(void*)((char*)expanders[LSUB].mem + extra);
		}
		if ( type < UCOL ) {
		    Glu->ucol = expanders[UCOL].mem =
			(void*)((char*)expanders[UCOL].mem + extra);
		}
		stack.top1 += extra;
		stack.used += extra;
		if ( type == UCOL ) {
		    stack.top1 += extra;   /* Add same amount for USUB */
		    stack.used += extra;
		}
		
	    } /* if ... */

	} /* else ... */
    }
#ifdef DEBUG
    printf("pxgstrf_expand[type %d]\n", type);
#endif
    expanders[type].size = new_len;
    *prev_len = new_len;
    if ( no_expand ) ++no_expand;
    
    return (void *) expanders[type].mem;
    
} /* expand */


/*
 * Compress the work[] array to remove fragmentation.
 */
void
pdgstrf_StackCompress(GlobalLU_t *Glu)
{
    register int iword, dword;
    char     *last, *fragment;
    int      *ifrom, *ito;
    double   *dfrom, *dto;
    int      *xlsub, *lsub, *xusub_end, *usub, *xlusup;
    double   *ucol, *lusup;
    
    iword = sizeof(int);
    dword = sizeof(double);

    xlsub  = Glu->xlsub;
    lsub   = Glu->lsub;
    xusub_end  = Glu->xusub_end;
    usub   = Glu->usub;
    xlusup = Glu->xlusup;
    ucol   = Glu->ucol;
    lusup  = Glu->lusup;
    
    dfrom = ucol;
    dto = (double *)((char*)lusup + xlusup[ndim] * dword);
    copy_mem_double(xusub_end[ndim-1], dfrom, dto);
    ucol = dto;

    ifrom = lsub;
    ito = (int *) ((char*)ucol + xusub_end[ndim-1] * iword);
    copy_mem_int(xlsub[ndim], ifrom, ito);
    lsub = ito;
    
    ifrom = usub;
    ito = (int *) ((char*)lsub + xlsub[ndim] * iword);
    copy_mem_int(xusub_end[ndim-1], ifrom, ito);
    usub = ito;
    
    last = (char*)usub + xusub_end[ndim-1] * iword;
    fragment = (char*) ((char*)stack.array + stack.top1 - last);
    stack.used -= (long int) fragment;
    stack.top1 -= (long int) fragment;

    Glu->ucol = ucol;
    Glu->lsub = lsub;
    Glu->usub = usub;
    
#ifdef CHK_EXPAND
    printf("pdgstrf_StackCompress: fragment %d\n", fragment);
    /* PrintStack("After compress", Glu);
    for (last = 0; last < ndim; ++last)
	print_lu_col("After compress:", last, 0);*/
#endif    
    
}

/*
 * Allocate storage for original matrix A
 */
void
dallocateA(int n, int nnz, double **a, int **asub, int **xa)
{
    *a    = (double *) doubleMalloc(nnz);
    *asub = (int *) intMalloc(nnz);
    *xa   = (int *) intMalloc(n+1);
}

double *doubleMalloc(int n)
{
    double *buf;
    buf = (double *) SUPERLU_MALLOC(n * sizeof(double)); 
    if ( !buf ) {
	fprintf(stderr, "SUPERLU_MALLOC failed for buf in doubleMalloc()");
	exit (1);
    }
    return (buf);
}

double *doubleCalloc(int n)
{
    double *buf;
    register int i;
    buf = (double *) SUPERLU_MALLOC(n * sizeof(double));
    if ( !buf ) {
	fprintf(stderr, "SUPERLU_MALLOC failed for buf in doubleCalloc()");
	exit (1);
    }
    for (i = 0; i < n; ++i) buf[i] = 0.;
    return (buf);
}



