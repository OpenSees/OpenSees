#include "pdsp_defs.h"

pdgstrf_threadarg_t *
pdgstrf_thread_init(SuperMatrix *A, SuperMatrix *L, SuperMatrix *U,
		    pdgstrf_options_t *pdgstrf_options, 
		    pxgstrf_shared_t *pxgstrf_shared,
		    Gstat_t *Gstat, int *info)
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
 * pdgstrf_thread_init() initializes the parallel data structures
 * for the multithreaded routine pdgstrf_thread().
 *
 * Arguments
 * =========
 *
 * A        (input) SuperMatrix*
 *	    Original matrix A, permutated by columns, of dimension
 *          (A->nrow, A->ncol). The type of A can be:
 *          Stype = NCP; Dtype = _D; Mtype = GE.
 *
 * L        (input) SuperMatrix*
 *          If pdgstrf_options->refact = YES, then use the existing
 *          storage in L to perform LU factorization;
 *          Otherwise, L is not accessed. L has types: 
 *          Stype = SCP, Dtype = _D, Mtype = TRLU.
 *
 * U        (input) SuperMatrix*
 *          If pdgstrf_options->refact = YES, then use the existing
 *          storage in U to perform LU factorization;
 *          Otherwise, U is not accessed. U has types:
 *          Stype = NCP, Dtype = _D, Mtype = TRU.
 *
 * pdgstrf_options (input) pdgstrf_options_t*
 *          The structure contains the parameters to control how the
 *          factorization is performed;
 *          See pdgstrf_options_t structure defined in pdsp_defs.h.
 *
 * pxgstrf_shared (output) pxgstrf_shared_t*
 *          The structure contains the shared task queue and the 
 *          synchronization variables for parallel factorization.
 *          See pxgstrf_shared_t structure defined in pdsp_defs.h.
 *
 * Gstat    (output) Gstat_t*
 *          Record all the statistics about the factorization; 
 *          See Gstat_t structure defined in util.h.
 *
 * info     (output) int*
 *          = 0: successful exit
 *          > 0: if pdgstrf_options->lwork = -1, info returns the estimated
 *               amount of memory (in bytes) required;
 *               Otherwise, it returns the number of bytes allocated when
 *               memory allocation failure occurred, plus A->ncol.
 *
 */
    static GlobalLU_t Glu; /* persistent to support repeated factors. */
    pdgstrf_threadarg_t *pdgstrf_threadarg;
    register int n, i, nprocs;
    NCPformat *Astore;
    int  *perm_c;
    int  *perm_r;
    int  *inv_perm_c; /* inverse of perm_c */
    int  *inv_perm_r; /* inverse of perm_r */
    int	 *xprune;  /* points to locations in subscript vector lsub[*].
			For column i, xprune[i] denotes the point where 
			structural pruning begins.
			I.e. only xlsub[i],..,xprune[i]-1 need to be
			traversed for symbolic factorization.     */
    int  *ispruned;/* flag to indicate whether column j is pruned */
    int   nzlumax;
    pxgstrf_relax_t *pxgstrf_relax;
    
    nprocs     = pdgstrf_options->nprocs;
    perm_c     = pdgstrf_options->perm_c;
    perm_r     = pdgstrf_options->perm_r;
    n          = A->ncol;
    Astore     = A->Store;
    inv_perm_r = (int *) intMalloc(n);
    inv_perm_c = (int *) intMalloc(n);
    xprune     = (int *) intMalloc(n);
    ispruned   = (int *) intCalloc(n);
    
    /* Pack shared data objects to each process. */
    pxgstrf_shared->inv_perm_r   = inv_perm_r;
    pxgstrf_shared->inv_perm_c   = inv_perm_c;
    pxgstrf_shared->xprune       = xprune;
    pxgstrf_shared->ispruned     = ispruned;
    pxgstrf_shared->A            = A;
    pxgstrf_shared->Glu          = &Glu;
    pxgstrf_shared->Gstat        = Gstat;
    pxgstrf_shared->info         = info;

    if ( pdgstrf_options->usepr ) {
	/* Compute the inverse of perm_r */
	for (i = 0; i < n; ++i) inv_perm_r[perm_r[i]] = i;
    }
    for (i = 0; i < n; ++i) inv_perm_c[perm_c[i]] = i;

    /* Initialization. */
    Glu.nsuper = -1;
    Glu.nextl  = 0;
    Glu.nextu  = 0;
    Glu.nextlu = 0;
    ifill(perm_r, n, EMPTY);

    /* Identify relaxed supernodes at the bottom of the etree. */
    pxgstrf_relax = (pxgstrf_relax_t *)
        SUPERLU_MALLOC((n+2) * sizeof(pxgstrf_relax_t));
    pxgstrf_relax_snode(n, pdgstrf_options, pxgstrf_relax);
    
    /* Initialize mutex variables, task queue, determine panels. */
    ParallelInit(n, pxgstrf_relax, pdgstrf_options, pxgstrf_shared);
    
    /* Set up memory image in lusup[*]. */
    nzlumax = PresetMap(n, A, pxgstrf_relax, pdgstrf_options, &Glu);
    if ( pdgstrf_options->refact == NO ) Glu.nzlumax = nzlumax;
    
    SUPERLU_FREE (pxgstrf_relax);

    /* Allocate global storage common to all the factor routines */
    *info = pdgstrf_MemInit(n, Astore->nnz, pdgstrf_options, L, U, &Glu);
    if ( *info ) return NULL;

    /* Prepare arguments to all threads. */
    pdgstrf_threadarg = (pdgstrf_threadarg_t *) 
        SUPERLU_MALLOC(nprocs * sizeof(pdgstrf_threadarg_t));
    for (i = 0; i < nprocs; ++i) {
        pdgstrf_threadarg[i].pnum = i;
        pdgstrf_threadarg[i].info = 0;
	pdgstrf_threadarg[i].pdgstrf_options = pdgstrf_options;
	pdgstrf_threadarg[i].pxgstrf_shared = pxgstrf_shared;
    }

#if ( DEBUGlevel==1 )
    printf("** pdgstrf_thread_init() called\n");
#endif

    return (pdgstrf_threadarg);
}
