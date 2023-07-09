/*! @file
 * \brief Implements parallel symbolic factorization
 *
 * <pre>
 * -- Parallel symbolic factorization routine  (version 2.3) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley - July 2003
 * INRIA France - January 2004
 * Laura Grigori
 *
 * November 1, 2007
 * Feburary 20, 2008
 * October 15, 2008
 *
 * The function symbfact_dist implements the parallel symbolic factorization
 * algorithm described in the paper:
 *
 * Parallel Symbolic Factorization for Sparse LU with Static Pivoting,
 * Laura Grigori, James W. Demmel and Xiaoye S. Li,
 * Pages 1289-1314, SIAM Journal on Scientific Computing, Volume 29, Issue 3.
 * </pre>
 */

/* limits.h:  the largest positive integer (INT_MAX) */
#include <limits.h>
#include <math.h>
#include "superlu_ddefs.h"
#include "psymbfact.h"

/*
 * Internal protypes
 */

static int_t *
intMalloc_symbfact(int_t );

static int_t *
intCalloc_symbfact(int_t );

static int_t
initParmsAndStats
(psymbfact_stat_t *PS);

static void
estimate_memUsage
(int_t, int, mem_usage_t *, float *, float *,
 Pslu_freeable_t *, Llu_symbfact_t *,
 vtcsInfo_symbfact_t *, comm_symbfact_t *, psymbfact_stat_t *);

static void
symbfact_free 
(int, int, Llu_symbfact_t *, vtcsInfo_symbfact_t *, comm_symbfact_t *);

static int_t
denseSep_symbfact 
(int , int_t, int, int, int, int_t *, int_t *, int, 
 int,  int, int_t, int_t, int_t *, int_t *, int_t *,
 int_t *, int_t *, MPI_Comm, MPI_Comm *, Llu_symbfact_t *,
 Pslu_freeable_t *_freeable, vtcsInfo_symbfact_t *, 
 comm_symbfact_t *, psymbfact_stat_t * );

static int_t
dnsUpSeps_symbfact
(int_t, int, int, int, int, int_t *, int_t *, int_t,
 Llu_symbfact_t *, Pslu_freeable_t *, vtcsInfo_symbfact_t *, 
 comm_symbfact_t *, psymbfact_stat_t *, int_t *, int_t *, int_t *);

static void
intraLvl_symbfact 
(SuperMatrix *, int, int, int, int, int, int_t *, int_t *, int, 
 int, int_t, int_t,  Pslu_freeable_t *, Llu_symbfact_t *, vtcsInfo_symbfact_t *, 
 comm_symbfact_t *, psymbfact_stat_t *, int_t *, int_t *, int_t *, int_t *, 
 int_t *, int_t *, int_t *, MPI_Comm, MPI_Comm *);

static void
initLvl_symbfact
(int_t, int, int_t, int_t, Pslu_freeable_t *, 
 Llu_symbfact_t *, vtcsInfo_symbfact_t *, psymbfact_stat_t *, MPI_Comm, 
 int_t *, int_t, int_t);

static void
createComm (int, int, MPI_Comm *, MPI_Comm *);

static void
freeComm (int, int, MPI_Comm *, MPI_Comm *);

static void
domain_symbfact
(SuperMatrix *, int, int, int,  int, int, int_t *, int_t *,
 int_t, int_t, Pslu_freeable_t *, Llu_symbfact_t *, vtcsInfo_symbfact_t *,
 comm_symbfact_t *, psymbfact_stat_t *, int_t *, int_t *, int_t *, int_t *, 
 int_t *, int_t *, int_t *);

static float
allocPrune_domain
(int_t, int_t, Llu_symbfact_t *, vtcsInfo_symbfact_t *, psymbfact_stat_t *);

static float
allocPrune_lvl
(Llu_symbfact_t *, vtcsInfo_symbfact_t *, psymbfact_stat_t *);

static int
symbfact_alloc
(int_t, int, Pslu_freeable_t *, Llu_symbfact_t *, 
 vtcsInfo_symbfact_t *, comm_symbfact_t *, psymbfact_stat_t *);

static float 
symbfact_mapVtcs
(int, int, int, SuperMatrix *, int_t *, int_t *, 
 Pslu_freeable_t *, vtcsInfo_symbfact_t *, int_t *, int_t, psymbfact_stat_t *);

static void 
symbfact_distributeMatrix 
(int, int, int, SuperMatrix *, int_t *, int_t *, matrix_symbfact_t *, 
 Pslu_freeable_t *, vtcsInfo_symbfact_t *, int_t *, MPI_Comm *);

static int_t
interLvl_symbfact
(SuperMatrix *, int, int, int, int, int, int, int, 
 int_t *, int_t *, int_t *, int_t *, int_t *, int_t *, int_t *,
 Llu_symbfact_t *, Pslu_freeable_t*, comm_symbfact_t *, vtcsInfo_symbfact_t *,
 psymbfact_stat_t *, MPI_Comm, MPI_Comm *);

static float
cntsVtcs 
(int_t, int, int, Pslu_freeable_t *, Llu_symbfact_t *, vtcsInfo_symbfact_t *, 
 int_t *, int_t *, int_t *, psymbfact_stat_t *, MPI_Comm *);

/************************************************************************/
float symbfact_dist
/************************************************************************/
(
 int         nprocs_num,  /* Input - no of processors */
 int         nprocs_symb, /* Input - no of processors for the symbolic
			     factorization */
 SuperMatrix *A,          /* Input - distributed input matrix */
 int_t       *perm_c,     /* Input - column permutation */
 int_t       *perm_r,     /* Input - row permutation */
 int_t       *sizes,      /* Input - sizes of each node in the separator tree */
 int_t       *fstVtxSep,  /* Input - first vertex of each node in the tree */
 Pslu_freeable_t *Pslu_freeable, /* Output - local L and U structure, 
				    global to local indexing information */
 MPI_Comm    *num_comm,   /* Input - communicator for numerical factorization */
 MPI_Comm    *symb_comm,  /* Input - communicator for symbolic factorization */
 mem_usage_t *symb_mem_usage
 )
{
/*! \brief
 *
 * <pre> 
 * Purpose
 * =======
 *   symbfact_dist() performs symbolic factorization of matrix A suitable
 *   for performing the supernodal Gaussian elimination with no pivoting (GEPP). 
 *   This routine computes the structure of one column of L and one row of U 
 *   at a time.  It uses:
 *        o distributed input matrix
 *        o supernodes
 *        o symmetric structure pruning
 *
 *
 * Arguments
 * =========
 *
 * nprocs_num (input) int
 *         Number of processors SuperLU_DIST is executed on, and the input 
 *         matrix is distributed on.
 *
 * nprocs_symb (input) int
 *         Number of processors on which the symbolic factorization is
 *         performed.  It is equal to the number of independent domains
 *         idenfied in the graph partitioning algorithm executed
 *         previously and has to be a power of 2.  It corresponds to
 *         number of leaves in the separator tree.
 *
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The
 *         number of the linear equations is A->nrow.  Matrix A is
 *         distributed in NRformat_loc format.
 *         Matrix A is not yet permuted by perm_c.
 *
 * perm_c  (input) int_t*
 *	   Column permutation vector of size A->ncol, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 *
 * perm_r  (input) int_t*
 *	   Row permutation vector of size A->nrow, which defines the 
 *         permutation matrix Pr; perm_r[i] = j means column i of A is 
 *         in position j in Pr*A.
 *
 * sizes   (input) int_t*
 *         Contains the number of vertices in each separator.
 *
 * fstVtxSep (input) int_t*
 *         Contains first vertex for each separator.
 *
 * Pslu_freeable (output) Pslu_freeable_t*
 *         Returns the local L and U structure, and global to local
 *         information on the indexing of the vertices.  Contains all
 *         the information necessary for performing the data
 *         distribution towards the numeric factorization.
 *				    
 * num_comm (input) MPI_Comm*
 *         Communicator for numerical factorization 
 *
 * symb_comm (input) MPI_Comm*
 *         Communicator for symbolic factorization 
 *
 * symb_mem_usage (input) mem_usage_t *
 *         Statistics on memory usage.
 *
 * Return value
 * ============
 *   < 0, number of bytes allocated on return from the symbolic factorization.
 *   > 0, number of bytes allocated when out of memory.
 *
 * Sketch of the algorithm
 * =======================
 *
 *  Distrbute the vertices on the processors using a subtree to
 *  subcube algorithm.
 *
 *  Redistribute the structure of the input matrix A according to the
 *  subtree to subcube computed previously for the symbolic
 *  factorization routine.  This implies in particular a distribution
 *  from nprocs_num processors to nprocs_symb processors.
 *
 *  Perform symbolic factorization guided by the separator tree provided by
 *  a graph partitioning algorithm.  The symbolic factorization uses a 
 *  combined left-looking, right-looking approach. 
 * </pre>
 */
  NRformat_loc *Astore;
  int iam, szSep, fstP, lstP, npNode, nlvls, lvl, p, iSep, jSep;
  int iinfo; /* return code */
  int_t m, n;
  int_t nextl, nextu, neltsZr, neltsTotal, nsuper_loc, szLGr, szUGr;
  int_t ind_blk, nsuper, vtx, min_mn, szsn;
  long long int nnzL, nnzU, nnzLU;
  float stat_loc[23], stat_glob[23], mem_glob[15];
  
  Llu_symbfact_t Llu_symbfact; /* local L and U and pruned L and U data structures */
  vtcsInfo_symbfact_t VInfo; /* local information on number of blocks,
				number of vertices in a block etc */
  matrix_symbfact_t   AS; /* temporary storage for the input matrix after redistribution */
  comm_symbfact_t CS;  /* information on communication */
  /* relaxation parameters (for future release) and 
     statistics collected during the symbolic factorization */
  psymbfact_stat_t PS; 
  /* temp array of size n, used as a marker by the subroutines */
  int_t *tempArray; 
  int_t i, j, k;
  int_t fstVtx, lstVtx, mark, fstVtx_lid, vtx_lid, maxNvtcsPProc;
  int_t nnz_asup_loc, nnz_ainf_loc, fill_rcmd;
  float totalMemLU, overestimMem;
  MPI_Comm *commLvls;  

  /* maximum block size */
  int_t  maxSzBlk;
  float flinfo;
#if ( PRNTlevel >= 1)
  float stat_msgs_l[10], stat_msgs_g[10]; 
#endif  
#if ( PROFlevel>=1 )
  double t, t_symbFact[3], t_symbFact_loc[3];
  double *time_lvlsT, *time_lvls, t1, t2, time_lvlsg[9];
#endif
  
  /* Initialization */
  MPI_Comm_rank ((*num_comm), &iam);
  commLvls = NULL;
#if ( DEBUGlevel>=1 )
  CHECK_MALLOC(iam, "Enter psymbfact()");
#endif
  initParmsAndStats (&PS);
  if (nprocs_symb != 1) {
    if (!(commLvls = (MPI_Comm *) SUPERLU_MALLOC(2*nprocs_symb*sizeof(MPI_Comm)))) {
      fprintf (stderr, "Malloc fails for commLvls[].");  
      return (PS.allocMem);
    }
    PS.allocMem += 2 * nprocs_symb * sizeof(MPI_Comm);
  }
  
  nlvls = (int) LOG2( nprocs_num ) + 1;
#if ( PROFlevel>=1 )
  time_lvlsT = (double *) SUPERLU_MALLOC(3*nprocs_symb*(nlvls+1) 
					 * sizeof(double));
  time_lvls  = (double *) SUPERLU_MALLOC(3*(nlvls+1) * sizeof(double));
  if (!time_lvls || !time_lvlsT) {
    fprintf (stderr, "Malloc fails for time_lvls[].");  
    return (PS.allocMem);
  }
  PS.allocMem += (3*nprocs_symb*(nlvls+1) + 3*(nlvls+1)) * sizeof(double);
#endif
  
  VInfo.xlsub_nextLvl  = 0;
  VInfo.xusub_nextLvl  = 0;
  VInfo.maxSzBlk = sp_ienv_dist(3);
  maxSzBlk = VInfo.maxSzBlk;
  
  mark = EMPTY;
  nsuper_loc = 0;
  nextl   = 0; nextu      = 0;
  neltsZr = 0; neltsTotal = 0;
  
  m = A->nrow;
  n = A->ncol;
  min_mn = SUPERLU_MIN( m, n );
  
  if (!(tempArray = intMalloc_symbfact(n))) {
    fprintf (stderr, "Malloc fails for tempArray[].\n");  
    return (PS.allocMem);
  }
  PS.allocMem += n * sizeof(int_t);
  
#if ( PROFlevel>=1 )  
  t = SuperLU_timer_();
#endif
  
  /* Distribute vertices on processors */
  if ((flinfo = 
       symbfact_mapVtcs (iam, nprocs_num, nprocs_symb, A, fstVtxSep, sizes, 
			 Pslu_freeable, &VInfo, tempArray, maxSzBlk, &PS)) > 0) 
    return (flinfo);

  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  
  /* Redistribute matrix A on processors following the distribution found
     in symbfact_mapVtcs.  Store the redistributed A temporarily into AS */
  symbfact_distributeMatrix (iam, nprocs_num, nprocs_symb,  A, 
			     perm_c, perm_r, &AS, 
			     Pslu_freeable, &VInfo, tempArray, num_comm);
  
  /* THE REST OF THE SYMBOLIC FACTORIZATION IS EXECUTED ONLY BY NPROCS_SYMB
     PROCESSORS */
  if ( iam < nprocs_symb ) {
    
#if ( PROFlevel>=1 )
    t_symbFact_loc[0] = SuperLU_timer_() - t;
    t = SuperLU_timer_();
    t_symbFact_loc[1] = t;
#endif

    /* Allocate storage common to the symbolic factor routines */
    if (iinfo = symbfact_alloc (n, nprocs_symb, Pslu_freeable, 
				&Llu_symbfact, &VInfo, &CS, &PS)) 
      return (PS.allocMem);
    /* Copy the redistributed input matrix AS at the end of the memory buffer
       allocated to store L and U.  That is, copy (AS.x_ainf, AS.ind_ainf) in
       (xlsub, lsub), (AS.x_asup, AS.ind_asup) in (xusub, usub).  Free the
       memory used to store the input matrix */
    nnz_ainf_loc = VInfo.nnz_ainf_loc;
    nnz_asup_loc = VInfo.nnz_asup_loc;
    j = Llu_symbfact.szUsub - VInfo.nnz_asup_loc;
    k = Llu_symbfact.szLsub - VInfo.nnz_ainf_loc;
    for (i = 0; i <= VInfo.nvtcs_loc; i++) {
      Llu_symbfact.xusub[i] = AS.x_asup[i] + j;
      Llu_symbfact.xlsub[i] = AS.x_ainf[i] + k;
    }
    
    for (i = 0; i < VInfo.nnz_asup_loc; i++, j++) 
      Llu_symbfact.usub[j] = AS.ind_asup[i];
    for (i = 0; i < VInfo.nnz_ainf_loc; i++, k++)   
      Llu_symbfact.lsub[k] = AS.ind_ainf[i];
    SUPERLU_FREE( AS.x_ainf );
    SUPERLU_FREE( AS.x_asup );  
    SUPERLU_FREE( AS.ind_ainf );  
    SUPERLU_FREE( AS.ind_asup );  

    if (nprocs_symb != 1) {
      createComm (iam, nprocs_symb, commLvls, symb_comm);    

#if ( PROFlevel>=1 )
      t_symbFact_loc[2] = SuperLU_timer_();
#endif
      if ((flinfo = cntsVtcs (n, iam, nprocs_symb, Pslu_freeable, &Llu_symbfact, 
			      &VInfo, tempArray, fstVtxSep, sizes, &PS, commLvls)) > 0) 
	return (flinfo);
			     
#if ( PROFlevel>=1 )
      t_symbFact_loc[2] = SuperLU_timer_() - t_symbFact_loc[2];
#endif
    }
    
    /* set to EMPTY marker[] array */
    for (i = 0; i < n; i++)
      tempArray[i] = EMPTY;
    
    szSep = nprocs_symb;
    iSep = 0;
    lvl = 0;
    while (szSep >= 1) {
      /* for each level in the separator tree */
      npNode = nprocs_symb / szSep; 
      fstP = 0; 
      /* for each node in the level */
      for (jSep = iSep; jSep < iSep + szSep; jSep++) {
	fstVtx = fstVtxSep[jSep];
	lstVtx  = fstVtx + sizes[jSep];
	/* if this is the first level */
	if (szSep == nprocs_symb) {
	  /* compute symbolic factorization for my domain */
	  if (fstP == iam) {
	    /* allocate storage for the pruned structures */
#if ( PROFlevel>=1 )
	    t1 = SuperLU_timer_();
#endif
	    if ((flinfo = allocPrune_domain (fstVtx, lstVtx, 
					     &Llu_symbfact, &VInfo, &PS)) > 0)
	      return (flinfo);
	    if (fstVtx < lstVtx)
	      VInfo.fstVtx_nextLvl = VInfo.begEndBlks_loc[2];
	    
	    domain_symbfact 
	      (A, iam, lvl, szSep, iSep, jSep, sizes, fstVtxSep, fstVtx, lstVtx, 
	       Pslu_freeable, &Llu_symbfact, &VInfo, &CS, &PS, tempArray, 
	       &mark, &nextl, &nextu, &neltsZr, &neltsTotal, &nsuper_loc);

	    PS.estimLSz = nextl;
	    PS.estimUSz = nextu;
	    if (nprocs_symb != 1) 
	      if((flinfo = allocPrune_lvl (&Llu_symbfact, &VInfo, &PS)) > 0)
		return (flinfo);
#if ( PROFlevel>=1 )
	    t2 = SuperLU_timer_();
	    time_lvls[lvl] = 0.; time_lvls[lvl+1] = 0.;
	    time_lvls[lvl + 2] = t2 - t1;
#endif
	  }
	}
	else {
	  lstP = fstP + npNode;
	  if (fstP <= iam && iam < lstP) {
#if ( PROFlevel>=1 )
	    t1 = SuperLU_timer_();	  
#endif
	    if (VInfo.filledSep != FILLED_SEPS)
	      initLvl_symbfact(n, iam, fstVtx, lstVtx,
			       Pslu_freeable, &Llu_symbfact, &VInfo, &PS, commLvls[jSep], 
			       tempArray, nextl, nextu);
#if ( PROFlevel>=1 )
	    t2 = SuperLU_timer_();
	    time_lvls[3*lvl] = t2 - t1;
#endif
	    interLvl_symbfact (A, iam, lvl, szSep, fstP, lstP,
			       iSep, jSep, sizes, fstVtxSep, 
			       &nextl, &nextu, &nsuper_loc, &mark, tempArray,
			       &Llu_symbfact, Pslu_freeable, &CS, &VInfo, &PS,
			       commLvls[jSep], symb_comm);
#if ( PROFlevel>=1 )
	    t1 = SuperLU_timer_();
	    time_lvls[3*lvl+1] = t1 - t2;
#endif
	    if (VInfo.filledSep != FILLED_SEPS)
	      intraLvl_symbfact 
		(A, iam, lvl, szSep, iSep, jSep, sizes, fstVtxSep, fstP, lstP, 
		 fstVtx, lstVtx, Pslu_freeable, &Llu_symbfact, &VInfo, &CS, &PS,
		 tempArray, &mark, &nextl, &nextu, &neltsZr, &neltsTotal, 
		 &nsuper_loc, commLvls[jSep], symb_comm);
#if ( PROFlevel>=1 )
	    t2 = SuperLU_timer_();
	    time_lvls[3*lvl+2] = t2 - t1;		 
#endif
	  }
	}
	fstP += npNode;
      }
      iSep += szSep;
      szSep = szSep / 2;
      lvl ++;
    }
  
    SUPERLU_FREE( tempArray );
    
    /* Set up global information and collect statistics */
    if (PS.maxSzLPr < Llu_symbfact.indLsubPr)
      PS.maxSzLPr = Llu_symbfact.indLsubPr;
    if (PS.maxSzUPr < Llu_symbfact.indUsubPr)
      PS.maxSzUPr = Llu_symbfact.indUsubPr;
    
    Llu_symbfact.xlsub[VInfo.nvtcs_loc] = nextl;
    Llu_symbfact.xusub[VInfo.nvtcs_loc] = nextu;
    fill_rcmd = SUPERLU_MAX( nextl / (nnz_ainf_loc+1), nextu / (nnz_asup_loc+1)) + 1;
    Pslu_freeable->xsup_beg_loc = intMalloc_dist (nsuper_loc+1);
    Pslu_freeable->xsup_end_loc = intMalloc_dist (nsuper_loc+1);
    if (!Pslu_freeable->xsup_beg_loc || !Pslu_freeable->xsup_end_loc) {
      fprintf (stderr, "Malloc fails for xsup_beg_loc, xsup_end_loc.");
      return (PS.allocMem);
    }
    PS.allocMem += 2 * (nsuper_loc+1) * sizeof(int_t);
    maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
    nnzL = 0; nnzU = 0;
    
    i = 0;
    nsuper = 0;
    ind_blk = 0;
    for (ind_blk = 0; ind_blk < VInfo.nblks_loc; ind_blk ++) {
      fstVtx = VInfo.begEndBlks_loc[2 * ind_blk];
      lstVtx = VInfo.begEndBlks_loc[2 * ind_blk + 1];
      fstVtx_lid = LOCAL_IND( Pslu_freeable->globToLoc[fstVtx] );
      nsuper = Pslu_freeable->supno_loc[fstVtx_lid];
      Pslu_freeable->xsup_beg_loc[nsuper] = fstVtx;
      szsn = 1;
      if (INT_MAX - nnzL <= Llu_symbfact.xlsub[fstVtx_lid + 1] - 
	  Llu_symbfact.xlsub[fstVtx_lid])
	printf ("PE[%d] ERR nnzL %ld\n", iam, nnzL); 
      if (INT_MAX - nnzU <= Llu_symbfact.xusub[fstVtx_lid + 1] - 
	  Llu_symbfact.xusub[fstVtx_lid])
	printf ("PE[%d] ERR nnzU %ld\n", iam, nnzU);
      
      j = Llu_symbfact.xlsub[fstVtx_lid + 1] - Llu_symbfact.xlsub[fstVtx_lid];
      k = Llu_symbfact.xusub[fstVtx_lid + 1] - Llu_symbfact.xusub[fstVtx_lid];
      nnzL += j;
      nnzU += k;

      for (vtx = fstVtx + 1, vtx_lid = fstVtx_lid + 1; 
	   vtx < lstVtx; vtx++, vtx_lid ++) {
	if (Pslu_freeable->supno_loc[vtx_lid] != nsuper) {
	  nsuper = Pslu_freeable->supno_loc[vtx_lid];
	  Pslu_freeable->xsup_end_loc[nsuper-1] = vtx;
	  Pslu_freeable->xsup_beg_loc[nsuper] = vtx;
	  szsn = 1;
	  j = Llu_symbfact.xlsub[vtx_lid + 1] - Llu_symbfact.xlsub[vtx_lid];
	  k = Llu_symbfact.xusub[vtx_lid + 1] - Llu_symbfact.xusub[vtx_lid];
	}
	else {
	  szsn ++;
	}
	nnzL += j - szsn + 1;
	nnzU += k - szsn + 1;
      }
      Pslu_freeable->xsup_end_loc[nsuper] = lstVtx;
    }
    Pslu_freeable->supno_loc[VInfo.nvtcs_loc] = nsuper_loc;
    Pslu_freeable->nvtcs_loc = VInfo.nvtcs_loc; 

    /* set up xsup data */
    Pslu_freeable->lsub = Llu_symbfact.lsub;
    Pslu_freeable->xlsub = Llu_symbfact.xlsub;
    Pslu_freeable->usub = Llu_symbfact.usub;
    Pslu_freeable->xusub = Llu_symbfact.xusub;
    Pslu_freeable->szLsub = Llu_symbfact.szLsub;
    Pslu_freeable->szUsub = Llu_symbfact.szUsub;
    
#if ( PROFlevel>=1 )
    t_symbFact_loc[1] = SuperLU_timer_() - t_symbFact_loc[1];
#endif  

#if ( PRNTlevel>=1 )
    estimate_memUsage (n, iam,  symb_mem_usage, 
		       &totalMemLU, &overestimMem, 
		       Pslu_freeable, &Llu_symbfact, &VInfo, &CS, &PS);
    stat_loc[0] = (float) nnzL;
    stat_loc[1] = (float) nnzU;  
    stat_loc[2] = (float) nsuper_loc;
    stat_loc[3] = (float) Pslu_freeable->xlsub[VInfo.nvtcs_loc];
    stat_loc[4] = (float) Pslu_freeable->xusub[VInfo.nvtcs_loc];
    stat_loc[5] = totalMemLU;
    stat_loc[6] = overestimMem;
    stat_loc[7] = totalMemLU - overestimMem;
    stat_loc[8] = (float) PS.maxSzBuf;
    stat_loc[9] = (float) PS.nDnsUpSeps;
    stat_loc[10] = (float) PS.nDnsCurSep;
    stat_loc[11] = (float) (Llu_symbfact.no_expand + Llu_symbfact.no_expcp +
			    Llu_symbfact.no_expand_pr);
    stat_loc[12] = (float) Llu_symbfact.no_expand;
    stat_loc[13] = (float) Llu_symbfact.no_expcp;
    stat_loc[14] = (float) Llu_symbfact.no_expand_pr;
    stat_loc[15] = (float) fill_rcmd;
    stat_loc[16] = PS.nops;
    stat_loc[17] = PS.fill_pelt[1];
    stat_loc[18] = PS.fill_pelt[4];
    stat_loc[19] = PS.fill_pelt[0];
    stat_loc[20] = PS.fill_pelt[2];
    stat_loc[21] = PS.fill_pelt[3];
    stat_loc[22] = PS.fill_pelt[5];
    
    MPI_Reduce (stat_loc, stat_glob, 23, MPI_FLOAT, 
		MPI_SUM, 0, (*symb_comm));
    MPI_Reduce (&(stat_loc[5]), mem_glob, 14, MPI_FLOAT, 
		MPI_MAX, 0, (*symb_comm));
    fill_rcmd = (int_t) mem_glob[10];
    PS.fill_pelt[0] = stat_glob[19];
    PS.fill_pelt[1] = mem_glob[12];
    PS.fill_pelt[2] = stat_glob[20];
    PS.fill_pelt[3] = stat_glob[21];
    PS.fill_pelt[4] = mem_glob[13];
    PS.fill_pelt[5] = stat_glob[22];
    if (PS.fill_pelt[2] == 0.) PS.fill_pelt[2] = 1.;
    if (PS.fill_pelt[5] == 0.) PS.fill_pelt[5] = 1.;
    
#if ( PROFlevel>=1 )
    MPI_Reduce (t_symbFact_loc, t_symbFact, 3, MPI_DOUBLE,
		MPI_MAX, 0, (*symb_comm));
    MPI_Gather (time_lvls, 3 * nlvls, MPI_DOUBLE,
		time_lvlsT, 3 * nlvls , MPI_DOUBLE,
		0, (*symb_comm));
#endif
    
    stat_msgs_l[0] = (float) PS.maxsz_msgSnd;
    stat_msgs_l[1] = (float) PS.maxsz_msgSnd;
    if (PS.maxsz_msgSnd < PS.maxsz_msgCol)
      stat_msgs_l[1] = PS.maxsz_msgCol;
    stat_msgs_l[2] = PS.no_shmSnd + PS.no_msgsSnd + 
      PS.no_shmRcvd + PS.no_msgsRcvd;
    stat_msgs_l[3] = stat_msgs_l[2] + PS.no_msgsCol;
    stat_msgs_l[4] = stat_msgs_l[2];
    stat_msgs_l[5] = stat_msgs_l[3]; 
    stat_msgs_l[6] = PS.no_msgsSnd;
    stat_msgs_l[7] = PS.no_msgsSnd + PS.no_msgsCol;  
    stat_msgs_l[8] = PS.sz_msgsSnd;
    stat_msgs_l[9] = PS.sz_msgsSnd + PS.sz_msgsCol;
    MPI_Reduce (stat_msgs_l, stat_msgs_g, 4, MPI_FLOAT,
		MPI_MAX, 0, (*symb_comm));
    MPI_Reduce (&(stat_msgs_l[4]), &(stat_msgs_g[4]), 6, MPI_FLOAT,
		MPI_SUM, 0, (*symb_comm));
    if (stat_msgs_g[6] == 0) stat_msgs_g[6] = 1;
    if (stat_msgs_g[7] == 0) stat_msgs_g[7] = 1;
    
    if (!iam) {
      nnzL   = (long long) stat_glob[0]; nnzU  = (long long) stat_glob[1];
      nsuper = (int_t) stat_glob[2];
      szLGr  = (int_t) stat_glob[3]; szUGr = (int_t) stat_glob[4];
      printf("\tMax szBlk          %ld\n", (long long) VInfo.maxSzBlk);
#if ( PRNTlevel>=2 )
      printf("\t relax_gen %.2f, relax_curSep %.2f, relax_seps %.2f\n",
	     PS.relax_gen, PS.relax_curSep, PS.relax_seps);
#endif
      printf("\tParameters: fill mem %ld fill pelt %ld\n",
	     (long long) sp_ienv_dist(6), (long long) PS.fill_par);
      printf("\tNonzeros in L       %ld\n", nnzL);
      printf("\tNonzeros in U       %ld\n", nnzU);
      nnzLU = nnzL + nnzU;
      printf("\tnonzeros in L+U-I   %ld\n", nnzLU);
      printf("\tNo of supers   %ld\n", (long long) nsuper);
      printf("\tSize of G(L)   %ld\n", (long long) szLGr);
      printf("\tSize of G(U)   %ld\n", (long long) szUGr);
      printf("\tSize of G(L+U) %ld\n", (long long) szLGr+szUGr);

      printf("\tParSYMBfact (MB)      :\tL\\U MAX %.2f\tAVG %.2f\n",
	     mem_glob[0]*1e-6, 
	     stat_glob[5]/nprocs_symb*1e-6);
#if ( PRNTlevel>=2 )
      printf("\tRL overestim (MB):\tL\\U MAX %.2f\tAVG %.2f\n",
	     mem_glob[1]*1e-6, 
	     stat_glob[6]/nprocs_symb*1e-6);
      printf("\tsnd/rcv buffers (MB):\tL\\U MAX %.2f\tAVG %.2f\n",
	     mem_glob[3]*1e-6, 
	     stat_glob[8]/nprocs_symb*1e-6);
      printf("\tSYMBfact 2*n+4*nvtcs_loc+2*maxNvtcsNds_loc:\tL\\U %.2f\n",
	     (float) (2 * n * sizeof(int_t)) *1e-6);
      printf("\tint_t %d, int %d, long int %d, short %d, float %d, double %d\n", 
	     sizeof(int_t), sizeof(int),  sizeof(long int), sizeof(short), sizeof(float),
	     sizeof(double));
      printf("\tDNS ALLSEPS:\t MAX %d\tAVG %.2f\n",
	     (int_t) mem_glob[4], stat_glob[9]/nprocs_symb);
      printf("\tDNS CURSEP:\t MAX %d\tAVG %.2f\n\n",
	     (int_t) mem_glob[5], stat_glob[10]/nprocs_symb);

      printf("\t MAX FILL Mem(L+U) / Mem(A) per processor %ld\n", fill_rcmd);    
      printf("\t      Per elt MAX %ld AVG %ld\n", 
	     (int_t) PS.fill_pelt[4], (int_t)(PS.fill_pelt[3]/PS.fill_pelt[5]));
      printf("\t      Per elt RL MAX %ld AVG %ld\n",
	     (int_t) PS.fill_pelt[1], (int_t)(PS.fill_pelt[0]/PS.fill_pelt[2]));
      printf("\tM Nops:\t MAX %.2f\tAVG %.2f\n",
	     mem_glob[11]*1e-6, (stat_glob[16]/nprocs_symb)*1e-6);
      
      
      printf("\tEXPANSIONS: MAX/AVG\n");
      printf("\tTOTAL: %d / %.2f\n",
	     (int_t) mem_glob[6], stat_glob[11]/nprocs_symb);
      printf("\tREALLOC: %.f / %.2f RL_CP %.f / %.2f PR_CP %.f / %.2f\n",
	     mem_glob[7], stat_glob[12]/nprocs_symb,
	     mem_glob[8], stat_glob[13]/nprocs_symb,
	     mem_glob[9], stat_glob[14]/nprocs_symb);
      
      printf ("\n\tDATA MSGS  noMsgs*10^3 %.3f/%.3f size (MB) %.3f/%.3f \n",
	      stat_msgs_g[2]*1e-3, stat_msgs_g[4]/nprocs_symb*1e-3,
	      stat_msgs_g[0]*1e-6, stat_msgs_g[8] / stat_msgs_g[6]*1e-6);
      printf ("\tTOTAL MSGS noMsgs*10^3 %.3f/%.3f size (MB) %.3f/%.3f \n",
	      stat_msgs_g[3]*1e-3, stat_msgs_g[5]/nprocs_symb*1e-3,
	      stat_msgs_g[1]*1e-6, stat_msgs_g[9]/stat_msgs_g[7]*1e-6);
#endif      

#if ( PROFlevel>=1 )
      printf("Distribute matrix time = %8.3f\n", t_symbFact[0]);
      printf("Count vertices time    = %8.3f\n", t_symbFact[2]);
      printf("Symbfact DIST time     = %8.3f\n", t_symbFact[1]);
      
      printf("\nLvl\t    Time\t    Init\t   Inter\t    Intra\n");
      time_lvlsg[0] = 0.;
      for (i = 0; i < nlvls; i++) {
	for (j = 1; j < 9; j++)
	  time_lvlsg[j] = 0.;
	for (p = 0; p < nprocs_symb; p++) {
	  k = p * 3 * nlvls;
	  t = time_lvlsT[i*3+k] + time_lvlsT[i*3+k+1] + time_lvlsT[i*3+k+2];
	  if (t > time_lvlsg[1]) {
	    time_lvlsg[1] = t; j = p;
	  }
	  time_lvlsg[2] += t;
	  if (time_lvlsT[i*3+k] > time_lvlsg[3])
	    time_lvlsg[3] = time_lvlsT[i*3+k];
	  time_lvlsg[4] += time_lvlsT[i*3+k];
	  if (time_lvlsT[i*3+k+1] > time_lvlsg[5])
	    time_lvlsg[5] = time_lvlsT[i*3+k+1];
	  time_lvlsg[6] += time_lvlsT[i*3+k+1];
	  if (time_lvlsT[i*3+k+2] > time_lvlsg[7])
	    time_lvlsg[7] = time_lvlsT[i*3+k+2];
	  time_lvlsg[8] += time_lvlsT[i*3+k+2];
	}
	time_lvlsg[0] += time_lvlsg[1];
	printf ("%d \t%.3f/%.3f\t%.3f/%.3f\t%.3f/%.3f\t%.3f/%.3f\n", i,
		time_lvlsg[1], time_lvlsg[2] / nprocs_symb,
		time_lvlsg[3], time_lvlsg[4] / nprocs_symb,
		time_lvlsg[5], time_lvlsg[6] /nprocs_symb,
		time_lvlsg[7], time_lvlsg[8] / nprocs_symb); 
      }
      printf("\t   %8.3f \n", time_lvlsg[0]);    
#endif
    }
#endif
#if ( PROFlevel>=1 )
    SUPERLU_FREE (time_lvls);
    SUPERLU_FREE (time_lvlsT);
#endif
    symbfact_free (iam, nprocs_symb, &Llu_symbfact, &VInfo, &CS);
  } /* if (iam < nprocs_symb) */  
  else {
    /* update Pslu_freeable before returning */
    Pslu_freeable->nvtcs_loc = 0; 
    Pslu_freeable->xlsub = NULL; Pslu_freeable->lsub = NULL; 
    Pslu_freeable->xusub = NULL; Pslu_freeable->usub = NULL; 
    Pslu_freeable->supno_loc = NULL;
    Pslu_freeable->xsup_beg_loc = NULL;     
    Pslu_freeable->xsup_end_loc = NULL;
    
    SUPERLU_FREE( tempArray );
    PS.allocMem -= n * sizeof(int_t);
  }

  if (iam < nprocs_symb && nprocs_symb != 1) 
    freeComm (iam, nprocs_symb, commLvls, symb_comm);     
  if (commLvls != NULL)
    SUPERLU_FREE( commLvls );
  
#if ( DEBUGlevel>=1 )
  CHECK_MALLOC(iam, "Exit psymbfact()");
#endif

  return (- PS.allocMem);
} /* SYMBFACT_DIST */


static int_t
initParmsAndStats
(
 psymbfact_stat_t *PS /* Output -statistics*/
)
/*! \brief
 * <pre> 
 * Purpose
 * =======
 * Initialize relaxation parameters and statistics variables
 * </pre>
 */
{
  int  i;

  PS->nDnsCurSep = 0;
  PS->nDnsUpSeps = 0;
  
  PS->relax_gen = 1.0;
  PS->relax_curSep = 1.0;
  PS->relax_seps = 1.0;
  PS->fill_par = sp_ienv_dist(6);
  PS->nops = 0.;
  PS->no_shmSnd = 0.;
  PS->no_msgsSnd = 0.;
  PS->maxsz_msgSnd = 0;
  PS->sz_msgsSnd = 0.;
  PS->no_shmRcvd = 0.;
  PS->no_msgsRcvd = 0.;
  PS->maxsz_msgRcvd = 0;
  PS->sz_msgsRcvd = 0.;
  PS->no_msgsCol = 0.;
  PS->maxsz_msgCol = 0;
  PS->sz_msgsCol = 0.;

  for (i = 0; i < 6; i++)
    PS->fill_pelt[i] = 0.;

  PS->estimUSz = 0;
  PS->estimLSz = 0;
  PS->maxSzLPr = 0;
  PS->maxSzUPr = 0;
  PS->maxSzBuf = 0;
  PS->szDnsSep = 0;  
  PS->allocMem = 0;

  return 0;
}

static float
cntsVtcs 
(
 int_t  n,           /* Input - order of the input matrix */
 int    iam,         /* Input - my processor number */
 int    nprocs_symb, /* Input - no of processors for symbolic factorization */
 Pslu_freeable_t *Pslu_freeable, /* Input -globToLoc and maxNvtcsPProc */
 Llu_symbfact_t  *Llu_symbfact, /* Input/Output -local L, U data structures */
 vtcsInfo_symbfact_t *VInfo,  /* Input - local info on vertices distribution */
 int_t            *tempArray, /* Input - temporary storage */
 int_t            *fstVtxSep, /* Input - first vertex of each node in the tree */
 int_t            *sizes,     /* Input - sizes of each node in the tree */
 psymbfact_stat_t *PS,  /* Input/Output -statistics */
 MPI_Comm         *commLvls
 )
/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 * 
 * Computes an estimation of the number of elements in columns of L
 * and rows of U.  Stores this information in cntelt_vtcs, and it will
 * be used in the right-looking symbolic factorization.
 * </pre>
 */
{
  int   fstP, lstP, szSep, npNode, i, j;
  int_t nvtcs_loc, ind_blk, vtx, vtx_lid, ii, jj, lv, vtx_elt, cur_blk;
  int_t fstVtx, lstVtx, fstVtx_blk, lstVtx_blk;
  int_t nelts, nelts_new_blk;
  int_t *xlsub, *lsub, *xusub, *usub, *globToLoc, maxNvtcsPProc;
  int_t *minElt_vtx, *cntelt_vtcs;
  
  /* Initialization */
  xlsub = Llu_symbfact->xlsub; lsub = Llu_symbfact->lsub;
  xusub = Llu_symbfact->xusub; usub = Llu_symbfact->usub;
  cntelt_vtcs = Llu_symbfact->cntelt_vtcs;
  globToLoc = Pslu_freeable->globToLoc;
  nvtcs_loc = VInfo->nvtcs_loc;
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  if (Llu_symbfact->szLsub - VInfo->nnz_ainf_loc > n)
    minElt_vtx = lsub;
  else { 
    /* allocate memory for minElt_vtx */
    if (!(minElt_vtx = intMalloc_dist(n))) {
      fprintf(stderr, "Malloc fails for minElt_vtx[].");
      return (PS->allocMem);
    }
    PS->allocMem += n * sizeof (int_t);
  } 
  
  for (ii = 0; ii < n; ii++) 
    tempArray[ii] = n;
  for (ii = 0; ii < nvtcs_loc; ii++)
    cntelt_vtcs[ii] = 0;

  szSep = nprocs_symb;
  i = 0;
  cur_blk = 0;
  vtx_lid = 0;
  while (szSep >= 1) {
    /* for each level in the separator tree */
    npNode = nprocs_symb / szSep; 
    fstP = 0; 
    /* for each node in the level */
    for (j = i; j < i + szSep; j++) {
      fstVtx = fstVtxSep[j];
      lstVtx  = fstVtx + sizes[j];
      lstP = fstP + npNode;

      if (fstP <= iam && iam < lstP) {      
	ind_blk = cur_blk;
	ii = vtx_lid;
	while (VInfo->begEndBlks_loc[ind_blk] < lstVtx && 
	       ind_blk < 2 * VInfo->nblks_loc) {	  
	  fstVtx_blk = VInfo->begEndBlks_loc[ind_blk];
	  lstVtx_blk = VInfo->begEndBlks_loc[ind_blk + 1];
	  ind_blk += 2;
	  for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, ii++) {
	    for (jj = xlsub[ii]; jj < xlsub[ii+1]; jj++) {
	      vtx_elt = lsub[jj];
	      if (tempArray[vtx_elt] == n) {
		tempArray[vtx_elt] = vtx;
	      }
	    }
	    for (jj = xusub[ii]; jj < xusub[ii+1]; jj++) {
	      vtx_elt = usub[jj];
	      if (tempArray[vtx_elt] == n) {
		tempArray[vtx_elt] = vtx;
	      }
	    }
	  }	  
	} 
	if (szSep == nprocs_symb) 
	  vtx_lid = ii;
	else {
	  MPI_Allreduce (&(tempArray[fstVtx]), &(minElt_vtx[fstVtx]), 
			 (int) (n - fstVtx), mpi_int_t, MPI_MIN, commLvls[j]);
#if ( PRNTlevel>=1 )
	  PS->no_msgsCol += (float) (2 * (int) LOG2( npNode ));
	  PS->sz_msgsCol += (float) (n - fstVtx);
	  if (PS->maxsz_msgCol < n - fstVtx) 
	    PS->maxsz_msgCol = n - fstVtx;      
#endif

	  nelts = 0;
	  for (ii = fstVtx; ii < lstVtx; ii++)
	    tempArray[ii] = 0;
	  for (ii = fstVtx; ii < n; ii++) {
	    if (minElt_vtx[ii] != n) {
	      if (minElt_vtx[ii] < fstVtx)
		nelts ++;
	      else
		tempArray[minElt_vtx[ii]] ++;
	      if (ii > lstVtx)
		tempArray[ii] = minElt_vtx[ii];
	    }
	  }
	
	  ind_blk = cur_blk;
	  lv = fstVtx;
	  while (VInfo->begEndBlks_loc[ind_blk] < lstVtx && 
		 ind_blk < 2 * VInfo->nblks_loc) {	  
	    fstVtx_blk = VInfo->begEndBlks_loc[ind_blk];
	    lstVtx_blk = VInfo->begEndBlks_loc[ind_blk + 1];
	    ind_blk += 2;
	    
	    for (ii = lv; ii < fstVtx_blk; ii++)
	      nelts += tempArray[ii];
	    lv = lstVtx_blk;

	    nelts_new_blk = 0;
	    for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid++) {
	      nelts_new_blk += tempArray[vtx];
	      cntelt_vtcs[vtx_lid] = nelts;
	    }
	    nelts += nelts_new_blk;
	  }
	} /* if (szSep != nprocs_symb) */
	cur_blk = ind_blk;
      }
      fstP += npNode;
    }
    i += szSep;
    szSep = szSep / 2;
  }
  /* free memory */
  if (minElt_vtx != lsub) {
    SUPERLU_FREE (minElt_vtx);
    PS->allocMem -= n * sizeof(int_t);
  }
  return (SUCCES_RET);
}

static float
symbfact_mapVtcs
(
 int iam,             /* Input -process number */
 int nprocs_num,      /* Input -number of processors */
 int nprocs_symb,     /* Input -number of procs for symbolic factorization */
 SuperMatrix *A,      /* Input -input distributed matrix A */
 int_t *fstVtxSep,    /* Input -first vertex in each separator */
 int_t *sizes,        /* Input -size of each separator in the separator tree */
 Pslu_freeable_t *Pslu_freeable, /* Output -globToLoc and maxNvtcsPProc 
				    computed */
 vtcsInfo_symbfact_t *VInfo, /* Output -local info on vertices distribution */
 int_t *tempArray,    /* Input -temp array of size n = order of the matrix */
 int_t  maxSzBlk,     /* Input -maximum number of vertices in a block */
 psymbfact_stat_t *PS /* Input/Output -statistics */
 ) 
{
/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 *  symbfact_mapVtcs maps the vertices of the graph of the input
 *  matrix A on nprocs_symb processors, using the separator tree
 *  returned by a graph partitioning algorithm from the previous step
 *  of the symbolic factorization.  The number of processors
 *  nprocs_symb must be a power of 2.
 *
 * Description of the algorithm
 * ============================
 *
 *  A subtree to subcube algorithm is used first to map the processors
 *  on the nodes of the separator tree.
 *
 *  For each node of the separator tree, its corresponding vertices
 *  are distributed on the processors affected to this node, using a
 *  block cyclic distribution.
 *
 *  After the distribution, fields of the VInfo structure are
 *  computed.  The array globToLoc and maxNvtcsPProc of Pslu_freeable
 *  are also computed.
 * </pre>
 */
  int szSep, npNode, firstP, p, iSep, jSep, ind_ap_s, ind_ap_d;
  int_t k, n, kk;
  int_t fstVtx, lstVtx;
  int_t fstVtxBlk, ind_blk;
  int_t noVtcsProc, noBlk;
  int_t nvtcs_loc; /* number of vertices owned by process iam */
  int_t nblks_loc; /* no of blocks owned by process iam */
  int_t *globToLoc;    /* global indexing to local indexing */
  int_t maxNvtcsPProc, maxNvtcsNds_loc, nvtcsNds_loc, maxNeltsVtx;
  int_t *begEndBlks_loc; /* begin and end vertex of each local block */
  int_t *vtcs_pe;  /* contains the number of vertices on each processor */
  int   *avail_pes; /* contains the processors to be used at each level */
  
  n = A->ncol;
  /* allocate memory */
  if (!(globToLoc = intMalloc_dist(n + 1))) {
    fprintf (stderr, "Malloc fails for globToLoc[].");
    return (PS->allocMem);
  }
  PS->allocMem += (n+1) * sizeof(int_t);
  if (!(avail_pes = (int *) SUPERLU_MALLOC(nprocs_symb*sizeof(int)))) {
    fprintf (stderr, "Malloc fails for avail_pes[].");
    return (PS->allocMem);
  }
  PS->allocMem += nprocs_symb*sizeof(int);
  if (!(vtcs_pe = (int_t *) SUPERLU_MALLOC(nprocs_symb*sizeof(int_t)))) {
    fprintf (stderr, "Malloc fails for vtcs_pe[].");
    return (PS->allocMem);
  }
  PS->allocMem += nprocs_symb*sizeof(int_t);
  
  /* Initialization */
  globToLoc[n] = n;  
  for (p = 0; p < nprocs_symb; p++) {
    vtcs_pe[p] = 0;
    avail_pes[p] = EMPTY;
  }
  nvtcs_loc = 0;
  nblks_loc = 0;
  maxNvtcsNds_loc = 0;
  maxNeltsVtx     = 0;
  
  /* distribute data among processors */
  szSep = nprocs_symb;
  iSep = 0;
  while (szSep >= 1) {
    /* for each level in the separator tree */
    npNode = nprocs_symb / szSep; 
    firstP = 0; 
    nvtcsNds_loc = 0;
    
    for (jSep = iSep; jSep < iSep + szSep; jSep++) {
      /* for each node in the level */
      fstVtx = fstVtxSep[jSep];
      lstVtx = fstVtx + sizes[jSep];
      if (firstP <= iam && iam < firstP + npNode)
	maxNeltsVtx += lstVtx - fstVtx;

      if (szSep == nprocs_symb) {
	/* leaves of the separator tree */
	for (k = fstVtx; k < lstVtx; k++) {
	  globToLoc[k] = (int_t) firstP;
	  vtcs_pe[firstP] ++;
	}
	if (firstP == iam) {	  
	  nvtcs_loc += lstVtx - fstVtx;
	  if (fstVtx != lstVtx)
	    nblks_loc ++;
	}
      }
      else {
	/* superior levels of the separator tree */
	k = fstVtx;
	noVtcsProc = maxSzBlk;
	fstVtxBlk = fstVtx;
	if ((jSep - iSep) % 2 == 0) ind_ap_d = (jSep - iSep) * npNode;
	/* first allocate processors from previous levels */	
	for (ind_ap_s = (jSep-iSep) * npNode; ind_ap_s < (jSep-iSep+1) * npNode; ind_ap_s ++) {
	  p = avail_pes[ind_ap_s];
	  if (p != EMPTY && k < lstVtx) {
	    /* for each column in the separator */	  
	    avail_pes[ind_ap_s] = EMPTY;
	    kk = 0;
	    while (kk < noVtcsProc && k < lstVtx) {
	      globToLoc[k] = p;
	      vtcs_pe[p] ++;
	      k ++;
	      kk ++;
	    }
	    if (p == iam) {
	      nvtcs_loc += kk;
	      nblks_loc ++;
	      nvtcsNds_loc += kk;
	    }
	  }
	  else {
	    if (p != EMPTY && k == lstVtx) {
	      avail_pes[ind_ap_s] = EMPTY;
	      avail_pes[ind_ap_d] = p; ind_ap_d ++;
	    }
	  }
	} 
	noBlk = 0;
	p = firstP + npNode;
	while (k < lstVtx) {
	  /* for each column in the separator */
	  kk = 0;
	  p = (int) (noBlk % (int_t) npNode) + firstP;
	  while (kk < noVtcsProc && k < lstVtx) {
	    globToLoc[k] = p;
	    vtcs_pe[p] ++;
	    k ++;
	    kk ++;
	  }
	  if (p == iam) {
	    nvtcs_loc += kk;
	    nblks_loc ++;
	    nvtcsNds_loc += kk;
	  }
	  noBlk ++;
	} /* while (k < lstVtx) */
	/* Add the unused processors to the avail_pes list of pes */
	for (p = p + 1; p < firstP + npNode; p ++) {
	  avail_pes[ind_ap_d] = p; ind_ap_d ++;
	}
      }
      firstP += npNode;
    }
    if (maxNvtcsNds_loc < nvtcsNds_loc && szSep != nprocs_symb)
      maxNvtcsNds_loc = nvtcsNds_loc;
    iSep += szSep;
    szSep = szSep / 2;
  }
  
#if ( PRNTlevel>=2 )
  if (!iam)
    PrintInt10 (" novtcs_pe", nprocs_symb, vtcs_pe);
#endif
  /* determine maximum number of vertices among processors */
  maxNvtcsPProc = vtcs_pe[0];
  vtcs_pe[0] = 0;
  for (p = 1; p < nprocs_symb; p++) {
    if (maxNvtcsPProc < vtcs_pe[p])
      maxNvtcsPProc = vtcs_pe[p];
    vtcs_pe[p] = 0;
  }
#if ( PRNTlevel>=2 )
  if (!iam)
    printf ("  MaxNvtcsPerProc %d MaxNvtcs/Avg %e\n\n", 
	    maxNvtcsPProc, ((float) maxNvtcsPProc * nprocs_symb)/(float)n);
#endif

  if (iam < nprocs_symb)
    if (!(begEndBlks_loc = intMalloc_symbfact(2 * nblks_loc + 1)))
      ABORT("Malloc fails for begEndBlks_loc[].");
  
  ind_blk = 0;
  k = 0;
  while (k < n) {
    p = globToLoc[k];
    if (p == iam) 
      begEndBlks_loc[ind_blk] = k;
    while (globToLoc[k] == p && k < n) {
      globToLoc[k] = globToLoc[k] * maxNvtcsPProc + vtcs_pe[p];
      vtcs_pe[p] ++;
      k ++;
    }
    if (p == iam) {
      begEndBlks_loc[ind_blk + 1] = k;
      ind_blk += 2;
    }
  }
  if (iam < nprocs_symb)
    begEndBlks_loc[2 * nblks_loc] = n;
 
  SUPERLU_FREE (avail_pes);
  SUPERLU_FREE (vtcs_pe);
  
  Pslu_freeable->maxNvtcsPProc   = maxNvtcsPProc;
  Pslu_freeable->globToLoc       = globToLoc;
  if (iam < nprocs_symb) {
    VInfo->maxNvtcsNds_loc = maxNvtcsNds_loc;
    VInfo->nblks_loc       = nblks_loc;
    VInfo->nvtcs_loc       = nvtcs_loc;
    VInfo->curblk_loc      = 0;
    VInfo->maxNeltsVtx     = maxNeltsVtx;
    VInfo->filledSep       = FALSE;
    VInfo->xlsub_nextLvl   = 0;
    VInfo->xusub_nextLvl   = 0;
    VInfo->begEndBlks_loc  = begEndBlks_loc;
    VInfo->fstVtx_nextLvl  = begEndBlks_loc[0];
  }
  return SUCCES_RET;
}

static void 
symbfact_distributeMatrix
(
 int   iam,             /* Input - my processor number */  
 int   nprocs_num,      /* Input - number of processors */
 int   nprocs_symb,     /* Input - number of processors for the
			   symbolic factorization */
 SuperMatrix *A,        /* Input - input matrix A */
 int_t *perm_c,         /* Input - column permutation */
 int_t *perm_r,         /* Input - row permutation */
 matrix_symbfact_t *AS, /* Output - temporary storage for the
			   redistributed matrix */
 Pslu_freeable_t *Pslu_freeable, /* Input - global to local information */
 vtcsInfo_symbfact_t *VInfo,  /* Input - local info on vertices
				 distribution */
 int_t  *tempArray,     /* Input/Output - temporary array of size n
			   (order of the matrix) */
 MPI_Comm    *num_comm  /* Input - communicator for nprocs_num procs */
 )
{
/*! \brief
 *
 * <pre>
 * Purpose 
 * =======
 *
 * Distribute input matrix A for the symbolic factorization routine.
 * Only structural information is distributed.  The redistributed
 * matrix has its rows and columns permuted according to perm_r and
 * perm_c. A is not modified during this routine.
 * </pre>
 */
/* Notations:
 * Ainf : inferior part of A, including diagonal.
 * Asup : superior part of A.
 */
  int p, p_irow, code_err, ainf_data;
  int_t n, m_loc, fst_row;
  int_t i, j, k, irow, jcol;
  NRformat_loc *Astore;
  int_t  nnz_loc, nnz_iam;    /* number of local nonzeros */
  int_t  nnz_remote; /* number of remote nonzeros to be sent */
  int_t  SendCnt; /* number of remote nonzeros to be sent */
  int_t  RecvCnt; /* number of remote nonzeros to be received */
  /* number of nonzeros to send/receive per processor */
  int_t  *nnzToSend, *nnzToRecv; 
  int_t *nnzAinf_toSnd; /* nnz in Ainf to send */
  /* VInfo data structures */
  int_t *globToLoc, *begEndBlks_loc, nblks_loc, nvtcs_loc, maxNvtcsPProc;
  
  int_t neltsRow, vtx, vtx_lid, nelts, ind;
  int_t *snd_aind, *rcv_aind;
  int_t *ptr_toSnd, *buf, *ptr_toRcv;
  /* matrix_symbfact_t *As data */
  int_t *x_ainf, *x_asup, *ind_ainf, *ind_asup;
  int  *intBuf1, *intBuf2, *intBuf3, *intBuf4;

  /* ------------------------------------------------------------
     INITIALIZATION.
     ------------------------------------------------------------*/
  Astore = (NRformat_loc *) A->Store;
  n = A->ncol;
  m_loc = Astore->m_loc;
  fst_row = Astore->fst_row;
  globToLoc = Pslu_freeable->globToLoc;
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  nnzToRecv = intCalloc_symbfact(3 * (int_t)nprocs_num);
  nnzToSend = nnzToRecv + nprocs_num;
  nnzAinf_toSnd = nnzToRecv + 2 * nprocs_num;

  /* --------------------------------------------------------------------- 
    COUNT THE NUMBER OF NONZEROS TO BE SENT TO EACH PROCESS, THEN ALLOCATE
    SPACE.  THIS ACCOUNTS FOR THE FIRST PASS OF A.
    ----------------------------------------------------------------------*/
  /* tempArray stores the number of nonzeros in each column of ainf */
  for (i = 0; i < n; i++)
    tempArray[i] = 0;
  for (i = 0; i < m_loc; i++) {
    irow   = perm_c[perm_r[i+fst_row]];  /* Row number in Pc*Pr*A */
    p_irow = OWNER(globToLoc[irow]);
    neltsRow = 0;

    for (j = Astore->rowptr[i]; j < Astore->rowptr[i+1]; j++) {
      jcol = perm_c[Astore->colind[j]];
      if (jcol <= irow) {
	p = OWNER(globToLoc[jcol]);
	if (tempArray[jcol] == 0) {
	  nnzToSend[p] += 2;
	  nnzAinf_toSnd[p] += 2;
	}
	tempArray[jcol] ++;
	nnzAinf_toSnd[p] ++;
      }
      else {
	p = p_irow;
	neltsRow ++;
      }
      nnzToSend[p] ++; 
    }
    if (neltsRow != 0) {
      nnzToSend[p_irow] += 2;
    }
  }
  
  /* add one entry which will separate columns of Ainf from rows
     of Asup */
  for (p = 0; p < nprocs_num; p++)
    if (nnzToSend[p] != 0)
      nnzToSend[p] ++;
  
  /* All-to-all communication */
  MPI_Alltoall (nnzToSend, 1, mpi_int_t, nnzToRecv, 1, mpi_int_t,
		(*num_comm));

  nnz_loc = SendCnt = RecvCnt = 0;
  for (p = 0; p < nprocs_num; p++) {
    if ( p != iam ) {
      SendCnt += nnzToSend[p];
      RecvCnt += nnzToRecv[p];
    } else {
      nnz_loc += nnzToRecv[p];
      nnzToSend[p] = 0;
    }
  }
  nnz_iam = nnz_loc + RecvCnt; /* Total nonzeros ended up in my process. */
  
  /* Allocate temporary storage for sending/receiving the A triplets. */
  if (!(snd_aind = intMalloc_symbfact(SendCnt)) && SendCnt != 0)
    ABORT("Malloc fails for snd_aind[].");
  if ( !(rcv_aind = intMalloc_symbfact(nnz_iam + 1)))
    ABORT("Malloc fails for rcv_aind[].");
  if ( !(ptr_toSnd = intCalloc_symbfact((int_t) nprocs_num)) )
    ABORT("Malloc fails for ptr_toSnd[].");
  if ( !(ptr_toRcv = intCalloc_symbfact((int_t) nprocs_num)) )
    ABORT("Malloc fails for ptr_toRcv[].");

  /* setup ptr_toSnd[p] to point to data in snd_aind to be send to 
   processor p */
  for (i = 0, j = 0, p = 0; p < nprocs_num; p++) {
    if ( p != iam ) 
      ptr_toSnd[p] = i;
    else
      ptr_toSnd[p] = j;
    i += nnzToSend[p]; 
    j += nnzToRecv[p];
  }

  for (i = 0; i < n; i++) {
    if (tempArray[i] != 0) {
      /* column i of Ainf will be send to a processor  */
      p = OWNER( globToLoc[i] );
      if (p == iam) {
	buf = &(rcv_aind[ptr_toSnd[p]]);
      }
      else {
	buf = &(snd_aind[ptr_toSnd[p]]);
      }
      buf[0] = tempArray[i];
      buf[1] = i;
      tempArray[i] = ptr_toSnd[p] + 2;
      ptr_toSnd[p] += 2 + buf[0];
    }
  }

  /* set ptr_toSnd to point to Asup data (stored by rows) */
  for (i = 0, j = 0, p = 0; p < nprocs_num; p++) {
    if ( p != iam ) {
      if (nnzToSend[p] != 0) { 
	snd_aind[i + nnzAinf_toSnd[p]] = EMPTY;
	ptr_toSnd[p] = i + nnzAinf_toSnd[p] + 1;
      }
    }
    else {
      if (nnzToRecv[p] != 0) {
	rcv_aind[j + nnzAinf_toSnd[p]] = EMPTY;
	ptr_toSnd[p] = j + nnzAinf_toSnd[p] + 1;
      }
    }
    i += nnzToSend[p]; 
    j += nnzToRecv[p];
  }

  /* ------------------------------------------------------------
     LOAD THE ENTRIES OF A INTO THE snd_aind STRUCTURE TO SEND.
     THIS ACCOUNTS FOR THE SECOND PASS OF A.
     For each processor, we store first the columns to be sent,
     and then the rows to be sent. For each row/column sent:
     entry 0            : x = number of elements in that row/column
     entry 1            : row/column number
     entries 2 .. x + 2 : row/column indices.
     ------------------------------------------------------------*/
  for (i = 0; i < m_loc; i++) {
    irow = perm_c[perm_r[i+fst_row]];  /* Row number in Pc*A */
    p_irow = OWNER( globToLoc[irow] );
    ptr_toSnd[p_irow] +=2;
    neltsRow = 0;
    for (j = Astore->rowptr[i]; j < Astore->rowptr[i+1]; j++) {
      jcol = perm_c[Astore->colind[j]];
      if (jcol <= irow) {
	p = OWNER( globToLoc[jcol] );
	k = tempArray[jcol];
	tempArray[jcol] ++;
	if (p == iam) { /* local */
	  rcv_aind[k] = irow; 
	}
	else {
	  snd_aind[k] = irow;
	}
      }
      else {
	p = p_irow;
	neltsRow ++;
	k = ptr_toSnd[p];
	ptr_toSnd[p] ++;
	if (p == iam) { /* local */
	  rcv_aind[k] = jcol;
	}
	else {
	  snd_aind[k] = jcol;
	}
      }
    }

    if (neltsRow == 0)
      ptr_toSnd[p_irow] -= 2;
    else {
      /* store entry 0 and entry 1 */
      if (p_irow == iam) { /* local */
	rcv_aind[ptr_toSnd[p_irow] - neltsRow - 2] = neltsRow;
	rcv_aind[ptr_toSnd[p_irow] - neltsRow - 1] = irow;
      }
      else { /* remote */
	snd_aind[ptr_toSnd[p_irow] - neltsRow - 2] = neltsRow;
	snd_aind[ptr_toSnd[p_irow] - neltsRow - 1] = irow;
      }
    }
  }
  
  /* reset ptr_toSnd to point to the beginning of the data for
     each processor (structure needed in MPI_Alltoallv */
  for (i = 0, j = 0, p = 0; p < nprocs_num; p++) {
    ptr_toSnd[p] = i;
    i += nnzToSend[p];
    ptr_toRcv[p] = j;
    j += nnzToRecv[p];
  }
  
  /* ------------------------------------------------------------
     PERFORM REDISTRIBUTION. THIS INVOLVES ALL-TO-ALL COMMUNICATION.
     Note: it uses MPI_Alltoallv.
     ------------------------------------------------------------*/
  if (nprocs_num > 1) {
#if defined (_LONGINT)
    intBuf1 = (int *) SUPERLU_MALLOC(4 * nprocs_num * sizeof(int));
    intBuf2 = intBuf1 + nprocs_num;
    intBuf3 = intBuf1 + 2 * nprocs_num;
    intBuf4 = intBuf1 + 3 * nprocs_num;
    
    for (p=0; p<nprocs_num; p++) {
      if (nnzToSend[p] > INT_MAX || ptr_toSnd[p] > INT_MAX ||
	  nnzToRecv[p] > INT_MAX || ptr_toRcv[p] > INT_MAX)
	ABORT("ERROR in symbfact_distributeMatrix size to send > INT_MAX\n");
      intBuf1[p] = (int) nnzToSend[p];
      intBuf2[p] = (int) ptr_toSnd[p];
      intBuf3[p] = (int) nnzToRecv[p];
      intBuf4[p] = (int) ptr_toRcv[p];
    }
    intBuf1[iam]=0; /* This corresponds to nnzToSend[iam] */
    intBuf3[iam]=0; /* This corresponds to nnzToRecv[iam] */
#else  /* Default */
    intBuf1 = nnzToSend;  intBuf2 = ptr_toSnd;
    intBuf3 = nnzToRecv;  intBuf4 = ptr_toRcv;
    i = nnzToRecv[iam]; 
    nnzToRecv[iam] = 0;
    nnzToSend[iam] = 0;
#endif

    MPI_Alltoallv (snd_aind, intBuf1, intBuf2, mpi_int_t, 
		   rcv_aind, intBuf3, intBuf4, mpi_int_t,
		   (*num_comm));

#if defined (_LONGINT)
    SUPERLU_FREE (intBuf1);
#else  /* Default */
    nnzToRecv[iam] = i;
#endif
  }
  
  /* ------------------------------------------------------------
     DEALLOCATE SEND STORAGE
     ------------------------------------------------------------*/
  if (snd_aind) SUPERLU_FREE( snd_aind );
  SUPERLU_FREE( ptr_toSnd );

  /* ------------------------------------------------------------
     CONVERT THE RECEIVED FORMAT INTO THE SYMBOLIC FORMAT.
     THIS IS PERFORMED ONLY BY NPROCS_SYMB PROCESSORS
     ------------------------------------------------------------*/
  if (iam < nprocs_symb) {
    nblks_loc = VInfo->nblks_loc;
    begEndBlks_loc = VInfo->begEndBlks_loc;
    nvtcs_loc = VInfo->nvtcs_loc;
    /* ------------------------------------------------------------
       Allocate space for storing indices of A after redistribution.
       ------------------------------------------------------------*/
    if (!(x_ainf = intCalloc_symbfact (nvtcs_loc + 1)))
      ABORT("Malloc fails for x_ainf[].");
    if (!(x_asup = intCalloc_symbfact (nvtcs_loc + 1)))
      ABORT("Malloc fails for x_asup[].");
    
    /* Initialize the array of columns/rows pointers */
    for (i = 0, p = 0; p < nprocs_num; p++) {
      ainf_data = TRUE;
      k = 0;
      while (k < nnzToRecv[p]) {
	j = rcv_aind[i + k];
	if (j == EMPTY) {
	  ainf_data = FALSE;
	  k ++;
	}
	else {
	  nelts = rcv_aind[i + k];
	  vtx = rcv_aind[i + k + 1];
	  vtx_lid = LOCAL_IND( globToLoc[vtx] );
	  k += nelts + 2;
	  if (ainf_data) 
	    x_ainf[vtx_lid] += nelts; 
	  else 
	    x_asup[vtx_lid] = nelts;
	}
      }
      i += nnzToRecv[p];
    }
    
    /* copy received information */
    vtx_lid = 0;
    for (i = 0, k = 0, j = 0; i < nblks_loc; i++) {
      for (vtx = begEndBlks_loc[2*i]; vtx < begEndBlks_loc[2*i+1]; vtx++, vtx_lid ++) {
	nelts = x_ainf[vtx_lid];
	x_ainf[vtx_lid] = k;
	k += nelts;
	nelts = x_asup[vtx_lid];
	x_asup[vtx_lid] = j;
	j += nelts;
	tempArray[vtx] = x_ainf[vtx_lid];
      }
    }
    x_ainf[nvtcs_loc] = k;
    x_asup[nvtcs_loc] = j;
    
    /* Allocate space for storing indices of A after conversion */
    if ( !(ind_ainf = intMalloc_symbfact(x_ainf[nvtcs_loc])) && x_ainf[nvtcs_loc] != 0 )
      ABORT("Malloc fails for ind_ainf[].");
    if ( !(ind_asup = intMalloc_symbfact(x_asup[nvtcs_loc])) && x_asup[nvtcs_loc] != 0)
      ABORT("Malloc fails for ind_asup[].");
    
    /* Copy the data into the row/column oriented storage */  
    for (i = 0, p = 0; p < nprocs_num; p++) {
      ainf_data = TRUE;
      k = 0;
      while (k < nnzToRecv[p]) {
	j = rcv_aind[i + k];
	if (ainf_data && j == EMPTY) {
	  ainf_data = FALSE;
	  k ++;
	}
	else {
	  nelts = rcv_aind[i + k];
	  vtx = rcv_aind[i + k + 1];
	  vtx_lid = LOCAL_IND( globToLoc[vtx] );
	  if (ainf_data) {
	    /* traverse ainf data */
	    ind = tempArray[vtx];
	    for (j = i + k + 2; j < i + k + 2 + nelts; j++, ind ++) 
	      ind_ainf[ind] = rcv_aind[j];
	    tempArray[vtx] = ind;
	  }
	  else {
	    /* traverse asup data */
	    ind = x_asup[vtx_lid];
	    for (j = i + k + 2; j < i + k + 2 + nelts; j++, ind ++) 
	      ind_asup[ind] = rcv_aind[j];
	  }
	  k += nelts + 2;
	}
      }
      i += nnzToRecv[p];
    }
    
    /* ------------------------------------------------------------
       DEALLOCATE TEMPORARY STORAGE
       ------------------------------------------------------------*/
    SUPERLU_FREE( ptr_toRcv );
    if (rcv_aind) SUPERLU_FREE( rcv_aind );
    if (nnzToRecv) SUPERLU_FREE( nnzToRecv );

    AS->x_ainf = x_ainf;
    AS->x_asup = x_asup;
    AS->ind_ainf = ind_ainf;
    AS->ind_asup = ind_asup;
    
    VInfo->nnz_asup_loc = x_asup[nvtcs_loc];
    VInfo->nnz_ainf_loc = x_ainf[nvtcs_loc];
  }
}

static
float allocPrune_lvl
(
 Llu_symbfact_t *Llu_symbfact, /* Input/Output - local L, U data
				  structures */
 vtcsInfo_symbfact_t *VInfo,   /* Input -local info on vertices
				  distribution */
 psymbfact_stat_t *PS          /* Input -statistics */
 )
/*! \brief
 *
 * <pre>
 * Allocate storage for data structures necessary for pruned graphs.
 * For those unpredictable size, make a guess as FILL * n.
 * Return value:
 *     0 if enough memory was available;
 *     otherwise, return the amount of space intended to allocate 
 *     when memory allocation failure occurred.
 * </pre>
 */
{
  int_t  lword;
  int_t  nzlmaxPr, nzumaxPr, *xlsubPr, *xusubPr, *lsubPr, *usubPr;
  int_t  nvtcs_loc, no_expand_pr, x_sz;
  float  alpha = 1.5;
  int_t  FILL = sp_ienv_dist(6);
  
  nvtcs_loc = VInfo->nvtcs_loc;
  
  no_expand_pr = 0;
  lword     = (int_t) sizeof(int_t);
  
  /* free memory allocated for the domain symbolic factorization */
  if (Llu_symbfact->szLsubPr)
    SUPERLU_FREE( Llu_symbfact->lsubPr );
  if (Llu_symbfact->szUsubPr)
    SUPERLU_FREE( Llu_symbfact->usubPr );
  if (Llu_symbfact->xlsubPr)
    SUPERLU_FREE( Llu_symbfact->xlsubPr );
  if (Llu_symbfact->xusubPr)
    SUPERLU_FREE( Llu_symbfact->xusubPr );
  
  Llu_symbfact->xlsub_rcvd = intMalloc_symbfact (VInfo->maxSzBlk + 1);
  Llu_symbfact->xusub_rcvd = intMalloc_symbfact (VInfo->maxSzBlk + 1);

  /* allocate memory to use during superior levels of sep_tree */
  x_sz = SUPERLU_MIN( VInfo->maxNvtcsNds_loc, VInfo->maxSzBlk);
  nzlmaxPr = 2 * FILL * VInfo->maxNvtcsNds_loc;
  nzumaxPr = 2 * FILL * VInfo->maxSzBlk;  

  /* Integer pointers for L\U factors */
  if (x_sz != 0) {
    xlsubPr   = intMalloc_symbfact(VInfo->maxNvtcsNds_loc + 1);
    xusubPr   = intMalloc_symbfact(VInfo->maxNvtcsNds_loc + 1);
    
    lsubPr = (int_t *) SUPERLU_MALLOC (nzlmaxPr * lword);
    usubPr = (int_t *) SUPERLU_MALLOC (nzumaxPr * lword);
    
    while ( !lsubPr || !usubPr ) {
      if ( lsubPr ) SUPERLU_FREE( lsubPr ); 
      if ( usubPr ) SUPERLU_FREE( usubPr );
      
      nzlmaxPr /= 2;     nzlmaxPr = alpha * (float) nzlmaxPr;
      nzumaxPr /= 2;     nzumaxPr = alpha * (float) nzumaxPr;
      
      if ( nzumaxPr < x_sz ) {
	fprintf(stderr, "Not enough memory to perform factorization.\n");
	return (PS->allocMem);
      }
      lsubPr  = (int_t *) SUPERLU_MALLOC(nzlmaxPr * lword);
      usubPr  = (int_t *) SUPERLU_MALLOC(nzumaxPr * lword);
      ++no_expand_pr;
    }
  }    
  else {
    xlsubPr = NULL; lsubPr = NULL;
    xusubPr = NULL; usubPr = NULL;
    nzlmaxPr = 0; nzumaxPr = 0;
  }
  
  if (VInfo->maxNvtcsNds_loc)
    Llu_symbfact->cntelt_vtcsA_lvl = 
      (int_t *) SUPERLU_MALLOC (VInfo->maxNvtcsNds_loc * lword);

  if (PS->maxSzLPr < Llu_symbfact->indLsubPr)
    PS->maxSzLPr = Llu_symbfact->indLsubPr;
  if (PS->maxSzUPr < Llu_symbfact->indUsubPr)
    PS->maxSzUPr = Llu_symbfact->indUsubPr;
  
  Llu_symbfact->lsubPr   = lsubPr;
  Llu_symbfact->xlsubPr  = xlsubPr;
  Llu_symbfact->usubPr   = usubPr;
  Llu_symbfact->xusubPr  = xusubPr;
  Llu_symbfact->szLsubPr = nzlmaxPr;
  Llu_symbfact->szUsubPr = nzumaxPr;
  Llu_symbfact->indLsubPr = 0;
  Llu_symbfact->indUsubPr = 0;

  Llu_symbfact->no_expand_pr += no_expand_pr;
  return 0;
}

static float 
allocPrune_domain
(
 int_t fstVtx,  /* Input - first vertex of current node */ 
 int_t lstVtx,  /* Input - last vertex of current node */
 Llu_symbfact_t *Llu_symbfact, /* Output - local L, U data
				  structures */
 vtcsInfo_symbfact_t *VInfo,   /* Input -local info on vertices
				  distribution */
 psymbfact_stat_t *PS           /* Input -statistics */
 )
/*! \brief
 *
 * <pre>
 * Allocate storage for data structures necessary for pruned graphs.
 * For those unpredictable size, make a guess as FILL * n.
 * Return value:
 *     0 if enough memory was available;
 *     otherwise, return the amount of space intended to allocate 
 *     when memory allocation failure occurred.
 * </pre>
 */
{
  int_t  lword;
  int_t  nzlmaxPr, nzumaxPr, *xlsubPr, *xusubPr, *lsubPr, *usubPr;
  int_t  nvtcs_loc, no_expand_pr, x_sz;
  float  alpha = 1.5;
  int_t  FILL = 2 * sp_ienv_dist(6);
  
  nvtcs_loc = VInfo->nvtcs_loc;
  
  no_expand_pr = 0;
  lword     = (int_t) sizeof(int_t);
  
  /* allocate memory to use during domain_symbolic routine */
  /* Guess for prune graph */
  x_sz = lstVtx - fstVtx;
  nzlmaxPr = nzumaxPr = 2*FILL * x_sz;
  
  /* Integer pointers for L\U factors */
  if (x_sz != 0) {
    xlsubPr   = intMalloc_symbfact(x_sz+1);
    xusubPr   = intMalloc_symbfact(x_sz+1);
    
    lsubPr = (int_t *) SUPERLU_MALLOC (nzlmaxPr * lword);
    usubPr = (int_t *) SUPERLU_MALLOC (nzumaxPr * lword);
    
    while ( !lsubPr || !usubPr ) {
      if ( lsubPr ) SUPERLU_FREE(lsubPr); 
      if ( usubPr ) SUPERLU_FREE(usubPr);
      
      nzlmaxPr /= 2;     nzlmaxPr = alpha * (float) nzlmaxPr;
      nzumaxPr /= 2;     nzumaxPr = alpha * (float) nzumaxPr;
      
      if ( nzumaxPr < x_sz ) {
	fprintf(stderr, "Not enough memory to perform factorization.\n");
	return (PS->allocMem);
      }
      lsubPr  = (void *) SUPERLU_MALLOC(nzlmaxPr * lword);
      usubPr  = (void *) SUPERLU_MALLOC(nzumaxPr * lword);
      ++no_expand_pr;
    }
  }    
  else {
    xlsubPr = NULL;
    xusubPr = NULL;
  }
  
  Llu_symbfact->lsubPr   = lsubPr;
  Llu_symbfact->xlsubPr  = xlsubPr;
  Llu_symbfact->usubPr   = usubPr;
  Llu_symbfact->xusubPr  = xusubPr;
  Llu_symbfact->szLsubPr = nzlmaxPr;
  Llu_symbfact->szUsubPr = nzumaxPr;
  Llu_symbfact->indLsubPr = 0;
  Llu_symbfact->indUsubPr = 0;
  Llu_symbfact->xlsub_rcvd = NULL;
  Llu_symbfact->xusub_rcvd = NULL;
  Llu_symbfact->cntelt_vtcsA_lvl = NULL;

  PS->maxSzLPr = Llu_symbfact->indLsubPr;
  PS->maxSzUPr = Llu_symbfact->indUsubPr;

  Llu_symbfact->no_expand_pr = no_expand_pr;
  Llu_symbfact->no_expcp = 0;
  return 0;
}

/************************************************************************/
static
int symbfact_alloc
/************************************************************************/
(
 int_t n,       /* Input - order of the matrix */
 int   nprocs,  /* Input - number of processors for the symbolic
		   factorization */  
 Pslu_freeable_t *Pslu_freeable, 
 Llu_symbfact_t *Llu_symbfact, /* Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo,   /* Input - local info on vertices
				  distribution */
 comm_symbfact_t *CS, /* Input -information on communication */
 psymbfact_stat_t *PS /* Input -statistics */
 )
/*! \brief
 *
 * <pre>
 * Allocate storage for the data structures common to symbolic factorization
 * routines. For those unpredictable size, make a guess as FILL * nnz(A).
 * Return value:
 *     0 if enough memory was available;
 *     otherwise, return the amount of space intended to allocate 
 *     when memory allocation failure occurred.
 * </pre>
 */
{
  int    nlvls, p;  /* no of levels in the separator tree */
  int_t  lword, no_expand;
  int_t  *xsup, *supno;
  int_t  *lsub, *xlsub;
  int_t  *usub, *xusub;
  int_t  nzlmax, nzumax, nnz_a_loc;
  int_t  nvtcs_loc, *cntelt_vtcs;
  float  alpha = 1.5;
  int_t  FILL = sp_ienv_dist(6);
  
  nvtcs_loc = VInfo->nvtcs_loc;
  nnz_a_loc = VInfo->nnz_ainf_loc + VInfo->nnz_asup_loc;
  nlvls = (int) LOG2( nprocs ) + 1;
  no_expand = 0;
  lword     = sizeof(int_t);
  
  /* Guess for L\U factors */
  nzlmax = nzumax = FILL * nnz_a_loc + 1;
  
  /* Integer pointers for L\U factors */
  supno  = intMalloc_symbfact(nvtcs_loc+1);
  xlsub  = intMalloc_symbfact(nvtcs_loc+1);
  xusub  = intMalloc_symbfact(nvtcs_loc+1);
  
  lsub = (void *) SUPERLU_MALLOC(nzlmax * lword);
  usub = (void *) SUPERLU_MALLOC(nzumax * lword);
  
  while ( !lsub || !usub ) {
    if (!lsub) SUPERLU_FREE(lsub); 
    if (!usub) SUPERLU_FREE(usub);
    
    nzlmax /= 2;     nzlmax = alpha * nzlmax;
    nzumax /= 2;     nzumax = alpha * nzumax;
    
    if ( nzumax < nnz_a_loc/2 ) {
      fprintf(stderr, "Not enough memory to perform factorization.\n");
      return (PS->allocMem);
    }
    lsub  = (void *) SUPERLU_MALLOC(nzlmax * lword);
    usub  = (void *) SUPERLU_MALLOC(nzumax * lword);
    ++no_expand;
  }
  
  if (nprocs == 1)
    cntelt_vtcs = NULL;
  else 
    cntelt_vtcs = intMalloc_symbfact (nvtcs_loc+1);
  
  /* allocate memory for communication data structures */
  CS->rcv_interLvl = intMalloc_symbfact (2 * (int_t) nprocs + 1);
  CS->snd_interLvl = intMalloc_symbfact (2 * (int_t) nprocs + 1);
  CS->ptr_rcvBuf   = intMalloc_symbfact (2 * (int_t) nprocs );
  CS->rcv_intraLvl = intMalloc_symbfact ((int_t) nprocs + 1);
  CS->snd_intraLvl = intMalloc_symbfact ((int_t) nprocs + 1);
  
  CS->snd_interSz  = intMalloc_symbfact ((int_t) nlvls + 1);
  CS->snd_LinterSz = intMalloc_symbfact ((int_t) nlvls + 1);  
  CS->snd_vtxinter = intMalloc_symbfact ((int_t) nlvls + 1);  
  CS->rcv_bufSz    = 0;
  CS->rcv_buf      = NULL;
  CS->snd_bufSz    = 0;
  CS->snd_buf      = NULL;

  for (p = 0; p < nprocs; p++) {
    CS->rcv_interLvl[p] = EMPTY;
    CS->snd_interLvl[p] = EMPTY;
    CS->rcv_intraLvl[p] = EMPTY;
    CS->snd_intraLvl[p] = EMPTY;
  }
  
  for (p = 0; p <= nlvls; p++) {
    CS->snd_vtxinter[p] = EMPTY;
    CS->snd_interSz[p]  = 0;
    CS->snd_LinterSz[p] = 0;
  }
  
  Pslu_freeable->supno_loc   = supno;
  Llu_symbfact->lsub   = lsub;
  Llu_symbfact->xlsub  = xlsub;
  Llu_symbfact->usub   = usub;
  Llu_symbfact->xusub  = xusub;
  Llu_symbfact->szLsub = nzlmax;
  Llu_symbfact->szUsub = nzumax;
  Llu_symbfact->cntelt_vtcs = cntelt_vtcs;
  
  Llu_symbfact->no_expand = no_expand;  
  
  return SUCCES_RET;
} /* SYMBFACT_ALLOC */

static int_t 
symbfact_vtx
(
 int_t n,         /* Input - order of the matrix */
 int   iam,       /* Input - my processor number */
 int_t vtx,       /* Input - vertex number to perform symbolic factorization */
 int_t vtx_lid,   /* Input - local vertex number */
 int_t vtx_prid,  /* Input - */
 int_t computeL,  /* Input - TRUE when compute column L(:,vtx)
		             otheriwse compute row U(vtx, :) */
 int   domain_symb,  /* Input - if TRUE, computation corresponds to the independent
			domain at the bottom of the separator tree */
 int_t fstVtx,       /* Input - first vertex of current node */ 
 int_t lstVtx,       /* Input - last vertex of current node */
 int_t snrep_lid,    /* local index of current supernode reprezentative */
 int_t szSn,         /* size of supernode with snrep_lid reprezentative */
 int_t *p_next,      /* next element in sub structure */
 int_t *marker,      
 int_t *sub_rcvd,    /* elements of node */
 int_t sub_rcvd_sz,  /* size of sub to be explored */
 Pslu_freeable_t *Pslu_freeable,
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo,    /* Input/Output - local info on vertices distribution */
 psymbfact_stat_t *PS,
 int_t *p_neltsVtxInit,
 int_t *p_neltsVtx,
 int_t *p_neltsVtx_CSep,
 int_t *p_neltsZrVtx,
 int_t *p_neltsMatched,
 int_t mark_vtx,
 int_t *p_prval_curvtx,
 int_t vtx_bel_othSn,
 int_t *p_vtx_bel_mySn
 )
{ 
  int_t x_aind_beg, x_aind_end;
  int_t k, vtx_elt, ind, pr, pr_lid, mem_error, ii, jj, compRcvd;
  int_t *xsub, *sub, *xsubPr, *subPr, *xsub_rcvd, *xsub_src, *sub_src;
  int_t pr_elt, next, prval_curvtx, maxNvtcsPProc;
  int_t  neltsVtx, neltsMatched, neltsZrVtx, neltsZrSn, neltsVtx_CSep;
  int_t  neltsVtxInit, kk;
  int   diagind, upd_lstSn;
  
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  upd_lstSn     = FALSE;
  diagind       = FALSE;
  prval_curvtx  = *p_prval_curvtx;
  neltsVtx_CSep = 0;
  next = *p_next;
  if (computeL) {
    xsub = Llu_symbfact->xlsub; sub = Llu_symbfact->lsub;
    xsub_rcvd = Llu_symbfact->xlsub_rcvd;
    xsubPr = Llu_symbfact->xusubPr; subPr = Llu_symbfact->usubPr;
  }
  else {
    xsub = Llu_symbfact->xusub; sub = Llu_symbfact->usub;
    xsub_rcvd = Llu_symbfact->xusub_rcvd;
    xsubPr = Llu_symbfact->xlsubPr; subPr = Llu_symbfact->lsubPr;
  }

  x_aind_beg = xsub[vtx_lid];
  x_aind_end = xsub[vtx_lid + 1];
  xsub[vtx_lid] = next;
  k = x_aind_beg;
  /* while (sub[k] != EMPTY && k < x_aind_end) { */
  while (k < x_aind_end) {
    if (sub[k] == EMPTY)
      k = x_aind_end;
    else {
      vtx_elt = sub[k];
      if (!computeL)
	if (marker[vtx_elt] == mark_vtx - 2)
	  if (vtx_elt < prval_curvtx)
	    prval_curvtx = vtx_elt;
      marker[vtx_elt] = mark_vtx;
      if (computeL && vtx_elt == vtx)
	diagind = TRUE;
      if (!computeL && vtx_elt == vtx)
	printf ("Pe[%d] ERROR diag elt in U part vtx " IFMT " dom_s %d fstV "
		IFMT " lstV " IFMT "\n", 
		iam, vtx, domain_symb, fstVtx, lstVtx);
      else {
	sub[next] = vtx_elt; 
	next ++;
      }
      if (vtx_elt < lstVtx) neltsVtx_CSep ++;
      k++;
    }
  }
  neltsVtxInit = k - x_aind_beg;
  PS->nops += neltsVtxInit;
  
  if (domain_symb) {
    if (computeL)
      VInfo->nnz_ainf_loc -= x_aind_end - x_aind_beg;
    else    
      VInfo->nnz_asup_loc -= x_aind_end - x_aind_beg;
  }

#ifdef TEST_SYMB
  printf ("compL %d vtx %d vtx_lid %d vtx_prid %d vtx_bel_othSn %d\n", 
	  computeL, vtx, vtx_lid, vtx_prid, vtx_bel_othSn);
  PrintInt10 ("A(:, v)", x_aind_end - x_aind_beg, &(sub[xsub[vtx_lid]]));
#endif

  ind = xsubPr[vtx_prid];
  if (vtx_bel_othSn == vtx)
    upd_lstSn = TRUE;

  while (ind != EMPTY || upd_lstSn) {
    if (upd_lstSn ) {
      upd_lstSn = FALSE;
      pr_lid = snrep_lid;
    }
    else {
      pr_lid = subPr[ind];
      ind = subPr[ind - 1];
    }
    
    if (!computeL)
      marker[vtx] = mark_vtx;
    if (pr_lid >= VInfo->nvtcs_loc) {
      compRcvd = TRUE;
      xsub_src = xsub_rcvd; sub_src = sub_rcvd;
      pr_lid -= VInfo->nvtcs_loc;
      k = xsub_src[pr_lid] + RCVD_IND;
    }
    else {
      compRcvd = FALSE;
      xsub_src = xsub; sub_src = sub;
      k = xsub_src[pr_lid];
    }

    PS->nops += xsub_src[pr_lid+1] - xsub_src[pr_lid];
    for (; k < xsub_src[pr_lid+1]; k++) {
      pr_elt = sub_src[k];
      if (pr_elt >= vtx && marker[pr_elt] != mark_vtx) {

	/* TEST available memory */
	if (next >= x_aind_end) {	
	  if (domain_symb) {
	    if (mem_error =
		psymbfact_LUXpandMem (iam, n, vtx, next, 0,
				      computeL, DOMAIN_SYMB, 1, 
				      Pslu_freeable, Llu_symbfact, VInfo, PS))
	      return (mem_error);
	  } else if (mem_error =
		     psymbfact_LUXpand (iam, n, EMPTY, vtx, &next, 0, 
					computeL, LL_SYMB, 1, 
					Pslu_freeable, Llu_symbfact, VInfo, PS))
	    return (mem_error);

	  x_aind_end = xsub[vtx_lid + 1];
	  if (computeL)   sub = Llu_symbfact->lsub; 
	  else   sub = Llu_symbfact->usub; 
	  if (!compRcvd) 
	    sub_src = sub;	  
	}

	sub[next] = pr_elt; next ++;

	if (pr_elt < lstVtx) neltsVtx_CSep ++;
	if (computeL && pr_elt == vtx)
	  diagind = TRUE;
	if (!computeL)
	  if (marker[pr_elt] == mark_vtx - 2)
	    if (pr_elt < prval_curvtx)
	      prval_curvtx = pr_elt;
	marker[pr_elt] = mark_vtx;	
      }
    }
  }

  /* Abort if the diagonal element is zero */
  if (computeL && diagind == FALSE) {
    printf("Pe[%d] At column " IFMT ", ", iam, vtx);
    ABORT("ParSymbFact() encounters zero diagonal");
  } 

  neltsVtx = next - xsub[vtx_lid];
  neltsZrVtx = 0; /* number of zero elements which would
		     be introduced in the vertex */
  neltsZrSn = 0; /* -"- in the supernode */
  neltsMatched = 0; 
  if (vtx != fstVtx) {
    for (k = xsub[snrep_lid]; k < xsub[snrep_lid+1]; k++) {
      vtx_elt = sub[k];
      if (vtx_elt >= vtx) {
	if ((vtx_elt > vtx && !computeL) || 
	    (vtx_elt >= vtx && computeL)) {
	  if (marker[vtx_elt] != mark_vtx)
	    neltsZrVtx ++;
	  else {
	    neltsMatched ++;
	  }
	}
	if (computeL && vtx_elt == vtx)
	  *p_vtx_bel_mySn = vtx;
	if (!computeL && vtx_elt == vtx + 1)
	  *p_vtx_bel_mySn = vtx + 1;
      }
    }
  }
  else {
    neltsMatched = neltsVtx;
    if (! computeL) 
      for (k = xsub[vtx_lid]; k < next; k++) {
	vtx_elt = sub[k];
	if (vtx_elt == vtx + 1)
	  *p_vtx_bel_mySn = vtx + 1;
      }
  }

  *p_neltsVtxInit  = neltsVtxInit;
  *p_neltsVtx      = neltsVtx;
  *p_neltsVtx_CSep = neltsVtx_CSep;
  *p_neltsZrVtx    = neltsZrVtx;
  *p_neltsMatched  = neltsMatched;
  *p_next          = next;
  *p_prval_curvtx  = prval_curvtx;
  return SUCCES_RET;
}

static int_t
updateRcvd_prGraph
(
 int_t n,         /* Input - order of the matrix */
 int   iam,       /* Input - my processor number */
 int_t *sub_rcvd,      /* elements of node */
 int_t sub_rcvd_sz,   /* Input - size of sub to be used in the update */
 int_t fstVtx_toUpd,  /* Input - first vertex to update */
 int_t lstVtx_toUpd,  /* Input - last vertex to update */
 int_t pr_offset,
 int   computeL,
 int_t *marker,
 Pslu_freeable_t *Pslu_freeable,
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo,   /* Input - local info on vertices distribution */
 psymbfact_stat_t *PS
 /*  marker: first elements of marker contain the nodes that will
     be used in the updates */
)
{
  int_t i, k, nelts, prVal, vtx_elt, vtx_elt_lid, ind;
  int_t vtx, vtx_lid, fstVtx_toUpd_lid, fstVtx_srcUpd_lid;
  int_t *xsub, *sub, *xsub_rcvd, *xsubPr, *subPr, szsubPr, *p_indsubPr;
  int_t maxNvtcsPProc, *globToLoc, mem_error;
  int_t nvtcs_toUpd, fstVtx_srcUpd, vtx_lid_p;
  
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  globToLoc     = Pslu_freeable->globToLoc;
  fstVtx_toUpd_lid = LOCAL_IND( globToLoc[fstVtx_toUpd] );
  nvtcs_toUpd = lstVtx_toUpd - fstVtx_toUpd;
  
  if (computeL) {
    xsub = Llu_symbfact->xlsub; sub = Llu_symbfact->lsub;
    xsub_rcvd = Llu_symbfact->xlsub_rcvd;
    xsubPr = Llu_symbfact->xlsubPr; subPr = Llu_symbfact->lsubPr;
    p_indsubPr = &(Llu_symbfact->indLsubPr);
    szsubPr = Llu_symbfact->szLsubPr;
  }
  else {
    xsub = Llu_symbfact->xusub; sub = Llu_symbfact->usub;
    xsub_rcvd = Llu_symbfact->xusub_rcvd;
    xsubPr = Llu_symbfact->xusubPr; subPr = Llu_symbfact->usubPr;
    p_indsubPr = &(Llu_symbfact->indUsubPr);
    szsubPr = Llu_symbfact->szUsubPr;
  }
  
  /* count number of elements in transpose representation of sub_rcvd */
  /* use marker to count those elements */
  for (i = 0; i < nvtcs_toUpd; i++)
    marker[i] = 0;
  for (i = 0; i <= VInfo->maxSzBlk; i++)
    xsub_rcvd[i] = 0;
  
  i = 0;
  fstVtx_srcUpd = EMPTY;
  while (i < sub_rcvd_sz) {
    vtx   = sub_rcvd[i + DIAG_IND];
    nelts = sub_rcvd[i + NELTS_IND];
    i += RCVD_IND;
    prVal = sub_rcvd[i];
    if (fstVtx_srcUpd == EMPTY) fstVtx_srcUpd = vtx;
    xsub_rcvd[vtx - fstVtx_srcUpd] = i - RCVD_IND;
    xsub_rcvd[vtx-fstVtx_srcUpd+1] = i + nelts;
    for (k = i; k < i + nelts; k++) {
      vtx_elt = sub_rcvd[k];
      if (vtx_elt > prVal)
	k = i + nelts;
      else {
	if (OWNER( globToLoc[vtx_elt] ) == iam) {
	  if (vtx_elt >= fstVtx_toUpd && vtx_elt < lstVtx_toUpd) {
	    vtx_elt_lid = LOCAL_IND( globToLoc[vtx_elt] ) - 
	      fstVtx_toUpd_lid;
	    marker[vtx_elt_lid] ++;
	  }
	}
      }
    }
    i += nelts;
  }

  vtx_lid = fstVtx_toUpd_lid - pr_offset;
  ind = 0;
  for (i = 0; i < nvtcs_toUpd; i++) {
    if (marker[i] != 0) {
      xsubPr[vtx_lid] = ind + 1;
      ind += 2* marker[i];
      marker[i] = xsubPr[vtx_lid] - 1;
    }
    vtx_lid ++;
  }
  
  if (ind == 0) 
    /* quick return if no update */
    return 0;

  /* test if enough memory in usubPr array */
  if (ind >= szsubPr) {
    if (mem_error = 
	psymbfact_prLUXpand (iam, ind, computeL, Llu_symbfact, PS))
      return (mem_error);
    if (computeL) 
      subPr = Llu_symbfact->lsubPr;  
    else 
      subPr = Llu_symbfact->usubPr;
  }
  *p_indsubPr = ind;
  
  i = 0;
  while (i < sub_rcvd_sz) {
    vtx   = sub_rcvd[i + DIAG_IND];
    nelts = sub_rcvd[i + NELTS_IND];
    i += RCVD_IND;
    prVal = sub_rcvd[i];
    for (k = i; k < i + nelts; k++) {
      vtx_elt = sub_rcvd[k];
      if (vtx_elt > prVal)
	k = i + nelts;
      else {
	if (OWNER( globToLoc[vtx_elt] ) == iam) {
	  if (vtx_elt >= fstVtx_toUpd && vtx_elt < lstVtx_toUpd) {
	    vtx_elt_lid = LOCAL_IND( globToLoc[vtx_elt] );
	    vtx_lid_p = vtx_elt_lid - pr_offset;
	    vtx_elt_lid -= fstVtx_toUpd_lid;
	    /* add vtx to structure of pruned graph */
	    if (marker[vtx_elt_lid] != xsubPr[vtx_lid_p] - 1) 
	      subPr[marker[vtx_elt_lid] - 2] = marker[vtx_elt_lid] + 1;
	    subPr[marker[vtx_elt_lid] + 1] = vtx - fstVtx_srcUpd + VInfo->nvtcs_loc;
	    subPr[marker[vtx_elt_lid]] = EMPTY;
	    marker[vtx_elt_lid] += 2;
	  }
	}
      }
    }
    i += nelts;
  }  
  
  for (i = fstVtx_toUpd; i < nvtcs_toUpd; i++)
    marker[i] = 0;
  return 0;
}

static int_t
update_prGraph 
(
 int   iam, 
 int_t n,           /* order of the matrix */
 int_t fstVtx_blk,  /* first vertex in block to factorize */
 int_t lstVtx_blk,  /* last vertex in block to factorize */
 int_t snrep_lid,   /* local index of current supernode reprezentative */
 int_t pr_offset,   /* offset in the indexing of prune structure */
 int_t prval_cursn, /* prune value of current supernode reprezentative */
 int_t xsub_snp1,   /* denotes xsub[snrep_lid + 1] */
 int   computeL,    /* Input - if 1, compute column L(:,vtx)
		               else compute row U(vtx, :) */
 Pslu_freeable_t *Pslu_freeable,
 Llu_symbfact_t *Llu_symbfact,   /* Input/Output - local L, U data structures */
 psymbfact_stat_t *PS
 )
{
  int_t k, mem_error;
  int_t kmin, kmax, ktemp, maxElt;
  int_t sn_elt, sn_elt_prid;
  int_t *globToLoc, maxNvtcsPProc;
  int_t *xsub, *sub, *xsubPr, *subPr;
  int_t *p_indsubPr, szsubPr;
  
  globToLoc     = Pslu_freeable->globToLoc;
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;

  if (computeL) {
    xsub = Llu_symbfact->xlsub; sub = Llu_symbfact->lsub;
    xsubPr = Llu_symbfact->xlsubPr; subPr = Llu_symbfact->lsubPr;
    p_indsubPr = &(Llu_symbfact->indLsubPr);
    szsubPr = Llu_symbfact->szLsubPr;
  }
  else {
    xsub = Llu_symbfact->xusub; sub = Llu_symbfact->usub;
    xsubPr = Llu_symbfact->xusubPr; subPr = Llu_symbfact->usubPr;
    p_indsubPr = &(Llu_symbfact->indUsubPr);
    szsubPr = Llu_symbfact->szUsubPr;
  }
  
  kmin = xsub[snrep_lid];
  kmax = xsub_snp1 - 1;
  if (prval_cursn != n)
    maxElt = prval_cursn;
  else
    maxElt = EMPTY;
  while (kmin <= kmax) {
    if (prval_cursn == n) {
      /* compute maximum element of L(:, vtx) */
      if (sub[kmin] > maxElt)
	maxElt = sub[kmin];
      kmin ++;
    }
    else {
      /* Do a quicksort-type partition. */    
      if (sub[kmax] > prval_cursn) 
	kmax--;
      else if (sub[kmin] <= prval_cursn)
	kmin++;
      else { /* kmin does'nt belong to G^s(L), and kmax belongs: 
	      * 	   interchange the two subscripts
	      */
	ktemp = sub[kmin];
	sub[kmin] = sub[kmax];
	sub[kmax] = ktemp;
	kmin ++;
	kmax --;
      }
    }
  }
  k = xsub[snrep_lid];
  while (sub[k] <= prval_cursn && k < xsub_snp1) {
    sn_elt = sub[k];
    if (sn_elt < lstVtx_blk) {
      sn_elt_prid = LOCAL_IND( globToLoc[sn_elt] ) - pr_offset;
      if ((*p_indsubPr) + 2 >= szsubPr) {
	if (mem_error = 
	    psymbfact_prLUXpand (iam, 0, computeL, Llu_symbfact, PS))
	  return (mem_error);
	if (computeL) {
	  subPr = Llu_symbfact->lsubPr;  szsubPr = Llu_symbfact->szLsubPr;
	}
	else {
	  subPr = Llu_symbfact->usubPr;  szsubPr = Llu_symbfact->szUsubPr;
	}
      }
      /* add krow to structure of pruned graph */
      subPr[(*p_indsubPr) + 1] = snrep_lid;
      subPr[(*p_indsubPr)] = xsubPr[sn_elt_prid];
      xsubPr[sn_elt_prid] = (*p_indsubPr) + 1;
      (*p_indsubPr) += 2;
    }
    if (sn_elt == maxElt) {
      /* move prune val in the first position */
      sub[k] = sub[xsub[snrep_lid]];
      sub[xsub[snrep_lid]] = sn_elt;
    }
    k ++; 
  }
  return SUCCES_RET;
}

static int_t
blk_symbfact
(SuperMatrix *A,
 int   iam,
 int   lvl, 
 int   szSep,
 int   ind_sizes1,
 int   ind_sizes2, 
 int_t *sizes,     /* Input - sizes of each node in the separator tree */
 int_t *fstVtxSep, /* Input - first vertex of each node in the tree */
 int_t fstVtx_loc, /* Input - first vertex local of the level */
 int_t fstVtx_blk,
 int_t lstVtx_blk,
 int_t *lsub_rcvd,      /* elements of node */
 int_t lsub_rcvd_sz,    /* size of sub to be explored */
 int_t *usub_rcvd,  
 int_t usub_rcvd_sz,
 Pslu_freeable_t *Pslu_freeable,   /* global LU data structures (modified) */
 Llu_symbfact_t *Llu_symbfact,   /* Input/Output - local L, U data structures */  
 vtcsInfo_symbfact_t *VInfo,  /* Input/Output - local info on vertices distribution */
 comm_symbfact_t *CS,
 psymbfact_stat_t *PS,
 int_t *marker,
 int_t *p_mark,    /* marker used to merge elements of vertices */
 int_t *p_nextl,   /* ptr to nextl in lsub structure */
 int_t *p_nextu,   /* ptr to nextu in usub structure */
 int_t *p_neltsZr, /* no of artificial zeros introduced so far */
 int_t *p_neltsTotal, /* no of nonzeros (including artificials) 
			 computed so far */
 int_t *p_nsuper_loc
 )
{
  int szSep_tmp, lvl_tmp, ii, jj;
  int_t  *xlsubPr, *xusubPr; 
  int_t  *xsup, *supno, *lsub, *xlsub, *usub, *xusub;
  int_t  vtx_lid, vtx_prid, vtx, vtx_super, vtx_elt, maxNvtcsPProc;
  int_t  ind, pr, pr_elt, newnext, k, vtx_elt_lid;
  int_t  nextl, nextu, nsuper_loc, nvtcs, n, mem_error;
  int_t  x_aind_beg, x_aind_end, i, szLp, xlsub_snp1, xusub_snp1;
  int_t  snrep, snrep_lid, szsn, vtxp1, *globToLoc, domain_symb;
  int_t lstVtx, neltsCurSep, maxNeltsVtx, fstVtx_loc_lid;
  /* supernode relaxation parameters */
  int_t  neltsVtx_L, neltsZrVtx_L, neltsMatched_L, neltsVtx_CSep_L;
  int_t  neltsVtx_U, neltsZrVtx_U, neltsMatched_U, neltsVtx_CSep_U;
  int_t  neltsZrSn_L, neltsZrSn_U, neltsZr, neltsTotal, 
    neltsZr_tmp, neltsTotal_tmp, neltsZrSn, neltsVtxInit_l, neltsVtxInit_u;
  /* next vertex belongs to current supernode pruned structure */
  int_t  vtx_bel_snL, vtx_bel_snU;
  /* marker variables */
  int_t  markl1_vtx, markl2_vtx, marku1_vtx, marku2_vtx;
  /* prune structure variables */
  int_t prval_cursn, prval_curvtx, pr_offset;
  /* variables for comms info */
  int_t neltSn_L, neltSn_U, lstVtx_tmp, stat;
  float relax_param, relax_seps;

  if (fstVtx_blk >= lstVtx_blk)
    return 0;
  
  /* Initializations */
  supno   = Pslu_freeable->supno_loc;
  lsub    = Llu_symbfact->lsub;   xlsub    = Llu_symbfact->xlsub;
  usub    = Llu_symbfact->usub;   xusub    = Llu_symbfact->xusub;
  xusubPr  = Llu_symbfact->xusubPr; 
  xlsubPr  = Llu_symbfact->xlsubPr;   
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  globToLoc     = Pslu_freeable->globToLoc;
  maxNeltsVtx   = VInfo->maxNeltsVtx;
  
  n          = A->ncol;
  nextl      = *p_nextl;
  nextu      = *p_nextu;
  neltsZr    = *p_neltsZr;
  neltsTotal = *p_neltsTotal;
  nsuper_loc = *p_nsuper_loc;
  marku2_vtx = *p_mark;
  lstVtx     = fstVtxSep[ind_sizes2] + sizes[ind_sizes2];

  snrep = fstVtx_blk; 
  snrep_lid = LOCAL_IND( globToLoc[fstVtx_blk] );
  szsn = 1;
  nvtcs = lstVtx_blk - fstVtx_blk;
  prval_cursn = n;
  vtx_bel_snL = EMPTY; vtx_bel_snU = EMPTY;
  
  /* set up to EMPTY xlsubPr[], xusubPr[] */
  if (PS->maxSzLPr < Llu_symbfact->indLsubPr)
    PS->maxSzLPr = Llu_symbfact->indLsubPr;
  if (PS->maxSzUPr < Llu_symbfact->indUsubPr)
    PS->maxSzUPr = Llu_symbfact->indUsubPr;
  for (i = 0; i < nvtcs; i++) {
    xlsubPr[i] = EMPTY;
    xusubPr[i] = EMPTY;
  }
  Llu_symbfact->indLsubPr = 0;
  Llu_symbfact->indUsubPr = 0;

  if (ind_sizes1 == 0) 
    domain_symb = TRUE;
  else {
    domain_symb = FALSE;
    fstVtx_loc_lid = LOCAL_IND( globToLoc[fstVtx_loc] );
  }
  
  vtx_prid = 0;
  vtx_lid = LOCAL_IND( globToLoc[fstVtx_blk] );
  pr_offset = vtx_lid;

  if (lsub_rcvd != NULL) {
    updateRcvd_prGraph (n, iam, lsub_rcvd, lsub_rcvd_sz,
			fstVtx_blk, lstVtx_blk, pr_offset, 1, marker,
			Pslu_freeable, Llu_symbfact, VInfo, PS);    
    updateRcvd_prGraph (n, iam, usub_rcvd, usub_rcvd_sz,
			fstVtx_blk, lstVtx_blk, pr_offset, 0, marker,
			Pslu_freeable, Llu_symbfact, VInfo, PS);
  }
  
  for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid ++, vtx_prid ++) {
    vtxp1 = vtx + 1;
    if (marku2_vtx +4 >= n) {
      /* reset to EMPTY marker array */
      for (i = 0; i < n; i++)
	marker[i] = EMPTY;
      marku2_vtx = EMPTY;
    }
    markl1_vtx = marku2_vtx + 1; markl2_vtx = markl1_vtx + 1;
    marku1_vtx = markl2_vtx + 1; marku2_vtx = marku1_vtx + 1;    

    prval_curvtx   = n;
    /* Compute nonzero structure L(:,vtx) */
    if (mem_error = 
	symbfact_vtx (n, iam, vtx, vtx_lid, vtx_prid, 1, domain_symb, 
		      fstVtx_blk,  lstVtx,
		      snrep_lid, szsn, &nextl,
		      marker, 
		      lsub_rcvd, lsub_rcvd_sz,
		      Pslu_freeable, Llu_symbfact, VInfo, PS, &neltsVtxInit_l,
		      &neltsVtx_L, &neltsVtx_CSep_L, &neltsZrVtx_L, 
		      &neltsMatched_L, markl1_vtx, &prval_curvtx, 
		      vtx_bel_snU, &vtx_bel_snL))
      return (mem_error);
    lsub = Llu_symbfact->lsub;

#ifdef TEST_SYMB
    PrintInt10 ("L(:, %d)", nextl - xlsub[vtx_lid], &(lsub[xlsub[vtx_lid]]));
#endif
    
    /* Compute nonzero structure of U(vtx,:) */
    if (mem_error = 
	symbfact_vtx (n, iam, vtx, vtx_lid, vtx_prid, 0, domain_symb, 
		      fstVtx_blk, lstVtx,
		      snrep_lid, szsn, &nextu,
		      marker, 
		      usub_rcvd, usub_rcvd_sz,
		      Pslu_freeable, Llu_symbfact, VInfo, PS, &neltsVtxInit_u,
		      &neltsVtx_U, &neltsVtx_CSep_U, &neltsZrVtx_U, 
		      &neltsMatched_U, marku1_vtx, &prval_curvtx,
		      vtx_bel_snL, &vtx_bel_snU))
      return (mem_error);
    usub = Llu_symbfact->usub;

#ifdef TEST_SYMB
    PrintInt10 ("U(%d, :)", nextu - xusub[vtx_lid], &(usub[xusub[vtx_lid]])); 
#endif
    
    /* update statistics on fill-in */
    if (!domain_symb) {
      stat = CEILING( (neltsVtxInit_l + neltsVtxInit_u), 2);
      if (Llu_symbfact->cntelt_vtcsA_lvl[vtx_lid - fstVtx_loc_lid] != stat) {
	stat = CEILING(stat, Llu_symbfact->cntelt_vtcsA_lvl[vtx_lid - fstVtx_loc_lid]);
	PS->fill_pelt[0] += (float) stat;
	if ((float) stat > PS->fill_pelt[1]) PS->fill_pelt[1] = (float) stat;
	PS->fill_pelt[2] += 1.;
      }
      stat = CEILING( (neltsVtx_L + neltsVtx_U), 2);
      stat = CEILING( stat, Llu_symbfact->cntelt_vtcsA_lvl[vtx_lid - fstVtx_loc_lid] );
      PS->fill_pelt[3] += (float) stat;
      if ((float) stat > PS->fill_pelt[4]) PS->fill_pelt[4] = (float) stat;
      PS->fill_pelt[5] += 1.;
    }     

    /* compute number of artificial zeros */
    neltsTotal = 0;
    neltsZr = 0;
    neltsZrSn_L    = neltsVtx_L - neltsMatched_L;
    neltsZrSn_U    = neltsVtx_U - neltsMatched_U;
    neltsZrSn      = neltsZrVtx_L + neltsZrVtx_U +
      (neltsZrSn_L + neltsZrSn_U) * szsn;
    neltsZr_tmp    = neltsZr + neltsZrSn;
    neltsTotal_tmp = neltsTotal + neltsZrSn + neltsVtx_L + neltsVtx_U;
    if (neltsTotal_tmp == 0)
      neltsTotal_tmp = 1;
    relax_param = (float) (neltsTotal_tmp - neltsZr_tmp) / neltsTotal_tmp;

#ifdef TEST_SYMB
    printf ("[%d] vtx %d pr %d szsn %d nVtx_L %d nZrSn_L %d nZrVtx_L %d\n",
	    iam, vtx, prval_curvtx, szsn,neltsVtx_L, neltsZrSn_L, neltsZrVtx_L);
    printf ("  [%d] nVtx_U %d, nZrSn_U %d nZrVtx_U %d nextl %d nextu %d\n",
	    iam, neltsVtx_U, neltsZrSn_U, neltsZrVtx_U, nextl, nextu);
    printf ("  [%d] nZr %d nZr_tmp %d nTot %d nTot_tmp %d rel %f test %d\n\n", 
	    iam, neltsZr, neltsZr_tmp, neltsTotal, neltsTotal_tmp,
	    relax_param, i);
#endif

    /* Check to see if vtx belongs in the same supernode as vtx-1 */
    supno[vtx_lid] = nsuper_loc;
    if (vtx == fstVtx_blk) {
      prval_cursn = prval_curvtx;
      neltsTotal += neltsVtx_L + neltsVtx_U;
    }
    else {
      if (maxNeltsVtx > 0) {
	relax_seps = (float) neltsVtx_L / (float) maxNeltsVtx;
	relax_seps *= (float) (neltsVtx_U+1) / (float) maxNeltsVtx;
      } 
      else
	relax_seps = 0.0;

      /* check if all upper separators are dense */
      if (relax_seps >= PS->relax_seps ) {
	VInfo->filledSep = FILLED_SEPS; 
	*p_nextl      = xlsub[vtx_lid];
	*p_nextu      = xusub[vtx_lid];
	nsuper_loc   += 1;
	*p_nsuper_loc = nsuper_loc;
	if (mem_error =
	    dnsUpSeps_symbfact (n, iam, szSep, ind_sizes1, ind_sizes2, 
				sizes, fstVtxSep, vtx,
				Llu_symbfact, Pslu_freeable, VInfo, CS, PS,
				p_nextl, p_nextu, p_nsuper_loc))
	  return (mem_error);
	/* set up neltsZr and neltsTotal */
	vtx = lstVtx_blk;
	return 0;
      } /* if all upper separators are dense */
      else {
	if (relax_param >= PS->relax_gen) {
	  /* vertex belongs to the same supernode */
	  if (prval_cursn > prval_curvtx || prval_cursn <= vtx)
	    prval_cursn = prval_curvtx;
	  neltsZr    = neltsZr_tmp;
	  neltsTotal = neltsTotal_tmp;
	  szsn ++;
	  /* add artificial zeros at the structure of current supernode */
	  newnext = xlsub[snrep_lid+1];
	  if (neltsZrSn_L != 0) {
	    for (k = xlsub[snrep_lid]; k < xlsub[snrep_lid+1]; k++) {
	      vtx_elt = lsub[k];
	      if (vtx_elt >= vtx) 
		marker[vtx_elt] = markl2_vtx;
	    }
	    for (k = xlsub[vtx_lid]; k < nextl; k++) {
	      vtx_elt = lsub[k];
	      if (marker[vtx_elt] != markl2_vtx) {
		/* add vtx_elt to the structure of snrep */
		lsub[newnext] = vtx_elt; newnext ++;
		marker[vtx_elt] = markl2_vtx;
	      }
	    }
	    xlsub[snrep_lid+1] = newnext;
	  }
	  xlsub[vtx_lid] = newnext;
	  nextl = newnext;
	  neltsVtx_L += neltsZrVtx_L;
	  
	  newnext = xusub[snrep_lid+1];
	  if (neltsZrSn_U != 0) {
	    for (k = xusub[snrep_lid]; k < xusub[snrep_lid+1]; k++) {
	      vtx_elt = usub[k];
	      if (vtx_elt >= vtx) {
		if (marker[vtx_elt] == markl2_vtx)
		  if (prval_cursn > vtx_elt && vtx_elt != vtx)
		    prval_cursn = vtx_elt;
		marker[vtx_elt] = marku2_vtx;
	      }
	    }
	    for (k = xusub[vtx_lid]; k < nextu; k++) {
	      vtx_elt = usub[k];
	      if (marker[vtx_elt] != marku2_vtx) {
		/* add vtx_elt to the structure of snrep */
		usub[newnext] = vtx_elt; newnext ++;
		if (marker[vtx_elt] == markl2_vtx)
		  if (prval_cursn > vtx_elt && vtx_elt != vtx)
		    prval_cursn = vtx_elt;
		marker[vtx_elt] = marku2_vtx;
	      }
	    }
	    if (marker[vtxp1] == marku2_vtx)
	      vtx_bel_snU = vtxp1;
	    xusub[snrep_lid+1] = newnext;
	  }
	  xusub[vtx_lid] = newnext;
	  nextu = newnext;
	  neltsVtx_U += neltsZrVtx_U;
	}  /* if ( relax_param >= PS->relax_param) */
      }  /* if (VInfo->filledSep != FILLED_SEPS) */
    } /* if (vtx != fstVtx_blk) */

    if ((relax_param < PS->relax_gen || vtx == lstVtx_blk-1) 
	&& VInfo->filledSep != FILLED_SEPS) {
      /* if a new supernode starts or is the last vertex */
      /* vtx starts a new supernode. Note we only store the
       * subscript set of the first column of a supernode.  */
      
      if (marker[vtxp1] == marku1_vtx)
	vtx_bel_snU = vtxp1;
      /* build the pruned structure */
      if (relax_param < PS->relax_gen
	  && vtx == lstVtx_blk - 1 && vtx != fstVtx_blk) 
	szLp = 2;
      else
	szLp = 1;
      if (vtx == fstVtx_blk) {
	xlsub_snp1 = nextl;
	xusub_snp1 = nextu;
      }
      else {
	xlsub_snp1 = xlsub[snrep_lid+1];
	xusub_snp1 = xusub[snrep_lid+1];	
      }
      while (szLp > 0) {
	szLp --;
#ifdef TEST_SYMB
	printf ("End sn %d szsn %d\n", nsuper_loc, szsn);
	printf ("BLD pr vtx %d snrep %d prval %d szLp %d\n",
		vtx, snrep, prval_cursn, szLp);
#endif
	
	update_prGraph (iam, n, fstVtx_blk, lstVtx_blk,
			snrep_lid, pr_offset, prval_cursn,
			xlsub_snp1, 1,
			Pslu_freeable, Llu_symbfact, PS);
	update_prGraph (iam, n, fstVtx_blk, lstVtx_blk,
			snrep_lid, pr_offset, prval_cursn,
			xusub_snp1, 0,
			Pslu_freeable, Llu_symbfact, PS);

#ifdef TEST_SYMB
	printf ("Adr lsub %p usub %p lsub %p pos %d usub %p pos %d\n", 
		&(lsub[xlsub[snrep_lid]]), &(usub[xusub[snrep_lid]]),
		lsub, xlsub[snrep_lid], usub, xusub[snrep_lid]);
	PrintInt10 ("Lsn", xlsub_snp1 - xlsub[snrep_lid],
		    &(lsub[xlsub[snrep_lid]]));
	PrintInt10 ("Usn", xusub_snp1 - xusub[snrep_lid],
		    &(usub[xusub[snrep_lid]]));
#endif

	if (prval_cursn >= lstVtx_blk) {
	  neltSn_L = xlsub_snp1 - xlsub[snrep_lid];
	  neltSn_U = xusub_snp1 - xusub[snrep_lid];
	  if (ind_sizes1 != 0) {
	    CS->snd_intraSz += neltSn_L + neltSn_U + 4;
	    CS->snd_LintraSz += neltSn_L + 2;
	  }
	  if (prval_cursn >= lstVtx) {
	    /* this supernode will be send to next layers of the tree */
	    lvl_tmp = lvl;
	    ii = ind_sizes1;
	    jj = ind_sizes2;
	    szSep_tmp = szSep;
	    lstVtx_tmp = lstVtx;
	    while (prval_cursn >= lstVtx_tmp && szSep_tmp != 1) {
	      jj = ii + szSep_tmp + (jj - ii) / 2;
	      ii += szSep_tmp;
	      lvl_tmp ++;
	      szSep_tmp = szSep_tmp / 2;
	      lstVtx_tmp = fstVtxSep[jj] + sizes[jj];
	      CS->snd_interSz[lvl_tmp] += neltSn_L + neltSn_U + 4;
	      CS->snd_LinterSz[lvl_tmp] += neltSn_L + 2;
	      if (CS->snd_vtxinter[lvl_tmp] == EMPTY)
		CS->snd_vtxinter[lvl_tmp] = snrep;
	    }
	  }
	}
	snrep = vtx;
	snrep_lid = vtx_lid;
	prval_cursn = prval_curvtx;
	szsn        = 1;
	xlsub_snp1  = nextl;
	xusub_snp1  = nextu;
      }
      if (relax_param < PS->relax_gen) {
	neltsTotal += neltsVtx_L + neltsVtx_U;
	nsuper_loc ++;	
	supno[vtx_lid] = nsuper_loc;
	if (marker[vtxp1] == marku1_vtx)
	  vtx_bel_snU = vtxp1;
	else
	  vtx_bel_snU = EMPTY;
      }
    }
    if (vtx == lstVtx_blk - 1)
      nsuper_loc ++;
    
    /* check if current separator is dense */
    if (!VInfo->filledSep) {
      relax_seps = (float) neltsVtx_CSep_L / (float) (lstVtx - vtx);
      relax_seps *= (float) (neltsVtx_CSep_U+1) / (float) (lstVtx - vtx);
      if (relax_seps >= PS->relax_curSep ) 
	VInfo->filledSep = FILLED_SEP;
    }
    maxNeltsVtx --;
  }
    
  *p_mark       = marku2_vtx + 1;
  *p_nextl      = nextl;
  *p_nextu      = nextu;
  *p_neltsZr    = neltsZr;
  *p_neltsTotal = neltsTotal;
  *p_nsuper_loc = nsuper_loc;

  return 0;
}

static void
domain_symbfact
(SuperMatrix *A,
 int   iam,        /* Input - my processor number */  
 int   lvl,        /* Input - current level in the separator tree */
 int   szSep,      /* Input - size of the current separator (node) */
 int   ind_sizes1,
 int   ind_sizes2, 
 int_t *sizes,     /* Input - sizes of each node in the separator tree */
 int_t *fstVtxSep, /* Input - first vertex of each node in the tree */
 int_t fstVtx,     /* Input - first vertex of current node */ 
 int_t lstVtx,     /* Input - last vertex of current node */ 
 Pslu_freeable_t *Pslu_freeable,   /* global LU data structures (modified) */
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo,  /* Input/Output - local info on vertices distribution */
 comm_symbfact_t *CS,
 psymbfact_stat_t *PS,
 int_t *marker,
 int_t *p_mark,    /* marker used to merge elements of vertices */
 int_t *p_nextl,   /* ptr to nextl in lsub structure */
 int_t *p_nextu,   /* ptr to nextu in usub structure */
 int_t *p_neltsZr, /* no of artificial zeros introduced so far */
 int_t *p_neltsTotal, /* no of nonzeros (including artificials) 
			 computed so far */
 int_t *p_nsuper_loc
 )
{
  int_t lstVtx_lid, maxNvtcsPProc; 

  /* call blk_symbfact */
  blk_symbfact (A, iam, lvl, 
		szSep, ind_sizes1, ind_sizes2, sizes, fstVtxSep,
		EMPTY, fstVtx, lstVtx, 
		NULL, EMPTY, NULL, EMPTY,
		Pslu_freeable, Llu_symbfact, VInfo, CS, PS,
		marker, p_mark,
		p_nextl, p_nextu, p_neltsZr, p_neltsTotal, 
		p_nsuper_loc);

  if (VInfo->filledSep != FILLED_SEPS) {
    maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
    if (fstVtx >= lstVtx)
      lstVtx_lid = 0;
    else 
      lstVtx_lid = LOCAL_IND( Pslu_freeable->globToLoc[lstVtx-1] ) + 1;
    VInfo->xlsub_nextLvl  = Llu_symbfact->xlsub[lstVtx_lid];
    Llu_symbfact->xlsub[lstVtx_lid] = *p_nextl;
    VInfo->xusub_nextLvl  = Llu_symbfact->xusub[lstVtx_lid];
    Llu_symbfact->xusub[lstVtx_lid] = *p_nextu;
  }
  VInfo->maxNeltsVtx -= lstVtx - fstVtx;
}


/*! \brief
 *
 * <pre>
 * Compute counts of rows/columns of current separator.
 * cntelt_vtcs[i] is 0 when i is nonzero before current separator
 * and n when i is zero before current separator.
 *
 * Set up nvtcsLvl_loc.
 * </pre>
 */
static void
initLvl_symbfact
(
 int_t n,       /* Input - order of the matrix */
 int   iam,     /* Input - my processor number */
 int_t fstVtx,  /* Input - first vertex of current node */   
 int_t lstVtx,  /* Input - last vertex of current node */   
 Pslu_freeable_t *Pslu_freeable,
 Llu_symbfact_t *Llu_symbfact, /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo, /* Input/Output - local info on vertices distribution */
 psymbfact_stat_t *PS,
 MPI_Comm ndComm,
 int_t  *marker,
 int_t  nextl,
 int_t  nextu
 ) 
{
  int_t *cntelt_vtcs, x_aind_beg, x_aind_end, x_aind_beg_l, x_aind_beg_u,
    nelts_asup, nelts_ainf;
  int_t nvtcsLvl_loc, fstVtx_loc, fstVtx_loc_lid, fstVtx_nextLvl;
  int_t curblk_loc, nblks_loc, ind_blk;
  int_t *lsub, *xlsub, *usub, *xusub;
  int_t *begEndBlks_loc, code_err, mem_error;
  int_t i, j, k, vtx, vtx_lid, fstVtx_blk, lstVtx_blk, vtx_elt, p, fill;
  int_t nelts, nelts_fill_l, nelts_fill_u, nelts_cnts, maxNvtcsPProc, *globToLoc;
  int_t use_fillcnts, cntelt_vtx_l, cntelt_vtx_u;
  MPI_Status status;
  
  fill = PS->fill_par;
  VInfo->filledSep = FALSE;
  
  /* Initializations */
  maxNvtcsPProc  = Pslu_freeable->maxNvtcsPProc;
  globToLoc      = Pslu_freeable->globToLoc;
  curblk_loc     = VInfo->curblk_loc;
  nblks_loc      = VInfo->nblks_loc;
  begEndBlks_loc = VInfo->begEndBlks_loc;
  cntelt_vtcs    = Llu_symbfact->cntelt_vtcs;
  lsub    = Llu_symbfact->lsub;   xlsub    = Llu_symbfact->xlsub;
  usub    = Llu_symbfact->usub;   xusub    = Llu_symbfact->xusub;
  
  /* compute nvtcsLvl_loc */
  nvtcsLvl_loc = 0;  
  ind_blk = curblk_loc;
  while (fstVtx > begEndBlks_loc[ind_blk] && ind_blk < 2 * nblks_loc) {
    ind_blk += 2;
  }
  curblk_loc = ind_blk;
  fstVtx_loc = begEndBlks_loc[ind_blk];
  while (begEndBlks_loc[ind_blk] < lstVtx && ind_blk < 2 * nblks_loc) {
    nvtcsLvl_loc += begEndBlks_loc[ind_blk + 1] - 
      begEndBlks_loc[ind_blk];
    ind_blk += 2;
  }
  fstVtx_nextLvl = begEndBlks_loc[ind_blk];
  VInfo->nvtcsLvl_loc = nvtcsLvl_loc;
  VInfo->curblk_loc = curblk_loc;
  
  fstVtx_loc_lid = LOCAL_IND( globToLoc[fstVtx_loc] );
  vtx_lid      = fstVtx_loc_lid;
  x_aind_beg_l = VInfo->xlsub_nextLvl;
  x_aind_beg_u = VInfo->xusub_nextLvl;
  nelts_cnts   = 0;
  nelts_fill_l = 0;
  nelts_fill_u = 0;
  ind_blk      = curblk_loc;
  
  while (begEndBlks_loc[ind_blk] < lstVtx && ind_blk < 2 * nblks_loc) {
    fstVtx_blk = begEndBlks_loc[ind_blk];
    lstVtx_blk = begEndBlks_loc[ind_blk + 1];
    ind_blk += 2;
    for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid ++) 
      nelts_cnts += cntelt_vtcs[vtx_lid];      
    nelts_fill_l += fill * (xlsub[vtx_lid] - x_aind_beg_l);
    nelts_fill_u += fill * (xusub[vtx_lid] - x_aind_beg_u);
    x_aind_beg_l = xlsub[vtx_lid];
    x_aind_beg_u = xusub[vtx_lid];
  }

  if (nvtcsLvl_loc != 0) {
    nelts_ainf = xlsub[vtx_lid] - VInfo->xlsub_nextLvl;
    nelts_asup = xusub[vtx_lid] - VInfo->xusub_nextLvl;
  }
  else {
    nelts_ainf = 0;
    nelts_asup = 0;
  }
  
  use_fillcnts = FALSE;
  if (nextl + nelts_cnts >= Llu_symbfact->szLsub - nelts_ainf ||
      nextu + nelts_cnts >= Llu_symbfact->szUsub - nelts_asup) { 
    use_fillcnts = TRUE;
  }
 
  use_fillcnts = TRUE; 
  
  if (use_fillcnts) {
    if (nextl + nelts_fill_l >= Llu_symbfact->szLsub - nelts_ainf)
      mem_error = 
	psymbfact_LUXpandMem (iam, n, fstVtx, nextl,
			      nextl + nelts_fill_l, LSUB,
			      RL_SYMB, 1, 
			      Pslu_freeable, Llu_symbfact, VInfo, PS);
    lsub = Llu_symbfact->lsub;
    if (nextu + nelts_fill_u >= Llu_symbfact->szUsub - nelts_asup) 
      mem_error = 
	psymbfact_LUXpandMem (iam, n, fstVtx, nextu,
			      nextu + nelts_fill_u, USUB,
			      RL_SYMB, 1, 
			      Pslu_freeable, Llu_symbfact, VInfo, PS);      
    usub = Llu_symbfact->usub;
  }

  /* init xlsub[fstVtx:lstVtx] and xusub[fstVtx:lstVtx] and
     copy elements of A[fstVtx:lstVtx, fstVtx:lstVtx] in lsub and usub */
  fstVtx_loc_lid = LOCAL_IND( globToLoc[fstVtx_loc] );
  x_aind_beg_l = VInfo->xlsub_nextLvl;
  x_aind_beg_u = VInfo->xusub_nextLvl;
  vtx_lid = fstVtx_loc_lid;
  ind_blk = curblk_loc;

  while (begEndBlks_loc[ind_blk] < lstVtx && ind_blk < 2 * nblks_loc) {
    fstVtx_blk = begEndBlks_loc[ind_blk];
    lstVtx_blk = begEndBlks_loc[ind_blk + 1];
    ind_blk += 2;

    for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid ++) {
      if (vtx_lid != fstVtx_loc_lid) {
	x_aind_beg_l = xlsub[vtx_lid];
	x_aind_beg_u = xusub[vtx_lid];
      }
      if (use_fillcnts) {
	cntelt_vtx_l = fill * (xlsub[vtx_lid+1] - x_aind_beg_l);
	cntelt_vtx_u = fill * (xusub[vtx_lid+1] - x_aind_beg_u);
      }
      else {
	cntelt_vtx_l = cntelt_vtcs[vtx_lid];
	cntelt_vtx_u = cntelt_vtcs[vtx_lid];
      }
      x_aind_end = xlsub[vtx_lid + 1];
      Llu_symbfact->cntelt_vtcsA_lvl[vtx_lid - fstVtx_loc_lid] = 
	CEILING( (xlsub[vtx_lid+1]-x_aind_beg_l + xusub[vtx_lid+1]-x_aind_beg_u), 2);
      
      xlsub[vtx_lid] = nextl;
      nelts = 0;
      for (k = x_aind_beg_l; k < x_aind_end; k++) {
	lsub[nextl] = lsub[k]; nextl ++;
	nelts ++;
      }
      if (nelts < cntelt_vtx_l) 
	lsub[nextl] = EMPTY; 
      nextl += cntelt_vtx_l - nelts;
      x_aind_end = xusub[vtx_lid + 1];
      xusub[vtx_lid] = nextu;
      nelts = 0;
      for (k = x_aind_beg_u; k < x_aind_end; k++) {
	usub[nextu] = usub[k]; nextu ++;
	nelts ++;
      }
      if (nelts < cntelt_vtx_u) 
	usub[nextu] = EMPTY; 
      nextu += cntelt_vtx_u - nelts;
    }
  }
 
  if (nvtcsLvl_loc == 0) {
    if (curblk_loc == 0)
      vtx_lid = 0;
    else {
      if (begEndBlks_loc[curblk_loc-1] == 0)
	vtx_lid = 0;
      else
	vtx_lid = LOCAL_IND( globToLoc[begEndBlks_loc[curblk_loc-1] - 1] ) + 1;
    }

    xlsub[vtx_lid] = nextl;
    xusub[vtx_lid] = nextu;
  }
  else {
    VInfo->xlsub_nextLvl   = xlsub[vtx_lid];
    xlsub[vtx_lid] = nextl;
    VInfo->xusub_nextLvl   = xusub[vtx_lid];
    xusub[vtx_lid] = nextu;
    if (PS->estimLSz < nextl)
      PS->estimLSz = nextl;
    if (PS->estimUSz < nextu)
      PS->estimUSz = nextu;
    
    VInfo->nnz_ainf_loc -= nelts_ainf;
    VInfo->nnz_asup_loc -= nelts_asup;
  }
  VInfo->fstVtx_nextLvl = fstVtx_nextLvl;
}


static int_t
expand_RL 
(
 int_t computeRcvd, /* if = 1, then update from receive buffer,
		       else update from own data */
 int_t n,
 int   iam,       /* process number */
 int_t *lsub_rcvd,      /* elements of node */
 int_t lsub_rcvd_sz,    /* size of sub to be explored */
 int_t *usub_rcvd,  
 int_t usub_rcvd_sz,
 int_t vtxXp,
 int_t vtx_upd_pr,    /* ind in pruned structure of upd vertex which 
			 doesn't fit into the alloc memory */
 int_t lstVtx_upd_pr, /* ind in pruned structure of lst vtx to update */
 int_t fstVtx_srcUpd, /* first vertex source of the updates */
 int_t lstVtx_srcUpd, /* last vertex source of the updates */
 int_t fstVtx_toUpd,  /* first vertex to update */
 int_t lstVtx_toUpd,  /* last vertex to update */
 int_t nvtcs_toUpd,   /* no of vertices to update */
 int   computeL,
 int_t *pmarkl,
 int_t *marker,
 Pslu_freeable_t *Pslu_freeable,
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo, /* Input/Output - local info on vertices distribution */
 psymbfact_stat_t *PS
 )
{
  int_t fstVtx_toUpd_lid, vtx_lid, vtx, vtx_elt, vtx_elt_lid, nextl, nelts_in;
  int_t i, ii, j, nelts, nelts_vtx, mpnelts, lvtx_lid, elt, vtxXp_lid;
  int_t *xusubPr, *usubPr, *xlsub, *lsub, *xusub, *usub;
  int_t markl, *globToLoc, maxNvtcsPProc;
  int_t mem_error, len_texp;
  
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  globToLoc     = Pslu_freeable->globToLoc;

  xusubPr = Llu_symbfact->xlsubPr; usubPr  = Llu_symbfact->lsubPr;
  if (computeL) {
    xlsub   = Llu_symbfact->xlsub;   lsub    = Llu_symbfact->lsub;
    xusub   = Llu_symbfact->xusub;   usub    = Llu_symbfact->usub;
  }
  else {
    xlsub   = Llu_symbfact->xusub;   lsub    = Llu_symbfact->usub;
    xusub   = Llu_symbfact->xlsub;   usub    = Llu_symbfact->lsub;
  }
  markl = *pmarkl + 1;
  fstVtx_toUpd_lid = LOCAL_IND( globToLoc[fstVtx_toUpd] );
  vtxXp_lid = LOCAL_IND( globToLoc[vtxXp] );
  nextl = xlsub[vtxXp_lid+1];
    
  lvtx_lid = EMPTY;
  if (lstVtx_srcUpd != EMPTY)
    lvtx_lid = LOCAL_IND( globToLoc[lstVtx_srcUpd - 1] );

  /* count the number of new elements, and update Llu_symbfact->cntelt_vtcs */
  vtx_lid = fstVtx_toUpd_lid;
  vtx_lid += vtx_upd_pr;
  len_texp = 0;
  for (i = vtx_upd_pr; i < lstVtx_upd_pr; i++, vtx_lid ++) { 
    nelts_vtx = xlsub[vtx_lid+1] - xlsub[vtx_lid];
    if (xusubPr[i] != xusubPr[i+1]) {
      j = xusubPr[i]; 
      vtx = usubPr[j];
      /* setup marker structure for already existing elements */
      ii = xlsub[vtx_lid];
      while (lsub[ii] != EMPTY && ii < xlsub[vtx_lid + 1]) {
	marker[lsub[ii]] = markl;
	ii ++;
      }
      nelts_vtx = ii - xlsub[vtx_lid];
      for (j = xusubPr[i] + 1; j < xusubPr[i+1]; j++) {
	vtx_elt = usubPr[j];
	ii = marker[vtx_elt];
	if (computeRcvd) {
	  nelts = lsub_rcvd[ii + NELTS_IND];
	  ii += RCVD_IND;
	  mpnelts = marker[vtx_elt] + nelts + RCVD_IND;
	}
	else {
	  vtx_elt_lid = LOCAL_IND( globToLoc[vtx_elt] );
	  if (vtx_elt_lid == lvtx_lid)
	    nelts = lsub_rcvd_sz - ii;
	  else
	    nelts = xlsub[vtx_elt_lid+1] - xlsub[vtx_elt_lid];
	  mpnelts = marker[vtx_elt] + nelts;
	}
	
	if (!computeL)
	  marker[vtx] = markl;
	for (ii; ii < mpnelts; ii++) {
	  elt = lsub_rcvd[ii];
	  if (elt >= vtx) {
	    if (marker[elt] != markl) {
	      /* add elt to structure of vtx */
	      marker[elt] = markl;
	      nelts_vtx ++;
	    }
	  }
	}
      }
      if (nelts_vtx != 0 && (nelts_vtx > xlsub[vtx_lid+1] - xlsub[vtx_lid])) {
	nelts_in = xlsub[vtx_lid+1] - xlsub[vtx_lid];
	if (nelts_in == 0) nelts_in = 1;
	j = nelts_vtx / nelts_in;
	if (nelts_vtx % nelts_in != 0) j++;
	nelts_vtx = j * nelts_in;
      }
      else
	nelts_vtx = xlsub[vtx_lid+1] - xlsub[vtx_lid];
      markl ++;
      if (markl == n) {
	/* reset marker array */
	for (j = fstVtx_toUpd; j < n; j++)
	  marker[j] = EMPTY;
	markl = 0;
      }
    }
    Llu_symbfact->cntelt_vtcs[vtx_lid] = nelts_vtx;
    len_texp += nelts_vtx;
  }
  for (; i < nvtcs_toUpd; i++, vtx_lid++) {
    nelts_vtx = xlsub[vtx_lid+1] - xlsub[vtx_lid];    
    Llu_symbfact->cntelt_vtcs[vtx_lid] = nelts_vtx;
    len_texp += nelts_vtx;
  }

  *pmarkl = markl;
  /* mark elements array */
  for (i = xlsub[vtxXp_lid]; i < nextl; i++) {
    marker[lsub[i]] = markl;
  }

  nextl = xlsub[vtxXp_lid+1];  
  if (mem_error = 
      psymbfact_LUXpand_RL (iam, n, vtxXp, nextl, len_texp, 
			    computeL, Pslu_freeable, Llu_symbfact, VInfo, PS))
    return (mem_error);		

  return 0;
}


static int_t
rl_update
(
 int   computeRcvd, /* if = 1, then update from receive buffer,
		       else update from own data */
 int_t n,
 int   iam,       /* process number */
 int_t *lsub_rcvd,      /* elements of node */
 int_t lsub_rcvd_sz,    /* size of sub to be explored */
 int_t *usub_rcvd,  
 int_t usub_rcvd_sz,
 int_t fstVtx_srcUpd, /* first vertex source of the updates */
 int_t lstVtx_srcUpd, /* last vertex source of the updates */
 int_t indBlk_srcUpd, /* block index of first vertex */
 int_t fstVtx_toUpd,  /* first vertex to update */
 int_t lstVtx_toUpd,  /* last vertex to update */
 int_t nvtcs_toUpd,   /* no of vertices to update */
 int   computeL,
 int_t *pmarkl,
 int_t *marker,
 Pslu_freeable_t *Pslu_freeable,
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo,  /* Input/Output - local info on vertices distribution */
 psymbfact_stat_t *PS
 /*  marker: first elements of marker contain the nodes that will
     be used in the updates */
 )
{
  int_t i, j, k, prVal, nelts, ind, nextl, ii, mpnelts, mem_error;
  int_t vtx, vtx_lid, vtx_elt, vtx_elt_lid, lvtx_lid;
  int_t fstVtx_toUpd_lid, markl, elt, vtx_loc, ind_blk;
  int_t *xusubPr, *usubPr, *xlsub, *lsub, *xusub, *usub;
  int_t fstVtx_upd, lstVtx_upd, maxNvtcsPProc, *globToLoc;
  int_t fstVtx_srcUpd_lid, nelts_vtx, expand;
  
  /* quick return */
  if (fstVtx_toUpd >= lstVtx_toUpd)
    return 0;

  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  globToLoc     = Pslu_freeable->globToLoc;
  
  fstVtx_upd = EMPTY;
  lstVtx_upd = EMPTY;
  xusubPr = Llu_symbfact->xlsubPr; usubPr  = Llu_symbfact->lsubPr;
  if (computeL) {
    xlsub   = Llu_symbfact->xlsub;   lsub    = Llu_symbfact->lsub;
    xusub   = Llu_symbfact->xusub;   usub    = Llu_symbfact->usub;
  }
  else {
    xlsub   = Llu_symbfact->xusub;   lsub    = Llu_symbfact->usub;
    xusub   = Llu_symbfact->xlsub;   usub    = Llu_symbfact->lsub;
  }
  markl = *pmarkl;
  fstVtx_toUpd_lid = LOCAL_IND( globToLoc[fstVtx_toUpd] );

  /* count number of elements in transpose representation of usub_rcvd */
  /* use marker to count those elements */
  for (i = 0; i < nvtcs_toUpd; i++)
    marker[i] = 0;
  
  i = 0;
  if (fstVtx_srcUpd != EMPTY) {
    fstVtx_srcUpd_lid = LOCAL_IND( globToLoc[fstVtx_srcUpd] );
    vtx_lid = fstVtx_srcUpd_lid;
  }
  lvtx_lid = EMPTY;
  if (lstVtx_srcUpd != EMPTY)
    lvtx_lid = LOCAL_IND( globToLoc[lstVtx_srcUpd - 1] );
  
  while (i < usub_rcvd_sz) {
    if (computeRcvd) {
      vtx   = usub_rcvd[i + DIAG_IND];
      nelts = usub_rcvd[i + NELTS_IND];
      i += RCVD_IND;
    }
    else {
      if (vtx_lid == lvtx_lid)
	nelts = usub_rcvd_sz - i;
      else
	nelts = xusub[vtx_lid + 1] - xusub[vtx_lid];
      vtx_lid ++;
    }
    prVal = usub_rcvd[i];
    for (k = i; k < i + nelts; k++) {
      vtx_elt = usub_rcvd[k];
      if (vtx_elt > prVal)
	k = i + nelts;
      else {
	if (OWNER( globToLoc[vtx_elt] ) == iam) {
	  if (vtx_elt >= fstVtx_toUpd && vtx_elt < lstVtx_toUpd) {
	    vtx_elt_lid = LOCAL_IND( globToLoc[vtx_elt] ) - 
	      fstVtx_toUpd_lid;
	    marker[vtx_elt_lid] ++;
	  }
	}
      }
    }
    i += nelts;
  }

  ind = 0;
  for (i = 0; i < nvtcs_toUpd; i++) {
    if (marker[i] != 0) {
      marker[i] ++;
      if (fstVtx_upd == EMPTY)
	fstVtx_upd = i;
      lstVtx_upd = i;
    }
    xusubPr[i] = ind;
    ind += marker[i];
    marker[i] = xusubPr[i];
  }
  xusubPr[i] = ind;
  lstVtx_upd ++;

  if (ind == 0) 
    /* quick return if no update */
    return 0;

  /* test if enough memory in usubPr array */
  if (ind > Llu_symbfact->szLsubPr) {
    if (mem_error = 
	psymbfact_prLUXpand (iam, ind, LSUB_PR, Llu_symbfact, PS))
      return (mem_error);
    usubPr  = Llu_symbfact->lsubPr;
  }
  
  i = 0;
  if (fstVtx_srcUpd != EMPTY) {
    vtx_loc = fstVtx_srcUpd;
    vtx_lid = LOCAL_IND( globToLoc[vtx_loc] );
    ind_blk = indBlk_srcUpd;
  }
  while (i < usub_rcvd_sz) {
    if (computeRcvd) {
      vtx   = usub_rcvd[i + DIAG_IND];
      nelts = usub_rcvd[i + NELTS_IND];
      i += RCVD_IND;
    }
    else {
      vtx = vtx_loc;
      if (vtx_lid == lvtx_lid)
	nelts = usub_rcvd_sz - i;
      else
	nelts = xusub[vtx_lid + 1] - xusub[vtx_lid];
      vtx_lid ++;
      vtx_loc ++;
      if (ind_blk != EMPTY)
	if (vtx_loc == VInfo->begEndBlks_loc[ind_blk+1]) {
	  ind_blk += 2;
	  vtx_loc = VInfo->begEndBlks_loc[ind_blk];
	}
    }

    prVal = usub_rcvd[i];
    for (k = i; k < i + nelts; k++) {
      vtx_elt = usub_rcvd[k];
      if (vtx_elt > prVal)
	k = i + nelts;
      else {
	if (OWNER( globToLoc[vtx_elt]) == iam) {
	  if (vtx_elt >= fstVtx_toUpd && vtx_elt < lstVtx_toUpd) {
	    vtx_elt_lid = LOCAL_IND( globToLoc[vtx_elt] ) - fstVtx_toUpd_lid;
	    /* add vtx_elt to the pruned structure */
	    if (marker[vtx_elt_lid] == xusubPr[vtx_elt_lid]) {
	      usubPr[marker[vtx_elt_lid]] = vtx_elt;
	      marker[vtx_elt_lid] ++;
	    }
	    usubPr[marker[vtx_elt_lid]] = vtx;
	    marker[vtx_elt_lid] ++;
	  }
	}
      }
    }
    i += nelts;
  }
  /* reset marker array */
  for (i = 0; i < nvtcs_toUpd; i++)
    marker[i] = EMPTY;
  if (fstVtx_srcUpd != EMPTY) {
    vtx_loc = fstVtx_srcUpd;
    vtx_lid = LOCAL_IND( globToLoc[vtx_loc] );
    ind_blk = indBlk_srcUpd;
  }
  i = 0;
  while (i < lsub_rcvd_sz) {
    if (computeRcvd) {
      vtx   = lsub_rcvd[i + DIAG_IND];
      nelts = lsub_rcvd[i + NELTS_IND];
      marker[vtx] = i;
      i += RCVD_IND;
    }
    else {
      vtx = vtx_loc;
      if (vtx_lid == lvtx_lid)
	nelts = lsub_rcvd_sz - i;
      else
	nelts = xlsub[vtx_lid + 1] - xlsub[vtx_lid];
      vtx_lid ++;
      marker[vtx] = i;
      vtx_loc ++;
      if (ind_blk != EMPTY)
	if (vtx_loc == VInfo->begEndBlks_loc[ind_blk+1]) {
	  ind_blk += 2;
	  vtx_loc = VInfo->begEndBlks_loc[ind_blk];
	}
    }
    i += nelts;
  }

  /* use the pruned structure to update symbolic factorization */
  vtx_lid = fstVtx_toUpd_lid;
  vtx_lid += fstVtx_upd;
  for (i = fstVtx_upd; i < lstVtx_upd; i++, vtx_lid ++) { 
    if (xusubPr[i] != xusubPr[i+1]) {
      j = xusubPr[i]; 
      vtx = usubPr[j];
      /* setup marker structure for already existing elements */
      ii = xlsub[vtx_lid];
      while (lsub[ii] != EMPTY && ii < xlsub[vtx_lid + 1]) {
	marker[lsub[ii]] = markl;
	ii ++;
      }
      PS->nops += ii - xlsub[vtx_lid];
      nextl = ii;
      for (j = xusubPr[i] + 1; j < xusubPr[i+1]; j++) {
	vtx_elt = usubPr[j];
	ii = marker[vtx_elt];
	if (computeRcvd) {
	  nelts = lsub_rcvd[ii + NELTS_IND];
	  ii += RCVD_IND;
	  mpnelts = marker[vtx_elt] + nelts + RCVD_IND;
	}
	else {
	  vtx_elt_lid = LOCAL_IND( globToLoc[vtx_elt] );
	  if (vtx_elt_lid == lvtx_lid)
	    nelts = lsub_rcvd_sz - ii;
	  else
	    nelts = xlsub[vtx_elt_lid+1] - xlsub[vtx_elt_lid];
	  mpnelts = marker[vtx_elt] + nelts;
	}
		
	if (!computeL)
	  marker[vtx] = markl;
	PS->nops += mpnelts - ii;
	for (ii; ii < mpnelts; ii++) {
	  elt = lsub_rcvd[ii];
	  if (elt >= vtx) {
	    if (marker[elt] != markl) {
	      /* add elt to structure of vtx */
	      if (nextl >= xlsub[vtx_lid + 1]) {
		if (mem_error = 
		    expand_RL (computeRcvd, n, iam, lsub_rcvd, lsub_rcvd_sz,
			       usub_rcvd, usub_rcvd_sz, vtx, i,
			       lstVtx_upd, fstVtx_srcUpd, lstVtx_srcUpd,
			       fstVtx_toUpd, lstVtx_toUpd, nvtcs_toUpd, computeL,
			       &markl, marker, Pslu_freeable, Llu_symbfact, VInfo, PS))
		    return (mem_error);
		if (computeL) {
		  lsub    = Llu_symbfact->lsub;
		  if (!computeRcvd) 
		    lsub_rcvd    = 
		      &(Llu_symbfact->lsub[Llu_symbfact->xlsub[fstVtx_srcUpd_lid]]);
		} else {
		  marker[vtx] = markl;
		  lsub    = Llu_symbfact->usub;
		  if (!computeRcvd) 
		    lsub_rcvd = 
		      &(Llu_symbfact->usub[Llu_symbfact->xusub[fstVtx_srcUpd_lid]]);
		}
	      }
	      lsub[nextl] = elt; nextl ++;
	      marker[elt] = markl;
	    }
	  }
	}
      }
      if (nextl < xlsub[vtx_lid+1])
	lsub[nextl] = EMPTY;
      markl ++;
      if (markl == n) {
	/* reset marker array */
	for (j = fstVtx_toUpd; j < n; j++)
	  marker[j] = EMPTY;
	markl = 0;
      }
    }
  }
  *pmarkl = markl;

  return 0;
}

static int_t
dnsUpSeps_symbfact
(
 int_t n,
 int   iam,      /* my processor number */
 int   szSep, 
 int   ind_sizes1,
 int   ind_sizes2, 
 int_t *sizes,     /* Input - sizes of each node in the separator tree */
 int_t *fstVtxSep, /* Input - first vertex of each node in the tree */
 int_t fstVtx_dns,
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 Pslu_freeable_t *Pslu_freeable,
 vtcsInfo_symbfact_t *VInfo,  /* Input/Output - local info on vertices distribution */
 comm_symbfact_t *CS,
 psymbfact_stat_t *PS,
 int_t *p_nextl,   /* ptr to nextl in lsub structure */
 int_t *p_nextu,   /* ptr to nextu in usub structure */
 int_t *p_nsuper_loc
 )
{
  int_t nextl, nextu, nsuper_loc, curblk_loc, mem_error;
  int_t vtx_elt, ind_blk, vtx, k;
  int_t *xlsub, *xusub, *lsub, *usub;
  int_t fstVtx_blk, fstVtx_blk_lid, vtx_lid, lstVtx_blk, fstVtx_lvl, lstVtx_lvl;
  int_t *globToLoc, maxNvtcsPProc;
  
  /* Initialization */
  xlsub = Llu_symbfact->xlsub; lsub = Llu_symbfact->lsub;
  xusub = Llu_symbfact->xusub; usub = Llu_symbfact->usub;

  globToLoc  = Pslu_freeable->globToLoc;
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  nextl      = *p_nextl;
  nextu      = *p_nextu;
  nsuper_loc = *p_nsuper_loc;
  curblk_loc = VInfo->curblk_loc;
  VInfo->nnz_ainf_loc = 0;
  VInfo->nnz_asup_loc = 0;

  if (fstVtx_dns == EMPTY)
    fstVtx_blk     = VInfo->begEndBlks_loc[curblk_loc];
  else 
    fstVtx_blk  = fstVtx_dns;
  if (fstVtx_blk == n)
    return 0;
  fstVtx_blk_lid = LOCAL_IND( globToLoc[fstVtx_blk] );
  vtx_lid        = fstVtx_blk_lid;
  xlsub[vtx_lid] = nextl;
  xusub[vtx_lid] = nextu;
  PS->nDnsUpSeps = 0; 
  
  while (szSep >= 1) {
    PS->nDnsUpSeps++; 
    fstVtx_lvl = fstVtxSep[ind_sizes2];
    lstVtx_lvl = fstVtxSep[ind_sizes2] + sizes[ind_sizes2];
    if (fstVtx_blk > fstVtx_lvl)
      vtx_elt = fstVtx_blk;
    else 
      vtx_elt = fstVtx_lvl;
    if (nextl + lstVtx_lvl - vtx_elt >= Llu_symbfact->szLsub) {
      if (mem_error =
	  psymbfact_LUXpandMem (iam, n, fstVtx_blk, nextl, 
				nextl + fstVtx_lvl - vtx_elt,
				LSUB, DNS_UPSEPS, 1,
				Pslu_freeable, Llu_symbfact, VInfo, PS))
	return (mem_error);
      lsub = Llu_symbfact->lsub;
    }
    if (nextu + lstVtx_lvl - vtx_elt >= Llu_symbfact->szUsub) {
      if (mem_error =
	  psymbfact_LUXpandMem (iam, n, fstVtx_blk, nextu, 
				nextu + fstVtx_lvl - vtx_elt,
				LSUB, DNS_UPSEPS, 1,
				Pslu_freeable, Llu_symbfact, VInfo, PS))
	return (mem_error);
      usub = Llu_symbfact->usub;
    }
    PS->nops += 2 * (lstVtx_lvl - vtx_elt);
    for (; vtx_elt < lstVtx_lvl; vtx_elt++) {
      lsub[nextl] = vtx_elt; nextl++;
      usub[nextu] = vtx_elt; nextu++;
    }
    ind_sizes2 = ind_sizes1 + szSep + (ind_sizes2 - ind_sizes1) / 2;
    ind_sizes1 += szSep;
    szSep = szSep / 2;
  }
  /* delete the diagonal element from the U structure */
  usub[xusub[fstVtx_blk_lid]] = usub[nextu - 1];
  nextu --;
  xlsub[fstVtx_blk_lid+1] = nextl;
  xusub[fstVtx_blk_lid+1] = nextu;

  vtx_lid = fstVtx_blk_lid;
  ind_blk = curblk_loc;
  while (ind_blk < 2 * VInfo->nblks_loc) {
    if (ind_blk != curblk_loc) {
      fstVtx_blk = VInfo->begEndBlks_loc[ind_blk];

      xlsub[vtx_lid] = nextl;
      xusub[vtx_lid] = nextu;

      for (k = xlsub[fstVtx_blk_lid]; k < xlsub[fstVtx_blk_lid+1]; k++) 
	if (lsub[k] >= fstVtx_blk) {
	  lsub[nextl] = lsub[k]; nextl ++;
	  if (nextl >= MEM_LSUB( Llu_symbfact, VInfo ))
 	    if (mem_error =
		psymbfact_LUXpandMem (iam, n, fstVtx_blk, nextl, 0,
				      LSUB, DNS_UPSEPS, 1,
				      Pslu_freeable, Llu_symbfact, VInfo, PS))
	      return (mem_error);
	  lsub = Llu_symbfact->lsub;
	}
      for (k = xusub[fstVtx_blk_lid]; k < xusub[fstVtx_blk_lid+1]; k++)
	if (usub[k] > fstVtx_blk) {
	  usub[nextu] = usub[k]; nextu ++;
	  if (nextu >= MEM_USUB( Llu_symbfact, VInfo ))
	    if (mem_error =
		psymbfact_LUXpandMem (iam, n, fstVtx_blk, nextu, 0,
				      USUB, DNS_UPSEPS, 1,
				      Pslu_freeable, Llu_symbfact, VInfo, PS))
	      return (mem_error);
	  usub = Llu_symbfact->usub;
	}
      PS->nops += xlsub[fstVtx_blk_lid+1] - xlsub[fstVtx_blk_lid];
      PS->nops += xusub[fstVtx_blk_lid+1] - xusub[fstVtx_blk_lid];
    }
    lstVtx_blk = VInfo->begEndBlks_loc[ind_blk + 1];
    for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid++) {
      Pslu_freeable->supno_loc[vtx_lid] = nsuper_loc;
      if (vtx > fstVtx_blk) {
	xlsub[vtx_lid] = nextl;
	xusub[vtx_lid] = nextu;
      }
    }
    ind_blk += 2;
    nsuper_loc ++;
  }
  
  *p_nextl = nextl;
  *p_nextu = nextu;
  *p_nsuper_loc = nsuper_loc;
/*   VInfo->curblk_loc = ind_blk; */
  
  return 0;
}

static int_t
dnsCurSep_symbfact
(
 int_t n,          /* Input - order of the matrix */
 int   iam,        /* Input - my processor number */
 int   ind_sizes1,
 int   ind_sizes2,
 int_t *sizes,     /* Input - sizes of each node in the separator tree */
 int_t *fstVtxSep, /* Input - first vertex of each node in the tree */
 int   szSep,
 int   npNode,
 int_t rcvd_dnsSep,
 int_t *p_nextl,     
 int_t *p_nextu, 
 int_t *p_mark,
 int_t *p_nsuper_loc,
 int_t *marker,   /* temporary array of size n */
 MPI_Comm ndCom,
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 Pslu_freeable_t *Pslu_freeable,
 vtcsInfo_symbfact_t *VInfo,  /* Input/Output - local info on vertices distribution */
 comm_symbfact_t *CS,
 psymbfact_stat_t *PS
 )
{
  int_t fstVtx_blk, fstVtx_dns, fstVtx_dns_lid, lstVtx_blk, 
    fstVtx, lstVtx, lstVtx_dns_lid;
  int_t ind_blk, i, vtx, vtx_lid, vtx_lid_x, nvtcs_upd, save_cnt, mem_error;
  int_t computeL, computeU, vtx_elt, j, cur_blk, snlid, snrep;
  int_t *sub, *xsub, *minElt_vtx, *cntelt_vtcs;
  int_t mark, next, *x_newelts, *x_newelts_L, *x_newelts_U;
  int_t *newelts_L, *newelts_U, *newelts;
  int_t *globToLoc, maxNvtcsPProc, lvl;
  int_t prval, kmin, kmax, maxElt, ktemp, prpos;
  float mem_dnsCS;

  if (!rcvd_dnsSep)
    VInfo->curblk_loc += 2;
  
  computeL = TRUE; computeU = TRUE;
  lstVtx_dns_lid = EMPTY;
  globToLoc = Pslu_freeable->globToLoc;
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  fstVtx = fstVtxSep[ind_sizes2];
  lstVtx = fstVtx + sizes[ind_sizes2];
  cur_blk = VInfo->curblk_loc;
  fstVtx_dns = VInfo->begEndBlks_loc[cur_blk];
  fstVtx_dns_lid = LOCAL_IND( globToLoc[fstVtx_dns] );
  lvl = (int_t) LOG2( npNode );
  x_newelts_U = NULL;
  newelts_L = NULL;
  newelts_U = NULL;
  mem_dnsCS = 0.;
  
  PS->nDnsCurSep ++;

  if (CS->rcv_bufSz > n - fstVtx_dns)
    minElt_vtx = CS->rcv_buf;
  else {
    if (!(minElt_vtx = intMalloc_symbfact(n - fstVtx_dns)))
      ABORT("Malloc fails for minElt_vtx[].");
    mem_dnsCS += n - fstVtx_dns;
  }
  
  while (computeL || computeU) {
    if (computeL) {
      sub = Llu_symbfact->lsub; xsub = Llu_symbfact->xlsub;
      x_newelts = Llu_symbfact->cntelt_vtcs;
      x_newelts_L = x_newelts;
    }
    else {
      sub = Llu_symbfact->usub; xsub = Llu_symbfact->xusub;
    }

    /* use minElt_vtx to determine starting vertex of each nonzero element */
    for (i = 0; i < n - fstVtx_dns; i++)
      minElt_vtx[i] = n;

    ind_blk = cur_blk;
    vtx_lid = fstVtx_dns_lid;
    nvtcs_upd = 0;
    while (VInfo->begEndBlks_loc[ind_blk] < lstVtx && 
	   ind_blk < 2 * VInfo->nblks_loc) {	  
      fstVtx_blk = VInfo->begEndBlks_loc[ind_blk];
      lstVtx_blk = VInfo->begEndBlks_loc[ind_blk + 1];
      ind_blk += 2;
      nvtcs_upd += lstVtx_blk - fstVtx_blk;
      for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid++) {
	j = xsub[vtx_lid];
	while (j < xsub[vtx_lid+1] && sub[j] != EMPTY) {
	  PS->nops ++;
	  vtx_elt = sub[j] - fstVtx_dns;
	  if (minElt_vtx[vtx_elt] == n) {
	    minElt_vtx[vtx_elt] = vtx;
	  }
	  j ++;
	}
      }	  
    }
    if (!computeL) {
      if (!(x_newelts_U = intMalloc_symbfact(nvtcs_upd + 1)))
	ABORT("Malloc fails for x_newelts_U[].");
      mem_dnsCS += nvtcs_upd + 1;
      x_newelts = x_newelts_U;
    }
    else {
      /* save the value in cntelt_vtcs[lstVtx_blk_lid] */
      save_cnt = x_newelts[vtx_lid];
      lstVtx_dns_lid = vtx_lid;
    }
    
    MPI_Allreduce (&(minElt_vtx[lstVtx - fstVtx_dns]), &(marker[lstVtx]), 
		   n - lstVtx, mpi_int_t, MPI_MIN, ndCom);

#if ( PRNTlevel>=1 )
    PS->no_msgsCol += (float) (2 * (int_t) LOG2( npNode ));
    PS->sz_msgsCol += (float) (n - lstVtx);
    if (PS->maxsz_msgCol < n - lstVtx) 
      PS->maxsz_msgCol = n - lstVtx;      
#endif
    
    /* use x_newelts to determine counts of elements starting in each vertex */
    for (vtx_lid = 0; vtx_lid < nvtcs_upd; vtx_lid++)
      x_newelts[vtx_lid] = 0;
    
    for (vtx = lstVtx; vtx < n; vtx++) {
      if (marker[vtx] != n) {
	vtx_elt = marker[vtx];
	if (OWNER( globToLoc[vtx_elt] ) == iam) {
	  x_newelts[ LOCAL_IND( globToLoc[vtx_elt] ) - fstVtx_dns_lid ] ++;
	}
	else {
	  /* find the first vertex > vtx_elt which belongs to iam */
	  ind_blk = cur_blk;
	  vtx_lid = 0;
	  while (vtx_elt > VInfo->begEndBlks_loc[ind_blk] &&
		 ind_blk < 2 * VInfo->nblks_loc) {
	    vtx_lid += VInfo->begEndBlks_loc[ind_blk+1] -
	      VInfo->begEndBlks_loc[ind_blk];
	    ind_blk += 2;
	  }
	  if (VInfo->begEndBlks_loc[ind_blk] < lstVtx) {
	    x_newelts[vtx_lid] ++;		
	    marker[vtx] = VInfo->begEndBlks_loc[ind_blk];
	  }			    
	  else
	    marker[vtx] = n;
	}
      }
    }
    
    /* set up beginning of new elements for each local vtx */
    i = 0;
    for (vtx_lid = 0; vtx_lid < nvtcs_upd; vtx_lid++) {
      j = x_newelts[vtx_lid];
      x_newelts[vtx_lid] = i;
      i += j;
    }
    x_newelts[vtx_lid] = i;
    newelts = NULL;
    if (i != 0) {
      if (!(newelts = intMalloc_symbfact(x_newelts[vtx_lid])))
	ABORT("Malloc fails for newelts[].");    
      mem_dnsCS += x_newelts[vtx_lid];
      
      for (vtx = lstVtx; vtx < n; vtx++) {
	if (marker[vtx] != n) {
	  vtx_elt = marker[vtx];
	  vtx_lid = LOCAL_IND( globToLoc[vtx_elt] ) - fstVtx_dns_lid;	  
	  newelts[x_newelts[vtx_lid]] = vtx;
	  x_newelts[vtx_lid] ++;
	}
      }
    }
    /* reset beginning of new elements for each local vertex */
    i = 0;
    for (vtx_lid = 0; vtx_lid < nvtcs_upd; vtx_lid++) {
      j = x_newelts[vtx_lid];
      x_newelts[vtx_lid] = i;
      i = j;
    }

    if (computeL == TRUE) {
      computeL = FALSE;
      newelts_L = newelts;
    }
    else {
      computeU = FALSE;
      newelts_U = newelts;
    }
  }
  
  for (i = fstVtx_dns; i < n; i++)
    marker[i] = EMPTY;
  mark = 0;
  
  /* update vertices */
  prval = n; 	    
  ind_blk = cur_blk;
  fstVtx_dns = VInfo->begEndBlks_loc[ind_blk];
  vtx_lid = LOCAL_IND( globToLoc[fstVtx_dns] );
  while (VInfo->begEndBlks_loc[ind_blk] < lstVtx && 
	 ind_blk < 2 * VInfo->nblks_loc) {	  
    fstVtx_blk = VInfo->begEndBlks_loc[ind_blk];
    lstVtx_blk = VInfo->begEndBlks_loc[ind_blk + 1];
    ind_blk += 2;
    for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid++) {
      vtx_lid_x = vtx_lid - fstVtx_dns_lid;
      Llu_symbfact->xlsub[vtx_lid] = *p_nextl;
      Llu_symbfact->xusub[vtx_lid] = *p_nextu;
      if (vtx == fstVtx_blk || x_newelts_L[vtx_lid_x+1] != x_newelts_L[vtx_lid_x] ||
	  x_newelts_U[vtx_lid_x+1] != x_newelts_U[vtx_lid_x]) {
	/* a new supernode starts */
	snlid = vtx_lid;
	snrep = vtx;
	if (mark + 2 > n) {
	  /* reset to EMPTY marker array */
	  for (i = 0; i < n; i++)
	    marker[i] = EMPTY;
	  mark = 0;
	}

	computeL = TRUE;
	computeU = FALSE;
	while (computeL || computeU) {
	  if (computeL) {
	    sub = Llu_symbfact->lsub; xsub = Llu_symbfact->xlsub;
	    x_newelts = x_newelts_L; newelts = newelts_L;
	    next = *p_nextl;
	  }
	  else {
	    sub = Llu_symbfact->usub; xsub = Llu_symbfact->xusub;
	    x_newelts = x_newelts_U; newelts = newelts_U;
	    next = *p_nextu;
	  }
	  xsub[vtx_lid] = next;

	  /* TEST available memory */
	  j = x_newelts[vtx_lid_x+1] + lstVtx - vtx;
	  if ((computeL && next+j >= MEM_LSUB(Llu_symbfact, VInfo)) ||
	      (computeU && next+j >= MEM_USUB(Llu_symbfact, VInfo))) {
	    if (mem_error =
		psymbfact_LUXpandMem (iam, n, vtx, next, next + j,
				      computeL, DNS_CURSEP, 1,
				      Pslu_freeable, Llu_symbfact, VInfo, PS))
	      return (mem_error);
	    if (computeL) sub = Llu_symbfact->lsub;
	    else sub = Llu_symbfact->usub; 
	  }
	  
	  if (computeL)  i = vtx;
	  else           i = vtx+1;
	  while (i < lstVtx) {
	    sub[next] = i; next ++;
	    i ++;
	  }
	  PS->nops += x_newelts[vtx_lid_x+1];
	  for (i = 0; i < x_newelts[vtx_lid_x+1]; i++) {
	    vtx_elt = newelts[i];
	    sub[next] = vtx_elt; next ++;
	    if (computeU && vtx_elt < prval 
		&& marker[vtx_elt] == mark-1)
	      prval = vtx_elt;
	    marker[vtx_elt] = mark;
	  }
	  if (computeL) {
	    computeL = FALSE; computeU = TRUE;
	    *p_nextl = next;
	  }
	  else {
	    computeU = FALSE;
	    *p_nextu = next;
	  }
	  mark ++;	  
	}	  
	if (vtx != fstVtx_blk)
	  (*p_nsuper_loc) ++;
      } /* a new supernode starts */
      /* vtx belongs to the curent supernode */
      Pslu_freeable->supno_loc[vtx_lid] = *p_nsuper_loc;
    } 
    (*p_nsuper_loc) ++;
  }
  
  if (ind_blk > 0) {
    /* if iam owns blocks of this level */
    i = *p_nextl - Llu_symbfact->xlsub[snlid];
    j = *p_nextu - Llu_symbfact->xusub[snlid];
    
    if (VInfo->begEndBlks_loc[ind_blk - 1] == lstVtx && i > 1 && j > 0) {
      /* if iam the last processor owning a block of this level */
      computeL = TRUE; computeU = FALSE;
      /* prune the structure */
      while (computeL || computeU) {	   
	if (computeL) {
	  sub = Llu_symbfact->lsub; xsub = Llu_symbfact->xlsub;
	  next = *p_nextl;
	  computeL = FALSE; computeU = TRUE;
	}
	else {
	  sub = Llu_symbfact->usub; xsub = Llu_symbfact->xusub;
	  next = *p_nextu;
	  computeU = FALSE;
	}
	
	kmin = xsub[snlid];
	kmax = next - 1;
	if (prval != n) {
	  maxElt = prval;
	  while (kmin <= kmax) {
	    /* Do a quicksort-type partition. */    
	    if (sub[kmax] > prval) 
	      kmax--;
	    else if (sub[kmin] <= prval) {
	      kmin++;
	    }
	    else { /* kmin does'nt belong to G^s(L), and kmax belongs: 
		    * 	   interchange the two subscripts
		    */
	      ktemp = sub[kmin];
	      sub[kmin] = sub[kmax];
	      sub[kmax] = ktemp;
	      kmin ++;
	      kmax --;
	    }
	    if (sub[kmin-1] == prval)
	      prpos = kmin - 1;
	  }
	}
	else {
	  maxElt = EMPTY;
	  while (kmin <= kmax) {
	    /* compute maximum element of L(:, vtx) */
	    if (sub[kmin] > maxElt) {
	      maxElt = sub[kmin];
	      prpos = kmin;
	    }
	    kmin ++;	      
	  }
	}
	ktemp = sub[xsub[snlid]];
	sub[xsub[snlid]] = maxElt;
	sub[prpos] = ktemp;
      }
      
      /* setup snd_interSz information */
      prval = Llu_symbfact->lsub[Llu_symbfact->xlsub[snlid]];
      if (prval >= lstVtx) {
	/* this supernode will be send to next layers of the tree */
	while (prval >= lstVtx && szSep != 1) {
	  ind_sizes2 = ind_sizes1 + szSep + (ind_sizes2 - ind_sizes1) / 2;
	  ind_sizes1 += szSep;
	  lvl ++;
	  szSep = szSep / 2;
	  lstVtx = fstVtxSep[ind_sizes2] + sizes[ind_sizes2];
	  CS->snd_interSz[lvl] += i + j + 4;
	  CS->snd_LinterSz[lvl] += i + 2;
	  if (CS->snd_vtxinter[lvl] == EMPTY)
	    CS->snd_vtxinter[lvl] = snrep;
	}
      }
    }
  }

  /* restore value in cntelt_vtcs */
  if (lstVtx_dns_lid != EMPTY)
    Llu_symbfact->cntelt_vtcs[lstVtx_dns_lid] = save_cnt;
  *p_mark = mark;
  if (minElt_vtx != CS->rcv_buf)
    SUPERLU_FREE (minElt_vtx);  
  SUPERLU_FREE (x_newelts_U);
  if (newelts_L) SUPERLU_FREE (newelts_L);
  if (newelts_U) SUPERLU_FREE (newelts_U);
  if (PS->szDnsSep < mem_dnsCS)
    PS->szDnsSep = mem_dnsCS;

  return 0;
}

/*! \brief

<pre>
   All processors affected to current node must call this routine
   when VInfo->filledSep == FILLED_SEP
   This is necessary since subsequent routines called from here use 
   MPI_allreduce among all processors affected to curent node
</pre>
*/

static int_t
denseSep_symbfact 
(
 int   rcvd_dnsSep, /* =1 if processor received info that the separator
		       became dense,
		       =0 if myPE determined that separator is full */
 int_t n,           /* Input - order of the matrix */
 int   iam,         /* Input - my processor number */
 int   ind_sizes1,
 int   ind_sizes2,
 int_t *sizes,     /* Input - sizes of each separator in the separator tree */
 int_t *fstVtxSep, /* Input - first vertex of each node in the tree */
 int   szSep, 
 int   fstP,        /* first pe affected current node */
 int   lstP,        /* last pe affected current node */
 int_t fstVtx_blkCyc, 
 int_t nblk_loc,    /* block number in the block cyclic distribution of current
		       supernode */
 int_t *p_nextl,
 int_t *p_nextu,
 int_t *p_mark,
 int_t *p_nsuper_loc,
 int_t *marker,
 MPI_Comm ndCom,
 MPI_Comm *symb_comm, /* Input - communicator for symbolic factorization */
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 Pslu_freeable_t *Pslu_freeable,
 vtcsInfo_symbfact_t *VInfo,  /* Input - local info on vertices distribution */
 comm_symbfact_t *CS,
 psymbfact_stat_t *PS
) 
{
  int   nprocsLvl, p, prvP, tag;
  int_t nmsgsToSnd, nmsgsToRcv;
  int_t ind_blk, mem_error;
  int_t *rcv_intraLvl;
  int_t fstVtx, lstVtx, cur_blk, lstVtx_blk, fstVtx_blk;
  int_t *globToLoc, maxNvtcsPProc;
  MPI_Status status;
  
  globToLoc = Pslu_freeable->globToLoc;
  maxNvtcsPProc = Pslu_freeable->maxNvtcsPProc;
  fstVtx = fstVtxSep[ind_sizes2];
  lstVtx = fstVtx + sizes[ind_sizes2];
  rcv_intraLvl = CS->rcv_intraLvl;
  cur_blk   = VInfo->curblk_loc;
  nprocsLvl = lstP - fstP;
  
  if (nblk_loc == 0) {
    nmsgsToSnd = 2; nmsgsToRcv = 1;
  }
  else {
    nmsgsToSnd = 1; nmsgsToRcv = 0;
    if (!rcvd_dnsSep) nmsgsToRcv ++;
  }
  if (iam == fstP && rcvd_dnsSep && nblk_loc == 1) 
    nmsgsToRcv ++;
  
  /* first exchange msgs with all processors affected to current node */
  ind_blk = cur_blk;
  while ((nmsgsToSnd || nmsgsToRcv) && VInfo->begEndBlks_loc[ind_blk] < lstVtx) {
    tag = (int) (tag_intraLvl + nblk_loc);
    if (nmsgsToSnd) {
      lstVtx_blk = VInfo->begEndBlks_loc[ind_blk + 1];
      if (lstVtx_blk != lstVtx) {
	p = OWNER( globToLoc[lstVtx_blk]);
	MPI_Send (&(rcv_intraLvl[fstP]), nprocsLvl, mpi_int_t, p,
		  tag, (*symb_comm));
#if ( PRNTlevel>=1 )
	PS->no_shmSnd += (float) 1;
#endif
      }
      nmsgsToSnd --;
    }
    ind_blk += 2;
    nblk_loc ++;
    tag = tag_intraLvl + nblk_loc;
    fstVtx_blk = VInfo->begEndBlks_loc[ind_blk];
    if (nmsgsToRcv && fstVtx_blk < lstVtx) {
      if (iam == fstP) tag --;
      prvP = OWNER( globToLoc[fstVtx_blk - 1]);
      MPI_Recv (&(rcv_intraLvl[fstP]), nprocsLvl, mpi_int_t, prvP,
		tag, (*symb_comm), &status);
#if ( PRNTlevel>=1 )
      PS->no_shmRcvd += (float) 1;
#endif
      nmsgsToRcv --;
    }
  }

  if (VInfo->filledSep == FILLED_SEP) {
    if (mem_error = 
	dnsCurSep_symbfact (n, iam, ind_sizes1, ind_sizes2, sizes, fstVtxSep, 
			    szSep, lstP - fstP, rcvd_dnsSep, p_nextl, 
			    p_nextu, p_mark, p_nsuper_loc, marker, ndCom,
			    Llu_symbfact, Pslu_freeable, VInfo, CS, PS))
      return (mem_error);
  }
  else if (rcvd_dnsSep) 
    if (mem_error = 
	dnsUpSeps_symbfact (n, iam, szSep, ind_sizes1, ind_sizes2, 
			    sizes, fstVtxSep, EMPTY,
			    Llu_symbfact, Pslu_freeable, VInfo, CS, PS,
			    p_nextl, p_nextu, p_nsuper_loc))
      return (mem_error);
  return 0;
}


static int_t
interLvl_symbfact
(
 SuperMatrix *A, /* Input - input matrix A */
 int   iam,      /* Input - my processor number */  
 int   lvl,      /* Input - current level in the separator tree */ 
 int   szSep,    /* Input - size of the current separator (node) */
 int   fstP,     /* Input - first processor assigned to current node */
 int   lstP,     /* Input - last processor assigned to current node */
 int   ind_sizes1,
 int   ind_sizes2, 
 int_t *sizes,     /* Input - sizes of each node in the separator tree */
 int_t *fstVtxSep, /* Input - first vertex of each node in the tree */
 int_t *p_nextl,
 int_t *p_nextu,
 int_t *p_nsuper_loc,
 int_t *pmark,   /* mark for symbfact */
 int_t *marker,  /* temp array used for marking */
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 Pslu_freeable_t *Pslu_freeable,
 comm_symbfact_t *CS,/* infos on communication data structures */
 vtcsInfo_symbfact_t *VInfo, /* Input/Output - local info on vertices distribution */
 psymbfact_stat_t *PS,
 MPI_Comm ndComm,
 MPI_Comm    *symb_comm /* Input - communicator for symbolic factorization */
 )
{
  MPI_Status  *status; 
  MPI_Request *request_snd, *request_rcv;
  
  int   nprocsLvl, rcvdP, p, filledSep_lvl;
  int   toSend, toSendL, toSendU;
  int_t *rcv_interLvl;
  int_t *snd_interLvl, *snd_interLvl1, *snd_interLvl2,
    snd_interLvlSz, snd_LinterLvlSz, snd_vtxLvl;
  int_t  vtx_elt, update_loc, code_err;
  int_t *lsub, *xlsub, *usub, *xusub;
  int_t *lsub_rcvd, lsub_rcvd_sz, *usub_rcvd, usub_rcvd_sz;
  int_t  n, mark, max_rcvSz; 
  int_t nextl, nextu, ind_blk, vtx_lid, k, count, nelts, 
    lstVtxLvl_loc, lstVtxLvl_loc_lid, mem_error;
  int_t fstVtx_blk, lstVtx_blk, i, j, vtx, prElt_L, prElt_U, 
    snd_indBlk, prElt_ind;
  int_t fstVtxLvl_loc, nvtcsLvl_loc, maxNvtcsPProc, *globToLoc, 
    fstVtx, lstVtx;
  int  ind1, nprocsToRcv, nprocsToSnd, ind2, ind_l, ind_u, ij, ik;
  int_t req_ind, sent_msgs, req_ind_snd;
  int_t initInfo_loc[2], initInfo_gl[2];

  /* Initialization */
  n = A->ncol;
  fstVtx          = fstVtxSep[ind_sizes2];
  lstVtx          = fstVtx + sizes[ind_sizes2];
  maxNvtcsPProc   = Pslu_freeable->maxNvtcsPProc;
  globToLoc       = Pslu_freeable->globToLoc;
  nprocsLvl       = lstP - fstP;
  rcv_interLvl    = CS->rcv_interLvl;
  snd_interLvl    = CS->snd_interLvl;
  snd_interLvlSz  = CS->snd_interSz[lvl];
  snd_LinterLvlSz = CS->snd_LinterSz[lvl];
  snd_vtxLvl      = CS->snd_vtxinter[lvl];
  fstVtxLvl_loc   = VInfo->begEndBlks_loc[VInfo->curblk_loc];
  nvtcsLvl_loc    = VInfo->nvtcsLvl_loc;
  request_snd = NULL;
  request_rcv = NULL;
  status = NULL;
  mark = *pmark;
  
  lsub    = Llu_symbfact->lsub;   xlsub    = Llu_symbfact->xlsub;
  usub    = Llu_symbfact->usub;   xusub    = Llu_symbfact->xusub;

  /* snd_vtxLvl denotes the first vertex from which iam needs
     to send data.  
     snd_interLvlSz denotes maximum size of the send data,
     snd_LinterLvlSz denotes send data corresponding to L part */

  /* determine maximum size of receive buffer and information
   on filled sep */
  if (snd_interLvlSz != 0) {
    if (snd_LinterLvlSz == 0) 
      snd_interLvlSz = 0;
    if (snd_interLvlSz - snd_LinterLvlSz == 0)
      snd_interLvlSz = 0;
  }
  
  initInfo_loc[0] = snd_interLvlSz;
  initInfo_loc[1] = (int_t) VInfo->filledSep;
  MPI_Allreduce (initInfo_loc, initInfo_gl, 2, 
		 mpi_int_t, MPI_MAX, ndComm);
#if ( PRNTlevel>=1 )
  PS->no_msgsCol += (float) (2 * (int_t) LOG2( nprocsLvl ));
  PS->sz_msgsCol += 2;
  if (PS->maxsz_msgCol < 2) 
    PS->maxsz_msgCol = 2;      
#endif  
  max_rcvSz = initInfo_gl[0];
  filledSep_lvl = (int) initInfo_gl[1];

  if (filledSep_lvl == FILLED_SEPS) {
    /* quick return if all upper separators are dense */
    if (VInfo->filledSep != FILLED_SEPS) {
      VInfo->filledSep = FILLED_SEPS;
      if (mem_error = 
	  dnsUpSeps_symbfact (n, iam, szSep, ind_sizes1, ind_sizes2, sizes, 
			      fstVtxSep,
			      EMPTY, Llu_symbfact, Pslu_freeable, VInfo, CS, PS,
			      p_nextl, p_nextu, p_nsuper_loc))
	return (mem_error);
    }
    return 0;
  }

  if (max_rcvSz == 0)
    /* quick return if no communication necessary */
    return 0; 
  
  /* allocate data for the send buffer */  
  if (snd_interLvlSz)
    if (CS->snd_bufSz < snd_interLvlSz) {
      PS->maxSzBuf += snd_interLvlSz - CS->snd_bufSz;
      if (CS->snd_bufSz != 0)
	/* not first time allocate memory */
	SUPERLU_FREE (CS->snd_buf);
      CS->snd_bufSz = snd_interLvlSz;
      if (!(CS->snd_buf = intMalloc_symbfact (snd_interLvlSz))) {
	ABORT("Malloc fails for snd_buf[].");
      }
    }
    
  /* snd_interLvl : to which processors the data need to be send 
   * information setup during the copy of data to be send in the buffer  
   * rcv_interLvl : from which processors iam receives update data  */
  for (p = 2*fstP; p < 2*lstP; p++)
    snd_interLvl[p] = EMPTY;

  if (snd_interLvlSz == 0 && nvtcsLvl_loc == 0) {
    code_err = MPI_Alltoall (&(snd_interLvl[2*fstP]), 2, mpi_int_t,
			     &(rcv_interLvl[2*fstP]), 2, mpi_int_t,
			     ndComm);
#if ( PRNTlevel>=1 )
    PS->no_msgsCol += (float) (2 * (int_t) LOG2( nprocsLvl ));
    PS->sz_msgsCol += 2;
    if (PS->maxsz_msgCol < 2) 
      PS->maxsz_msgCol = 2;      
#endif  
    return 0;
  }
  
  /* in interLvlInfos, 
   * obtain from which processors iam receives update information */
  update_loc = FALSE;
  nextl = 0; 
  nextu = snd_LinterLvlSz;
  if (snd_interLvlSz != 0) {
    /* copy data to be send */
    /* find index block from where to send data */
    ind_blk = VInfo->curblk_loc;
    while (snd_vtxLvl < VInfo->begEndBlks_loc[ind_blk]) {
      ind_blk -= 2;
    }
    snd_indBlk = ind_blk;
    vtx_lid = LOCAL_IND( globToLoc[snd_vtxLvl] );
    for (; ind_blk < VInfo->curblk_loc; ind_blk += 2) {
      fstVtx_blk = VInfo->begEndBlks_loc[ind_blk];
      if (ind_blk == snd_indBlk)
	fstVtx_blk = snd_vtxLvl;
      lstVtx_blk = VInfo->begEndBlks_loc[ind_blk + 1];
      for (vtx = fstVtx_blk; vtx < lstVtx_blk; vtx++, vtx_lid ++) {
	toSendL = FALSE; toSendU = FALSE;
	if (xlsub[vtx_lid] != xlsub[vtx_lid+1] && 
	    xusub[vtx_lid] != xusub[vtx_lid+1]) {
	  k = xlsub[vtx_lid];
	  prElt_L = lsub[k];
	  j = xusub[vtx_lid];
	  prElt_U = usub[j];
	  if (prElt_L >= fstVtx || prElt_U >= fstVtx) {
	    if (prElt_L >= fstVtx)
	      while (lsub[k] <= prElt_L && k < xlsub[vtx_lid + 1]) {
		vtx_elt = lsub[k];
		if (vtx_elt >= fstVtx && vtx_elt < lstVtx) {
		  p = OWNER( globToLoc[vtx_elt] );
		  if (p != iam) {
		    /* vtx will be send to another processor */
		    snd_interLvl[2*p] = TRUE;
		    toSendL = TRUE;
		  }
		  else
		    update_loc = TRUE;
		}
		k++;
	      }
	    if (prElt_U >= fstVtx)
	      while (usub[j] <= prElt_U && j < xusub[vtx_lid + 1]) {
		vtx_elt = usub[j];
		if (vtx_elt >= fstVtx && vtx_elt < lstVtx) {
		  p = OWNER( globToLoc[vtx_elt] );
		  if (p != iam) {
		    /* vtx will be send to another processor */
		    snd_interLvl[2*p+1] = TRUE;
		    toSendU = TRUE;
		  }
		  else
		    update_loc = TRUE;
		}
		j ++;
	      }
	    if (toSendL || toSendU) {
	      /* L(:, vtx) and U(vtx, :) will be send to processors */
	      CS->snd_buf[nextu + DIAG_IND]  = vtx;
	      nelts = xusub[vtx_lid+1] - xusub[vtx_lid];
	      CS->snd_buf[nextu + NELTS_IND] = nelts;
	      nextu += 2;
	      for (j = xusub[vtx_lid]; j < xusub[vtx_lid+1]; j++, nextu ++) {
		CS->snd_buf[nextu] = usub[j]; 
	      }
	      CS->snd_buf[nextl + DIAG_IND] = vtx;
	      nelts = xlsub[vtx_lid+1] - xlsub[vtx_lid];
	      CS->snd_buf[nextl + NELTS_IND] = nelts; 
	      nextl += 2;
	      for (j = xlsub[vtx_lid]; j < xlsub[vtx_lid+1]; j++, nextl ++) {
		CS->snd_buf[nextl] = lsub[j];
	      }
	    }
	  }
	}
      }
    }
    lstVtxLvl_loc = vtx;
    lstVtxLvl_loc_lid = vtx_lid;
  }
  
  if (nextl == 0 || nextu - snd_LinterLvlSz == 0) {
    for (p = 2*fstP; p < 2*lstP; p++)
      snd_interLvl[p] = EMPTY;
  }
  
  nprocsToSnd = 0;
  for (p = 2*fstP; p < 2*lstP; p +=2) {
    if (snd_interLvl[p] != EMPTY || snd_interLvl[p+1] != EMPTY) {
      snd_interLvl[p] = nextl;
      snd_interLvl[p+1] = nextu - snd_LinterLvlSz;
      nprocsToSnd ++;
    }
  }
  
  MPI_Alltoall (&(snd_interLvl[2*fstP]), 2, mpi_int_t,
		&(rcv_interLvl[2*fstP]), 2, mpi_int_t, ndComm);    
#if ( PRNTlevel>=1 )
  PS->no_msgsCol += (float) (2 * (int_t) LOG2( nprocsLvl ));
  PS->sz_msgsCol += 2 * nprocsLvl;
  if (PS->maxsz_msgCol < 2 * nprocsLvl) 
    PS->maxsz_msgCol = 2 * nprocsLvl;      
#endif    

  max_rcvSz = 0;
  nprocsToRcv = 0;
  for (p = 2*fstP; p < 2*lstP; p +=2) {
    CS->ptr_rcvBuf[p] = max_rcvSz;
    if (rcv_interLvl[p] != EMPTY) 
      max_rcvSz += rcv_interLvl[p];
    CS->ptr_rcvBuf[p+1] = max_rcvSz;
    if (rcv_interLvl[p+1] != EMPTY) 
      max_rcvSz += rcv_interLvl[p+1];
    if (rcv_interLvl[p] != EMPTY || rcv_interLvl[p+1] != EMPTY) 
      nprocsToRcv ++;
  }

  /* allocate data for the receive buffer */  
  if (CS->rcv_bufSz < max_rcvSz) {
    PS->maxSzBuf += max_rcvSz - CS->rcv_bufSz;
    if (CS->rcv_bufSz != 0) /* not first time allocate memory */
      SUPERLU_FREE (CS->rcv_buf);
    CS->rcv_bufSz = max_rcvSz;
    if (!(CS->rcv_buf = intMalloc_symbfact (max_rcvSz))) {
      ABORT("Malloc fails for rcv_buf[].");
    }
  }
  
  /* allocate memory for status arrays */
  if (nprocsToSnd)
    if ( !(request_snd = (MPI_Request*) 
	   SUPERLU_MALLOC(2 * nprocsToSnd * sizeof(MPI_Request))))
      ABORT("Not enough memory when allocating MPI_Request");
  if (nprocsToRcv)
    if ( !(request_rcv = (MPI_Request*) 
	   SUPERLU_MALLOC(2 * nprocsToRcv * sizeof(MPI_Request))))
      ABORT("Not enough memory when allocating MPI_Request");
  if (nprocsToRcv || nprocsToSnd)
    if ( !(status = (MPI_Status*) 
	   SUPERLU_MALLOC(2 * (lstP-fstP) * sizeof(MPI_Status))))
      ABORT("Not enough memory when allocating MPI_Request");
  
  /* determine if we have to send data */
  i = 0;
  for (toSend = fstP, p = 2*fstP; p < 2*lstP; toSend++, p+=2) 
    if (snd_interLvl[p] != EMPTY && toSend != iam) {
      MPI_Isend (CS->snd_buf, nextl, mpi_int_t, toSend,
		 tag_interLvl_LData, (*symb_comm), &(request_snd[2*i]));
      MPI_Isend (&(CS->snd_buf[snd_LinterLvlSz]), 
		 nextu - snd_LinterLvlSz, mpi_int_t, toSend,
		 tag_interLvl_UData, (*symb_comm), &(request_snd[2*i+1]));
      i++;
#if ( PRNTlevel>=1 )
      PS->no_msgsSnd += (float) 2;
      PS->sz_msgsSnd += (float) (nextl + nextu - snd_LinterLvlSz);
      if (PS->maxsz_msgSnd < nextl) PS->maxsz_msgSnd = nextl;
      if (PS->maxsz_msgSnd < nextu - snd_LinterLvlSz) 
	PS->maxsz_msgSnd = nextu - snd_LinterLvlSz;      
#endif
    }
  
  if (update_loc) {
    /* use own data to update symbolic factorization */
    vtx_lid = LOCAL_IND( globToLoc[snd_vtxLvl] );
    lsub_rcvd    = &(lsub[xlsub[vtx_lid]]);
    lsub_rcvd_sz = xlsub[lstVtxLvl_loc_lid] - xlsub[vtx_lid];
    usub_rcvd    = &(usub[xusub[vtx_lid]]);
    usub_rcvd_sz = xusub[lstVtxLvl_loc_lid] - xusub[vtx_lid];
    
    mem_error = 
      rl_update (0, n, iam, lsub_rcvd, lsub_rcvd_sz,
		 usub_rcvd, usub_rcvd_sz, snd_vtxLvl, EMPTY, snd_indBlk,
		 fstVtxLvl_loc, lstVtx, nvtcsLvl_loc,
		 1, &mark, marker, Pslu_freeable, Llu_symbfact, VInfo, PS);

    lsub_rcvd    = &(Llu_symbfact->lsub[xlsub[vtx_lid]]);
    lsub_rcvd_sz = xlsub[lstVtxLvl_loc_lid] - xlsub[vtx_lid];
    usub_rcvd    = &(Llu_symbfact->usub[xusub[vtx_lid]]);
    usub_rcvd_sz = xusub[lstVtxLvl_loc_lid] - xusub[vtx_lid];
    lsub = Llu_symbfact->lsub; usub = Llu_symbfact->usub;
    mem_error = 
      rl_update (0, n, iam, usub_rcvd, usub_rcvd_sz,
		 lsub_rcvd, lsub_rcvd_sz, snd_vtxLvl, EMPTY, snd_indBlk,
		 fstVtxLvl_loc, lstVtx, nvtcsLvl_loc,
		 0, &mark, marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
    lsub = Llu_symbfact->lsub; usub = Llu_symbfact->usub;
  }

  /* post non-blocking receives for all the incoming messages */
  i = 0;
  for (rcvdP = fstP, p = 2*fstP; p < 2*lstP; rcvdP++, p += 2) 
    if (rcv_interLvl[p] != EMPTY) {
      lsub_rcvd    = &(CS->rcv_buf[CS->ptr_rcvBuf[p]]);
      MPI_Irecv (lsub_rcvd, rcv_interLvl[p], mpi_int_t, rcvdP,
		 tag_interLvl_LData, (*symb_comm), &(request_rcv[i]));
      usub_rcvd    = &(CS->rcv_buf[CS->ptr_rcvBuf[p+1]]);
      MPI_Irecv (usub_rcvd, rcv_interLvl[p+1], mpi_int_t, rcvdP,
		 tag_interLvl_UData, (*symb_comm), &(request_rcv[i+1]));
      i += 2;
#if ( PRNTlevel>=1 )
      PS->no_msgsRcvd += (float) 2;
      PS->sz_msgsRcvd += (float) (rcv_interLvl[p] + rcv_interLvl[p+1]);
      if (PS->maxsz_msgRcvd < rcv_interLvl[p])
	PS->maxsz_msgRcvd = rcv_interLvl[p];
      if (PS->maxsz_msgRcvd < rcv_interLvl[p+1])
	PS->maxsz_msgRcvd = rcv_interLvl[p+1];
#endif
    }
  
  /* wait until messages are received and update local data */
  for (i = 0; i < nprocsToRcv; i++) {
    MPI_Waitany (2*nprocsToRcv, request_rcv, &ind1, status);
    ij = 0;
    for (p = fstP; p < lstP; p++)
      if (rcv_interLvl[2*p] != EMPTY) {
	if (ij <= ind1 && ind1 < ij+2) {
	  rcvdP = p; p = lstP;
	  if (ind1 == ij) ind2 = ij+1;
	  else ind2 = ind1 - 1;
	  ind_l = ij; ind_u = ij+1;
	}
	ij += 2;
      }
    MPI_Get_count (status, mpi_int_t, &ij);
    MPI_Wait (&(request_rcv[ind2]), status);
    MPI_Get_count (status, mpi_int_t, &ik);    
    if (ind1 == ind_l) {
      lsub_rcvd_sz = ij;
      usub_rcvd_sz = ik;
    } else {
      lsub_rcvd_sz = ik;
      usub_rcvd_sz = ij;
    }
    lsub_rcvd    = &(CS->rcv_buf[CS->ptr_rcvBuf[2*rcvdP]]);
    usub_rcvd    = &(CS->rcv_buf[CS->ptr_rcvBuf[2*rcvdP+1]]);
    
    /* use received data to update symbolic factorization information */
    mem_error = 
      rl_update (1, n, iam, lsub_rcvd, lsub_rcvd_sz,
		 usub_rcvd, usub_rcvd_sz, EMPTY, EMPTY, EMPTY,
		 fstVtxLvl_loc, lstVtx, nvtcsLvl_loc,
		 1, &mark, marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
    lsub = Llu_symbfact->lsub;
    mem_error = 
      rl_update (1, n, iam, usub_rcvd, usub_rcvd_sz,
		 lsub_rcvd, lsub_rcvd_sz, EMPTY, EMPTY, EMPTY,
		 fstVtxLvl_loc, lstVtx, nvtcsLvl_loc,
		 0, &mark, marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
    usub = Llu_symbfact->usub;      
  }
  
  if (nprocsToSnd)
    MPI_Waitall (2*nprocsToSnd, request_snd, status);

  *pmark = mark;
  if (request_snd != NULL) SUPERLU_FREE (request_snd);
  if (request_rcv != NULL) SUPERLU_FREE (request_rcv);
  if (status != NULL) SUPERLU_FREE (status);

  return 0;
}

static void
freeComm
(
 int   iam,          /* Input -my processor number */
 int   nprocs,       /* Input -number of procs for the symbolic fact. */
 MPI_Comm *commLvls, /* Input -communicators for the nodes in the sep tree */
 MPI_Comm *symb_comm /* Input - communicator for symbolic factorization */
 )
{
  int szSep, i, j, k;
  int np, npNode, fstP, lstP, ind;

  i = 2 * nprocs - 2;
  MPI_Comm_free (&(commLvls[i]));
  
  szSep = 2;
  i -= szSep;
  
  while (i > 0) {
    /* for each level in the separator tree */
    npNode = nprocs / szSep; 
    fstP = 0; 
    /* for each node in the level */
    for (j = i; j < i + szSep; j++) {
      lstP = fstP + npNode;
      if (fstP <= iam && iam < lstP) {
	ind = j;
      }
      fstP += npNode;
    }
    MPI_Comm_free ( &(commLvls[ind]) );
    szSep *= 2;
    i -= szSep;
  }  
}

static void
createComm 
(
 int   iam,          /* Input -my processor number */
 int   nprocs,       /* Input -number of procs for the symbolic factorization */
 MPI_Comm *commLvls, /* Output -communicators for the nodes in the sep tree */
 MPI_Comm *symb_comm
 )
{
  int szSep, i, j, jj, k, *pranks;
  int np, npNode, fstP, lstP, p, code_err, ind, col, key;
  
  for (i=0; i < 2*nprocs; i++)
    commLvls[i] = MPI_COMM_NULL;

  /* Make a list of the processes in the new communicator. */
  pranks = (int *) SUPERLU_MALLOC( nprocs * sizeof(int) );
  
  i = 2 * nprocs - 2;
  MPI_Comm_dup ((*symb_comm), &(commLvls[i]));
  szSep = 2;
  i -= szSep;

  while (i > 0) {
    /* for each level in the separator tree */
    npNode = nprocs / szSep; 
    fstP = 0; 
    /* for each node in the level */
    for (j = i; j < i + szSep; j++) {
      lstP = fstP + npNode;
      if (fstP <= iam && iam < lstP) {
	ind = j;
	key = iam - fstP;
	col = fstP;
      }
      fstP += npNode;
    }
    MPI_Comm_split ((*symb_comm), col, key, &(commLvls[ind]) );
    
    szSep *= 2;
    i -= szSep;
  }
  
  SUPERLU_FREE (pranks);
}

static void
intraLvl_symbfact 
(
 SuperMatrix *A, /* Input - original matrix A  */
 int   iam,      /* Input - my processor number */
 int   lvl,      /* Input - current level in the separator tree */
 int   szSep,    /* Input - size of the current separator(node) */
 int   ind_sizes1,
 int   ind_sizes2, 
 int_t *sizes,     /* Input - sizes of each node in the separator tree */
 int_t *fstVtxSep, /* Input - first vertex of each node in the tree */
 int   fstP,     /* Input - first processor assigned to current node */
 int   lstP,     /* Input - last processor assigned to current node */
 int_t fstVtx,   /* Input - first vertex of current node */
 int_t lstVtx,   /* Input - last vertex of current node */
 Pslu_freeable_t *Pslu_freeable,   /* global LU data structures (modified) */
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo, /* Input/Output - local info on vertices distribution */
 comm_symbfact_t *CS,
 psymbfact_stat_t *PS,
 int_t *marker,
 int_t *p_mark,    /* marker used to merge elements of vertices */
 int_t *p_nextl,   /* ptr to nextl in lsub structure */
 int_t *p_nextu,   /* ptr to nextu in usub structure */
 int_t *p_neltsZr, /* no of artificial zeros introduced so far */
 int_t *p_neltsTotal, /* no of nonzeros (including artificials) 
			 computed so far */
 int_t *p_nsuper_loc,
 MPI_Comm ndComm,
 MPI_Comm    *symb_comm /* Input - communicator for symbolic factorization */
 )
{
  int nprocsLvl, p, prvP, rcvP;
  int toSend, rcvd_prvP, index_req[2];
  int_t fstVtx_loc_lid, fstVtx_loc, vtx, vtxLvl, curblk_loc, denseSep;
  int_t fstVtx_blk, fstVtx_blk_lid, lstVtx_blk, lstVtx_blk_lid, tag;
  int_t nvtcs_blk, xusub_end, xlsub_end, prv_fstVtx_blk;
  int_t n;
  int_t *rcv_intraLvl, *snd_intraLvl;
  int_t *lsub_rcvd, lsub_rcvd_sz, *usub_rcvd, usub_rcvd_sz;
  int_t nmsgsRcvd, nmsgsTRcv, sz_msg;
  int_t nvtcsLvl_loc, nextl, nextu, ind_blk, snd_vtxLvl, maxNeltsVtx_in;
  int_t count, vtx_loc, mem_error, lstBlkRcvd;
  int_t fstVtx_blk_loc, fstBlk, vtx_lid, prElt, nelts, j, nvtcs_toUpd;
  int_t snd_LinterLvlSz, fstVtx_blk_loc_lid, prElt_ind, maxNmsgsToRcv;
  int_t *xlsub, *xusub, *lsub, *usub;
  int_t *globToLoc, maxNvtcsPProc, nblk_loc, upd_myD, r, fstVtx_blkCyc;
  int_t k, prElt_L, prElt_U, vtx_elt, fstVtx_toUpd;
  int intSzMsg;

  MPI_Status status[4];
  MPI_Request request[4];
  
  /* Initializations */
  lsub    = Llu_symbfact->lsub;   xlsub    = Llu_symbfact->xlsub;
  usub    = Llu_symbfact->usub;   xusub    = Llu_symbfact->xusub;
  
  /* max number of msgs this processor can receive during 
     intraLvl_symbfact routine */
  maxNmsgsToRcv  = (lstVtx - fstVtx) / VInfo->maxSzBlk + 1;
  maxNeltsVtx_in = VInfo->maxNeltsVtx;
  globToLoc      = Pslu_freeable->globToLoc;
  maxNvtcsPProc  = Pslu_freeable->maxNvtcsPProc;
  n = A->ncol;
  nprocsLvl       = lstP - fstP;
  rcv_intraLvl    = CS->rcv_intraLvl;
  snd_intraLvl    = CS->snd_intraLvl;
  nvtcsLvl_loc    = VInfo->nvtcsLvl_loc;
  nmsgsTRcv       = 0;
  nmsgsRcvd       = 0;
  nblk_loc        = 0;
  nvtcs_toUpd     = nvtcsLvl_loc;
  fstVtx_blk      = fstVtx;
  denseSep        = FALSE;

  /* determine first vertex that belongs to fstP */
  k = fstVtx;
  fstVtx_blkCyc = n;
  while (k < lstVtx && fstVtx_blkCyc == n) {
    p = OWNER( globToLoc[k] );
    if (p == fstP)
      fstVtx_blkCyc = k;
    k += VInfo->maxSzBlk;
  }

  for (p = fstP; p < lstP; p++)
    rcv_intraLvl[p] = 0;

  for (r = 0; r < 3; r++) 
    request[r] = MPI_REQUEST_NULL;

  fstVtx_loc = VInfo->begEndBlks_loc[VInfo->curblk_loc];
  fstVtx_loc_lid = LOCAL_IND( globToLoc[fstVtx_loc] ); 
  vtx = fstVtx_loc;
  if (fstVtx_loc >= fstVtx_blkCyc)
    nblk_loc = 1;
  while (VInfo->begEndBlks_loc[VInfo->curblk_loc] < lstVtx && !VInfo->filledSep) {
    CS->snd_intraSz  = 0;
    CS->snd_LintraSz = 0;

    lstBlkRcvd     = FALSE;
    prv_fstVtx_blk = fstVtx_blk;
    fstVtx_blk     = VInfo->begEndBlks_loc[VInfo->curblk_loc];
    lstVtx_blk     = VInfo->begEndBlks_loc[VInfo->curblk_loc + 1];
    fstVtx_toUpd   = VInfo->begEndBlks_loc[VInfo->curblk_loc + 2];
    fstVtx_blk_lid = LOCAL_IND( globToLoc[fstVtx_blk] );
    lstVtx_blk_lid = LOCAL_IND( globToLoc[lstVtx_blk - 1] + 1);
    nvtcs_blk      = lstVtx_blk - fstVtx_blk;
    nvtcs_toUpd   -= nvtcs_blk;
    nmsgsTRcv      = n;
    VInfo->maxNeltsVtx -= fstVtx_blk - prv_fstVtx_blk;

    index_req[0] = EMPTY;
    for (r = 0; r < 3; r++) 
      request[r] = MPI_REQUEST_NULL;
    if (fstVtx_blk != fstVtx) {
      /* if not the first vertex of the level */
      prvP           = OWNER( globToLoc[fstVtx_blk - 1] );
      rcvd_prvP      = FALSE;
      /* receive info on number messages to receive */
      tag = tag_intraLvl + nblk_loc;
      if (iam == fstP)  tag --;
      
      MPI_Irecv (&(rcv_intraLvl[fstP]), nprocsLvl, mpi_int_t, prvP,
		 tag, (*symb_comm), &(request[1]));

      while (!rcvd_prvP || nmsgsRcvd < nmsgsTRcv) {
	if (index_req[0] != 1) {
	  MPI_Irecv (&sz_msg, 1, mpi_int_t, 
		     MPI_ANY_SOURCE, tag_intraLvl_szMsg, 
		     (*symb_comm), &(request[0]));  
	  if (sz_msg > INT_MAX)
	    ABORT("ERROR in intraLvl_symbfact size to send > INT_MAX\n");
	}
	MPI_Waitany (2, request, index_req, status);
	if (index_req[0] == 1) {
	  /* receive information on no msgs to receive */
#if ( PRNTlevel>=1 )
	  PS->no_shmRcvd ++;
#endif
	  rcvd_prvP = TRUE;
	  nmsgsTRcv = rcv_intraLvl[iam];
	  /* if dense separator was detected by one of the 
	     previous processors ... */
	  if (nmsgsTRcv > maxNmsgsToRcv) {
	    VInfo->filledSep = (int) nmsgsTRcv / maxNmsgsToRcv;
	    nmsgsTRcv = nmsgsTRcv % maxNmsgsToRcv;
	  }
	  
	  if (nmsgsTRcv == nmsgsRcvd) {
	    /* MPI_Cancel (&(request[0])); */
	    MPI_Send (&r, 1, mpi_int_t, iam, 
		      tag_intraLvl_szMsg, (*symb_comm));
	    MPI_Wait (&(request[0]), status);	    
	  }
	}
	if (index_req[0] == 0) {
	  nmsgsRcvd ++;
	  if (nmsgsTRcv == nmsgsRcvd)  lstBlkRcvd = TRUE; 
	  rcvP = status->MPI_SOURCE;

	  /* allocate enough space to receive data */
	  if (CS->rcv_bufSz < sz_msg) {
	    PS->maxSzBuf += sz_msg - CS->rcv_bufSz;
	    if (CS->rcv_bufSz != 0)
	      /* not first time allocate memory */
	      SUPERLU_FREE (CS->rcv_buf);
	    CS->rcv_bufSz = sz_msg;
	    if (!(CS->rcv_buf = intMalloc_symbfact (sz_msg))) {
	      ABORT("Malloc fails for rcv_buf[].");
	    }
	  }
	  
	  /* use received data to update symbolic factorization */
	  lsub_rcvd = CS->rcv_buf;
	  MPI_Recv (lsub_rcvd, sz_msg, mpi_int_t, 
		    rcvP, tag_intraLvl_LData, (*symb_comm), status);
	  MPI_Get_count (status, mpi_int_t, &intSzMsg);
	  lsub_rcvd_sz = intSzMsg;
	  usub_rcvd    = &(CS->rcv_buf[lsub_rcvd_sz]);
	  MPI_Recv (usub_rcvd, sz_msg - lsub_rcvd_sz, 
		    mpi_int_t, rcvP,
		    tag_intraLvl_UData, (*symb_comm), status);
	  MPI_Get_count (status, mpi_int_t, &intSzMsg);
	  usub_rcvd_sz = intSzMsg;
#if ( PRNTlevel>=1 )
	  PS->no_shmRcvd ++;
	  PS->no_msgsRcvd += (float) 2;
	  PS->sz_msgsRcvd += (float) sz_msg;
	  if (PS->maxsz_msgRcvd < lsub_rcvd_sz) PS->maxsz_msgRcvd = lsub_rcvd_sz;
	  if (PS->maxsz_msgRcvd < usub_rcvd_sz) PS->maxsz_msgRcvd = usub_rcvd_sz;
#endif

	  if (!lstBlkRcvd) {
	    mem_error = 
	      rl_update (1, n, iam, lsub_rcvd, lsub_rcvd_sz,
			 usub_rcvd, usub_rcvd_sz, EMPTY, EMPTY, EMPTY,
			 fstVtx_blk, lstVtx, nvtcs_blk + nvtcs_toUpd,
			 1, p_mark,
			 marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
	    lsub = Llu_symbfact->lsub;
	    mem_error = 
	      rl_update (1, n, iam, usub_rcvd, usub_rcvd_sz,
			 lsub_rcvd, lsub_rcvd_sz, EMPTY, EMPTY, EMPTY, 
			 fstVtx_blk, lstVtx, nvtcs_blk + nvtcs_toUpd,
			 0, p_mark,
			 marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
	    usub = Llu_symbfact->usub;
	  }
	}
      }
    }
  
    if (VInfo->filledSep) {
      mem_error = 
	denseSep_symbfact (1, n, iam, ind_sizes1, ind_sizes2, sizes, fstVtxSep,
			   szSep, fstP, lstP, fstVtx_blkCyc, nblk_loc,
			   p_nextl, p_nextu, p_mark, p_nsuper_loc, marker,
			   ndComm, symb_comm, Llu_symbfact, Pslu_freeable, VInfo, CS, PS);
    }
    else {
      /* compute symbolic factorization for this block */
      if (!lstBlkRcvd) {
	lsub_rcvd = NULL; usub_rcvd = NULL;
      }

      blk_symbfact (A, iam, lvl, 
		    szSep, ind_sizes1, ind_sizes2, sizes, fstVtxSep,
		    fstVtx_loc, fstVtx_blk, lstVtx_blk, 
		    lsub_rcvd, lsub_rcvd_sz, usub_rcvd, usub_rcvd_sz,
		    Pslu_freeable, Llu_symbfact, VInfo, CS, PS,
		    marker, p_mark,
		    p_nextl, p_nextu, p_neltsZr, p_neltsTotal, 
		    p_nsuper_loc);
      lsub = Llu_symbfact->lsub;
      usub = Llu_symbfact->usub; 	 
      
      if (lstVtx_blk != lstVtx) {
	/* if this is not the last block of the level */
	if (VInfo->filledSep == FILLED_SEPS ||
	    ( VInfo->filledSep == FILLED_SEP && 
	      ((lstVtx - lstVtx_blk > VInfo->maxSzBlk * nprocsLvl && nblk_loc > 0) ||
	       (lstVtx - fstVtx_blkCyc > VInfo->maxSzBlk * nprocsLvl && nblk_loc == 0))))
	  /* if current separator is dense and this is not the last block, 
	     then ... */
	  denseSep = TRUE;
	else
	  /* separator dense but not enough uncomputed blocks 
	     in the separator to take advantage of it */
	  VInfo->filledSep = FALSE;
	
	if (VInfo->filledSep == FILLED_SEPS) {
	  for (p = fstP; p < lstP; p++)
	    rcv_intraLvl[p] = maxNmsgsToRcv * VInfo->filledSep + rcv_intraLvl[p];
	  denseSep_symbfact (0, n, iam, ind_sizes1, ind_sizes2, sizes, fstVtxSep,
			     szSep, fstP, lstP, fstVtx_blkCyc, nblk_loc,
			     p_nextl, p_nextu, p_mark, p_nsuper_loc, marker, ndComm, 
			     symb_comm, Llu_symbfact, Pslu_freeable, VInfo, CS, PS);
	}
	else {
	  /* send blk to next procs and update the rest of my own blocks */
	  if (lstBlkRcvd) {
	    mem_error = 
	      rl_update (1, n, iam, lsub_rcvd, lsub_rcvd_sz,
			 usub_rcvd, usub_rcvd_sz, EMPTY, EMPTY, EMPTY,
			 fstVtx_toUpd, lstVtx, nvtcs_toUpd,
			 1, p_mark,
			 marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
	    lsub = Llu_symbfact->lsub;
	    mem_error = 
	      rl_update (1, n, iam, usub_rcvd, usub_rcvd_sz,
			 lsub_rcvd, lsub_rcvd_sz, EMPTY, EMPTY, EMPTY, 
			 fstVtx_toUpd, lstVtx, nvtcs_toUpd,
			 0, p_mark,
			 marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
	    usub = Llu_symbfact->usub;
	  }

	  upd_myD = FALSE;
	  /* determine processors to which send this block
	     and copy data to be sent */
	  for (p = fstP; p < lstP; p++)
	    snd_intraLvl[p] = FALSE;
	  nextl = 0; 
	  nextu = nextl + CS->snd_LintraSz;
	  
	  /* allocate enough space to receive data */
	  if (CS->rcv_bufSz < CS->snd_intraSz) {
	    PS->maxSzBuf += CS->snd_intraSz - CS->rcv_bufSz;
	    if (CS->rcv_bufSz != 0)
	      /* not first time allocate memory */
	      SUPERLU_FREE (CS->rcv_buf);
	    CS->rcv_bufSz = CS->snd_intraSz;
	    if (!(CS->rcv_buf = intMalloc_symbfact (CS->snd_intraSz))) {
	      ABORT("Malloc fails for rcv_buf[].");
	    }
	  }

	  for (vtx = fstVtx_blk, vtx_lid = fstVtx_blk_lid; 
	       vtx < lstVtx_blk; vtx++, vtx_lid ++) {
	    toSend = FALSE;
	    k = xlsub[vtx_lid];
	    prElt_L = lsub[k];
	    j = xusub[vtx_lid];
	    prElt_U = usub[j];

	    if (prElt_L >= lstVtx_blk || prElt_U >= lstVtx_blk) {
	      if (vtx == lstVtx_blk - 1) {
		xlsub_end = *p_nextl;
		xusub_end = *p_nextu;
	      }
	      else {
		xlsub_end = xlsub[vtx_lid + 1];
		xusub_end = xusub[vtx_lid + 1];
	      }
	      if (prElt_L >= lstVtx_blk) {
		while (lsub[k] <= prElt_L && k < xlsub_end) {
		  vtx_elt = lsub[k];
		  if (vtx_elt >= lstVtx_blk && vtx_elt < lstVtx) {
		    p = OWNER( globToLoc[vtx_elt] );
		    if (p != iam) {
		      /* vtx will be send to another processor */
		      snd_intraLvl[p] = TRUE;
		      toSend = TRUE;
		    }
		    else {
		      upd_myD = TRUE;
		    }
		  }
		  k++;
		}
	      }
	      if (prElt_U >= lstVtx_blk) {
		while (usub[j] <= prElt_U && j < xusub_end) {
		  vtx_elt = usub[j];
		  if (vtx_elt >= lstVtx_blk && vtx_elt < lstVtx) {
		    p = OWNER( globToLoc[vtx_elt] );
		    if (p != iam) {
		      /* vtx will be send to another processor */
		      snd_intraLvl[p] = TRUE;
		      toSend = TRUE;
		    }
		    else {
		      upd_myD = TRUE;
		    }
		  }
		  j ++;
		}
	      }
	      if (toSend) {
		/* L(:, vtx) and U(vtx, :) will be send to processors */
		nelts = xusub_end - xusub[vtx_lid];
		CS->rcv_buf[nextu + DIAG_IND]  = vtx;
		CS->rcv_buf[nextu + NELTS_IND] = nelts;
		nextu += 2;
		for (j = xusub[vtx_lid]; j < xusub_end; j++) {
		  CS->rcv_buf[nextu] = usub[j]; nextu ++;
		}
		
		nelts = xlsub_end - xlsub[vtx_lid];
		CS->rcv_buf[nextl + DIAG_IND] = vtx;
		CS->rcv_buf[nextl + NELTS_IND] = nelts; 
		nextl += 2;
		for (j = xlsub[vtx_lid]; j < xlsub_end; j++) {
		  CS->rcv_buf[nextl] = lsub[j]; nextl ++;
		}
	      }
	    }
	  }
	  for (p = fstP; p < lstP; p++) 
	    if (snd_intraLvl[p])
	      rcv_intraLvl[p] ++;

	  if (VInfo->filledSep == FILLED_SEP) {
	    for (p = fstP; p < lstP; p++)
	      rcv_intraLvl[p] = maxNmsgsToRcv * VInfo->filledSep + 
		rcv_intraLvl[p];
	  }
	  else {
	    /* send to the owner of the next block info on no of messages */
	    p = OWNER( globToLoc[lstVtx_blk] );
	    tag = tag_intraLvl + nblk_loc;
	    
	    MPI_Isend (&(rcv_intraLvl[fstP]), nprocsLvl, mpi_int_t, p,
		       tag, (*symb_comm), request);
#if ( PRNTlevel>=1 )
	    PS->no_shmSnd ++;
#endif
	  }

	  /* there is data to be send */
	  sz_msg = nextl + nextu - CS->snd_LintraSz;
	  for (p = fstP; p < lstP; p++) {
	    if (p != iam && snd_intraLvl[p]) {
	      MPI_Isend (&sz_msg, 1, mpi_int_t, p,
			 tag_intraLvl_szMsg, (*symb_comm), &(request[1]));
	      MPI_Isend (CS->rcv_buf, nextl, mpi_int_t, p,
			 tag_intraLvl_LData, (*symb_comm), &(request[2]));
	      MPI_Isend (&(CS->rcv_buf[CS->snd_LintraSz]), 
			 nextu - CS->snd_LintraSz, mpi_int_t, p,
			 tag_intraLvl_UData, (*symb_comm), &(request[3]));
	      MPI_Waitall(3, &(request[1]), &(status[1]));
#if ( PRNTlevel>=1 )
	      PS->no_shmSnd ++;
	      PS->no_msgsSnd += (float) 2;
	      PS->sz_msgsSnd += (float) sz_msg;
	      if (PS->maxsz_msgSnd < nextl) PS->maxsz_msgSnd = nextl;
	      if (PS->maxsz_msgSnd < nextu - CS->snd_LintraSz) 
		PS->maxsz_msgSnd = nextu - CS->snd_LintraSz;
#endif
	    }
	  }
	  if (VInfo->filledSep != FILLED_SEP) {
	    MPI_Wait (request, status);      
	  }

	  /* update rest of vertices */
	  if (upd_myD) {
	    lsub_rcvd_sz = (*p_nextl) - xlsub[fstVtx_blk_lid];
	    lsub_rcvd    = &(lsub[xlsub[fstVtx_blk_lid]]);
	    usub_rcvd_sz = (*p_nextu) - xusub[fstVtx_blk_lid];
	    usub_rcvd    = &(usub[xusub[fstVtx_blk_lid]]);
	    
	    mem_error =
	      rl_update (0, n, iam, lsub_rcvd, lsub_rcvd_sz,
			 usub_rcvd, usub_rcvd_sz, fstVtx_blk, lstVtx_blk,
			 EMPTY,
			 fstVtx_toUpd, lstVtx, nvtcs_toUpd,
			 1, p_mark,
			 marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
	    lsub = Llu_symbfact->lsub;
	    lsub_rcvd    = &(lsub[xlsub[fstVtx_blk_lid]]);
	    mem_error =
	      rl_update (0, n, iam, usub_rcvd, usub_rcvd_sz,
			 lsub_rcvd, lsub_rcvd_sz, fstVtx_blk, lstVtx_blk,
			 EMPTY,
			 fstVtx_toUpd, lstVtx, nvtcs_toUpd,
			 0, p_mark,
			 marker, Pslu_freeable, Llu_symbfact, VInfo, PS);
	    usub = Llu_symbfact->usub;
	  }
	  if (VInfo->filledSep == FILLED_SEP)
	    denseSep_symbfact (0, n, iam, ind_sizes1, ind_sizes2, sizes, fstVtxSep,
			       szSep, fstP, lstP, fstVtx_blkCyc, nblk_loc,
			       p_nextl, p_nextu, p_mark, p_nsuper_loc, marker, ndComm, 
			       symb_comm, Llu_symbfact, Pslu_freeable, VInfo, CS, PS);
	}
      }
    }
    VInfo->curblk_loc += 2;
    nblk_loc ++;
  }
  
  /* update maxNeltsVtx */
  VInfo->maxNeltsVtx = maxNeltsVtx_in - lstVtx + fstVtx;
  
  /* if current separator dense, then reset value of filledSep */
  if (VInfo->filledSep == FILLED_SEP)
    VInfo->filledSep = FALSE;
}

static void
symbfact_free 
(
 int   iam,    /* Input - my processor number */
 int   nprocs, /* Input - number of processors for the symbolic factorization */
 Llu_symbfact_t *Llu_symbfact,  /* Input/Output - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo, /* Input/Output - local info on vertices distribution */
 comm_symbfact_t *CS
 )
{
  /* free memory corresponding to prune structure */
  if (Llu_symbfact->szLsubPr != 0)
    SUPERLU_FREE( Llu_symbfact->lsubPr );
  if (Llu_symbfact->szUsubPr != 0)
    SUPERLU_FREE( Llu_symbfact->usubPr );
  if (Llu_symbfact->xlsubPr != NULL)
    SUPERLU_FREE( Llu_symbfact->xlsubPr );
  if (Llu_symbfact->xusubPr != NULL)
    SUPERLU_FREE( Llu_symbfact->xusubPr );
  
  if (Llu_symbfact->xlsub_rcvd != NULL)
    SUPERLU_FREE( Llu_symbfact->xlsub_rcvd);
  if (Llu_symbfact->xusub_rcvd != NULL)
    SUPERLU_FREE( Llu_symbfact->xusub_rcvd);
  
  if (Llu_symbfact->cntelt_vtcs != NULL)
    SUPERLU_FREE( Llu_symbfact->cntelt_vtcs);
  if (Llu_symbfact->cntelt_vtcsA_lvl != NULL)
    SUPERLU_FREE( Llu_symbfact->cntelt_vtcsA_lvl);
  
  if (CS->rcv_bufSz != 0)
    SUPERLU_FREE( CS->rcv_buf );
  if (CS->snd_bufSz != 0)
    SUPERLU_FREE( CS->snd_buf );
  
  SUPERLU_FREE( VInfo->begEndBlks_loc);
  SUPERLU_FREE( CS->rcv_interLvl);
  SUPERLU_FREE( CS->snd_interLvl);
  SUPERLU_FREE( CS->ptr_rcvBuf);
  SUPERLU_FREE( CS->rcv_intraLvl);
  SUPERLU_FREE( CS->snd_intraLvl);
  SUPERLU_FREE( CS->snd_interSz);
  SUPERLU_FREE( CS->snd_LinterSz);
  SUPERLU_FREE( CS->snd_vtxinter);  
}

static void
estimate_memUsage
(
 int_t n,  /* Input - order of the matrix */
 int iam,  /* Input - my processor number */
 mem_usage_t *symb_mem_usage,
 float *p_totalMemLU,   /* Output -memory used for symbolic factorization */
 float *p_overestimMem, /* Output -memory allocated during to right looking 
			   overestimation memory usage */
 Pslu_freeable_t *Pslu_freeable,   /* global LU data structures (modified) */
 Llu_symbfact_t *Llu_symbfact,  /* Input - local L, U data structures */
 vtcsInfo_symbfact_t *VInfo, /* Input - local info on vertices distribution */
 comm_symbfact_t *CS,
 psymbfact_stat_t *PS
 )
{
  int_t nvtcs_loc, lword, nsuper_loc;
  float lu_mem, other_mem, overestimMem;
  
  nvtcs_loc = VInfo->nvtcs_loc;
  nsuper_loc = Pslu_freeable->supno_loc[nvtcs_loc];
  lword     = sizeof(int_t);
  
  /* memory for xlsub, xusub, supno_loc, cntelt_vtcs */
  lu_mem = 4.0 * (float) nvtcs_loc * (float) lword;
  /* memory for xlsubPr, xusubPr */
  lu_mem += 2.0 * (float) VInfo->maxNvtcsNds_loc * (float) lword;
  
  if (PS->estimLSz < Llu_symbfact->xlsub[nvtcs_loc])
    PS->estimLSz = Llu_symbfact->xlsub[nvtcs_loc];
  if (PS->estimUSz < Llu_symbfact->xusub[nvtcs_loc])
    PS->estimUSz = Llu_symbfact->xusub[nvtcs_loc];
  
  lu_mem += (float) PS->estimLSz * lword;
  lu_mem += (float) PS->estimUSz * lword;
  lu_mem += (float) PS->maxSzLPr * lword;
  lu_mem += (float) PS->maxSzUPr * lword;
  lu_mem += (float) PS->szDnsSep * lword;
  /* memory for globToLoc, tempArray */
  lu_mem += (float) 2* (float) n * lword;
  lu_mem += (float) PS->maxSzBuf * lword;
  
  overestimMem  = (float) (PS->estimLSz - Llu_symbfact->xlsub[nvtcs_loc]) * lword;
  overestimMem += (float) (PS->estimUSz - Llu_symbfact->xusub[nvtcs_loc]) * lword;
  
  *p_totalMemLU = lu_mem;  
  *p_overestimMem = overestimMem;
  
  symb_mem_usage->for_lu = (float) ((3 * nvtcs_loc + 2 * nsuper_loc) * lword);
  symb_mem_usage->for_lu += (float) (Llu_symbfact->xlsub[nvtcs_loc] * lword); 
  symb_mem_usage->for_lu += (float) (Llu_symbfact->xusub[nvtcs_loc] * lword);   
  symb_mem_usage->total = lu_mem;
}


static int_t *
intMalloc_symbfact(int_t n)
{
  int_t *buf;
  if (n == 0)
    buf = NULL;
  else
    buf = (int_t *) SUPERLU_MALLOC(n * sizeof(int_t));
  return buf;
}

static int_t *
intCalloc_symbfact(int_t n)
{
  int_t *buf;
  register int_t i;

  if (n == 0)
    buf = NULL;
  else
    buf = (int_t *) SUPERLU_MALLOC(n * sizeof(int_t));
  if ( buf )
    for (i = 0; i < n; i++) buf[i] = 0;
  return (buf);
}

