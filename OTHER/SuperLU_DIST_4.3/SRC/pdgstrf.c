

/*! @file
 * \brief Performs LU factorization in parallel
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.3) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 *
 * Modified:
 *     September 1, 1999
 *     Feburary 7, 2001  use MPI_Isend/MPI_Irecv
 *     October 15, 2008  latency-reducing panel factorization
 *     July    12, 2011  static scheduling and arbitrary look-ahead
 *     March   13, 2013  change NTAGS to MPI_TAG_UB value
 *     September 24, 2015 replace xLAMCH by xMACH, using C99 standard.
 *     December 31, 2015 rename xMACH to xMACH_DIST
 *
 * Sketch of the algorithm 
 *
 * ======================= 
 *    
 * The following relations hold:
 *     * A_kk = L_kk * U_kk
 *     * L_ik = Aik * U_kk^(-1)
 *     * U_kj = L_kk^(-1) * A_kj
 *
 *              ----------------------------------
 *              |   |                            |
 *              ----|-----------------------------
 *              |   | \ U_kk|                    |
 *              |   |   \   |        U_kj        |
 *              |   |L_kk \ |         ||         |
 *              ----|-------|---------||----------
 *              |   |       |         \/         |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   | L_ik ==>       A_ij        |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   |       |                    |
 *              ----------------------------------
 *
 * Handle the first block of columns separately.
 *     * Factor diagonal and subdiagonal blocks and test for exact
 *       singularity. ( pdgstrf2(0), one column at a time )
 *     * Compute block row of U
 *     * Update trailing matrix
 *
 * Loop over the remaining blocks of columns.
 *   mycol = MYCOL( iam, grid );
 *   myrow = MYROW( iam, grid );
 *   N = nsupers;
 *   For (k = 1; k < N; ++k) {
 *       krow = PROW( k, grid );
 *       kcol = PCOL( k, grid );
 *       Pkk = PNUM( krow, kcol, grid );
 *
 *     * Factor diagonal and subdiagonal blocks and test for exact
 *       singularity.
 *       if ( mycol == kcol ) {
 *           pdgstrf2(k), one column at a time
 *       }
 *
 *     * Parallel triangular solve
 *       if ( iam == Pkk ) multicast L_k,k to this process row;
 *       if ( myrow == krow && mycol != kcol ) {
 *          Recv L_k,k from process Pkk;
 *          for (j = k+1; j < N; ++j)
 *              if ( PCOL( j, grid ) == mycol && A_k,j != 0 )
 *                 U_k,j = L_k,k \ A_k,j;
 *       }
 *
 *     * Parallel rank-k update
 *       if ( myrow == krow ) multicast U_k,k+1:N to this process column;
 *       if ( mycol == kcol ) multicast L_k+1:N,k to this process row;
 *       if ( myrow != krow ) {
 *          Pkj = PNUM( krow, mycol, grid );
 *          Recv U_k,k+1:N from process Pkj;
 *       }
 *       if ( mycol != kcol ) {
 *          Pik = PNUM( myrow, kcol, grid );
 *          Recv L_k+1:N,k from process Pik;
 *       }
 *       for (j = k+1; k < N; ++k) {
 *          for (i = k+1; i < N; ++i)
 *              if ( myrow == PROW( i, grid ) && mycol == PCOL( j, grid )
 *                   && L_i,k != 0 && U_k,j != 0 )
 *                 A_i,j = A_i,j - L_i,k * U_k,j;
 *       }
 *  }
 *
 * </pre>
 */

#include <math.h>
/*#include "mkl.h"*/
#include "superlu_ddefs.h"

#ifdef GPU_ACC
#include "cublas_utils.h"
/*#include "cublas_dgemm.h"*/
// #define NUM_CUDA_STREAMS 16
// #define NUM_CUDA_STREAMS 16
#endif 

/* Various defininations     */
/* 
    Name    : SUPERNODE_PROFILE  
    Purpose : For SuperNode Level profiling of various measurements such as gigaflop/sec
    obtained,bandwidth achived:
    Overhead : Low 
*/
// #define SUPERNODE_PROFILE   

/* 
    Name    :   BAELINE
    Purpose : baseline to compare performance against
    Overhead : NA : this wont be used for running experiments
*/
// #define BASELINE

/* 
    Name    :   PHI_FRAMEWORK
    Purpose : To simulate and test algorithm used for offloading Phi
    Overhead : NA : this wont be used for running experiments
*/
#define PHI_FRAMEWORK

#define PDGSTRF2 pdgstrf2_trsm
#define PDGSTRS2 pdgstrs2_omp

extern void PDGSTRF2 (superlu_options_t *, int_t, int_t, double,
                        Glu_persist_t *, gridinfo_t *, LocalLU_t *,
                        MPI_Request *, int, SuperLUStat_t *, int *);
#ifdef _CRAY
extern void PDGSTRS2 (int_t, int_t, Glu_persist_t *, gridinfo_t *,
                      LocalLU_t *, SuperLUStat_t *, _fcd, _fcd, _fcd);
#else
extern void PDGSTRS2 (int_t, int_t, Glu_persist_t *, gridinfo_t *,
                      LocalLU_t *, SuperLUStat_t *);
#endif

#define ISORT                   /* Note: qsort() has bug on Mac */

#ifdef ISORT
extern void isort (int_t N, int_t * ARRAY1, int_t * ARRAY2);
extern void isort1 (int_t N, int_t * ARRAY);

#else

int
superlu_sort_perm (const void *arg1, const void *arg2)
{
    const int_t *val1 = (const int_t *) arg1;
    const int_t *val2 = (const int_t *) arg2;
    return (*val2 < *val1);
}
#endif


int get_thread_per_process()
{   
    char* ttemp; 
    ttemp = getenv("THREAD_PER_PROCESS");

    if(ttemp) return atoi(ttemp);
    else return 1;
}

int
get_mic_offload ()
{
    char *ttemp;
    ttemp = getenv ("SUPERLU_MIC_OFFLOAD");

    if (ttemp)
        return atoi (ttemp);
    else
        return 0;
}

int_t
get_max_buffer_size ()
{
    char *ttemp;
    ttemp = getenv ("MAX_BUFFER_SIZE");
    if (ttemp)
        return atoi (ttemp);
    else
        return 5000000;
}

int_t
get_cublas_nb ()
{
    char *ttemp;
    ttemp = getenv ("CUBLAS_NB");
    if (ttemp)
        return atoi (ttemp);
    else
        return 64;
}

int_t
get_num_cuda_streams ()
{
    char *ttemp;
    ttemp = getenv ("NUM_CUDA_STREAMS");
    if (ttemp)
        return atoi (ttemp);
    else
        return 8;
}

/*int omp_get_num_threads (void);
  int omp_get_thread_num (void);*/

int AssignMic(int my_rank)
{
    return (my_rank+1)%2;
}

/************************************************************************/

#include "dscatter.c"

/************************************************************************/


/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * PDGSTRF performs the LU factorization in parallel.
 *
 * Arguments
 * =========
 *
 * options (input) superlu_options_t*
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *         The following field should be defined:
 *         o ReplaceTinyPivot (yes_no_t)
 *           Specifies whether to replace the tiny diagonals by
 *           sqrt(epsilon)*norm(A) during LU factorization.
 *
 * m      (input) int
 *        Number of rows in the matrix.
 *
 * n      (input) int
 *        Number of columns in the matrix.
 *
 * anorm  (input) double
 *        The norm of the original matrix A, or the scaled A if
 *        equilibration was done.
 *
 * LUstruct (input/output) LUstruct_t*
 *         The data structures to store the distributed L and U factors.
 *         The following fields should be defined:
 *
 *         o Glu_persist (input) Glu_persist_t*
 *           Global data structure (xsup, supno) replicated on all processes,
 *           describing the supernode partition in the factored matrices
 *           L and U:
 *         xsup[s] is the leading column of the s-th supernode,
 *             supno[i] is the supernode number to which column i belongs.
 *
 *         o Llu (input/output) LocalLU_t*
 *           The distributed data structures to store L and U factors.
 *           See superlu_ddefs.h for the definition of 'LocalLU_t'.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh. It contains the MPI communicator, the number
 *        of process rows (NPROW), the number of process columns (NPCOL),
 *        and my process rank. It is an input argument to all the
 *        parallel routines.
 *        Grid can be initialized by subroutine SUPERLU_GRIDINIT.
 *        See superlu_ddefs.h for the definition of 'gridinfo_t'.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics on runtime and floating-point operation count.
 *        See util.h for the definition of 'SuperLUStat_t'.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 * </pre>
 */
int_t
pdgstrf(superlu_options_t * options, int m, int n, double anorm,
       LUstruct_t * LUstruct, gridinfo_t * grid, SuperLUStat_t * stat, int *info)
{
#ifdef _CRAY
    _fcd ftcs = _cptofcd ("N", strlen ("N"));
    _fcd ftcs1 = _cptofcd ("L", strlen ("L"));
    _fcd ftcs2 = _cptofcd ("N", strlen ("N"));
    _fcd ftcs3 = _cptofcd ("U", strlen ("U"));
#endif
    double zero = 0.0, alpha = 1.0, beta = 0.0;
    int_t *xsup;
    int_t *lsub, *lsub1, *usub, *Usub_buf;
    int_t **Lsub_buf_2, **Usub_buf_2;
    double **Lval_buf_2, **Uval_buf_2;          /* pointers to starts of bufs */
    double *lusup, *lusup1, *uval, *Uval_buf;   /* pointer to current buf     */
    int_t fnz, i, ib, ijb, ilst, it, iukp, jb, jj, klst, knsupc,
        lb, lib, ldv, ljb, lptr, lptr0, lptrj, luptr, luptr0, luptrj,
        nlb, nub, nsupc, rel, rukp, il, iu;
    int_t Pc, Pr;
    int iam, kcol, krow, yourcol, mycol, myrow, pi, pj;
    int j, k, lk, nsupers;  /* k - current panel to work on */
    int k0;        /* counter of the next supernode to be factored */
    int kk, kk0, kk1, kk2, jj0; /* panels in the look-ahead window */
    int iukp0, rukp0, flag0, flag1;
    int nsupr, nbrow, segsize;
    int msg0, msg2;
    int_t **Ufstnz_br_ptr, **Lrowind_bc_ptr;
    double **Unzval_br_ptr, **Lnzval_bc_ptr;
    int_t *index;
    double *nzval;
    int_t *iuip, *ruip; /* Pointers to U index/nzval; size ceil(NSUPERS/Pr). */
    double *ucol;
    int *indirect, *indirect2;
    double *tempv, *tempv2d;
    int iinfo;
    int *ToRecv, *ToSendD, **ToSendR;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;
    superlu_scope_t *scp;
    float s_eps;
    double thresh;
    double *tempU2d, *tempu;
    int full, ldt, ldu, lead_zero, ncols, ncb, nrb, p, pr, pc, nblocks;
    int_t *etree_supno_l, *etree_supno, *blocks, *blockr, *Ublock, *Urows,
        *Lblock, *Lrows, *perm_u, *sf_block, *sf_block_l, *nnodes_l,
        *nnodes_u, *edag_supno_l, *recvbuf, **edag_supno;
    float edag_supno_l_bytes;
#ifdef ISORT
    int_t *iperm_u;
#endif
    int *msgcnt;   /* Count the size of the message xfer'd in each buffer:
		    *     0 : transferred in Lsub_buf[]
		    *     1 : transferred in Lval_buf[]
		    *     2 : transferred in Usub_buf[]
		    *     3 : transferred in Uval_buf[]
		    */
    int **msgcnts, **msgcntsU;
    int *factored, *factoredU, nnodes, *sendcnts, *sdispls, *recvcnts,
        *rdispls, *srows, *rrows;
    etree_node *head, *tail, *ptr;
    int *num_child;
    int num_look_aheads, look_id, *look_ahead;
    int_t *perm_c_supno, *iperm_c_supno;
    MPI_Request *recv_req, **recv_reqs, **send_reqs, **send_reqs_u,
        **recv_reqs_u;
    MPI_Request *send_req, *U_diag_blk_send_req = NULL;
    MPI_Status status;
    void *attr_val;
    int flag;

    int iword = sizeof (int_t);
    int dword = sizeof (double);

    double scatter_timer = 0;
    double gemm_timer = 0;

    /* For measuring load imbalence in omp threads*/
    double omp_load_imblc = 0.0;
    double *omp_loop_time;

    double CPUOffloadTimer      = 0;
    double CPUOffloadFlop       = 0;
    double CPUOffloadMop        = 0;
    double schur_flop_timer     = 0.0;
    double pdgstrf2_timer       = 0.0;
    double pdgstrs2_timer       = 0.0;
    double lookaheadupdatetimer = 0.0;

#if !defined( GPU_ACC )
    /* Counter for couting memory operations */
    double scatter_mem_op_counter  = 0.0;
    double scatter_mem_op_timer    = 0.0;
    double scatterL_mem_op_counter = 0.0;
    double scatterL_mem_op_timer   = 0.0;
    double scatterU_mem_op_counter = 0.0;
    double scatterU_mem_op_timer   = 0.0;

    double LookAheadRowSepTimer    = 0.0;
    double LookAheadRowSepMOP      = 0.0;
    double GatherTimer             = 0.0;
    double GatherMOP               = 0.0;
    double LookAheadGEMMTimer      = 0.0;
    double LookAheadGEMMFlOp       = 0.0;
    double LookAheadScatterTimer   = 0.0;
    double LookAheadScatterMOP     = 0.0;
    double schur_flop_counter      = 0.0;
#endif

#if ( DEBUGlevel>=2 )
    int_t num_copy = 0, num_update = 0;
#endif
#if ( PRNTlevel==3 )
    int zero_msg = 0, total_msg = 0;
#endif
#if ( PROFlevel>=1 )
    double t1, t2;
    float msg_vol = 0, msg_cnt = 0;
#endif

    /* Test the input parameters. */
    *info = 0;
    if (m < 0)
        *info = -2;
    else if (n < 0)
        *info = -3;
    if (*info) {
        pxerr_dist ("pdgstrf", grid, -*info);
        return (-1);
    }

    /* Quick return if possible. */
    if (m == 0 || n == 0) return 0;
 
    /* 
     * Initialization.  
     */
    iam = grid->iam;
    Pc = grid->npcol; 
    Pr = grid->nprow;
    myrow = MYROW (iam, grid);
    mycol = MYCOL (iam, grid);
    nsupers = Glu_persist->supno[n - 1] + 1;
    xsup = Glu_persist->xsup;
    s_eps = smach_dist("Epsilon");
    thresh = s_eps * anorm;

    MPI_Attr_get (MPI_COMM_WORLD, MPI_TAG_UB, &attr_val, &flag);
    if (!flag) {
        fprintf (stderr, "Could not get TAG_UB\n");
        return (-1);
    }
    int tag_ub = *(int *) attr_val;

#if ( PRNTlevel>=1 )
    if (!iam)
        printf ("MPI tag upper bound = %d\n", tag_ub);
#endif

#if ( DEBUGlevel>=1 )
    if (s_eps == 0.0)
        printf (" ***** warning s_eps = %e *****\n", s_eps);
    CHECK_MALLOC (iam, "Enter pdgstrf()");
#endif

    stat->ops[FACT]      = 0.0;
    stat->current_buffer = 0.0;
    stat->peak_buffer    = 0.0;
    stat->gpu_buffer     = 0.0;

    /* make sure the range of look-ahead window [0, MAX_LOOKAHEADS-1] */
    num_look_aheads = SUPERLU_MAX(0, SUPERLU_MIN(options->num_lookaheads, MAX_LOOKAHEADS - 1));

    if (Pr * Pc > 1) {
        if (!(U_diag_blk_send_req =
              (MPI_Request *) SUPERLU_MALLOC (Pr * sizeof (MPI_Request))))
            ABORT ("Malloc fails for U_diag_blk_send_req[].");
	/* flag no outstanding Isend */
        U_diag_blk_send_req[myrow] = MPI_REQUEST_NULL; /* used 0 before */

        /* allocating buffers for look-ahead */
        i = Llu->bufmax[0];
        if (i != 0) {
            if ( !(Llu->Lsub_buf_2[0] = intMalloc_dist ((num_look_aheads + 1) * ((size_t) i))) )
                ABORT ("Malloc fails for Lsub_buf.");
            for (jj = 0; jj < num_look_aheads; jj++)
                Llu->Lsub_buf_2[jj + 1] = Llu->Lsub_buf_2[jj] + i;
        }
        i = Llu->bufmax[1];
        if (i != 0) {
            if (!(Llu->Lval_buf_2[0] = doubleMalloc_dist ((num_look_aheads + 1) * ((size_t) i))))
                ABORT ("Malloc fails for Lval_buf[].");
            for (jj = 0; jj < num_look_aheads; jj++)
                Llu->Lval_buf_2[jj + 1] = Llu->Lval_buf_2[jj] + i;
        }
        i = Llu->bufmax[2];
        if (i != 0) {
            if (!(Llu->Usub_buf_2[0] = intMalloc_dist ((num_look_aheads + 1) * i)))
                ABORT ("Malloc fails for Usub_buf_2[].");
            for (jj = 0; jj < num_look_aheads; jj++)
                Llu->Usub_buf_2[jj + 1] = Llu->Usub_buf_2[jj] + i;
        }
        i = Llu->bufmax[3];
        if (i != 0) {
            if (!(Llu->Uval_buf_2[0] = doubleMalloc_dist ((num_look_aheads + 1) * i)))
                ABORT ("Malloc fails for Uval_buf_2[].");
            for (jj = 0; jj < num_look_aheads; jj++)
                Llu->Uval_buf_2[jj + 1] = Llu->Uval_buf_2[jj] + i;
        }
    }

    log_memory( (Llu->bufmax[0] + Llu->bufmax[2]) * (num_look_aheads + 1) 
		* iword +
		(Llu->bufmax[1] + Llu->bufmax[3]) * (num_look_aheads + 1) 
		* dword, stat );

    /* creating pointers to the look-ahead buffers */
    if (! (Lsub_buf_2 = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (int_t *))))
        ABORT ("Malloc fails for Lsub_buf_2[].");
    if (! (Lval_buf_2 = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (double *))))
        ABORT ("Malloc fails for Lval_buf_2[].");
    if (! (Usub_buf_2 = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (int_t *))))
        ABORT ("Malloc fails for Uval_buf_2[].");
    if (! (Uval_buf_2 = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (double *))))
        ABORT ("Malloc fails for buf_2[].");
    for (i = 0; i <= num_look_aheads; i++) {
        Lval_buf_2[i] = Llu->Lval_buf_2[i];
        Lsub_buf_2[i] = Llu->Lsub_buf_2[i];
        Uval_buf_2[i] = Llu->Uval_buf_2[i];
        Usub_buf_2[i] = Llu->Usub_buf_2[i];
    }

    if (!(msgcnts = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (int *))))
        ABORT ("Malloc fails for msgcnts[].");
    if (!(msgcntsU = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (int *))))
        ABORT ("Malloc fails for msgcntsU[].");
    for (i = 0; i <= num_look_aheads; i++) {
        if (!(msgcnts[i] = SUPERLU_MALLOC (4 * sizeof (int))))
            ABORT ("Malloc fails for msgcnts[].");
        if (!(msgcntsU[i] = SUPERLU_MALLOC (4 * sizeof (int))))
            ABORT ("Malloc fails for msgcntsU[].");
    }

    if (! (recv_reqs_u = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (MPI_Request *))))
        ABORT ("Malloc fails for recv_reqs_u[].");
    if (! (send_reqs_u = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (MPI_Request *))))
        ABORT ("Malloc fails for send_reqs_u[].");
    if (! (send_reqs = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (MPI_Request *))))
        ABORT ("Malloc fails for send_reqs_u[].");
    if (! (recv_reqs = SUPERLU_MALLOC ((1 + num_look_aheads) * sizeof (MPI_Request *))))
        ABORT ("Malloc fails for recv_reqs[].");
    for (i = 0; i <= num_look_aheads; i++) {
        if (!(recv_reqs_u[i] = (MPI_Request *) SUPERLU_MALLOC (2 * sizeof (MPI_Request))))
            ABORT ("Malloc fails for recv_req_u[i].");
        if (!(send_reqs_u[i] = (MPI_Request *) SUPERLU_MALLOC (2 * Pr * sizeof (MPI_Request))))
            ABORT ("Malloc fails for send_req_u[i].");
        if (!(send_reqs[i] = (MPI_Request *) SUPERLU_MALLOC (2 * Pc * sizeof (MPI_Request))))
            ABORT ("Malloc fails for send_reqs[i].");
        if (!(recv_reqs[i] = (MPI_Request *) SUPERLU_MALLOC (4 * sizeof (MPI_Request))))
            ABORT ("Malloc fails for recv_req[].");
        send_reqs[i][0] = send_reqs[i][1] = MPI_REQUEST_NULL;
        recv_reqs[i][0] = recv_reqs[i][1] = MPI_REQUEST_NULL;
    }

    if (!(factored = SUPERLU_MALLOC (nsupers * sizeof (int_t))))
        ABORT ("Malloc fails for factored[].");
    if (!(factoredU = SUPERLU_MALLOC (nsupers * sizeof (int_t))))
        ABORT ("Malloc fails for factoredU[].");
    for (i = 0; i < nsupers; i++) factored[i] = factoredU[i] = -1;
    log_memory(2 * nsupers * iword, stat);

    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel default(shared)
    {
        if (omp_get_thread_num () == 0) {
            num_threads = omp_get_num_threads ();
        }
    }
#endif

#if 0
    omp_loop_time = (double *) _mm_malloc (sizeof (double) * num_threads,64);
#else
    omp_loop_time = (double *) doubleMalloc_dist(num_threads);
#endif

#if ( PRNTlevel>=1 )
    if(!iam) printf(".. Starting with %d OpenMP threads \n", num_threads );
    double tt1 = SuperLU_timer_ ();
#endif
    nblocks = 0;
    ncb = nsupers / Pc;
    nrb = nsupers / Pr;

    int nstreams = get_num_cuda_streams ();
    /* int nstreams = NUM_CUDA_STREAMS;    */
    /* in order to have dynamic scheduling */
    int *full_u_cols;
    int *blk_ldu;
#if 0
    full_u_cols = (int_t *) _mm_malloc (sizeof (int_t) * ncb,64);
    blk_ldu = (int_t *) _mm_malloc (sizeof (int_t) * ncb,64);
#else
    full_u_cols = SUPERLU_MALLOC(ncb * sizeof(int));
    blk_ldu = SUPERLU_MALLOC(ncb * sizeof(int));
#endif
    log_memory(2 * ncb * iword, stat);

    /* array holding last column blk for each partition,
       used in SchCompUdt--CUDA.c         */
#if 0
    int *stream_end_col = (int_t *) _mm_malloc (sizeof (int_t) * nstreams,64);
#else
    int *stream_end_col = SUPERLU_MALLOC( nstreams * sizeof(int) );
#endif

    /* insert a check condition here */

#if 0  /* Sherry: not used? */
    /* This bunch is used for static scheduling */
    pair *full_col_count = (pair *) _mm_malloc (sizeof (pair) * ncb,64);
    int_t *count_cols, *sum_cols, *partition;
    count_cols = (int_t *) _mm_malloc (sizeof (int_t) * num_threads,64);
    sum_cols = (int_t *) _mm_malloc (sizeof (int_t) * num_threads,64);
    partition = (int_t *) _mm_malloc (sizeof (int_t) * num_threads * ncb,64);
    int_t ldp = ncb;
#endif

    /* ##################################################################
     *  Compute a good static schedule based on the factorization task graph.
     * ################################################################## */
    perm_c_supno = SUPERLU_MALLOC (2 * nsupers * sizeof (int_t));
    iperm_c_supno = perm_c_supno + nsupers;

    static_schedule(options, m, n, LUstruct, grid, stat,
		    perm_c_supno, iperm_c_supno, info);

#if ( DEBUGlevel >= 2 )
    PrintInt10("schedule:perm_c_supno", nsupers, perm_c_supno);
    
    printf("[%d] .. Turn off static schedule for debugging ..\n", iam);
    /* Turn off static schedule */
    for (i = 0; i < nsupers; ++i) perm_c_supno[i] = iperm_c_supno[i] = i;
#endif
     /* ################################################################## */

    /* constructing look-ahead table to indicate the last dependency */
    int *look_ahead_l;
    stat->num_look_aheads = num_look_aheads;

    look_ahead_l = SUPERLU_MALLOC (nsupers * sizeof (int));
    look_ahead = SUPERLU_MALLOC (nsupers * sizeof (int));
    for (lb = 0; lb < nsupers; lb++) look_ahead_l[lb] = -1;
    log_memory(3 * nsupers * iword, stat);

    /* go through U-factor */
    for (lb = 0; lb < nrb; ++lb) {
        ib = lb * Pr + myrow;
        index = Llu->Ufstnz_br_ptr[lb];
        if (index) { /* Not an empty row */
            k = BR_HEADER;
            for (j = 0; j < index[0]; ++j) {
                jb = index[k];
                if (jb != ib)
                    look_ahead_l[jb] =
                        SUPERLU_MAX (iperm_c_supno[ib], look_ahead_l[jb]);
                k += UB_DESCRIPTOR + SuperSize (index[k]);
            }
        }
    }
    if (myrow < nsupers % grid->nprow) {
        ib = nrb * Pr + myrow;
        index = Llu->Ufstnz_br_ptr[nrb];
        if (index) {             /* Not an empty row */
            k = BR_HEADER;
            for (j = 0; j < index[0]; ++j) {
                jb = index[k];
                if (jb != ib)
                    look_ahead_l[jb] =
                        SUPERLU_MAX (iperm_c_supno[ib], look_ahead_l[jb]);
                k += UB_DESCRIPTOR + SuperSize (index[k]);
            }
        }
    }

    if (options->SymPattern == NO) {
        /* go through L-factor */
        for (lb = 0; lb < ncb; lb++) {
            ib = lb * Pc + mycol;
            index = Llu->Lrowind_bc_ptr[lb];
            if (index) {
                k = BC_HEADER;
                for (j = 0; j < index[0]; j++) {
                    jb = index[k];
                    if (jb != ib)
                        look_ahead_l[jb] =
                            SUPERLU_MAX (iperm_c_supno[ib], look_ahead_l[jb]);
                    k += LB_DESCRIPTOR + index[k + 1];
                }
            }
        }
        if (mycol < nsupers % grid->npcol) {
            ib = ncb * Pc + mycol;
            index = Llu->Lrowind_bc_ptr[ncb];
            if (index) {
                k = BC_HEADER;
                for (j = 0; j < index[0]; j++) {
                    jb = index[k];
                    if (jb != ib)
                        look_ahead_l[jb] =
                            SUPERLU_MAX (iperm_c_supno[ib], look_ahead_l[jb]);
                    k += LB_DESCRIPTOR + index[k + 1];
                }
            }
        }
    }
    MPI_Allreduce (look_ahead_l, look_ahead, nsupers, MPI_INT, MPI_MAX, grid->comm);
    SUPERLU_FREE (look_ahead_l);

#ifdef ISORT
    iperm_u = SUPERLU_MALLOC (nsupers * sizeof (int_t));
    perm_u = SUPERLU_MALLOC (nsupers * sizeof (int_t));
#else
    perm_u = SUPERLU_MALLOC (2 * nsupers * sizeof (int_t));
#endif
    log_memory(nsupers * iword, stat);

#if ( PRNTlevel>=1 )
    if (grid->iam == 0)
        printf (" * init: %e seconds\n", SuperLU_timer_ () - tt1);
#endif

    k = sp_ienv_dist (3);       /* max supernode size */
#if 0
    if ( !(Llu->ujrow = doubleMalloc_dist(k*(k+1)/2)) )
         ABORT("Malloc fails for ujrow[].");
#else
    /* Instead of half storage, we'll do full storage */
    if (!(Llu->ujrow = doubleMalloc_dist (k * k)))
        ABORT ("Malloc fails for ujrow[].");
    log_memory(k * k * iword, stat);
#endif

#if ( PRNTlevel>=1 )
    if (!iam) {
        printf (".. thresh = s_eps %e * anorm %e = %e\n", s_eps, anorm,
                thresh);
        printf
            (".. Buffer size: Lsub %ld\tLval %ld\tUsub %ld\tUval %ld\tLDA %ld\n",
             (long int) Llu->bufmax[0], (long int) Llu->bufmax[1],
             (long int) Llu->bufmax[2], (long int) Llu->bufmax[3],
             (long int) Llu->bufmax[4]);
    }
#endif
   
    Lrowind_bc_ptr = Llu->Lrowind_bc_ptr;
    Lnzval_bc_ptr = Llu->Lnzval_bc_ptr;
    Ufstnz_br_ptr = Llu->Ufstnz_br_ptr;
    Unzval_br_ptr = Llu->Unzval_br_ptr;
    ToRecv = Llu->ToRecv; 
    ToSendD = Llu->ToSendD;
    ToSendR = Llu->ToSendR;

    ldt = sp_ienv_dist (3);     /* Size of maximum supernode */
    k = CEILING (nsupers, Pr);  /* Number of local block rows */

    /* Following circuit is for finding maximum block size */
    int local_max_row_size = 0;
    int max_row_size;

    for (int i = 0; i < nsupers; ++i) {
        int tpc = PCOL (i, grid);
        if (mycol == tpc) {
            lk = LBj (i, grid);
            lsub = Lrowind_bc_ptr[lk];
            if (lsub != NULL) {
                local_max_row_size = SUPERLU_MAX (local_max_row_size, lsub[1]);
            }
        }

    }

    /* Max row size is global reduction of within A row */
    MPI_Allreduce (&local_max_row_size, &max_row_size, 1, MPI_INT, MPI_MAX, (grid->rscp.comm));

    /* Buffer size is max of of look ahead window */
    /* int_t buffer_size =
         SUPERLU_MAX (max_row_size * num_threads * ldt,
                      get_max_buffer_size ());           */
            
    int cublas_nb = get_cublas_nb();
#ifdef GPU_ACC
    int buffer_size  = SUPERLU_MAX(max_row_size*nstreams*cublas_nb,get_max_buffer_size());
#else 
    int Threads_per_process = get_thread_per_process();
    int buffer_size  = SUPERLU_MAX(max_row_size*Threads_per_process*ldt,get_max_buffer_size());
#endif 
 
    /* symmetric assumption */
    /* Note that in following expression 8 can be anything
       as long as its not too big */
    int bigu_size = 8 * sp_ienv_dist (3) * (max_row_size);

#if ( PRNTlevel>=1 )
    if(!iam) printf("[%d] .. BIG U size %d \n", iam, bigu_size);
#endif

#ifdef GPU_ACC

    // printf("hello 1\n");
    double* bigU;
    if ( checkCuda(cudaHostAlloc((void**)&bigU,  bigu_size * sizeof(double), cudaHostAllocDefault)) )
        ABORT("Malloc fails for dgemm buffer U ");

    int bigv_size = buffer_size;
#if ( PRNTlevel>=1 )
    if (!iam) printf("[%d] .. BIG V size %d\n", iam, bigv_size);
#endif
    double* bigV;
    if ( checkCuda(cudaHostAlloc((void**)&bigV, bigv_size * sizeof(double) ,cudaHostAllocDefault)) )
        ABORT("Malloc fails for dgemm buffer V");
 
    DisplayHeader();

#if ( PRNTlevel>=1 )
    printf(" Starting with %d Cuda Streams \n",nstreams );
#endif

    cublasHandle_t *handle;
    handle = (cublasHandle_t *) SUPERLU_MALLOC(sizeof(cublasHandle_t)*nstreams);
    for(int i = 0; i < nstreams; i++) handle[i] = create_handle();

    // creating streams 
    cudaStream_t *streams;
    streams = (cudaStream_t *) SUPERLU_MALLOC(sizeof(cudaStream_t)*nstreams);
    for (int i = 0; i < nstreams; ++i)
        checkCuda( cudaStreamCreate(&streams[i]) );
    
    // allocating data in device 
    double *dA, *dB, *dC;
    cudaError_t cudaStat;
#if 0
    // cudaStat = cudaMalloc( (void**)&dA, m*k*sizeof(double));
    // HOw much should be the size of dA?
    // for time being just making it 
    // cudaStat = cudaMalloc( (void**)&dA, ((max_row_size*sp_ienv_dist(3)))* sizeof(double));
#endif

    cudaStat = cudaMalloc( (void**)&dA, max_row_size*sp_ienv_dist(3)* sizeof(double));
    if (cudaStat!= cudaSuccess) {
        fprintf(stderr, "!!!! Error in allocating A in the device %ld \n",m*k*sizeof(double) );
        return 1;
    }

    // size of B should be max_supernode_size*buffer

    cudaStat = cudaMalloc((void**)&dB, bigu_size * sizeof(double));
    if (cudaStat!= cudaSuccess) {
        fprintf(stderr, "!!!! Error in allocating B in the device %ld \n",n*k*sizeof(double));
        return 1;
    }

    cudaStat = cudaMalloc((void**)&dC, buffer_size* sizeof(double) );
    if (cudaStat!= cudaSuccess) {
        fprintf(stderr, "!!!! Error in allocating C in the device \n" );
        return 1;
    }

    stat->gpu_buffer += ( max_row_size * sp_ienv_dist(3) 
			  + bigu_size + buffer_size ) * dword;

#else  /* not CUDA */
    
    double* bigU;
    if ( !(bigU = doubleMalloc_dist(bigu_size)) )
        ABORT ("Malloc fails for dgemm u buff U"); 
          //Maximum size of of bigU= sqrt(buffsize) ?

    int bigv_size = 8 * ldt * ldt * num_threads;
#if ( PRNTlevel>=1 )
    if (!iam) printf("[%d] .. BIG V size %d\n", iam, bigv_size);
#endif
    double *bigV;
    if ( !(bigV = doubleMalloc_dist(bigv_size)) )
        ABORT ("Malloc failed for dgemm buffer V");

#endif

    log_memory((bigv_size + bigu_size) * dword, stat);

    // mlock(bigU,(bigu_size) * sizeof (double));   

#if ( PRNTlevel>=1 )
    if(!iam) {
	printf ("  Max row size is %d \n", max_row_size);
        printf ("  Using buffer_size of %d \n", buffer_size);
        printf ("  Threads per process %d \n", num_threads);
    }
#endif

    if (!(tempv2d = doubleCalloc_dist (2 * ((size_t) ldt) * ldt)))
        ABORT ("Calloc fails for tempv2d[].");
    tempU2d = tempv2d + ldt * ldt;
    if (!(indirect = SUPERLU_MALLOC (ldt * num_threads * sizeof(int))))
        ABORT ("Malloc fails for indirect[].");
    if (!(indirect2 = SUPERLU_MALLOC (ldt * num_threads * sizeof(int))))
        ABORT ("Malloc fails for indirect[].");
    if (!(iuip = intMalloc_dist (k)))  ABORT ("Malloc fails for iuip[].");
    if (!(ruip = intMalloc_dist (k)))  ABORT ("Malloc fails for ruip[].");

    log_memory(2 * ldt *ldt * dword + 2 * ldt * num_threads * iword
	       + 2 * k * iword, stat);

    int_t *lookAheadFullRow,*lookAheadStRow,*lookAhead_lptr,*lookAhead_ib,
        *RemainFullRow,*RemainStRow,*Remain_lptr,*Remain_ib;

    lookAheadFullRow   = intMalloc_dist( (num_look_aheads+1) );
    lookAheadStRow     = intMalloc_dist( (num_look_aheads+1) );
    lookAhead_lptr     = intMalloc_dist( (num_look_aheads+1) );
    lookAhead_ib       = intMalloc_dist( (num_look_aheads+1) );

    int_t mrb=    (nsupers+Pr-1) / Pr;
    int_t mcb=    (nsupers+Pc-1) / Pc;
    
    RemainFullRow   = intMalloc_dist(mrb); 
    RemainStRow     = intMalloc_dist(mrb);
#if 0
    Remain_lptr     = (int *) _mm_malloc(sizeof(int)*mrb,1);
#else
    Remain_lptr     = intMalloc_dist(mrb);
#endif
    // mlock(Remain_lptr, sizeof(int)*mrb );
    Remain_ib       = intMalloc_dist(mrb);
    
    Remain_info_t *Remain_info;
#if 0
    Remain_info = (Remain_info_t *) _mm_malloc(mrb*sizeof(Remain_info_t),64);
#else
    Remain_info = (Remain_info_t *) SUPERLU_MALLOC(mrb*sizeof(Remain_info_t));
#endif
    log_memory(4 * mrb * iword + mrb * sizeof(Remain_info_t), stat);

    double *lookAhead_L_buff, *Remain_L_buff;
    Ublock_info_t *Ublock_info;
    ldt = sp_ienv_dist (3);       /* max supernode size */
    lookAhead_L_buff = doubleMalloc_dist(ldt*ldt* (num_look_aheads+1) );
    log_memory(ldt * ldt * (num_look_aheads+1) * dword, stat);

#if 0
    Remain_L_buff = (double *) _mm_malloc( sizeof(double)*(Llu->bufmax[1]),64);
    Ublock_info = (Ublock_info_t *) _mm_malloc(mcb*sizeof(Ublock_info_t),64);
    int * Ublock_info_iukp = (int *) _mm_malloc(mcb*sizeof(int),64);
    int * Ublock_info_rukp = (int *) _mm_malloc(mcb*sizeof(int),64);
    int * Ublock_info_jb = (int *) _mm_malloc(mcb*sizeof(int),64);
#else
    Remain_L_buff = doubleMalloc_dist(Llu->bufmax[1]);
    Ublock_info = (Ublock_info_t *) SUPERLU_MALLOC(mcb*sizeof(Ublock_info_t));
    int *Ublock_info_iukp = (int *) SUPERLU_MALLOC(mcb*sizeof(int));
    int *Ublock_info_rukp = (int *) SUPERLU_MALLOC(mcb*sizeof(int));
    int *Ublock_info_jb = (int *) SUPERLU_MALLOC(mcb*sizeof(int));
#endif
    log_memory(Llu->bufmax[1] * dword, stat);

    double NetSchurUpTimer = 0;
    double pdgstrfTimer= SuperLU_timer_();

    /* ##################################################################
       ** Handle first block column separately to start the pipeline. **
       ################################################################## */
    look_id = 0;
    msgcnt = msgcnts[0];
    send_req = send_reqs[0];
    recv_req = recv_reqs[0];

    k0 = 0;
    k = perm_c_supno[0];
    kcol = PCOL (k, grid);
    krow = PROW (k, grid);
    if (mycol == kcol) {
        double ttt1 = SuperLU_timer_();

	/* panel factorization */
        PDGSTRF2 (options, k0, k, thresh, Glu_persist, grid, Llu,
                  U_diag_blk_send_req, tag_ub, stat, info);

        pdgstrf2_timer += SuperLU_timer_()-ttt1; 

        scp = &grid->rscp;      /* The scope of process row. */

        /* Multicasts numeric values of L(:,0) to process rows. */
        lk = LBj (k, grid);     /* Local block number. */
        lsub = Lrowind_bc_ptr[lk];
        lusup = Lnzval_bc_ptr[lk];
        if (lsub) {
            msgcnt[0] = lsub[1] + BC_HEADER + lsub[0] * LB_DESCRIPTOR;
            msgcnt[1] = lsub[1] * SuperSize (k);
        } else {
            msgcnt[0] = msgcnt[1] = 0;
        }

        for (pj = 0; pj < Pc; ++pj) {
            if (ToSendR[lk][pj] != EMPTY) {
#if ( PROFlevel>=1 )
                TIC (t1);
#endif

                MPI_Isend (lsub, msgcnt[0], mpi_int_t, pj, SLU_MPI_TAG (0, 0) /* 0 */ ,
                           scp->comm, &send_req[pj]);
                MPI_Isend (lusup, msgcnt[1], MPI_DOUBLE, pj, SLU_MPI_TAG (1, 0) /* 1 */ ,
                           scp->comm, &send_req[pj + Pc]);
#if ( DEBUGlevel>=2 )
                printf ("[%d] first block cloumn Send L(:,%4d): lsub %4d, lusup %4d to Pc %2d\n",
                        iam, 0, msgcnt[0], msgcnt[1], pj);
#endif

#if ( PROFlevel>=1 )
                TOC (t2, t1);
                stat->utime[COMM] += t2;
                msg_cnt += 2;
                msg_vol += msgcnt[0] * iword + msgcnt[1] * dword;
#endif
            } /* end if */
        }  /* end for pj ... */
    } else {  /* Post immediate receives. */
        if (ToRecv[k] >= 1) {   /* Recv block column L(:,0). */
            scp = &grid->rscp;  /* The scope of process row. */
            MPI_Irecv (Lsub_buf_2[0], Llu->bufmax[0], mpi_int_t, kcol,
                       SLU_MPI_TAG (0, 0) /* 0 */ ,
                       scp->comm, &recv_req[0]);
            MPI_Irecv (Lval_buf_2[0], Llu->bufmax[1], MPI_DOUBLE, kcol,
                       SLU_MPI_TAG (1, 0) /* 1 */ ,
                       scp->comm, &recv_req[1]);
        }
    } /* end if mycol == 0 */

    factored[k] = 0; /* flag column k as factored. */

    /* post receive of first U-row */
    if (myrow != krow) {
        if (ToRecv[k] == 2) {   /* Recv block row U(k,:). */
            scp = &grid->cscp;  /* The scope of process column. */
            Usub_buf = Llu->Usub_buf_2[0];
            Uval_buf = Llu->Uval_buf_2[0];
            MPI_Irecv (Usub_buf, Llu->bufmax[2], mpi_int_t, krow,
                       SLU_MPI_TAG (2, 0) /* 2%tag_ub */ ,
                       scp->comm, &recv_reqs_u[0][0]);
            MPI_Irecv (Uval_buf, Llu->bufmax[3], MPI_DOUBLE, krow,
                       SLU_MPI_TAG (3, 0) /* 3%tag_ub */ ,
                       scp->comm, &recv_reqs_u[0][1]);
        }
    }

    /* ##################################################################
       **** MAIN LOOP ****
       ################################################################## */
    for (k0 = 0; k0 < nsupers; ++k0) {
        k = perm_c_supno[k0];

        /* ============================================ *
         * ======== look-ahead the new columns ======== *
         * ============================================ */
        /* tt1 = SuperLU_timer_(); */
        if (k0 == 0) { /* look-ahead all the columns in the window */
            kk1 = k0 + 1;
            kk2 = SUPERLU_MIN (k0 + num_look_aheads, nsupers - 1);
        } else {  /* look-ahead one new column after the current window */
            kk1 = k0 + num_look_aheads;
            kk2 = SUPERLU_MIN (kk1, nsupers - 1);
        }

        for (kk0 = kk1; kk0 <= kk2; kk0++) {
	    /* loop through look-ahead window */

            kk = perm_c_supno[kk0]; /* use the ordering from static schedule */
            look_id = kk0 % (1 + num_look_aheads); /* which column in window */

            if (look_ahead[kk] < k0) { /* does not depend on current column */
                kcol = PCOL (kk, grid);
                if (mycol == kcol) {

                    /* Panel factorization -- Factor diagonal and subdiagonal
                       L blocks and test for exact singularity.  */
                    factored[kk] = 0; /* flag column kk as factored */
                    double ttt1 = SuperLU_timer_();

                    PDGSTRF2 (options, kk0, kk, thresh, Glu_persist,
                              grid, Llu, U_diag_blk_send_req, tag_ub, stat, info);
                     pdgstrf2_timer += SuperLU_timer_() - ttt1; 

                    /* Multicasts numeric values of L(:,kk) to process rows. */
                    /* ttt1 = SuperLU_timer_(); */
                    msgcnt = msgcnts[look_id];  /* point to the proper count array */
                    send_req = send_reqs[look_id];

                    lk = LBj (kk, grid);    /* Local block number. */
                    lsub1 = Lrowind_bc_ptr[lk];
                    if (lsub1) {
                        msgcnt[0] = lsub1[1] + BC_HEADER + lsub1[0] * LB_DESCRIPTOR;
                        msgcnt[1] = lsub1[1] * SuperSize (kk);
                    } else {
                        msgcnt[0] = 0;
                        msgcnt[1] = 0;
                    }
                    scp = &grid->rscp;  /* The scope of process row. */
                    for (pj = 0; pj < Pc; ++pj) {
                        if (ToSendR[lk][pj] != EMPTY) {
                            lusup1 = Lnzval_bc_ptr[lk];
                            MPI_Isend (lsub1, msgcnt[0], mpi_int_t, pj,
                                       SLU_MPI_TAG (0, kk0),  /* (4*kk0)%tag_ub */
                                       scp->comm, &send_req[pj]);
                            MPI_Isend (lusup1, msgcnt[1], MPI_DOUBLE, pj,
                                       SLU_MPI_TAG (1, kk0),  /* (4*kk0+1)%tag_ub */
                                       scp->comm, &send_req[pj + Pc]);
#if ( DEBUGlevel>=2 )
			    printf ("[%d] -1- Send L(:,%4d): #lsub1 %4d, #lusup1 %4d right to Pj %2d\n",
				    iam, kk, msgcnt[0], msgcnt[1], pj);
#endif
                        }
                    }
                    /* stat->time9 += SuperLU_timer_() - ttt1; */
                } else {     /* Post Recv of block column L(:,kk). */
                    /* double ttt1 = SuperLU_timer_(); */
                    if (ToRecv[kk] >= 1) {
                        scp = &grid->rscp;  /* The scope of process row. */
                        recv_req = recv_reqs[look_id];

                        MPI_Irecv (Lsub_buf_2[look_id], Llu->bufmax[0],
                                   mpi_int_t, kcol, SLU_MPI_TAG (0, kk0), /* (4*kk0)%tag_ub */
                                   scp->comm, &recv_req[0]);
                        MPI_Irecv (Lval_buf_2[look_id], Llu->bufmax[1],
                                   MPI_DOUBLE, kcol,
                                   SLU_MPI_TAG (1, kk0), /* (4*kk0+1)%tag_ub */
                                   scp->comm, &recv_req[1]);
                    }
                    /* stat->time10 += SuperLU_timer_() - ttt1; */
                }  /* end if mycol == Pc(kk) */
            }  /* end if look-ahead */

            /* post irecv for U-row look-ahead */
            krow = PROW (kk, grid);
            if (myrow != krow) {
                if (ToRecv[kk] == 2) { /* post iRecv block row U(k,:). */
                    scp = &grid->cscp;  /* The scope of process column. */
                    Usub_buf = Llu->Usub_buf_2[look_id];
                    Uval_buf = Llu->Uval_buf_2[look_id];

                    MPI_Irecv (Usub_buf, Llu->bufmax[2], mpi_int_t, krow,
                               SLU_MPI_TAG (2, kk0) /* (4*kk0+2)%tag_ub */ ,
                               scp->comm, &recv_reqs_u[look_id][0]);
                    MPI_Irecv (Uval_buf, Llu->bufmax[3], MPI_DOUBLE, krow,
                               SLU_MPI_TAG (3, kk0) /* (4*kk0+3)%tag_ub */ ,
                               scp->comm, &recv_reqs_u[look_id][1]);
                }
            }

        }  /* end for each column in look-ahead window */

        /* stat->time4 += SuperLU_timer_()-tt1; */

        /* ================================= *
         * == looking-ahead the U columns == *
         * ================================= */
        kk1 = k0;
        kk2 = SUPERLU_MIN (k0 + num_look_aheads, nsupers - 1);
        for (kk0 = kk1; kk0 < kk2; kk0++) {
            kk = perm_c_supno[kk0];
            if (factoredU[kk0] != 1 && look_ahead[kk] < k0) {
                kcol = PCOL (kk, grid);
                krow = PROW (kk, grid);
                lk = LBj (kk, grid);    /* Local block number. */

                look_id = kk0 % (1 + num_look_aheads);
                msgcnt = msgcntsU[look_id];
                recv_req = recv_reqs[look_id];

                /* ================================================= *
                 * checking if diagonal block has been received      *
                 * for panel factorization of U in look-ahead window *
                 * ================================================= */

                if (mycol == kcol) {
                    flag0 = flag1 = 1;
                    msgcnt[0] = msgcnt[1] = -1;
                } else {
                    flag0 = flag1 = 0;
                    if (ToRecv[kk] >= 1) {
                        if (recv_req[0] != MPI_REQUEST_NULL) {
                            MPI_Test (&recv_req[0], &flag0, &status);
                            if (flag0) {
                                MPI_Get_count (&status, mpi_int_t, &msgcnt[0]);
                                recv_req[0] = MPI_REQUEST_NULL;
                            }
                        } else flag0 = 1;

                        if (recv_req[1] != MPI_REQUEST_NULL) {
                            MPI_Test (&recv_req[1], &flag1, &status);
                            if (flag1) {
                                MPI_Get_count (&status, mpi_int_t, &msgcnt[1]);
                                recv_req[1] = MPI_REQUEST_NULL;
                            }
                        } else flag1 = 1;
                    } else msgcnt[0] = 0;
                }

                if (flag0 && flag1) {
                    /* tt1 = SuperLU_timer_(); */
                    scp = &grid->cscp;  /* The scope of process column. */
                    if (myrow == krow) {
                        factoredU[kk0] = 1;
                        /* Parallel triangular solve across process row *krow* --
                           U(k,j) = L(k,k) \ A(k,j).  */
                        /* double ttt2 = SuperLU_timer_(); */
                        double ttt2 = SuperLU_timer_();
#ifdef _OPENMP
#pragma omp parallel
#endif
			{
                            PDGSTRS2 (kk0, kk, Glu_persist, grid, Llu,
                                      stat);
                        }
    
                        pdgstrs2_timer += SuperLU_timer_()-ttt2;
                        /* stat->time8 += SuperLU_timer_()-ttt2; */

                        /* Multicasts U(k,:) to process columns. */
                        lk = LBi (kk, grid);
                        usub = Ufstnz_br_ptr[lk];
                        uval = Unzval_br_ptr[lk];
                        if (usub) {
                            msgcnt[2] = usub[2];
                            msgcnt[3] = usub[1];
                        } else {
                            msgcnt[2] = msgcnt[3] = 0;
                        }

                        if (ToSendD[lk] == YES) {
                            for (pi = 0; pi < Pr; ++pi) {
                                if (pi != myrow) {
#if ( PROFlevel>=1 )
                                    TIC (t1);
#endif

                                    MPI_Isend (usub, msgcnt[2], mpi_int_t, pi,
                                               SLU_MPI_TAG (2, kk0), /* (4*kk0+2)%tag_ub */
                                               scp->comm, &send_reqs_u[look_id][pi]);
                                    MPI_Isend (uval, msgcnt[3], MPI_DOUBLE,
                                               pi, SLU_MPI_TAG (3, kk0), /* (4*kk0+3)%tag_ub */
                                               scp->comm, &send_reqs_u[look_id][pi + Pr]);

#if ( PROFlevel>=1 )
                                    TOC (t2, t1);
                                    stat->utime[COMM] += t2;
                                    msg_cnt += 2;
                                    msg_vol += msgcnt[2] * iword + msgcnt[3] * dword;
#endif
#if ( DEBUGlevel>=2 )
                                    printf ("[%d] Send U(%4d,:) to Pr %2d\n",
                                            iam, k, pi);
#endif
                                }   /* if pi ... */
                            }   /* for pi ... */
                        }       /* if ToSendD ... */

                        /* stat->time2 += SuperLU_timer_()-tt1; */

                    } /* end if myrow == krow */
                } /* end if flag0 ... */
            } /* end if factoredU[] ... */
        } /* end for kk0 ... */

        /* ============================================== *
         * == start processing the current row of U == *
         * ============================================== */
        knsupc = SuperSize (k);
        krow = PROW (k, grid);
        kcol = PCOL (k, grid);

        /* tt1 = SuperLU_timer_(); */
        look_id = k0 % (1 + num_look_aheads);
        recv_req = recv_reqs[look_id];
        send_req = send_reqs[look_id];
        msgcnt = msgcnts[look_id];
        Usub_buf = Llu->Usub_buf_2[look_id];
        Uval_buf = Llu->Uval_buf_2[look_id];

        if (mycol == kcol) {
            lk = LBj (k, grid); /* Local block number. */

            for (pj = 0; pj < Pc; ++pj) {
                /* Wait for Isend to complete before using lsub/lusup. */
                if (ToSendR[lk][pj] != EMPTY) {
                    MPI_Wait (&send_req[pj], &status);
                    MPI_Wait (&send_req[pj + Pc], &status);
                }
            }
            lsub = Lrowind_bc_ptr[lk];
            lusup = Lnzval_bc_ptr[lk];
        } else {
            if (ToRecv[k] >= 1) { /* Recv block column L(:,k). */
                scp = &grid->rscp;  /* The scope of process row. */
                /* ============================================ *
                 * waiting for L(:,kk) for outer-product uptate *
                 * if iam in U(kk,:) then                       *
                 *  the diagonal block did not reach in time    *
                 *  for panel factorization of U(k,:)           *
                 * ============================================ */
#if ( PROFlevel>=1 )
                TIC (t1);
#endif
                if (recv_req[0] != MPI_REQUEST_NULL) {
                    MPI_Wait (&recv_req[0], &status);
                    MPI_Get_count (&status, mpi_int_t, &msgcnt[0]);
                    recv_req[0] = MPI_REQUEST_NULL;
                } else {
                    msgcnt[0] = msgcntsU[look_id][0];
#if (DEBUGlevel>=2)
		    printf("\t[%d] k=%d, look_id=%d, recv_req[0] == MPI_REQUEST_NULL, msgcnt[0] = %d\n", 
			   iam, k, look_id, msgcnt[0]);
#endif
                }

                if (recv_req[1] != MPI_REQUEST_NULL) {
                    MPI_Wait (&recv_req[1], &status);
                    MPI_Get_count (&status, MPI_DOUBLE, &msgcnt[1]);
                    recv_req[1] = MPI_REQUEST_NULL;
                } else {
                    msgcnt[1] = msgcntsU[look_id][1];
#if (DEBUGlevel>=2)
		    printf("\t[%d] k=%d, look_id=%d, recv_req[1] == MPI_REQUEST_NULL, msgcnt[1] = %d\n", 
			   iam, k, look_id, msgcnt[1]);
#endif
                }

#if ( PROFlevel>=1 )
                TOC (t2, t1);
                stat->utime[COMM] += t2;
#endif
#if ( DEBUGlevel>=2 )
                printf("[%d] Recv L(:,%4d): #lsub %4d, #lusup %4d from Pc %2d\n",
                     iam, k, msgcnt[0], msgcnt[1], kcol);
                fflush (stdout);
#endif

#if ( PRNTlevel==3 )
                ++total_msg;
                if (!msgcnt[0])  ++zero_msg;
#endif
            } else {
                msgcnt[0] = 0;
	    }

            lsub = Lsub_buf_2[look_id];
            lusup = Lval_buf_2[look_id];
        }                       /* if mycol = Pc(k) */
        /* stat->time1 += SuperLU_timer_()-tt1; */

        scp = &grid->cscp;      /* The scope of process column. */

        /* tt1 = SuperLU_timer_(); */
        if (myrow == krow) {
            lk = LBi (k, grid);
            usub = Ufstnz_br_ptr[lk];
            uval = Unzval_br_ptr[lk];

            if (factoredU[k0] == -1) {
                /* Parallel triangular solve across process row *krow* --
                   U(k,j) = L(k,k) \ A(k,j).  */
                 double ttt2 = SuperLU_timer_(); 
#ifdef _OPENMP
#pragma omp parallel
#endif
                {
                    PDGSTRS2 (k0, k, Glu_persist, grid, Llu, stat);
                }
                pdgstrs2_timer += SuperLU_timer_() - ttt2; 

                /* Multicasts U(k,:) along process columns. */
                if (usub) {
                    msgcnt[2] = usub[2];
                    msgcnt[3] = usub[1];
                } else {
                    msgcnt[2] = msgcnt[3] = 0;
                }

                if (ToSendD[lk] == YES) {
                    for (pi = 0; pi < Pr; ++pi) {
                        if (pi != myrow) {
#if ( PROFlevel>=1 )
                            TIC (t1);
#endif
                            MPI_Send (usub, msgcnt[2], mpi_int_t, pi,
                                      SLU_MPI_TAG (2, k0), /* (4*k0+2)%tag_ub */
                                      scp->comm);
                            MPI_Send (uval, msgcnt[3], MPI_DOUBLE, pi,
                                      SLU_MPI_TAG (3, k0), /* (4*k0+3)%tag_ub */ 
                                      scp->comm);
#if ( PROFlevel>=1 )
                            TOC (t2, t1);
                            stat->utime[COMM] += t2;
                            msg_cnt += 2;
                            msg_vol += msgcnt[2] * iword + msgcnt[3] * dword;
#endif
#if ( DEBUGlevel>=2 )
                            printf ("[%d] Send U(%4d,:) down to Pr %2d\n", iam, k, pi);
#endif
                        }       /* if pi ... */
                    }           /* for pi ... */
                }               /* if ToSendD ... */
            } else {

                /* =========================================== *
                 * waiting for U(k,:) for outer-product update *
                 * =========================================== */
                if (ToSendD[lk] == YES) {
                    for (pi = 0; pi < Pr; ++pi) {
                        if (pi != myrow) {
                            MPI_Wait (&send_reqs_u[look_id][pi], &status);
                            MPI_Wait (&send_reqs_u[look_id][pi + Pr], &status);
                        }
                    }
                }
                msgcnt[2] = msgcntsU[look_id][2];
                msgcnt[3] = msgcntsU[look_id][3];
            }
            /* stat->time2 += SuperLU_timer_()-tt1; */
        } else {    /* myrow != krow */

            /* ========================================= *
             * wait for U(k,:) for outer-product updates *
             * ========================================= */

            if (ToRecv[k] == 2) { /* Recv block row U(k,:). */
#if ( PROFlevel>=1 )
                TIC (t1);
#endif
                MPI_Wait (&recv_reqs_u[look_id][0], &status);
                MPI_Get_count (&status, mpi_int_t, &msgcnt[2]);
                MPI_Wait (&recv_reqs_u[look_id][1], &status);
                MPI_Get_count (&status, MPI_DOUBLE, &msgcnt[3]);

#if ( PROFlevel>=1 )
                TOC (t2, t1);
                stat->utime[COMM] += t2;
#endif
                usub = Usub_buf;
                uval = Uval_buf;
#if ( DEBUGlevel>=2 )
                printf ("[%d] Recv U(%4d,:) from Pr %2d\n", iam, k, krow);
#endif
#if ( PRNTlevel==3 )
                ++total_msg;
                if (!msgcnt[2])  ++zero_msg;
#endif
            } else {
                msgcnt[2] = 0;
	    }
            /* stat->time6 += SuperLU_timer_()-tt1; */
        }                       /* if myrow == Pr(k) */

        /*
         * Parallel rank-k update; pair up blocks L(i,k) and U(k,j).
         *  for (j = k+1; k < N; ++k) {
         *     for (i = k+1; i < N; ++i)
         *         if ( myrow == PROW( i, grid ) && mycol == PCOL( j, grid )
         *              && L(i,k) != 0 && U(k,j) != 0 )
         *             A(i,j) = A(i,j) - L(i,k) * U(k,j);
         */
        msg0 = msgcnt[0];
        msg2 = msgcnt[2];
        /* tt1 = SuperLU_timer_(); */
        if (msg0 && msg2) {     /* L(:,k) and U(k,:) are not empty. */
            nsupr = lsub[1];    /* LDA of lusup. */
            if (myrow == krow) { /* Skip diagonal block L(k,k). */
                lptr0 = BC_HEADER + LB_DESCRIPTOR + lsub[BC_HEADER + 1];
                luptr0 = knsupc;
                nlb = lsub[0] - 1;
            } else {
                lptr0 = BC_HEADER;
                luptr0 = 0;
                nlb = lsub[0];
            }
            iukp = BR_HEADER;   /* Skip header; Pointer to index[] of U(k,:) */
            rukp = 0;           /* Pointer to nzval[] of U(k,:) */
            nub = usub[0];      /* Number of blocks in the block row U(k,:) */
            klst = FstBlockC (k + 1);

            /* -------------------------------------------------------------
               Update the look-ahead block columns A(:,k+1:k+num_look_ahead)
               ------------------------------------------------------------- */
            iukp0 = iukp;
            rukp0 = rukp;
            /* reorder the remaining columns in bottome-up */
            /* TAU_STATIC_TIMER_START("LOOK_AHEAD_UPDATE"); */
            for (jj = 0; jj < nub; jj++) {
#ifdef ISORT
                iperm_u[jj] = iperm_c_supno[usub[iukp]];    /* Global block number of block U(k,j). */
                perm_u[jj] = jj;
#else
                perm_u[2 * jj] = iperm_c_supno[usub[iukp]]; /* Global block number of block U(k,j). */
                perm_u[2 * jj + 1] = jj;
#endif
                jb = usub[iukp];    /* Global block number of block U(k,j). */
                nsupc = SuperSize (jb);
                iukp += UB_DESCRIPTOR;  /* Start fstnz of block U(k,j). */
                iukp += nsupc;
            }
            iukp = iukp0;
#ifdef ISORT
            isort (nub, iperm_u, perm_u);
#else
            qsort (perm_u, (size_t) nub, 2 * sizeof (int_t),
                   &superlu_sort_perm);
#endif
            j = jj0 = 0;

/************************************************************************/
            double ttx =SuperLU_timer_();

#include "dlook_ahead_update.c"

            lookaheadupdatetimer += SuperLU_timer_() - ttx;
/************************************************************************/

            /*ifdef OMP_LOOK_AHEAD */
            /* TAU_STATIC_TIMER_STOP("LOOK_AHEAD_UPDATE"); */
        }                       /* if L(:,k) and U(k,:) not empty */

        /* stat->time3 += SuperLU_timer_()-tt1; */

        /* ================== */
        /* == post receive == */
        /* ================== */
        kk1 = SUPERLU_MIN (k0 + num_look_aheads, nsupers - 1);
        for (kk0 = k0 + 1; kk0 <= kk1; kk0++) {
            kk = perm_c_supno[kk0];
            kcol = PCOL (kk, grid);

            if (look_ahead[kk] == k0) {
                if (mycol != kcol) {
                    if (ToRecv[kk] >= 1) {
                        scp = &grid->rscp;  /* The scope of process row. */

                        look_id = kk0 % (1 + num_look_aheads);
                        recv_req = recv_reqs[look_id];
                        MPI_Irecv (Lsub_buf_2[look_id], Llu->bufmax[0],
                                   mpi_int_t, kcol, SLU_MPI_TAG (0, kk0), /* (4*kk0)%tag_ub */
                                   scp->comm, &recv_req[0]);
                        MPI_Irecv (Lval_buf_2[look_id], Llu->bufmax[1],
                                   MPI_DOUBLE, kcol,
                                   SLU_MPI_TAG (1, kk0), /* (4*kk0+1)%tag_ub */
                                   scp->comm, &recv_req[1]);
                    }
                } else {
                    lk = LBj (kk, grid);    /* Local block number. */
                    lsub1 = Lrowind_bc_ptr[lk];
                    lusup1 = Lnzval_bc_ptr[lk];
                    if (factored[kk] == -1) {
                        /* Factor diagonal and subdiagonal blocks and
			   test for exact singularity.  */
                        factored[kk] = 0; /* flag column kk as factored */
                        double ttt1 = SuperLU_timer_(); 
                        PDGSTRF2 (options, kk0, kk, thresh,
                                  Glu_persist, grid, Llu, U_diag_blk_send_req,
                                  tag_ub, stat, info);
                        pdgstrf2_timer += SuperLU_timer_() - ttt1; 

                        /* Process column *kcol+1* multicasts numeric
			   values of L(:,k+1) to process rows. */
                        look_id = kk0 % (1 + num_look_aheads);
                        send_req = send_reqs[look_id];
                        msgcnt = msgcnts[look_id];

                        if (lsub1) {
                            msgcnt[0] = lsub1[1] + BC_HEADER + lsub1[0] * LB_DESCRIPTOR;
                            msgcnt[1] = lsub1[1] * SuperSize (kk);
                        } else {
                            msgcnt[0] = 0;
                            msgcnt[1] = 0;
                        }

                        scp = &grid->rscp;  /* The scope of process row. */
                        for (pj = 0; pj < Pc; ++pj) {
                            if (ToSendR[lk][pj] != EMPTY) {
                                MPI_Isend (lsub1, msgcnt[0], mpi_int_t, pj,
                                           SLU_MPI_TAG (0, kk0), /* (4*kk0)%tag_ub */
                                           scp->comm, &send_req[pj]);
                                MPI_Isend (lusup1, msgcnt[1], MPI_DOUBLE, pj,
                                           SLU_MPI_TAG (1, kk0), /* (4*kk0+1)%tag_ub */
                                           scp->comm, &send_req[pj + Pc]);
                            }
                        }
                    }           /* for pj ... */
                }
            }
        }

        double tsch = SuperLU_timer_();

/************************************************************************/

#ifdef GPU_ACC

#include "dSchCompUdt-cuda.c"

#else 

/*#include "SchCompUdt--Phi-2Ddynamic-alt.c"*/
#include "dSchCompUdt-2Ddynamic.c"

#endif 
	/*uncomment following to compare against SuperLU 3.3 baseline*/
        /* #include "SchCompUdt--baseline.c"  */
/************************************************************************/
        
        NetSchurUpTimer += SuperLU_timer_()-tsch;

    }  /* for k0 = 0, ... */

    /* ##################################################################
       ** END MAIN LOOP: for k0 = ...
       ################################################################## */
    
    pdgstrfTimer= SuperLU_timer_()-pdgstrfTimer;

    /* updating total flops */
#if ( PRNTlevel>=1 )
    if (!iam) {
        printf("Time in scattering %lf \n",scatter_timer );
        printf("Time in dgemm %lf \n", gemm_timer );
        printf("Total time spent in schur update is    \t\t: %5.2lf seconds,\n",NetSchurUpTimer );
        printf("Total Time in Factorization            \t\t: %5.2lf seconds, \n", pdgstrfTimer);
        printf("Time (other GEMM and Scatter)          \t\t: %5.2lf seconds, \n", pdgstrfTimer-schur_flop_timer);
        printf("Total time spent in schur update when offload    \t\t: %5.2lf seconds,\n",CPUOffloadTimer );
    }
#endif
    
#if ( DEBUGlevel>=2 )
    for (i = 0; i < Pr * Pc; ++i) {
        if (iam == i) {
            dPrintLblocks(iam, nsupers, grid, Glu_persist, Llu);
            dPrintUblocks(iam, nsupers, grid, Glu_persist, Llu);
            printf ("(%d)\n", iam);
            PrintInt10 ("Recv", nsupers, Llu->ToRecv);
        }
        MPI_Barrier (grid->comm);
    }
#endif

    // printf("Debug : MPI buffers 1\n");

    /********************************************************
     * Free memory                                          *
     ********************************************************/

    if (Pr * Pc > 1) {
        SUPERLU_FREE (Lsub_buf_2[0]);   /* also free Lsub_buf_2[1] */
        SUPERLU_FREE (Lval_buf_2[0]);   /* also free Lval_buf_2[1] */
        if (Llu->bufmax[2] != 0)
            SUPERLU_FREE (Usub_buf_2[0]);
        if (Llu->bufmax[3] != 0)
            SUPERLU_FREE (Uval_buf_2[0]);
        if (U_diag_blk_send_req[myrow] != MPI_REQUEST_NULL) {
            /* wait for last Isend requests to complete, deallocate objects */
            for (krow = 0; krow < Pr; ++krow) {
                if (krow != myrow)
                    MPI_Wait (U_diag_blk_send_req + krow, &status);
            }
        }
        SUPERLU_FREE (U_diag_blk_send_req);
    }

    log_memory( -((Llu->bufmax[0] + Llu->bufmax[2]) * (num_look_aheads + 1) * iword +
		  (Llu->bufmax[1] + Llu->bufmax[3]) * (num_look_aheads + 1) * dword),
		stat );
    
    SUPERLU_FREE (Lsub_buf_2);
    SUPERLU_FREE (Lval_buf_2);
    SUPERLU_FREE (Usub_buf_2);
    SUPERLU_FREE (Uval_buf_2);
    SUPERLU_FREE (perm_c_supno);
    SUPERLU_FREE (perm_u);
#ifdef ISORT
    SUPERLU_FREE (iperm_u);
#endif
    SUPERLU_FREE (look_ahead);
    SUPERLU_FREE (factoredU);
    SUPERLU_FREE (factored);
    log_memory(-(6 * nsupers * iword), stat);


    for (i = 0; i <= num_look_aheads; i++) {
        SUPERLU_FREE (msgcnts[i]);
        SUPERLU_FREE (msgcntsU[i]);
    }
    SUPERLU_FREE (msgcnts);
    SUPERLU_FREE (msgcntsU);

    for (i = 0; i <= num_look_aheads; i++) {
        SUPERLU_FREE (send_reqs_u[i]);
        SUPERLU_FREE (recv_reqs_u[i]);
        SUPERLU_FREE (send_reqs[i]);
        SUPERLU_FREE (recv_reqs[i]);
    }

    SUPERLU_FREE (recv_reqs_u);
    SUPERLU_FREE (send_reqs_u);
    SUPERLU_FREE (recv_reqs);
    SUPERLU_FREE (send_reqs);

    // printf("Debug : MPI buffers 3\n");

#ifdef GPU_ACC
    checkCuda (cudaFreeHost (bigV));
    checkCuda (cudaFreeHost (bigU));
    cudaFree( (void*)dA ); /* Sherry added */
    cudaFree( (void*)dB );
    cudaFree( (void*)dC );
    SUPERLU_FREE( handle );
    SUPERLU_FREE( streams );
#else
    SUPERLU_FREE (bigV);
    SUPERLU_FREE (bigU);
#endif

    log_memory(-(bigv_size + bigu_size) * dword, stat);
    // printf("Debug : MPI buffers 5\n");

    SUPERLU_FREE (Llu->ujrow);
    SUPERLU_FREE (tempv2d);
    SUPERLU_FREE (indirect);
    SUPERLU_FREE (indirect2); /* Sherry added */
    SUPERLU_FREE (iuip);
    SUPERLU_FREE (ruip);

    ldt = sp_ienv_dist(3);
    log_memory( -(3 * ldt *ldt * dword + 2 * ldt * num_threads * iword
		  + 2 * k * iword), stat );

    /* Sherry added */
    SUPERLU_FREE(omp_loop_time);
    SUPERLU_FREE(full_u_cols);
    SUPERLU_FREE(blk_ldu);
    log_memory(-2 * ncb * dword, stat);

    SUPERLU_FREE(stream_end_col);
    SUPERLU_FREE(lookAheadFullRow);
    SUPERLU_FREE(lookAheadStRow);
    SUPERLU_FREE(lookAhead_lptr);
    SUPERLU_FREE(lookAhead_ib);

    SUPERLU_FREE(RemainFullRow);
    SUPERLU_FREE(RemainStRow);
    SUPERLU_FREE(Remain_lptr);
    SUPERLU_FREE(Remain_ib);
    SUPERLU_FREE(Remain_info);
    SUPERLU_FREE(lookAhead_L_buff);
    SUPERLU_FREE(Remain_L_buff);
    log_memory( -(4 * mrb * iword + mrb * sizeof(Remain_info_t) + 
		  ldt * ldt * (num_look_aheads + 1) * dword +
		  Llu->bufmax[1] * dword), stat );

    SUPERLU_FREE(Ublock_info);
    SUPERLU_FREE(Ublock_info_iukp);
    SUPERLU_FREE(Ublock_info_rukp);
    SUPERLU_FREE(Ublock_info_jb);


#if ( PROFlevel>=1 )
    TIC (t1);
#endif

    /* Prepare error message - find the smallesr index i that U(i,i)==0 */
    if ( *info == 0 ) *info = n + 1;
    MPI_Allreduce (info, &iinfo, 1, MPI_INT, MPI_MIN, grid->comm);
    if ( iinfo == n + 1 ) *info = 0;
    else *info = iinfo;

    // printf("test out\n");

#if ( PROFlevel>=1 )
    TOC (t2, t1);
    stat->utime[COMM] += t2;
    {
        float msg_vol_max, msg_vol_sum, msg_cnt_max, msg_cnt_sum;

        MPI_Reduce (&msg_cnt, &msg_cnt_sum,
                    1, MPI_FLOAT, MPI_SUM, 0, grid->comm);
        MPI_Reduce (&msg_cnt, &msg_cnt_max,
                    1, MPI_FLOAT, MPI_MAX, 0, grid->comm);
        MPI_Reduce (&msg_vol, &msg_vol_sum,
                    1, MPI_FLOAT, MPI_SUM, 0, grid->comm);
        MPI_Reduce (&msg_vol, &msg_vol_max,
                    1, MPI_FLOAT, MPI_MAX, 0, grid->comm);
        if (!iam) {
            printf ("\tPDGSTRF comm stat:"
                    "\tAvg\tMax\t\tAvg\tMax\n"
                    "\t\t\tCount:\t%.0f\t%.0f\tVol(MB)\t%.2f\t%.2f\n",
                    msg_cnt_sum / Pr / Pc, msg_cnt_max,
                    msg_vol_sum / Pr / Pc * 1e-6, msg_vol_max * 1e-6);
        }
    }
#endif

#if ( PRNTlevel==3 )
    MPI_Allreduce (&zero_msg, &iinfo, 1, MPI_INT, MPI_SUM, grid->comm);
    if (!iam)
        printf (".. # msg of zero size\t%d\n", iinfo);
    MPI_Allreduce (&total_msg, &iinfo, 1, MPI_INT, MPI_SUM, grid->comm);
    if (!iam)
        printf (".. # total msg\t%d\n", iinfo);
#endif

#if ( DEBUGlevel>=2 )
    for (i = 0; i < Pr * Pc; ++i) {
        if (iam == i) {
            dPrintLblocks (iam, nsupers, grid, Glu_persist, Llu);
            dPrintUblocks (iam, nsupers, grid, Glu_persist, Llu);
            printf ("(%d)\n", iam);
            PrintInt10 ("Recv", nsupers, Llu->ToRecv);
        }
        MPI_Barrier (grid->comm);
    }
#endif

#if ( DEBUGlevel>=3 )
    printf ("(%d) num_copy=%d, num_update=%d\n", iam, num_copy, num_update);
#endif
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC (iam, "Exit pdgstrf()");
#endif

    return 0;
} /* PDGSTRF */

