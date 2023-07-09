/*! @file
 * \brief Utilities functions
 *
 * <pre>
 * -- Distributed SuperLU routine (version 2.3) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * February 1, 2003
 * Modified: March 31, 2013
 * </pre>
 */

#include <math.h>
#include "superlu_ddefs.h"

/*! \brief Deallocate the structure pointing to the actual storage of the matrix. */
void
Destroy_SuperMatrix_Store_dist(SuperMatrix *A)
{
    SUPERLU_FREE ( A->Store );
}

void
Destroy_CompCol_Matrix_dist(SuperMatrix *A)
{
    NCformat *Astore = A->Store;
    SUPERLU_FREE( Astore->rowind );
    SUPERLU_FREE( Astore->colptr );
    if ( Astore->nzval ) SUPERLU_FREE( Astore->nzval );
    SUPERLU_FREE( Astore );
}

void
Destroy_CompRowLoc_Matrix_dist(SuperMatrix *A)
{
    NRformat_loc *Astore = A->Store;
    SUPERLU_FREE( Astore->rowptr );
    SUPERLU_FREE( Astore->colind );
    SUPERLU_FREE( Astore->nzval );
    SUPERLU_FREE( Astore );
}

void
Destroy_CompRow_Matrix_dist(SuperMatrix *A)
{
    SUPERLU_FREE( ((NRformat *)A->Store)->rowptr );
    SUPERLU_FREE( ((NRformat *)A->Store)->colind );
    SUPERLU_FREE( ((NRformat *)A->Store)->nzval );
    SUPERLU_FREE( A->Store );
}

void
Destroy_SuperNode_Matrix_dist(SuperMatrix *A)
{
    SUPERLU_FREE ( ((SCformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((SCformat *)A->Store)->rowind_colptr );
    SUPERLU_FREE ( ((SCformat *)A->Store)->nzval );
    SUPERLU_FREE ( ((SCformat *)A->Store)->nzval_colptr );
    SUPERLU_FREE ( ((SCformat *)A->Store)->col_to_sup );
    SUPERLU_FREE ( ((SCformat *)A->Store)->sup_to_col );
    SUPERLU_FREE ( A->Store );
}

/*! \brief A is of type Stype==NCP */
void
Destroy_CompCol_Permuted_dist(SuperMatrix *A)
{
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colbeg );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colend );
    SUPERLU_FREE ( A->Store );
}

/*! \brief A is of type Stype==DN */
void
Destroy_Dense_Matrix_dist(SuperMatrix *A)
{
    DNformat* Astore = A->Store;
    SUPERLU_FREE (Astore->nzval);
    SUPERLU_FREE ( A->Store );
}

/*! \brief Destroy distributed L & U matrices. */
void
Destroy_LU(int_t n, gridinfo_t *grid, LUstruct_t *LUstruct)
{
    int_t i, nb, nsupers;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;

#if ( DEBUGlevel>=1 )
    int iam;
    MPI_Comm_rank( MPI_COMM_WORLD, &iam );
    CHECK_MALLOC(iam, "Enter Destroy_LU()");
#endif

    nsupers = Glu_persist->supno[n-1] + 1;

    nb = CEILING(nsupers, grid->npcol);
    for (i = 0; i < nb; ++i) 
	if ( Llu->Lrowind_bc_ptr[i] ) {
	    SUPERLU_FREE (Llu->Lrowind_bc_ptr[i]);
#ifdef GPU_ACC
	    checkCuda(cudaFreeHost(Llu->Lnzval_bc_ptr[i]));
#else
	    SUPERLU_FREE (Llu->Lnzval_bc_ptr[i]);
#endif
	}
    SUPERLU_FREE (Llu->Lrowind_bc_ptr);
    SUPERLU_FREE (Llu->Lnzval_bc_ptr);

    nb = CEILING(nsupers, grid->nprow);
    for (i = 0; i < nb; ++i)
	if ( Llu->Ufstnz_br_ptr[i] ) {
	    SUPERLU_FREE (Llu->Ufstnz_br_ptr[i]);
	    SUPERLU_FREE (Llu->Unzval_br_ptr[i]);
	}
    SUPERLU_FREE (Llu->Ufstnz_br_ptr);
    SUPERLU_FREE (Llu->Unzval_br_ptr);

    /* The following can be freed after factorization. */
    SUPERLU_FREE(Llu->ToRecv);
    SUPERLU_FREE(Llu->ToSendD);
    SUPERLU_FREE(Llu->ToSendR[0]);
    SUPERLU_FREE(Llu->ToSendR);

    /* The following can be freed only after iterative refinement. */
    SUPERLU_FREE(Llu->ilsum);
    SUPERLU_FREE(Llu->fmod);
    SUPERLU_FREE(Llu->fsendx_plist[0]);
    SUPERLU_FREE(Llu->fsendx_plist);
    SUPERLU_FREE(Llu->bmod);
    SUPERLU_FREE(Llu->bsendx_plist[0]);
    SUPERLU_FREE(Llu->bsendx_plist);
    SUPERLU_FREE(Llu->mod_bit);

    SUPERLU_FREE(Glu_persist->xsup);
    SUPERLU_FREE(Glu_persist->supno);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit Destroy_LU()");
#endif
}

/*! \brief Allocate storage in ScalePermstruct */
void ScalePermstructInit(const int_t m, const int_t n,
			 ScalePermstruct_t *ScalePermstruct)
{
    ScalePermstruct->DiagScale = NOEQUIL;
    if ( !(ScalePermstruct->perm_r = intMalloc_dist(m)) )
	ABORT("Malloc fails for perm_r[].");
    if ( !(ScalePermstruct->perm_c = intMalloc_dist(n)) )
	ABORT("Malloc fails for perm_c[].");
}

/*! \brief Deallocate ScalePermstruct */
void ScalePermstructFree(ScalePermstruct_t *ScalePermstruct)
{
    SUPERLU_FREE(ScalePermstruct->perm_r);
    SUPERLU_FREE(ScalePermstruct->perm_c);
    switch ( ScalePermstruct->DiagScale ) {
      case ROW:
	SUPERLU_FREE(ScalePermstruct->R);
	break;
      case COL:
	SUPERLU_FREE(ScalePermstruct->C);
	break;
      case BOTH:
	SUPERLU_FREE(ScalePermstruct->R);
	SUPERLU_FREE(ScalePermstruct->C);
	break;
    }
}

/*! \brief Allocate storage in LUstruct */
void LUstructInit(const int_t n, LUstruct_t *LUstruct)
{
    if ( !(LUstruct->etree = intMalloc_dist(n)) )
	ABORT("Malloc fails for etree[].");
    if ( !(LUstruct->Glu_persist = (Glu_persist_t *)
	   SUPERLU_MALLOC(sizeof(Glu_persist_t))) )
	ABORT("Malloc fails for Glu_persist_t.");
    if ( !(LUstruct->Llu = (LocalLU_t *)
	   SUPERLU_MALLOC(sizeof(LocalLU_t))) )
	ABORT("Malloc fails for LocalLU_t.");
}

/*! \brief Deallocate LUstruct */
void LUstructFree(LUstruct_t *LUstruct)
{
#if ( DEBUGlevel>=1 )
    int iam;
    MPI_Comm_rank( MPI_COMM_WORLD, &iam );
    CHECK_MALLOC(iam, "Enter LUstructFree()");
#endif

    SUPERLU_FREE(LUstruct->etree);
    SUPERLU_FREE(LUstruct->Glu_persist);
    SUPERLU_FREE(LUstruct->Llu);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit LUstructFree()");
#endif
}

/*! \brief
 *
 * <pre>
 * Count the total number of nonzeros in factors L and U,  and in the 
 * symmetrically reduced L. 
 * </pre>
 */
void
countnz_dist(const int_t n, int_t *xprune,
	     long long int *nnzL, long long int *nnzU, 
	     Glu_persist_t *Glu_persist, Glu_freeable_t *Glu_freeable)
{
    int_t  fnz, fsupc, i, j, nsuper;
    int_t  jlen, irep;
    long long int nnzL0;
    int_t  *supno, *xsup, *xlsub, *xusub, *usub;

    supno  = Glu_persist->supno;
    xsup   = Glu_persist->xsup;
    xlsub  = Glu_freeable->xlsub;
    xusub  = Glu_freeable->xusub;
    usub   = Glu_freeable->usub;
    *nnzL  = 0;
    *nnzU  = 0;
    nnzL0  = 0;
    nsuper = supno[n];

    if ( n <= 0 ) return;

    /* 
     * For each supernode in L.
     */
    for (i = 0; i <= nsuper; i++) {
	fsupc = xsup[i];
	jlen = xlsub[fsupc+1] - xlsub[fsupc];

	for (j = fsupc; j < xsup[i+1]; j++) {
	    *nnzL += jlen;
	    *nnzU += j - fsupc + 1;
	    jlen--;
	}
	irep = xsup[i+1] - 1;
	nnzL0 += xprune[irep] - xlsub[irep];
    }
    
    /* printf("\tNo of nonzeros in symm-reduced L = %ld\n", nnzL0);*/
    
    /* For each column in U. */
    for (j = 0; j < n; ++j) {
	for (i = xusub[j]; i < xusub[j+1]; ++i) {
	    fnz = usub[i];
	    fsupc = xsup[supno[fnz]+1];
	    *nnzU += fsupc - fnz;
	}
    }
}



/*! \brief
 *
 * <pre>
 * Fix up the data storage lsub for L-subscripts. It removes the subscript
 * sets for structural pruning,	and applies permuation to the remaining
 * subscripts.
 * </pre>
 */
long long int
fixupL_dist(const int_t n, const int_t *perm_r, 
	    Glu_persist_t *Glu_persist, Glu_freeable_t *Glu_freeable)
{
    register int_t nsuper, fsupc, nextl, i, j, k, jstrt;
    register long long int lsub_size;
    int_t          *xsup, *lsub, *xlsub;

    if ( n <= 1 ) return 0;

    xsup   = Glu_persist->xsup;
    lsub   = Glu_freeable->lsub;
    xlsub  = Glu_freeable->xlsub;
    nextl  = 0;
    nsuper = (Glu_persist->supno)[n];
    lsub_size = xlsub[n];
    
    /* 
     * For each supernode ...
     */
    for (i = 0; i <= nsuper; i++) {
	fsupc = xsup[i];
	jstrt = xlsub[fsupc];
	xlsub[fsupc] = nextl;
	for (j = jstrt; j < xlsub[fsupc+1]; j++) {
	    lsub[nextl] = perm_r[lsub[j]]; /* Now indexed into P*A */
	    nextl++;
  	}
	for (k = fsupc+1; k < xsup[i+1]; k++) 
	    	xlsub[k] = nextl;	/* Other columns in supernode i */

    }

    xlsub[n] = nextl;
    return lsub_size;
}

/*! \brief Set the default values for the options argument.
 */
void set_default_options_dist(superlu_options_t *options)
{
    options->Fact              = DOFACT;
    options->Equil             = YES;
    options->ParSymbFact       = NO;
    options->ColPerm           = METIS_AT_PLUS_A;
    options->RowPerm           = LargeDiag;
    options->ReplaceTinyPivot  = YES;
    options->IterRefine        = SLU_DOUBLE;
    options->Trans             = NOTRANS;
    options->SolveInitialized  = NO;
    options->RefineInitialized = NO;
    options->PrintStat         = YES;
    options->num_lookaheads    = 10;
    options->lookahead_etree   = NO;
    options->SymPattern        = NO;
}

/*! \brief Print the options setting.
 */
void print_options_dist(superlu_options_t *options)
{
    if ( options->PrintStat == NO ) return;

    printf("**************************************************\n");
    printf(".. options:\n");
    printf("**    Fact             : %4d\n", options->Fact);
    printf("**    Equil            : %4d\n", options->Equil);
    printf("**    ParSymbFact      : %4d\n", options->ParSymbFact);
    printf("**    ColPerm          : %4d\n", options->ColPerm);
    printf("**    RowPerm          : %4d\n", options->RowPerm);
    printf("**    ReplaceTinyPivot : %4d\n", options->ReplaceTinyPivot);
    printf("**    IterRefine       : %4d\n", options->IterRefine);
    printf("**    Trans            : %4d\n", options->Trans);
    printf("**    num_lookaheads   : %4d\n", options->num_lookaheads);
    printf("**    SymPattern       : %4d\n", options->SymPattern);
    printf("**    lookahead_etree  : %4d\n", options->lookahead_etree);
    printf("**************************************************\n");
}

/*! \brief Print the blocking parameters.
 */
void print_sp_ienv_dist(superlu_options_t *options)
{
    if ( options->PrintStat == NO ) return;

    printf("**************************************************\n");
    printf(".. blocking parameters from sp_ienv():\n");
    printf("**    relaxation           : " IFMT "\n", sp_ienv_dist(2));
    printf("**    max supernode        : " IFMT "\n", sp_ienv_dist(3));
    printf("**    estimated fill ratio : " IFMT "\n", sp_ienv_dist(6));
    printf("**************************************************\n");
}


/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *   Set up the communication pattern for the triangular solution.
 * 
 * Arguments
 * =========
 *
 * n      (input) int (global)
 *        The dimension of the linear system.
 *
 * m_loc  (input) int (local)
 *        The local row dimension of the distributed input matrix.
 *
 * nrhs   (input) int (global)
 *        Number of right-hand sides.
 *
 * fst_row (input) int (global)
 *        The row number of matrix B's first row in the global matrix.
 *
 * perm_r (input) int* (global)
 *        The row permutation vector.
 *
 * perm_c (input) int* (global)
 *        The column permutation vector.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 * </pre>
 */
int_t
pxgstrs_init(int_t n, int_t m_loc, int_t nrhs, int_t fst_row,
	     int_t perm_r[], int_t perm_c[], gridinfo_t *grid,
	     Glu_persist_t *Glu_persist, SOLVEstruct_t *SOLVEstruct)
{

    int *SendCnt, *SendCnt_nrhs, *RecvCnt, *RecvCnt_nrhs;
    int *sdispls, *sdispls_nrhs, *rdispls, *rdispls_nrhs;
    int *itemp, *ptr_to_ibuf, *ptr_to_dbuf;
    int_t *row_to_proc;
    int_t i, gbi, k, l, num_diag_procs, *diag_procs;
    int_t irow, q, knsupc, nsupers, *xsup, *supno;
    int   iam, p, pkk, procs;
    pxgstrs_comm_t *gstrs_comm;

    procs = grid->nprow * grid->npcol;
    iam = grid->iam;
    gstrs_comm = SOLVEstruct->gstrs_comm;
    xsup = Glu_persist->xsup;
    supno = Glu_persist->supno;
    nsupers = Glu_persist->supno[n-1] + 1;
    row_to_proc = SOLVEstruct->row_to_proc;

    /* ------------------------------------------------------------
       SET UP COMMUNICATION PATTERN FOR ReDistribute_B_to_X.
       ------------------------------------------------------------*/
    if ( !(itemp = SUPERLU_MALLOC(8*procs * sizeof(int))) )
        ABORT("Malloc fails for B_to_X_itemp[].");
    SendCnt      = itemp;
    SendCnt_nrhs = itemp +   procs;
    RecvCnt      = itemp + 2*procs;
    RecvCnt_nrhs = itemp + 3*procs;
    sdispls      = itemp + 4*procs;
    sdispls_nrhs = itemp + 5*procs;
    rdispls      = itemp + 6*procs;
    rdispls_nrhs = itemp + 7*procs;

    /* Count the number of elements to be sent to each diagonal process.*/
    for (p = 0; p < procs; ++p) SendCnt[p] = 0;
    for (i = 0, l = fst_row; i < m_loc; ++i, ++l) {
        irow = perm_c[perm_r[l]]; /* Row number in Pc*Pr*B */
	gbi = BlockNum( irow );
	p = PNUM( PROW(gbi,grid), PCOL(gbi,grid), grid ); /* Diagonal process */
	++SendCnt[p];
    }
  
    /* Set up the displacements for alltoall. */
    MPI_Alltoall(SendCnt, 1, MPI_INT, RecvCnt, 1, MPI_INT, grid->comm);
    sdispls[0] = rdispls[0] = 0;
    for (p = 1; p < procs; ++p) {
        sdispls[p] = sdispls[p-1] + SendCnt[p-1];
        rdispls[p] = rdispls[p-1] + RecvCnt[p-1];
    }
    for (p = 0; p < procs; ++p) {
        SendCnt_nrhs[p] = SendCnt[p] * nrhs;
	sdispls_nrhs[p] = sdispls[p] * nrhs;
        RecvCnt_nrhs[p] = RecvCnt[p] * nrhs;
	rdispls_nrhs[p] = rdispls[p] * nrhs;
    }

    /* This is saved for repeated solves, and is freed in pxgstrs_finalize().*/
    gstrs_comm->B_to_X_SendCnt = SendCnt;

    /* ------------------------------------------------------------
       SET UP COMMUNICATION PATTERN FOR ReDistribute_X_to_B.
       ------------------------------------------------------------*/
    /* This is freed in pxgstrs_finalize(). */
    if ( !(itemp = SUPERLU_MALLOC(8*procs * sizeof(int))) )
        ABORT("Malloc fails for X_to_B_itemp[].");
    SendCnt      = itemp;
    SendCnt_nrhs = itemp +   procs;
    RecvCnt      = itemp + 2*procs;
    RecvCnt_nrhs = itemp + 3*procs;
    sdispls      = itemp + 4*procs;
    sdispls_nrhs = itemp + 5*procs;
    rdispls      = itemp + 6*procs;
    rdispls_nrhs = itemp + 7*procs;

    /* Count the number of X entries to be sent to each process.*/
    for (p = 0; p < procs; ++p) SendCnt[p] = 0;
    num_diag_procs = SOLVEstruct->num_diag_procs;
    diag_procs = SOLVEstruct->diag_procs;

    for (p = 0; p < num_diag_procs; ++p) { /* for all diagonal processes */
	pkk = diag_procs[p];
	if ( iam == pkk ) {
	    for (k = p; k < nsupers; k += num_diag_procs) {
		knsupc = SuperSize( k );
		irow = FstBlockC( k );
		for (i = 0; i < knsupc; ++i) {
#if 0
		    q = row_to_proc[inv_perm_c[irow]];
#else
		    q = row_to_proc[irow];
#endif
		    ++SendCnt[q];
		    ++irow;
		}
	    }
	}
    }

    MPI_Alltoall(SendCnt, 1, MPI_INT, RecvCnt, 1, MPI_INT, grid->comm);
    sdispls[0] = rdispls[0] = 0;
    sdispls_nrhs[0] = rdispls_nrhs[0] = 0;
    SendCnt_nrhs[0] = SendCnt[0] * nrhs;
    RecvCnt_nrhs[0] = RecvCnt[0] * nrhs;
    for (p = 1; p < procs; ++p) {
        sdispls[p] = sdispls[p-1] + SendCnt[p-1];
        rdispls[p] = rdispls[p-1] + RecvCnt[p-1];
        sdispls_nrhs[p] = sdispls[p] * nrhs;
        rdispls_nrhs[p] = rdispls[p] * nrhs;
	SendCnt_nrhs[p] = SendCnt[p] * nrhs;
	RecvCnt_nrhs[p] = RecvCnt[p] * nrhs;
    }

    /* This is saved for repeated solves, and is freed in pxgstrs_finalize().*/
    gstrs_comm->X_to_B_SendCnt = SendCnt;

    if ( !(ptr_to_ibuf = SUPERLU_MALLOC(2*procs * sizeof(int))) )
        ABORT("Malloc fails for ptr_to_ibuf[].");
    gstrs_comm->ptr_to_ibuf = ptr_to_ibuf;
    gstrs_comm->ptr_to_dbuf = ptr_to_ibuf + procs;

    return 0;
} /* PXGSTRS_INIT */


void pxgstrs_finalize(pxgstrs_comm_t *gstrs_comm)
{
    SUPERLU_FREE(gstrs_comm->B_to_X_SendCnt);
    SUPERLU_FREE(gstrs_comm->X_to_B_SendCnt);
    SUPERLU_FREE(gstrs_comm->ptr_to_ibuf);
    SUPERLU_FREE(gstrs_comm);
}


/*! \brief Diagnostic print of segment info after panel_dfs().
 */
void print_panel_seg_dist(int_t n, int_t w, int_t jcol, int_t nseg, 
			  int_t *segrep, int_t *repfnz)
{
    int_t j, k;
    
    for (j = jcol; j < jcol+w; j++) {
	printf("\tcol " IFMT ":\n", j);
	for (k = 0; k < nseg; k++)
	    printf("\t\tseg " IFMT ", segrep " IFMT ", repfnz " IFMT "\n", k, 
			segrep[k], repfnz[(j-jcol)*n + segrep[k]]);
    }

}

void
PStatInit(SuperLUStat_t *stat)
{
    register int_t i;

    if ( !(stat->utime = SUPERLU_MALLOC(NPHASES*sizeof(double))) )
	ABORT("Malloc fails for stat->utime[]");
    if ( !(stat->ops = (flops_t *) SUPERLU_MALLOC(NPHASES * sizeof(flops_t))) )
	ABORT("SUPERLU_MALLOC fails for stat->ops[]");
    for (i = 0; i < NPHASES; ++i) {
        stat->utime[i] = 0.;
        stat->ops[i] = 0.;
    }
    stat->TinyPivots = stat->RefineSteps = 0;
}

void
PStatPrint(superlu_options_t *options, SuperLUStat_t *stat, gridinfo_t *grid)
{
    double  *utime = stat->utime;
    flops_t *ops = stat->ops;
    int_t   iam = grid->iam;
    flops_t flopcnt, factflop, solveflop;

    if ( options->PrintStat == NO ) return;

    if ( !iam && options->Fact != FACTORED ) {
	printf("**************************************************\n");
	printf("**** Time (seconds) ****\n");

        if ( options->Equil != NO )
	    printf("\tEQUIL time         %8.2f\n", utime[EQUIL]);
	if ( options->RowPerm != NOROWPERM )
	    printf("\tROWPERM time       %8.2f\n", utime[ROWPERM]);
	if ( options->ColPerm != NATURAL )
	    printf("\tCOLPERM time       %8.2f\n", utime[COLPERM]);
        printf("\tSYMBFACT time      %8.2f\n", utime[SYMBFAC]);
	printf("\tDISTRIBUTE time    %8.2f\n", utime[DIST]);

    }

    MPI_Reduce(&ops[FACT], &flopcnt, 1, MPI_FLOAT, MPI_SUM,
	       0, grid->comm);
    factflop = flopcnt;
    if ( !iam && options->Fact != FACTORED ) {
	printf("\tFACTOR time        %8.2f\n", utime[FACT]);
	if ( utime[FACT] != 0.0 )
	    printf("\tFactor flops\t%e\tMflops \t%8.2f\n",
		   flopcnt,
		   flopcnt*1e-6/utime[FACT]);
    }
	
    MPI_Reduce(&ops[SOLVE], &flopcnt, 1, MPI_FLOAT, MPI_SUM, 
	       0, grid->comm);
    solveflop = flopcnt;
    if ( !iam ) {
	printf("\tSOLVE time         %8.2f\n", utime[SOLVE]);
	if ( utime[SOLVE] != 0.0 )
	    printf("\tSolve flops\t%e\tMflops \t%8.2f\n",
		   flopcnt,
		   flopcnt*1e-6/utime[SOLVE]);
	if ( options->IterRefine != NOREFINE ) {
	    printf("\tREFINEMENT time    %8.2f\tSteps%8d\n\n",
		   utime[REFINE], stat->RefineSteps);
	}
	printf("**************************************************\n");
    }

#if ( PROFlevel>=1 )
    fflush(stdout);
    MPI_Barrier( grid->comm );

    {
	int_t i, P = grid->nprow*grid->npcol;
	flops_t b, maxflop;
	if ( !iam ) printf("\n.. FACT time breakdown:\tcomm\ttotal\n");
	for (i = 0; i < P; ++i) {
	    if ( iam == i) {
		printf("\t\t(%d)%8.2f%8.2f\n", iam, utime[COMM], utime[FACT]);
		fflush(stdout);
	    }
	    MPI_Barrier( grid->comm );
	}
	if ( !iam ) printf("\n.. FACT ops distribution:\n");
	for (i = 0; i < P; ++i) {
	    if ( iam == i ) {
		printf("\t\t(%d)\t%e\n", iam, ops[FACT]);
		fflush(stdout);
	    }
	    MPI_Barrier( grid->comm );
	}
	MPI_Reduce(&ops[FACT], &maxflop, 1, MPI_FLOAT, MPI_MAX, 0, grid->comm);
	if ( !iam ) {
	    b = factflop/P/maxflop;
	    printf("\tFACT load balance: %.2f\n", b);
	}
	if ( !iam ) printf("\n.. SOLVE ops distribution:\n");
	for (i = 0; i < P; ++i) {
	    if ( iam == i ) {
		printf("\t\t%d\t%e\n", iam, ops[SOLVE]);
		fflush(stdout);
	    }
	    MPI_Barrier( grid->comm );
	}
	MPI_Reduce(&ops[SOLVE], &maxflop, 1, MPI_FLOAT, MPI_MAX, 0,grid->comm);
	if ( !iam ) {
	    b = solveflop/P/maxflop;
	    printf("\tSOLVE load balance: %.2f\n", b);
	}
    }
#endif

/*  if ( !iam ) fflush(stdout);  CRASH THE SYSTEM pierre.  */
}

void
PStatFree(SuperLUStat_t *stat)
{
    SUPERLU_FREE(stat->utime);
    SUPERLU_FREE(stat->ops);
}

/*! \brief Fills an integer array with a given value.
 */
void ifill_dist(int_t *a, int_t alen, int_t ival)
{
    register int_t i;
    for (i = 0; i < alen; i++) a[i] = ival;
}


void
get_diag_procs(int_t n, Glu_persist_t *Glu_persist, gridinfo_t *grid,
	       int_t *num_diag_procs, int_t **diag_procs, int_t **diag_len)
{
    int_t i, j, k, knsupc, nprow, npcol, nsupers, pkk;
    int_t *xsup;

    i = j = *num_diag_procs = pkk = 0;
    nprow = grid->nprow;
    npcol = grid->npcol;
    nsupers = Glu_persist->supno[n-1] + 1;
    xsup = Glu_persist->xsup;

    do {
	++(*num_diag_procs);
	i = (++i) % nprow;
	j = (++j) % npcol;
	pkk = PNUM( i, j, grid );
    } while ( pkk != 0 ); /* Until wrap back to process 0 */
    if ( !(*diag_procs = intMalloc_dist(*num_diag_procs)) )
	ABORT("Malloc fails for diag_procs[]");
    if ( !(*diag_len = intCalloc_dist(*num_diag_procs)) )
	ABORT("Calloc fails for diag_len[]");
    for (i = j = k = 0; k < *num_diag_procs; ++k) {
	pkk = PNUM( i, j, grid );
	(*diag_procs)[k] = pkk;
	i = (++i) % nprow;
	j = (++j) % npcol;
    }
    for (k = 0; k < nsupers; ++k) {
	knsupc = SuperSize( k );
	i = k % *num_diag_procs;
	(*diag_len)[i] += knsupc;
    }
}


/*! \brief Get the statistics of the supernodes 
 */
#define NBUCKS 10
static 	int_t	max_sup_size;

void super_stats_dist(int_t nsuper, int_t *xsup)
{
    register int_t nsup1 = 0;
    int_t          i, isize, whichb, bl, bh;
    int_t          bucket[NBUCKS];

    max_sup_size = 0;

    for (i = 0; i <= nsuper; i++) {
	isize = xsup[i+1] - xsup[i];
	if ( isize == 1 ) nsup1++;
	if ( max_sup_size < isize ) max_sup_size = isize;	
    }

    printf("    Supernode statistics:\n\tno of super = " IFMT "\n", nsuper+1);
    printf("\tmax supernode size = " IFMT "\n", max_sup_size);
    printf("\tno of size 1 supernodes = " IFMT "\n", nsup1);

    /* Histogram of the supernode sizes */
    ifill_dist (bucket, NBUCKS, 0);

    for (i = 0; i <= nsuper; i++) {
        isize = xsup[i+1] - xsup[i];
        whichb = (float) isize / max_sup_size * NBUCKS;
        if (whichb >= NBUCKS) whichb = NBUCKS - 1;
        bucket[whichb]++;
    }
    
    printf("\tHistogram of supernode sizes:\n");
    for (i = 0; i < NBUCKS; i++) {
        bl = (float) i * max_sup_size / NBUCKS;
        bh = (float) (i+1) * max_sup_size / NBUCKS;
        printf("\tsnode: " IFMT "-" IFMT "\t\t" IFMT "\n", bl+1, bh, bucket[i]);
    }

}


/*! \brief Check whether repfnz[] == EMPTY after reset.
 */
void check_repfnz_dist(int_t n, int_t w, int_t jcol, int_t *repfnz)
{
    int_t jj, k;

    for (jj = jcol; jj < jcol+w; jj++) 
	for (k = 0; k < n; k++)
	    if ( repfnz[(jj-jcol)*n + k] != EMPTY ) {
		fprintf(stderr, "col " IFMT ", repfnz_col[" IFMT "] = " IFMT "\n",
			jj, k, repfnz[(jj-jcol)*n + k]);
		ABORT("check_repfnz_dist");
	    }
}

void PrintInt10(char *name, int_t len, int_t *x)
{
    register int_t i;
    
    printf("%10s:", name);
    for (i = 0; i < len; ++i) {
	if ( i % 10 == 0 ) printf("\n\t[" IFMT "-" IFMT "]", i, i+9);
	printf(IFMT, x[i]);
    }
    printf("\n");
}

void PrintInt32(char *name, int len, int *x)
{
    register int i;
    
    printf("%10s:", name);
    for (i = 0; i < len; ++i) {
	if ( i % 10 == 0 ) printf("\n\t[%2d-%2d]", i, i+9);
	printf("%6d", x[i]);
    }
    printf("\n");
}

int file_PrintInt10(FILE *fp, char *name, int_t len, int_t *x)
{
    register int_t i;
    
    fprintf(fp, "%10s:", name);
    for (i = 0; i < len; ++i) {
	if ( i % 10 == 0 ) fprintf(fp, "\n\t[" IFMT "-" IFMT "]", i, i+9);
	fprintf(fp, IFMT, x[i]);
    }
    fprintf(fp, "\n");
    return 0;
}

int file_PrintInt32(FILE *fp, char *name, int len, int *x)
{
    register int i;
    
    fprintf(fp, "%10s:", name);
    for (i = 0; i < len; ++i) {
	if ( i % 10 == 0 ) fprintf(fp, "\n\t[%2d-%2d]", i, i+9);
	fprintf(fp, "%6d", x[i]);
    }
    fprintf(fp, "\n");
    return 0;
}

int_t
CheckZeroDiagonal(int_t n, int_t *rowind, int_t *colbeg, int_t *colcnt)
{
    register int_t i, j, zd, numzd = 0;

    for (j = 0; j < n; ++j) {
	zd = 0;
	for (i = colbeg[j]; i < colbeg[j]+colcnt[j]; ++i) {
	    /*if ( iperm[rowind[i]] == j ) zd = 1;*/
	    if ( rowind[i] == j ) { zd = 1; break; }
	}
	if ( zd == 0 ) {
#if ( PRNTlevel>=2 )
	    printf(".. Diagonal of column %d is zero.\n", j);
#endif
	    ++numzd;
	}
    }

    return numzd;
}


/* --------------------------------------------------------------------------- */
void isort(int_t N, int_t *ARRAY1, int_t *ARRAY2)
{
/*
 * Purpose
 * =======
 * Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing
 * order of ARRAY1.
 *
 * Arguments
 * =========
 * N       (input) INTEGER
 *          On entry, specifies the size of the arrays.
 *
 * ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
 *          On entry, contains the array to be sorted.
 *          On exit, contains the sorted array.
 *
 * ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
 *          On entry, contains the array to be sorted.
 *          On exit, contains the sorted array.
 */
  int_t IGAP, I, J;
  int_t TEMP;
  IGAP = N / 2;
  while (IGAP > 0) {
    for (I = IGAP; I < N; I++) {
    J = I - IGAP;
    while (J >= 0) {
      if (ARRAY1[J] > ARRAY1[J + IGAP]) {
        TEMP = ARRAY1[J];
        ARRAY1[J] = ARRAY1[J + IGAP];
        ARRAY1[J + IGAP] = TEMP;
        TEMP = ARRAY2[J];
        ARRAY2[J] = ARRAY2[J + IGAP];
        ARRAY2[J + IGAP] = TEMP;
        J = J - IGAP;
      } else {
        break;
      }
    }
  }
    IGAP = IGAP / 2;
  }
}


void isort1(int_t N, int_t *ARRAY)
{
/*
 * Purpose
 * =======
 * Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing
 * order of ARRAY1.
 *
 * Arguments
 * =========
 * N       (input) INTEGER
 *          On entry, specifies the size of the arrays.
 *
 * ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
 *          On entry, contains the array to be sorted.
 *          On exit, contains the sorted array.
 *
 * ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
 *          On entry, contains the array to be sorted.
 *          On exit, contains the sorted array.
 */
  int_t IGAP, I, J;
  int_t TEMP;
  IGAP = N / 2;
  while (IGAP > 0) {
  for (I = IGAP; I < N; I++) {
    J = I - IGAP;
    while (J >= 0) {
      if (ARRAY[J] > ARRAY[J + IGAP]) {
        TEMP = ARRAY[J];
        ARRAY[J] = ARRAY[J + IGAP];
        ARRAY[J + IGAP] = TEMP;
        J = J - IGAP;
      } else {
        break;
      }
    }
  }
    IGAP = IGAP / 2;
  }
}

void log_memory(long long cur_bytes, SuperLUStat_t *stat) {
    stat->current_buffer += (float) cur_bytes;
    if (cur_bytes > 0) {
	stat->peak_buffer = 
	    SUPERLU_MAX(stat->peak_buffer, stat->current_buffer);
    }
}

void print_memorylog(SuperLUStat_t *stat, char *msg) {
    printf("__ %s (MB):\n\tcurrent_buffer : %8.2f\tpeak_buffer : %8.2f\n",
	   msg, stat->current_buffer, stat->peak_buffer);
}
