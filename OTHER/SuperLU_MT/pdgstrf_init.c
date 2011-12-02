#include "pdsp_defs.h"


void
pdgstrf_init(int nprocs, yes_no_t refact, int panel_size, int relax,
	     double diag_pivot_thresh, yes_no_t usepr, double drop_tol,
	     int *perm_c, int *perm_r, void *work, int lwork,
	     SuperMatrix *A, SuperMatrix *AC, 
	     pdgstrf_options_t *pdgstrf_options, Gstat_t *Gstat)
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
 * pdgstrf_init() initializes the option structure pdgstrf_options, using 
 * the user-input parameters. These options control how the factorization
 * will be performed by routine pdgstf().
 * 
 * In addition, it calls sp_colorder() to compute a postordered etree[], 
 * colcnt_h[] and super_part_h, and permute the columns of A using the 
 * permutation vector perm_c[]. See sp_colorder.c for details.
 *
 * Arguments
 * =========
 *
 * nprocs (input) int
 *        Number of processes used to perform LU factorization by pdgstrf().
 *
 * refact (input) yes_no_t
 *        Specifies whether we want to use perm_r from a previous factor.
 *
 * panel_size (input) int
 *        A panel consists of at most panel_size consecutive columns.
 *
 * relax  (input) int
 *        To control degree of relaxing supernodes. If the number
 *        of nodes (columns) in a subtree of the elimination tree is less
 *        than relax, this subtree is considered as one supernode,
 *        regardless of the row structures of those columns.
 *
 * diag_pivot_thresh (input) double
 *        Diagonal pivoting threshold. At step j of the Gaussian elimination,
 *        if abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
 *        use A_jj as pivot. 0 <= diag_pivot_thresh <= 1. The default
 *        value is 1, corresponding to partial pivoting.
 *
 * drop_tol (input) double (NOT IMPLEMENTED)
 *	  Drop tolerance parameter. At step j of the Gaussian elimination,
 *        if abs(A_ij)/(max_i abs(A_ij)) < drop_tol, drop entry A_ij.
 *        0 <= drop_tol <= 1. The default value of drop_tol is 0, 
 *        corresponding to not dropping any entry.
 *
 * perm_c (input) int*, dimension A->ncol
 *	  Column permutation vector, which defines the 
 *        permutation matrix Pc; perm_c[i] = j means column i of A is 
 *        in position j in A*Pc.
 *        When search for diagonal, perm_c[*] is applied to the
 *        row subscripts of A, so that diagonal threshold pivoting
 *        can find the diagonal of A, instead of that of A*Pc.
 *
 * perm_r (input/output) int*, dimension A->nrow
 *        Row permutation vector which defines the permutation matrix Pr,
 *        perm_r[i] = j means row i of A is in position j in Pr*A.
 *        If usepr = NO, perm_r is output argument;
 *        If usepr = YES, the pivoting routine will try to use the input
 *           perm_r, unless a certain threshold criterion is violated.
 *           In that case, perm_r is overwritten by a new permutation
 *           determined by partial pivoting or diagonal threshold pivoting.
 *
 * work   (input) void* of size lwork
 *        User-supplied work space and space for the output data structures.
 *        Not referenced if lwork = 0;
 *
 * lwork  (input) int
 *        Specifies the length of work array.
 *        = 0:  allocate space internally by system malloc;
 *        > 0:  use user-supplied work array of length lwork in bytes,
 *              returns error if space runs out.
 *        = -1: the routine guesses the amount of space needed without
 *              performing the factorization, and returns it in
 *              superlu_memusage->total_needed; no other side effects.
 *
 * A      (input) SuperMatrix*
 *        Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *        of linear equations is A->nrow. Currently, the type of A can be:
 *        Stype = NC or NCP; Dtype = _D; Mtype = GE. In the future, more
 *        general A can be handled.
 *
 * AC     (output) SuperMatrix*
 *        The resulting matrix after applied the column permutation
 *        perm_c[] to matrix A. The type of AC can be:
 *        Stype = NCP; Dtype = _D; Mtype = GE.
 *
 * pdgstrf_options (output) pdgstrf_options_t*
 *        The structure defines the parameters to control how the sparse
 *        LU factorization is performed, and will be input to pdgstrf().
 *
 * Gstat  (output) Gstat_t*
 *        Record the time used in sp_colorder phase.
 *
 */
    double t;

    pdgstrf_options->nprocs = nprocs;
    pdgstrf_options->refact = refact;
    pdgstrf_options->panel_size = panel_size;
    pdgstrf_options->relax = relax;
    pdgstrf_options->diag_pivot_thresh = diag_pivot_thresh;
    pdgstrf_options->usepr = usepr;
    pdgstrf_options->drop_tol = drop_tol;
    /* 
     * The following should be retained for repeated factorizations.
     */
    pdgstrf_options->perm_c = perm_c;
    pdgstrf_options->perm_r = perm_r;
    pdgstrf_options->work = work;
    pdgstrf_options->lwork = lwork;

    t = SuperLU_timer_();
    sp_colorder(A, perm_c, pdgstrf_options, AC);
    Gstat->utime[ETREE] = SuperLU_timer_() - t;

#if ( DEBUGlevel==1 )
    printf("** pdgstrf_init() called\n");
#endif
}

