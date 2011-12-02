/*! @file
 * \brief Definitions which are precision-neutral
 *
 * <pre>
 * -- Distributed SuperLU routine (version 2.2) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * November 1, 2007
 * Feburary 20, 2008
 * </pre>
 */

#ifndef __SUPERLU_DEFS /* allow multiple inclusions */
#define __SUPERLU_DEFS

/*
 * File name:	superlu_defs.h
 * Purpose:     Definitions which are precision-neutral
 */
#ifdef _CRAY
#include <fortran.h>
#include <string.h>
#endif
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>


/* Define my integer size int_t */
#ifdef _CRAY
typedef short int_t;
/*#undef int   Revert back to int of default size. */
#define mpi_int_t   MPI_SHORT
#elif defined (_LONGINT)
typedef long long int int_t;
#define mpi_int_t   MPI_LONG
#else /* Default */
typedef int int_t;
#define mpi_int_t   MPI_INT
#endif

#include "superlu_enum_consts.h"
#include "Cnames.h"
#include "supermatrix.h"
#include "util_dist.h"
#include "psymbfact.h"


/***********************************************************************
 * Constants
 ***********************************************************************/
/* 
 * For each block column of L, the index[] array contains both the row 
 * subscripts and the integers describing the size of the blocks.
 * The organization of index[] looks like:
 *
 *     [ BLOCK COLUMN HEADER (size BC_HEADER)
 *           number of blocks 
 *           number of row subscripts, i.e., LDA of nzval[]
 *       BLOCK 0                                        <----
 *           BLOCK DESCRIPTOR (of size LB_DESCRIPTOR)  |
 *               block number (global)                      |
 *               number of full rows in the block           |
 *           actual row subscripts                          |
 *       BLOCK 1                                            | Repeat ...
 *           BLOCK DESCRIPTOR                               | number of blocks
 *               block number (global)                      | 
 *               number of full rows in the block           |
 *           actual row subscripts                          |
 *       .                                                  |
 *       .                                                  |
 *       .                                              <----
 *     ]
 *
 * For each block row of U, the organization of index[] looks like:
 *
 *     [ BLOCK ROW HEADER (of size BR_HEADER)
 *           number of blocks 
 *           number of entries in nzval[]
 *           number of entries in index[]
 *       BLOCK 0                                        <----
 *           BLOCK DESCRIPTOR (of size UB_DESCRIPTOR)  |
 *               block number (global)                      |
 *               number of nonzeros in the block            |
 *           actual fstnz subscripts                        |
 *       BLOCK 1                                            | Repeat ...
 *           BLOCK DESCRIPTOR                               | number of blocks
 *               block number (global)                      |
 *               number of nonzeros in the block            |
 *           actual fstnz subscripts                        |
 *       .                                                  |
 *       .                                                  |
 *       .                                              <----
 *     ]
 *
 */
#define BC_HEADER      2
#define LB_DESCRIPTOR  2
#define BR_HEADER      3
#define UB_DESCRIPTOR  2
#define NBUFFERS       5

/*
 * Communication tags
 */
    /* For numeric factorization. */
#define NTAGS    10000
#define UjROW    10
#define UkSUB    11
#define UkVAL    12
#define LkSUB    13
#define LkVAL    14
#define LkkDIAG  15
    /* For triangular solves. */
#define XK_H     2  /* The header preceeding each X block. */
#define LSUM_H   2  /* The header preceeding each MOD block. */
#define GSUM     20 
#define Xk       21
#define Yk       22
#define LSUM     23

/* 
 * Communication scopes
 */
#define COMM_ALL      100
#define COMM_COLUMN   101
#define COMM_ROW      102

/*
 * Matrix distribution for sparse matrix-vector multiplication
 */
#define SUPER_LINEAR     11
#define SUPER_BLOCK      12

/*
 * No of marker arrays used in the symbolic factorization, each of size n
 */
#define NO_MARKER     3



/***********************************************************************
 * Macros
 ***********************************************************************/
#define IAM(comm)    { int rank; MPI_Comm_rank ( comm, &rank ); rank};
#define MYROW(iam,grid) ( (iam) / grid->npcol )
#define MYCOL(iam,grid) ( (iam) % grid->npcol )
#define BlockNum(i)     ( supno[i] )
#define FstBlockC(bnum) ( xsup[bnum] )
#define SuperSize(bnum) ( xsup[bnum+1]-xsup[bnum] )
#define LBi(bnum,grid)  ( (bnum)/grid->nprow )/* Global to local block rowwise */
#define LBj(bnum,grid)  ( (bnum)/grid->npcol )/* Global to local block columnwise*/
#define PROW(bnum,grid) ( (bnum) % grid->nprow )
#define PCOL(bnum,grid) ( (bnum) % grid->npcol )
#define PNUM(i,j,grid)  ( (i)*grid->npcol + j ) /* Process number at coord(i,j) */
#define CEILING(a,b)    ( ((a)%(b)) ? ((a)/(b) + 1) : ((a)/(b)) )
    /* For triangular solves */
#define RHS_ITERATE(i)                    \
        for (i = 0; i < nrhs; ++i)
#define X_BLK(i)                          \
        ilsum[i] * nrhs + (i+1) * XK_H
#define LSUM_BLK(i)                       \
        ilsum[i] * nrhs + (i+1) * LSUM_H

#define SuperLU_timer_  SuperLU_timer_dist_
#define LOG2(x)   (log10((double) x) / log10(2.0))


#if ( VAMPIR>=1 ) 
#define VT_TRACEON    VT_traceon()
#define VT_TRACEOFF   VT_traceoff()
#else
#define VT_TRACEON 
#define VT_TRACEOFF
#endif


/***********************************************************************
 * New data types
 ***********************************************************************/

/* 
 *   Define the 2D mapping of matrix blocks to process grid.
 *
 *   Process grid:
 *     Processes are numbered (0 : P-1).
 *     P = Pr x Pc, where Pr, Pc are the number of process rows and columns.
 *     (pr,pc) is the coordinate of IAM; 0 <= pr < Pr, 0 <= pc < Pc.
 *
 *   Matrix blocks:
 *     Matrix is partitioned according to supernode partitions, both
 *     column and row-wise. 
 *     The k-th block columns (rows) contains columns (rows) (s:t), where
 *             s=xsup[k], t=xsup[k+1]-1.
 *     Block A(I,J) contains
 *             rows from (xsup[I]:xsup[I+1]-1) and
 *             columns from (xsup[J]:xsup[J+1]-1)
 *
 *  Mapping of matrix entry (i,j) to matrix block (I,J):
 *     (I,J) = ( supno[i], supno[j] )
 *
 *  Mapping of matrix block (I,J) to process grid (pr,pc):
 *     (pr,pc) = ( MOD(I,NPROW), MOD(J,NPCOL) )
 *  
 *  (xsup[nsupers],supno[n]) are replicated on all processors.
 *
 */

/*-- Communication subgroup */
typedef struct {
    MPI_Comm comm;        /* MPI communicator */
    int Np;               /* number of processes */
    int Iam;              /* my process number */
} superlu_scope_t;

/*-- Process grid definition */
typedef struct {
    MPI_Comm comm;        /* MPI communicator */
    superlu_scope_t rscp; /* row scope */
    superlu_scope_t cscp; /* column scope */
    int iam;              /* my process number in this scope */
    int_t nprow;          /* number of process rows */
    int_t npcol;          /* number of process columns */
} gridinfo_t;


/*
 *-- The structures are determined by SYMBFACT and used thereafter.
 *
 * (xsup,supno) describes mapping between supernode and column:
 *	xsup[s] is the leading column of the s-th supernode.
 *      supno[i] is the supernode no to which column i belongs;
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 * This is allocated during symbolic factorization SYMBFACT.
 */
typedef struct {
    int_t     *xsup;
    int_t     *supno;
} Glu_persist_t;

/*
 *-- The structures are determined by SYMBFACT and used by DDISTRIBUTE.
 * 
 * (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 * (xusub,usub): xusub[i] points to the starting location of column i
 *      in usub[]. For each U-segment, only the row index of first nonzero
 *      is stored in usub[].
 *
 *      Each U column consists of a number of full segments. Each full segment
 *      starts from a leading nonzero, running up to the supernode (block)
 *      boundary. (Recall that the column-wise supernode partition is also
 *      imposed on the rows.) Because the segment is full, we don't store all
 *      the row indices. Instead, only the leading nonzero index is stored.
 *      The rest can be found together with xsup/supno pair.
 *      For example, 
 *          usub[xsub[j+1]] - usub[xsub[j]] = number of segments in column j.
 *          for any i in usub[], 
 *              supno[i]         = block number in which i belongs to
 *  	        xsup[supno[i]+1] = first row of the next block
 *              The nonzeros of this segment are: 
 *                  i, i+1 ... xsup[supno[i]+1]-1 (only i is stored in usub[])
 *
 */
typedef struct {
    int_t     *lsub;     /* compressed L subscripts */
    int_t     *xlsub;
    int_t     *usub;     /* compressed U subscripts */
    int_t     *xusub;
    int_t     nzlmax;    /* current max size of lsub */
    int_t     nzumax;    /*    "    "    "      usub */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
    int_t     *llvl;     /* keep track of level in L for level-based ILU */
    int_t     *ulvl;     /* keep track of level in U for level-based ILU */
} Glu_freeable_t;


/* 
 *-- The structure used to store matrix A of the linear system and
 *   several vectors describing the transformations done to matrix A.
 *
 * A      (SuperMatrix*)
 *        Matrix A in A*X=B, of dimension (A->nrow, A->ncol).
 *        The number of linear equations is A->nrow. The type of A can be:
 *        Stype = SLU_NC; Dtype = SLU_D; Mtype = SLU_GE.
 *         
 * DiagScale  (DiagScale_t)
 *        Specifies the form of equilibration that was done.
 *        = NOEQUIL: No equilibration.
 *        = ROW:  Row equilibration, i.e., A was premultiplied by diag(R).
 *        = COL:  Column equilibration, i.e., A was postmultiplied by diag(C).
 *        = BOTH: Both row and column equilibration, i.e., A was replaced 
 *                 by diag(R)*A*diag(C).
 *
 * R      double*, dimension (A->nrow)
 *        The row scale factors for A.
 *        If DiagScale = ROW or BOTH, A is multiplied on the left by diag(R).
 *        If DiagScale = NOEQUIL or COL, R is not defined.
 *
 * C      double*, dimension (A->ncol)
 *        The column scale factors for A.
 *        If DiagScale = COL or BOTH, A is multiplied on the right by diag(C).
 *        If DiagScale = NOEQUIL or ROW, C is not defined.
 *         
 * perm_r (int*) dimension (A->nrow)
 *        Row permutation vector which defines the permutation matrix Pr,
 *        perm_r[i] = j means row i of A is in position j in Pr*A.
 *
 * perm_c (int*) dimension (A->ncol)
 *	  Column permutation vector, which defines the 
 *        permutation matrix Pc; perm_c[i] = j means column i of A is 
 *        in position j in A*Pc.
 *
 */
typedef struct {
    DiagScale_t DiagScale;
    double *R;
    double *C; 
    int_t  *perm_r;
    int_t  *perm_c;
} ScalePermstruct_t;

/* 
 *-- This contains the options used to control the solution process.
 *
 * Fact   (fact_t)
 *        Specifies whether or not the factored form of the matrix
 *        A is supplied on entry, and if not, how the matrix A should
 *        be factorizaed.
 *        = DOFACT: The matrix A will be factorized from scratch, and the
 *             factors will be stored in L and U.
 *        = SamePattern: The matrix A will be factorized assuming
 *             that a factorization of a matrix with the same sparsity
 *             pattern was performed prior to this one. Therefore, this
 *             factorization will reuse column permutation vector 
 *             ScalePermstruct->perm_c and the column elimination tree
 *             LUstruct->etree.
 *        = SamePattern_SameRowPerm: The matrix A will be factorized
 *             assuming that a factorization of a matrix with the same
 *             sparsity	pattern and similar numerical values was performed
 *             prior to this one. Therefore, this factorization will reuse
 *             both row and column scaling factors R and C, both row and
 *             column permutation vectors perm_r and perm_c, and the
 *             data structure set up from the previous symbolic factorization.
 *        = FACTORED: On entry, L, U, perm_r and perm_c contain the 
 *              factored form of A. If DiagScale is not NOEQUIL, the matrix
 *              A has been equilibrated with scaling factors R and C.
 *
 * Equil  (yes_no_t)
 *        Specifies whether to equilibrate the system (scale A's row and
 *        columns to have unit norm).
 *
 * ColPerm (colperm_t)
 *        Specifies what type of column permutation to use to reduce fill.
 *        = NATURAL: use the natural ordering 
 *        = MMD_ATA: use minimum degree ordering on structure of A'*A
 *        = MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
 *        = COLAMD: use approximate minimum degree column ordering
 *        = MY_PERMC: use the ordering specified by the user
 *         
 * Trans  (trans_t)
 *        Specifies the form of the system of equations:
 *        = NOTRANS: A * X = B        (No transpose)
 *        = TRANS:   A**T * X = B     (Transpose)
 *        = CONJ:    A**H * X = B     (Transpose)
 *
 * IterRefine (IterRefine_t)
 *        Specifies whether to perform iterative refinement.
 *        = NO: no iterative refinement
 *        = SINGLE: perform iterative refinement in single precision
 *        = DOUBLE: perform iterative refinement in double precision
 *        = EXTRA: perform iterative refinement in extra precision
 *
 * DiagPivotThresh (double, in [0.0, 1.0]) (only for sequential SuperLU)
 *        Specifies the threshold used for a diagonal entry to be an
 *        acceptable pivot.
 *
 * SymmetricMode (yest_no_t)
 *        Specifies whether to use symmetric mode. Symmetric mode gives 
 *        preference to diagonal pivots, and uses an (A'+A)-based column
 *        permutation algorithm.
 *
 * PivotGrowth (yes_no_t)
 *        Specifies whether to compute the reciprocal pivot growth.
 *
 * ConditionNumber (ues_no_t)
 *        Specifies whether to compute the reciprocal condition number.
 *
 * RowPerm (rowperm_t) (only for SuperLU_DIST or ILU)
 *        Specifies whether to permute rows of the original matrix.
 *        = NO: not to permute the rows
 *        = LargeDiag: make the diagonal large relative to the off-diagonal
 *        = MY_PERMR: use the permutation given by the user
 *
 * ILU_DropRule (int)
 *        Specifies the dropping rule:
 *	  = DROP_BASIC:   Basic dropping rule, supernodal based ILUTP(tau).
 *	  = DROP_PROWS:   Supernodal based ILUTP(p,tau), p = gamma * nnz(A)/n.
 *	  = DROP_COLUMN:  Variant of ILUTP(p,tau), for j-th column,
 *			      p = gamma * nnz(A(:,j)).
 *	  = DROP_AREA:    Variation of ILUTP, for j-th column, use
 *			      nnz(F(:,1:j)) / nnz(A(:,1:j)) to control memory.
 *	  = DROP_DYNAMIC: Modify the threshold tau during factorizaion:
 *			  If nnz(L(:,1:j)) / nnz(A(:,1:j)) > gamma
 *				  tau_L(j) := MIN(tau_0, tau_L(j-1) * 2);
 *			  Otherwise
 *				  tau_L(j) := MAX(tau_0, tau_L(j-1) / 2);
 *			  tau_U(j) uses the similar rule.
 *			  NOTE: the thresholds used by L and U are separate.
 *	  = DROP_INTERP:  Compute the second dropping threshold by
 *	                  interpolation instead of sorting (default).
 *  		          In this case, the actual fill ratio is not
 *			  guaranteed to be smaller than gamma.
 *   	  Note: DROP_PROWS, DROP_COLUMN and DROP_AREA are mutually exclusive.
 *	  ( Default: DROP_BASIC | DROP_AREA )
 *
 * ILU_DropTol (double)
 *        numerical threshold for dropping.
 *
 * ILU_FillFactor (double) 
 *        Gamma in the secondary dropping.
 *
 * ILU_Norm (norm_t)
 *        Specify which norm to use to measure the row size in a
 *        supernode: infinity-norm, 1-norm, or 2-norm.
 *
 * ILU_FillTol (double)
 *        numerical threshold for zero pivot perturbation.
 *
 * ILU_MILU (milu_t)
 *        Specifies which version of MILU to use.
 *
 * ILU_MILU_Dim (double) 
 *        Dimension of the PDE if available.
 *
 * ReplaceTinyPivot (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether to replace the tiny diagonals by
 *        sqrt(epsilon)*||A|| during LU factorization.
 *
 * SolveInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        triangular solve.
 *
 * RefineInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        sparse matrix-vector multiplication routine needed in iterative
 *        refinement.
 *
 * PrintStat (yes_no_t)
 *        Specifies whether to print the solver's statistics.
 */
typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    double        DiagPivotThresh;
    yes_no_t      SymmetricMode;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    int 	  ILU_DropRule;
    double	  ILU_DropTol;    /* threshold for dropping */
    double	  ILU_FillFactor; /* gamma in the secondary dropping */
    norm_t	  ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
    double	  ILU_FillTol;    /* threshold for zero pivot perturbation */
    milu_t	  ILU_MILU;
    double	  ILU_MILU_Dim;   /* Dimension of PDE (if available) */
    yes_no_t      ParSymbFact;
    yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
    yes_no_t      PrintStat;
} superlu_options_t;

typedef struct {
    float for_lu;
    float total;
    int_t expansions;
    int_t nnzL, nnzU;
} mem_usage_t;


/***********************************************************************
 * Function prototypes
 ***********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

extern void    set_default_options_dist(superlu_options_t *);
extern void    print_options_dist(superlu_options_t *);
extern void    Destroy_CompCol_Matrix_dist(SuperMatrix *);
extern void    Destroy_SuperNode_Matrix_dist(SuperMatrix *);
extern void    Destroy_SuperMatrix_Store_dist(SuperMatrix *);
extern void    Destroy_CompCol_Permuted_dist(SuperMatrix *);
extern void    Destroy_CompRowLoc_Matrix_dist(SuperMatrix *);
extern void    Destroy_CompRow_Matrix_dist(SuperMatrix *);
extern void    sp_colorder (superlu_options_t*, SuperMatrix*, int_t*, int_t*,
			    SuperMatrix*);
extern int_t   sp_coletree_dist (int_t *, int_t *, int_t *, int_t, int_t,
				 int_t *);
extern void    countnz_dist (const int_t, int_t *, int_t *, int_t *,
			     Glu_persist_t *, Glu_freeable_t *);
extern int_t   fixupL_dist (const int_t, const int_t *, Glu_persist_t *,
			    Glu_freeable_t *);
extern int_t   *TreePostorder_dist (int_t, int_t *);
extern float   slamch_(char *);
extern double  dlamch_(char *);
extern void    *superlu_malloc_dist (size_t);
extern void    superlu_free_dist (void*);
extern int_t   *intMalloc_dist (int_t);
extern int_t   *intCalloc_dist (int_t);

/* Auxiliary routines */
extern double  SuperLU_timer_ ();
extern void    superlu_abort_and_exit_dist(char *);
extern int_t   sp_ienv_dist (int_t);
extern int     lsame_ (char *, char *);
extern int     xerbla_ (char *, int *);
extern void    ifill_dist (int_t *, int_t, int_t);
extern void    super_stats_dist (int_t, int_t *);
extern void    ScalePermstructInit(const int_t, const int_t, 
				   ScalePermstruct_t *);
extern void    ScalePermstructFree(ScalePermstruct_t *);
extern void  superlu_gridinit(MPI_Comm, int_t, int_t, gridinfo_t *);
extern void  superlu_gridmap(MPI_Comm, int_t, int_t, int_t [], int_t,
			     gridinfo_t *);
extern void  superlu_gridexit(gridinfo_t *);
extern void  get_perm_c_dist(int_t, int_t, SuperMatrix *, int_t *);
extern void  a_plus_at_dist(const int_t, const int_t, int_t *, int_t *,
			    int_t *, int_t **, int_t **);
extern void  bcast_tree(void *, int, MPI_Datatype, int, int,
			gridinfo_t *, int, int *);
extern int_t symbfact(superlu_options_t *, int, SuperMatrix *, int_t *,
                      int_t *, Glu_persist_t *, Glu_freeable_t *);
extern int_t symbfact_SubInit(fact_t, void *, int_t, int_t, int_t, int_t,
			      Glu_persist_t *, Glu_freeable_t *);
extern int_t symbfact_SubXpand(int_t, int_t, int_t, MemType, int_t *,
			       Glu_freeable_t *);
extern int_t symbfact_SubFree(Glu_freeable_t *);
extern void  get_diag_procs(int_t, Glu_persist_t *, gridinfo_t *, int_t *,
			    int_t **, int_t **);
extern int_t QuerySpace_dist(int_t, int_t, Glu_freeable_t *, mem_usage_t *);
extern int   xerbla_ (char *, int *);
extern void  pxerbla (char *, gridinfo_t *, int_t);
extern void  PStatInit(SuperLUStat_t *);
extern void  PStatFree(SuperLUStat_t *);
extern void  PStatPrint(superlu_options_t *, SuperLUStat_t *, gridinfo_t *);


/* Prototypes for parallel symbolic factorization */
extern float symbfact_dist
(int,  int, SuperMatrix *, int_t *, int_t *,  int_t *, int_t *,
 Pslu_freeable_t *, MPI_Comm *, MPI_Comm *,  mem_usage_t *);

/* Get the column permutation using parmetis */
extern float get_perm_c_parmetis 
(SuperMatrix *, int_t *, int_t *, int, int, 
 int_t **, int_t **, gridinfo_t *, MPI_Comm *);

/* Auxiliary routines for memory expansions used during
   the parallel symbolic factorization routine */

extern int_t psymbfact_LUXpandMem
(int_t, int_t, int_t, int_t, int_t, int_t, int_t, int_t, 
 Pslu_freeable_t *, Llu_symbfact_t *,  vtcsInfo_symbfact_t *, psymbfact_stat_t *);

extern int_t psymbfact_LUXpand
(int_t, int_t, int_t, int_t, int_t *, int_t, int_t, int_t, int_t, 
 Pslu_freeable_t *, Llu_symbfact_t *,  vtcsInfo_symbfact_t *, psymbfact_stat_t *);

extern int_t psymbfact_LUXpand_RL
(int_t, int_t, int_t, int_t, int_t, int_t,
 Pslu_freeable_t *, Llu_symbfact_t *, vtcsInfo_symbfact_t *, psymbfact_stat_t *);

extern int_t psymbfact_prLUXpand
(int_t,  int_t,
 MemType, Llu_symbfact_t *, psymbfact_stat_t *);

/* Routines for debugging */
extern void  print_panel_seg_dist(int_t, int_t, int_t, int_t, int_t *, int_t *);
extern void  check_repfnz_dist(int_t, int_t, int_t, int_t *);
extern int_t CheckZeroDiagonal(int_t, int_t *, int_t *, int_t *);
extern void  PrintDouble5(char *, int_t, double *);
extern void  PrintInt10(char *, int_t, int_t *);
extern int   file_PrintInt10(FILE *, char *, int_t, int_t *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_DEFS */
