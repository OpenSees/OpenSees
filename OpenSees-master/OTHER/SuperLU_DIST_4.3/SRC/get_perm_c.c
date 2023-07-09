/*! @file
 * \brief Gets matrix permutation
 *
 * <pre>
 * -- Distributed SuperLU routine (version 2.1) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * November 1, 2007
 * Feburary 20, 2008
 * </pre>
 *
 * Last update: 7/27/2011  fix a bug with metis ordering on empty graph.
 *
 */

#include "superlu_ddefs.h"


void
get_metis(
	  int_t n,         /* dimension of matrix B */
	  int_t bnz,       /* number of nonzeros in matrix A. */
	  int_t *b_colptr, /* column pointer of size n+1 for matrix B. */
	  int_t *b_rowind, /* row indices of size bnz for matrix B. */
	  int_t *perm_c    /* out - the column permutation vector. */
	  )
{
    /*#define METISOPTIONS 8*/
#define METISOPTIONS 40
    int_t metis_options[METISOPTIONS];
    int_t i, nm, numflag = 0; /* C-Style ordering */
    int_t *perm, *iperm;
    int_t *b_colptr_int, *b_rowind_int;
    extern int check_perm_dist(char *what, int_t n, int_t *perm);

    extern int METIS_NodeND(int_t*, int_t*, int_t*, int_t*, int_t*,
			    int_t*, int_t*);

    metis_options[0] = 0; /* Use Defaults for now */

    perm = (int_t*) SUPERLU_MALLOC(2*n * sizeof(int_t));
    if (!perm) ABORT("SUPERLU_MALLOC fails for perm.");
    iperm = perm + n;
    nm = n;

#if 0
#if defined(_LONGINT)
    /* Metis can only take 32-bit integers */

    if ( !(b_colptr_int = (int*) SUPERLU_MALLOC((n+1) * sizeof(int))) )
	 ABORT("SUPERLU_MALLOC fails for b_colptr_int.");
    for (i = 0; i < n+1; ++i) b_colptr_int[i] = b_colptr[i];
    SUPERLU_FREE(b_colptr);
    
    if ( !(b_rowind_int = (int*) SUPERLU_MALLOC(bnz * sizeof(int))) )
	ABORT("SUPERLU_MALLOC fails for b_rowind_int.");

    for (i = 0; i < bnz; ++i) b_rowind_int[i] = b_rowind[i];
    SUPERLU_FREE(b_rowind);
#else
    b_colptr_int = b_colptr;
    b_rowind_int = b_rowind;
#endif
#endif

    /* Call metis */
#undef USEEND
#ifdef USEEND
    METIS_EdgeND(&nm, b_colptr_int, b_rowind_int, &numflag, metis_options,
		 perm, iperm);
#else

    /* Earlier version 3.x.x */
    /* METIS_NodeND(&nm, b_colptr, b_rowind, &numflag, metis_options,
       perm, iperm);*/

    /* Latest version 4.x.x */
    METIS_NodeND(&nm, b_colptr, b_rowind, NULL, NULL, perm, iperm);

    /*check_perm_dist("metis perm",  n, perm);*/
#endif

    /* Copy the permutation vector into SuperLU data structure. */
    for (i = 0; i < n; ++i) perm_c[i] = iperm[i];

#if 0
    SUPERLU_FREE(b_colptr_int);
    SUPERLU_FREE(b_rowind_int);
#else
    SUPERLU_FREE(b_colptr);
    SUPERLU_FREE(b_rowind);
#endif
    SUPERLU_FREE(perm);
}

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * Form the structure of A'*A. A is an m-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'*A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (ata_colptr, ata_rowind).
 *
 * This routine is modified from GETATA routine by Tim Davis.
 * The complexity of this algorithm is: SUM_{i=1,m} r(i)^2,
 * i.e., the sum of the square of the row counts.
 *
 * Questions
 * =========
 *     o  Do I need to withhold the *dense* rows?
 *     o  How do I know the number of nonzeros in A'*A?
 * </pre>
 */
void
getata_dist(
	    const int_t m,    /* number of rows in matrix A. */
	    const int_t n,    /* number of columns in matrix A. */
	    const int_t nz,   /* number of nonzeros in matrix A */
	    int_t *colptr,    /* column pointer of size n+1 for matrix A. */
	    int_t *rowind,    /* row indices of size nz for matrix A. */
	    int_t *atanz,     /* out - on exit, returns the actual number of
				 nonzeros in matrix A'*A. */
	    int_t **ata_colptr, /* out - size n+1 */
	    int_t **ata_rowind  /* out - size *atanz */
	    )
{

    register int_t i, j, k, col, num_nz, ti, trow;
    int_t *marker, *b_colptr, *b_rowind;
    int_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */

    if ( !(marker = (int_t*) SUPERLU_MALLOC( (SUPERLU_MAX(m,n)+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int_t*) SUPERLU_MALLOC( (m+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC t_colptr[]");
    if ( !(t_rowind = (int_t*) SUPERLU_MALLOC( nz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for t_rowind[]");

    
    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < m; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < m; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j)
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}

    
    /* ----------------------------------------------------------------
       compute B = T * A, where column j of B is:

       Struct (B_*j) =    UNION   ( Struct (T_*k) )
                        A_kj != 0

       do not include the diagonal entry
   
       ( Partition A as: A = (A_*1, ..., A_*n)
         Then B = T * A = (T * A_*1, ..., T * A_*n), where
         T * A_*j = (T_*1, ..., T_*m) * A_*j.  )
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    num_nz++;
		}
	    }
	}
    }
    *atanz = num_nz;
    
    /* Allocate storage for A'*A */
    if ( !(*ata_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for ata_colptr[]");
    if ( *atanz ) {
	if ( !(*ata_rowind = (int_t*)SUPERLU_MALLOC(*atanz*sizeof(int_t)) ) ) {
	    fprintf(stderr, ".. atanz = %ld\n", (long long) *atanz);
	    ABORT("SUPERLU_MALLOC fails for ata_rowind[]");
	}
    }
    b_colptr = *ata_colptr; /* aliasing */
    b_rowind = *ata_rowind;
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	b_colptr[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    b_rowind[num_nz++] = trow;
		}
	    }
	}
    }
    b_colptr[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
}

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * Form the structure of A'+A. A is an n-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'+A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (b_colptr, b_rowind).
 * </pre>
 */
void
at_plus_a_dist(
	       const int_t n,    /* number of columns in matrix A. */
	       const int_t nz,   /* number of nonzeros in matrix A */
	       int_t *colptr,    /* column pointer of size n+1 for matrix A. */
	       int_t *rowind,    /* row indices of size nz for matrix A. */
	       int_t *bnz,       /* out - on exit, returns the actual number of
				    nonzeros in matrix A'+A. */
	       int_t **b_colptr, /* out - size n+1 */
	       int_t **b_rowind  /* out - size *bnz */
	       )
{

    register int_t i, j, k, col, num_nz;
    int_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */
    int_t *marker;

    if ( !(marker = (int_t*) SUPERLU_MALLOC( n * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for t_colptr[]");
    if ( !(t_rowind = (int_t*) SUPERLU_MALLOC( nz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails t_rowind[]");

    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }

    t_colptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}
    }


    /* ----------------------------------------------------------------
       compute B = A + T, where column j of B is:

       Struct (B_*j) = Struct (A_*k) UNION Struct (T_*k)

       do not include the diagonal entry
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}
    }
    *bnz = num_nz;


    /* Allocate storage for A+A' */
    if ( !(*b_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for b_colptr[]");
    if ( *bnz ) {
	if ( !(*b_rowind = (int_t*) SUPERLU_MALLOC( *bnz * sizeof(int_t)) ) )
	    ABORT("SUPERLU_MALLOC fails for b_rowind[]");
    }
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	(*b_colptr)[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}
    }
    (*b_colptr)[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
} /* at_plus_a_dist */

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * GET_PERM_C_DIST obtains a permutation matrix Pc, by applying the multiple
 * minimum degree ordering code by Joseph Liu to matrix A'*A or A+A',
 * or using approximate minimum degree column ordering by Davis et. al.
 * The LU factorization of A*Pc tends to have less fill than the LU 
 * factorization of A.
 *
 * Arguments
 * =========
 *
 * ispec   (input) colperm_t
 *         Specifies what type of column permutation to use to reduce fill.
 *         = NATURAL: natural ordering (i.e., Pc = I)
 *         = MMD_AT_PLUS_A: minimum degree ordering on structure of A'+A
 *         = MMD_ATA: minimum degree ordering on structure of A'*A
 *         = METIS_AT_PLUS_A: MeTis on A'+A
 * 
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of the linear equations is A->nrow. Currently, the type of A 
 *         can be: Stype = SLU_NC; Dtype = SLU_D; Mtype = SLU_GE.
 *         In the future, more general A can be handled.
 *
 * perm_c  (output) int*
 *	   Column permutation vector of size A->ncol, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 * </pre>
 */
void
get_perm_c_dist(int_t pnum, int_t ispec, SuperMatrix *A, int_t *perm_c)

{
    NCformat *Astore = A->Store;
    int_t m, n, bnz = 0, *b_colptr, *b_rowind, i;
    int_t delta, maxint, nofsub, *invp;
    int_t *dhead, *qsize, *llist, *marker;
    double t, SuperLU_timer_();

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC((int)pnum, "Enter get_perm_c_dist()");
#endif

    m = A->nrow;
    n = A->ncol;

    t = SuperLU_timer_();

    switch ( ispec ) {

        case NATURAL: /* Natural ordering */
	      for (i = 0; i < n; ++i) perm_c[i] = i;
#if ( PRNTlevel>=1 )
	      if ( !pnum ) printf(".. Use natural column ordering\n");
#endif
	      return;

        case MMD_AT_PLUS_A: /* Minimum degree ordering on A'+A */
	      if ( m != n ) ABORT("Matrix is not square");
	      at_plus_a_dist(n, Astore->nnz, Astore->colptr, Astore->rowind,
			     &bnz, &b_colptr, &b_rowind);
	      t = SuperLU_timer_() - t;
	      /*printf("Form A'+A time = %8.3f\n", t);*/
#if ( PRNTlevel>=1 )
	      if ( !pnum ) printf(".. Use minimum degree ordering on A'+A.\n");
#endif
	      break;

        case MMD_ATA: /* Minimum degree ordering on A'*A */
	      getata_dist(m, n, Astore->nnz, Astore->colptr, Astore->rowind,
			  &bnz, &b_colptr, &b_rowind);
	      t = SuperLU_timer_() - t;
	      /*printf("Form A'*A time = %8.3f\n", t);*/
#if ( PRNTlevel>=1 )
	      if ( !pnum ) printf(".. Use minimum degree ordering on A'*A\n");
#endif
	      break;

        case METIS_AT_PLUS_A: /* METIS ordering on A'+A */
	      if ( m != n ) ABORT("Matrix is not square");
	      at_plus_a_dist(n, Astore->nnz, Astore->colptr, Astore->rowind,
			     &bnz, &b_colptr, &b_rowind);

	      if ( bnz ) { /* non-empty adjacency structure */
		  get_metis(n, bnz, b_colptr, b_rowind, perm_c);
	      } else { /* e.g., diagonal matrix */
		  for (i = 0; i < n; ++i) perm_c[i] = i;
		  SUPERLU_FREE(b_colptr);
		  /* b_rowind is not allocated in this case */
	      }

#if ( PRNTlevel>=1 )
	      if ( !pnum ) printf(".. Use METIS ordering on A'+A\n");
#endif
	      return;

        default:
	      ABORT("Invalid ISPEC");
    }

    if ( bnz ) {
	t = SuperLU_timer_();

	/* Initialize and allocate storage for GENMMD. */
	delta = 0; /* DELTA is a parameter to allow the choice of nodes
		      whose degree <= min-degree + DELTA. */
	maxint = 2147483647; /* 2**31 - 1 */
	invp = (int_t *) SUPERLU_MALLOC((n+delta)*sizeof(int_t));
	if ( !invp ) ABORT("SUPERLU_MALLOC fails for invp.");
	dhead = (int_t *) SUPERLU_MALLOC((n+delta)*sizeof(int_t));
	if ( !dhead ) ABORT("SUPERLU_MALLOC fails for dhead.");
	qsize = (int_t *) SUPERLU_MALLOC((n+delta)*sizeof(int_t));
	if ( !qsize ) ABORT("SUPERLU_MALLOC fails for qsize.");
	llist = (int_t *) SUPERLU_MALLOC(n*sizeof(int_t));
	if ( !llist ) ABORT("SUPERLU_MALLOC fails for llist.");
	marker = (int_t *) SUPERLU_MALLOC(n*sizeof(int_t));
	if ( !marker ) ABORT("SUPERLU_MALLOC fails for marker.");
	
	/* Transform adjacency list into 1-based indexing required by GENMMD.*/
	for (i = 0; i <= n; ++i) ++b_colptr[i];
	for (i = 0; i < bnz; ++i) ++b_rowind[i];
	
	genmmd_dist_(&n, b_colptr, b_rowind, perm_c, invp, &delta, dhead, 
		     qsize, llist, marker, &maxint, &nofsub);

	/* Transform perm_c into 0-based indexing. */
	for (i = 0; i < n; ++i) --perm_c[i];

	SUPERLU_FREE(invp);
	SUPERLU_FREE(dhead);
	SUPERLU_FREE(qsize);
	SUPERLU_FREE(llist);
	SUPERLU_FREE(marker);
	SUPERLU_FREE(b_rowind);

	t = SuperLU_timer_() - t;
	/*    printf("call GENMMD time = %8.3f\n", t);*/

    } else { /* Empty adjacency structure */
	for (i = 0; i < n; ++i) perm_c[i] = i;
    }

    SUPERLU_FREE(b_colptr);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC((int) pnum, "Exit get_perm_c_dist()");
#endif
} /* get_perm_c_dist */
