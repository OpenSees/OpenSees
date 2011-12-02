/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */

#ifndef __SUPERLU_SUPERMATRIX /* allow multiple inclusions */
#define __SUPERLU_SUPERMATRIX

/********************************************
 * The matrix types are defined as follows. *
 ********************************************/
typedef enum {
    NC,        /* column-wise, no supernode */
    NCP,       /* column-wise, column-permuted, no supernode 
                  (The consecutive columns of nonzeros, after permutation,
		   may not be stored  contiguously.) */
    NR,        /* row-wize, no supernode */
    SC,        /* column-wise, supernode */
    SCP,       /* supernode, column-wise, permuted */    
    SR,        /* row-wise, supernode */
    DN         /* Fortran style column-wise storage for dense matrix */
} Stype_t;

typedef enum {
    _S,         /* single */
    _D,         /* double */
    _C,         /* single complex */
    _Z          /* double complex */
} Dtype_t;

typedef enum {
    GE,        /* general */
    TRLU,      /* lower triangular, unit diagonal */
    TRUU,      /* upper triangular, unit diagonal */
    TRL,       /* lower triangular */
    TRU,       /* upper triangular */
    SYL,       /* symmetric, store lower half */
    SYU,       /* symmetric, store upper half */
    HEL,       /* Hermitian, store lower half */
    HEU        /* Hermitian, store upper half */
} Mtype_t;

typedef struct {
	Stype_t Stype; /* Storage type: interprets the storage structure 
		   	  pointed to by *Store. */
	Dtype_t Dtype; /* Data type. */
	Mtype_t Mtype; /* Matrix type: describes the mathematical property of 
			  the matrix. */
	int  nrow;     /* number of rows */
	int  ncol;     /* number of columns */
	void *Store;   /* pointer to the actual storage of the matrix */
} SuperMatrix;

/***********************************************
 * The storage schemes are defined as follows. *
 ***********************************************/

/* Stype == NC (Also known as Harwell-Boeing sparse matrix format) */
typedef struct {
    int  nnz;	  /* number of nonzeros in the matrix */
    void *nzval;  /* pointer to array of nonzero values, packed by column */
    int  *rowind; /* pointer to array of row indices of the nonzeros */
    int  *colptr; /* pointer to array of beginning of columns in nzval[] 
                     and rowind[]  */
                  /* Note:
		     Zero-based indexing is used;
		     colptr[] has ncol+1 entries, the last one pointing
		         beyond the last column, so that colptr[ncol] = nnz. */
} NCformat;

/* Stype == NR (Also known as row compressed storage (RCS). */
typedef struct {
    int  nnz;	  /* number of nonzeros in the matrix */
    void *nzval;  /* pointer to array of nonzero values, packed by column */
    int  *colind; /* pointer to array of column indices of the nonzeros */
    int  *rowptr; /* pointer to array of beginning of rows in nzval[] 
                     and colind[]  */
                  /* Note:
		     Zero-based indexing is used;
		     rowptr[] has nrow+1 entries, the last one pointing
		         beyond the last column, so that rowptr[nrow] = nnz. */
} NRformat;

/* Stype == SC */
typedef struct {
  int  nnz;	     /* number of nonzeros in the matrix */
  int  nsuper;       /* number of supernodes, minus 1 */
  void *nzval;       /* pointer to array of nonzero values, packed by column */
  int  *nzval_colptr;/* pointer to array of beginning of columns in nzval[] */
  int  *rowind;      /* pointer to array of compressed row indices of 
			rectangular supernodes */
  int *rowind_colptr;/* pointer to array of beginning of columns in rowind[] */
  int *col_to_sup;   /* col_to_sup[j] is the supernode number to which column 
			j belongs; mapping from column to supernode number. */
  int *sup_to_col;   /* sup_to_col[s] points to the start of the s-th 
			supernode; mapping from supernode number to column.
		        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
		              sup_to_col: 0 1 2 4 7 12           (nsuper=4) */
                     /* Note:
		        Zero-based indexing is used;
		        nzval_colptr[], rowind_colptr[], col_to_sup and
		        sup_to_col[] have ncol+1 entries, the last one
		        pointing beyond the last column.         */
} SCformat;

/* Stype == SCP */
typedef struct {
  int  nnz;	     /* number of nonzeros in the matrix */
  int  nsuper;       /* number of supernodes */
  void *nzval;       /* pointer to array of nonzero values, packed by column */
  int  *nzval_colbeg;/* nzval_colbeg[j] points to beginning of column j
			in nzval[] */
  int  *nzval_colend;/* nzval_colend[j] points to one past the last element
			of column j in nzval[] */
  int  *rowind;      /* pointer to array of compressed row indices of 
			rectangular supernodes */
  int *rowind_colbeg;/* rowind_colbeg[j] points to beginning of column j
			in rowind[] */
  int *rowind_colend;/* rowind_colend[j] points to one past the last element
			of column j in rowind[] */
  int *col_to_sup;   /* col_to_sup[j] is the supernode number to which column 
			j belongs; mapping from column to supernode number. */
  int *sup_to_colbeg; /* sup_to_colbeg[s] points to the start of the s-th 
			supernode; mapping from supernode number to column.*/
  int *sup_to_colend; /* sup_to_colend[s] points to one past the end of the
			  s-th supernode; mapping from supernode number to
			  column.
		        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
		              sup_to_colbeg: 0 1 2 4 7              (nsuper=4)
			      sup_to_colend: 1 2 4 7 12                    */
                     /* Note:
		        Zero-based indexing is used;
		        nzval_colptr[], rowind_colptr[], col_to_sup and
		        sup_to_col[] have ncol+1 entries, the last one
		        pointing beyond the last column.         */
} SCPformat;

/* Stype == NCP */
typedef struct {
    int nnz;	  /* number of nonzeros in the matrix */
    void *nzval;  /* pointer to array of nonzero values, packed by column */
    int *rowind;  /* pointer to array of row indices of the nonzeros */
		  /* Note: nzval[]/rowind[] always have the same length */
    int *colbeg;  /* colbeg[j] points to the beginning of column j in nzval[] 
                     and rowind[]  */
    int *colend;  /* colend[j] points to one past the last element of column
		     j in nzval[] and rowind[]  */
		  /* Note:
		     Zero-based indexing is used;
		     The consecutive columns of the nonzeros may not be 
		     contiguous in storage, because the matrix has been 
		     postmultiplied by a column permutation matrix. */
} NCPformat;

/* Stype == DN */
typedef struct {
    int lda;      /* leading dimension */
    void *nzval;  /* array of size lda*ncol to represent a dense matrix */
} DNformat;


/*********************************************************
 * Macros used for easy access of sparse matrix entries. *
 *********************************************************/
#define L_SUB_START(col)     ( Lstore->rowind_colbeg[col] )
#define L_SUB_END(col)       ( Lstore->rowind_colend[col] )
#define L_SUB(ptr)           ( Lstore->rowind[ptr] )
#define L_NZ_START(col)      ( Lstore->nzval_colbeg[col] )
#define L_NZ_END(col)        ( Lstore->nzval_colend[col] )
#define L_FST_SUPC(superno)  ( Lstore->sup_to_colbeg[superno] )
#define L_LAST_SUPC(superno) ( Lstore->sup_to_colend[superno] )
#define U_NZ_START(col)      ( Ustore->colbeg[col] )
#define U_NZ_END(col)        ( Ustore->colend[col] )
#define U_SUB(ptr)           ( Ustore->rowind[ptr] )

#endif  /* __SUPERLU_SUPERMATRIX */
