/*! @file
 * \brief Header for utilities
 */

#ifndef __SUPERLU_UTIL /* allow multiple inclusions */
#define __SUPERLU_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "superlu_enum_consts.h"

/*
 * Macros
 */
#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit_dist(msg)
#endif

#define ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#ifndef USER_MALLOC
#define USER_MALLOC(size) superlu_malloc_dist(size)
#endif

#define SUPERLU_MALLOC(size) USER_MALLOC(size)

#ifndef USER_FREE
#define USER_FREE(addr) superlu_free_dist(addr)
#endif

#define SUPERLU_FREE(addr) USER_FREE(addr)

#define CHECK_MALLOC(pnum, where) {                 \
    extern long int superlu_malloc_total;        \
    printf("(%d) %s: superlu_malloc_total (MB) %.6f\n", \
	   pnum, where, superlu_malloc_total*1e-6); \
}

#define SUPERLU_MAX(x, y) 	( (x) > (y) ? (x) : (y) )
#define SUPERLU_MIN(x, y) 	( (x) < (y) ? (x) : (y) )

    
/* 
 * Constants 
 */
#define EMPTY	(-1)
#ifndef FALSE
#define FALSE	(0)
#endif
#ifndef TRUE
#define TRUE	(1)
#endif

/*
 * Type definitions
 */
typedef float    flops_t;
typedef unsigned char Logical;

/*
#ifdef _CRAY
#define int short
#endif
*/

typedef struct {
    int     *panel_histo; /* histogram of panel size distribution */
    double  *utime;       /* running time at various phases */
    flops_t *ops;         /* operation count at various phases */
    int     TinyPivots;   /* number of tiny pivots */
    int     RefineSteps;  /* number of iterative refinement steps */
    int     num_look_aheads; /* number of look ahead */
    /*-- new --*/
    float   current_buffer; /* bytes allocated for buffer in numerical factorization */
    float   peak_buffer;    /* monitor the peak buffer size (bytes) */
    float   gpu_buffer;     /* monitor the buffer allocated on GPU (bytes) */
} SuperLUStat_t;

/* Headers for 2 types of dynamatically managed memory */
typedef struct e_node {
    int size;      /* length of the memory that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

typedef struct {
    int  size;
    int  used;
    int  top1;  /* grow upward, relative to &array[0] */
    int  top2;  /* grow downward */
    void *array;
} LU_stack_t;

/* Constants */
#define GluIntArray(n)   (5 * (n) + 5)
#define NO_MEMTYPE  6      /* 0: lusup;
			      1: ucol;
			      2: lsub;
			      3: usub
			      4: llvl; level number in L for ILU(k)
			      5: ulvl; level number in U for ILU(k)
                           */

/* Macros to manipulate stack */
#define StackFull(x)         ( x + stack.used >= stack.size )
#define NotDoubleAlign(addr) ( (long)addr & 7 )
#define DoubleAlign(addr)    ( ((long)addr + 7) & ~7L )
#define TempSpace(n, w)      ( (2*w + 4 + NO_MARKER)*m*sizeof(int) + \
			      (w + 1)*n*sizeof(double) )
#define Reduce(alpha)        ((alpha + 1) / 2)  /* i.e. (alpha-1)/2 + 1 */

#define FIRSTCOL_OF_SNODE(i)	(xsup[i])

#if ( PROFlevel>=1 )
#define TIC(t)          t = SuperLU_timer_()
#define TOC(t2, t1)     t2 = SuperLU_timer_() - t1
#else
#define TIC(t)
#define TOC(t2, t1)
#endif

/*********************************************************
 * Macros used for easy access of sparse matrix entries. *
 *********************************************************/
#define L_SUB_START(col)     ( Lstore->rowind_colptr[col] )
#define L_SUB(ptr)           ( Lstore->rowind[ptr] )
#define L_NZ_START(col)      ( Lstore->nzval_colptr[col] )
#define L_FST_SUPC(superno)  ( Lstore->sup_to_col[superno] )
#define U_NZ_START(col)      ( Ustore->colptr[col] )
#define U_SUB(ptr)           ( Ustore->rowind[ptr] )

#endif /* __SUPERLU_UTIL */
