/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SuperLU_MT_util.h,v $
                                                                        
                                                                        
/*
 * -- SuperLU MT routine (alpha version) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */

#ifndef __SUPERLU_UTIL /* allow multiple inclusions */
#define __SUPERLU_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

/* Macros */
#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit(msg)
#endif

#define ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#ifndef USER_MALLOC
#define USER_MALLOC(size) superlu_malloc(size)
#endif

#define SUPERLU_MALLOC(size) USER_MALLOC(size)

#ifndef USER_FREE
#define USER_FREE(addr) superlu_free(addr)
#endif

#define SUPERLU_FREE(addr) USER_FREE(addr)


#define MAX(x, y) 	( (x) > (y) ? (x) : (y) )
#define MIN(x, y) 	( (x) < (y) ? (x) : (y) )

/* 
 * Constants 
 */
#define EMPTY	(-1)
#define FALSE	0
#define TRUE	1

/*
 * Type definitions
 */
typedef float    flops_t;
typedef unsigned char Logical;
typedef enum {
    RELAX,
    ETREE,
    EQUIL,
    FINDDOMAIN,
    FACT,
    DFS,
    FLOAT,
    TRSV,
    GEMV,
    RCOND,
    TRISOLVE,
    SOLVE,
    REFINE,
    FERR,
    NPHASES
} PhaseType;

/* ----------------------------------------------
    The definitions below are used for profiling.
   ---------------------------------------------- */

/* The statistics to be kept by each processor. */
typedef struct {
    int	    panels;    /* number of panels taken */
    float   fcops;     /* factor floating-point operations */
    double  fctime;    /* factor time */
    int     skedwaits; /* how many times the processor fails to get a task */
    double  skedtime;  /* time spent in the scheduler */
    double  cs_time;   /* time spent in the critical sections */
    double  spintime;  /* spin-wait time */
    int     pruned;
    int     unpruned;
} procstat_t;


/* Statistics about each panel. */

typedef struct {
    int    size;      /* size of the panel */
    int    pnum;      /* which processor grabs this panel */
    double starttime; /* at what time this panel is assigned to a proc */
    double fctime;    /* factorization time */
    float  flopcnt;   /* floating-point operations */
    int    pipewaits; /* how many times the panel waited during pipelining */
    double spintime;  /* spin waiting time during pipelining */
} panstat_t;

/* How was a panel selected by the scheduler */
typedef enum {NOPIPE, DADPAN, PIPE} how_selected_t;


typedef struct {
     flops_t flops;
     int     nzs;
     double  fctime;
} stat_relax_t;

typedef struct {
     flops_t flops;
     int nzs;
     double fctime;
} stat_col_t;

typedef struct {
     int ncols;
     flops_t flops;
     int nzs;
     double fctime;
} stat_snode_t;

/* -------------------------------------------------------------------
   The definitions below are used to simulate parallel execution time.
   ------------------------------------------------------------------- */
typedef struct {
    float est;  /* earliest (possible) start time of the panel */
    float pdiv; /* time in flops spent in the (inner) panel factorization */
} cp_panel_t;

typedef struct {
    float eft;  /* earliest finishing time */
    float pmod; /* pmod update to the ancestor panel */
} desc_eft_t;
		   
/* All statistics. */
typedef struct {
    int     	*panel_histo;	/* Panel size distribution */
    double  	*utime;
    flops_t 	*ops;
    procstat_t 	*procstat;
    panstat_t	*panstat;
    int      	num_panels;
    float     	dom_flopcnt;
    float     	flops_last_P_panels;
    /**/
    stat_relax_t *stat_relax;
    stat_col_t *stat_col; 
    stat_snode_t *stat_snode; 
    int *panhows;
    cp_panel_t *cp_panel; /* panels on the critical path */
    desc_eft_t *desc_eft; /* all we need to know from descendants */
    int        *cp_firstkid, *cp_nextkid; /* linked list of children */
    int        *height;
    float      *flops_by_height;
} Gstat_t;

struct Branch {
    int root, first_desc, which_bin;
    struct Branch *next;
};


#if 0

/* Statistics for supernode and panel size */
int 	no_panels;
float   sum_w;          /* Sum (Wi) */
float 	sum_np_w;       /* Sum (Npi*Wi) */
int 	max_np;          
int     no_sups;
float   sum_sup;        /* Sum (Supi) */
int     max_sup;     
flops_t reuse_flops;    /* Triangular solve and matrix vector multiply */
float   reuse_data;     /* Doubles in updating supernode */

/* Statistics for blas operations */
int     num_blas;       /* no of BLAS2 operations, including trsv/gemv */
int     max_blas_n;     /* max dimension n in tri-solve and mat-vec */
int     min_blas_n;     /* min dimension n in tri-solve and mat-vec */
float   sum_blas_n;     /* sum of "        "        " */
int     max_gemv_m;     /* max dimension m in mat-vec */
int     min_gemv_m;     /* max dimension m in mat-vec */
float   sum_gemv_m;     /* sum of "        "        " */
int     lda_blas_m;
int     lda_blas_n;
flops_t *gemv_ops;      /* flops distribution on (m,n) */
flops_t *trsv_ops;      /* flops distribution on n */

#define i_trsv_ops(i)      trsv_ops[i]
#define ij_gemv_ops(i,j)   gemv_ops[j*lda_blas_m + i]

#endif

#endif /* __SUPERLU_UTIL */

