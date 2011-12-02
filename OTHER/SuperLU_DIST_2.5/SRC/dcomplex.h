/*! @file
 * \brief Header for dcomplex.c
 *
 * <pre>
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 * </pre>
 */

/* 
 * This header file is to be included in source files z*.c
 */
#ifndef __SUPERLU_DCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_DCOMPLEX

#include <mpi.h>

typedef struct { double r, i; } doublecomplex;

/*
 * These variables will be defined to be MPI datatypes for complex
 * and double complex. I'm too lazy to declare
 * these guys external in every file that needs them.
 */
extern MPI_Datatype SuperLU_MPI_DOUBLE_COMPLEX;


/* Macro definitions */

/*! \brief Complex Addition c = a + b */
#define z_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/*! \brief Complex Subtraction c = a - b */
#define z_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/*! \brief Complex-Double Multiplication */
#define zd_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/*! \brief Complex-Complex Multiplication */
#define zz_mult(c, a, b) { \
	double cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

/*! \brief Complex equality testing */
#define z_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in dcomplex.c */
void   slud_z_div(doublecomplex *, doublecomplex *, doublecomplex *);
double slud_z_abs(doublecomplex *);     /* exact */
double slud_z_abs1(doublecomplex *);    /* approximate */


#ifdef __cplusplus
  }
#endif


#endif  /* __SUPERLU_DCOMPLEX */
