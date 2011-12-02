/*
 * GKlib.h
 * 
 * George's library of most frequently used routines
 *
 * Started 4/10/95
 * $Id: GKlib.h,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#ifndef _GKLIB_H_
#define _GKLIB_H_

#include <rename.h>

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define TMR_CLEAR   10
#define TMR_START   20
#define TMR_STOP    30


/*************************************************************************
* Data structures
**************************************************************************/
struct timer {
  double time;   /* The time in seconds */
  long count;    /* Count how many times the counter was started */
  int status;    /* Keeps track what the timer is doing */
};
typedef struct timer timer;



/*************************************************************************
* Nice Macros
**************************************************************************/
#define log2(x) ((int) (log((double)(x)) / log(2.0)))
#define sign(a, b) ((b) >= 0 ? ((a) >= 0.0 ? a : -a) : ((a) >= 0.0 ? -a : a))
#define amax(a, b) ((a) >= (b) ? (a) : (b))
#define amin(a, b) ((a) >= (b) ? (b) : (a))


#ifdef GK_ASSERT
#   define ASSERT(expr)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define ASSERT(expr) ;
#endif 

#ifdef GK_ASSERT
#   define ASSERTP(expr,msg)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        printf msg ; \
        printf("\n"); \
        abort();                                                \
    }
#else
#   define ASSERTP(expr,msg) ;
#endif 


/*************************************************************************
* Function prototypes
**************************************************************************/
/* blas.c */
int  iasum(int, int *);
void icopy(int, int *, int *);
void iset(int, int, int *);
int  iamax(int, int *);
int  iamin(int, int *);

float sasum(int, float *);

double dasum(int, double *);
void   daxpy(int, double, double *, double *);
void   dcopy(int, double *, double *);
double ddot(int, double *, double *);
double dnrm2(int, double *);
void   dscal(int, double, double *);
void   dswap(int, double *, double *);


/* file.c */
FILE *GKfopen(char *, char *, char *);
void GKfclose(FILE *);

/* memory.c */
int *imalloc(int, char *);
int *ismalloc(int, int, char *);
void *GKmalloc(int, char *);
void GKfree(void *,...);


/* sort.c */
void iincsort(int, int *);
void idecsort(int, int *);

/* timer.c */
#ifndef METISLIB
void cleartimer(timer *);
void starttimer(timer *);
void stoptimer(timer *);
void printtimer(timer *, char *);
double gettimer(timer *);
double seconds(void);
#endif

/* util.c */
void *errexit(char *,...);


#endif
