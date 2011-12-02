/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * blas.c
 *
 * This file contains functions that implement psudo BLAS functions
 * The functions that are included are:
 *   *asum, *set, *axpy, *copy, *scale, *nrm2, *dot
 *
 * Started 9/15/94
 * George
 *
 * $Id: GKlib.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "GKlib.h"

/*************************************************************************
* These functions add the elements of a vector
**************************************************************************/
int iasum(int n, int *x)
{
  int sum = 0;
  int i;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

float sasum(int n, float *x)
{
  float sum = 0.0;
  int i;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

double dasum(int n, double *x)
{
  double sum = 0.0;
  int i;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function computes the saxpy operation  y = y +ax
**************************************************************************/
void daxpy(int n, double a, double *x, double *y)
{
  int i;

  for (i=0; i<n; i++)
    y[i] += a*x[i];
}


/*************************************************************************
* These functions copies a vector x to y
**************************************************************************/
void icopy(int n, int *x, int *y)
{
  memcpy((void *)y, (void *)x, sizeof(int)*n);
}

void dcopy(int n, double *x, double *y)
{
  memcpy((void *)y, (void *)x, sizeof(double)*n);
}


/*************************************************************************
* This function performs the dot product of two vectors
**************************************************************************/
double ddot(int n, double *x, double *y)
{
  int i;
  double dot=0.0;

  for (i=0; i<n; i++)
    dot += x[i]*y[i];

  return dot;
}


/*************************************************************************
* This function computes the 2-norm
**************************************************************************/
double dnrm2(int n, double *x)
{
  int i;
  double norm=0.0;

  for (i=0; i<n; i++) 
    norm += x[i]*x[i];

  return sqrt(norm);
}


/*************************************************************************
* This function scales a vector by a consant
**************************************************************************/
void dscal(int n, double a, double *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = a*x[i];
}


/*************************************************************************
* This function swaps two vectors
**************************************************************************/
void dswap(int n, double *x, double *y)
{
  int i;
  double temp;

  for (i=0; i<n; i++) {
    temp = x[i];
    x[i] = y[i];
    y[i] = temp;
  }

}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
void iset(int n, int val, int *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int iamax(int n, int *x)
{
  int max = 0;
  int i;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int iamin(int n, int *x)
{
  int min = 0;
  int i;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}



/*
 * file.c
 *
 * This file contains some simple io functions
 *
 * Started 4/10/95
 * $Id: GKlib.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

/*************************************************************************
* This function opens a file
**************************************************************************/
FILE *GKfopen(char *fname, char *mode, char *msg)
{
  FILE *fp;
  char errmsg[256];

  fp = fopen(fname, mode);
  if (fp != NULL)
    return fp;

  sprintf(errmsg,"file: %s, mode: %s, [%s]", fname, mode, msg);
  perror(msg);
  exit(0);
}


void GKfclose(FILE *fp)
{
  fclose(fp);
}
/*
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 8/27/94
 * George
 *
 * $Id: GKlib.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *imalloc(int n, char *msg)
{
  int *ptr;

  if (n == 0)
    return NULL;

  ptr = (int *)GKmalloc(sizeof(int)*n, msg);

  return ptr;

}

/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *ismalloc(int n, int ival, char *msg)
{
  int *ptr;

  if (n == 0)
    return NULL;

  ptr = (int *)GKmalloc(sizeof(int)*n, msg);
  iset(n, ival, ptr);

  return ptr;
}


/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *GKmalloc(int nbytes, char *msg)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL) 
    errexit("Memory allocation failed for %s. Requested size: %d bytes\n", msg, nbytes);

  return ptr;
}


/*************************************************************************
* This function is my wrapper around free, allows multiple pointers
**************************************************************************/
void GKfree(void *ptr1,...)
{
  va_list plist;
  void *ptr;

  if (ptr1 != NULL)
    free(ptr1);
  ptr1 = NULL;

  va_start(plist, ptr1);

  while ((int)(ptr = va_arg(plist, void *)) != -1) {
    if (ptr != NULL)
      free(ptr);
    ptr = NULL;
  }

  va_end(plist);

}

/*
 * sort.c
 *
 * This function implements various sorting algorithms
 *
 * Started 10/27/94
 * George
 *
 */

static int incint(const void *, const void *);
static int decint(const void *, const void *);


/*************************************************************************
* These functions sorts an array of XXX
**************************************************************************/
void iincsort(int n, int *a)
{
  qsort((void *)a, (size_t)n, (size_t)sizeof(int), incint);
}


void idecsort(int n, int *a)
{
  qsort((void *)a, (size_t)n, (size_t)sizeof(int), decint);
}




/*************************************************************************
* This function compares 2 ints for sorting in inc order
**************************************************************************/
static int incint(const void *v1, const void *v2)
{
  return (*((int *)v1) - *((int *)v2));
}

/*************************************************************************
* This function compares 2 ints for sorting in dec order
**************************************************************************/
static int decint(const void *v1, const void *v2)
{
  return (*((int *)v2) - *((int *)v1));
}


#ifndef METISLIB
/*
 * timer.c
 *
 * This file contains function that implements timers on any UNIX system. 
 * This timers report the actual user and system time spent in the program
 *
 * Started 1/21/94
 * George
 */

/*************************************************************************
* This function clears/initializes a timer
**************************************************************************/
void cleartimer(timer *tmr)
{

  tmr->time = 0.0;
  tmr->count = 0L;
  tmr->status = TMR_CLEAR;
}


/*************************************************************************
* This function starts a timer
**************************************************************************/
void starttimer(timer *tmr)
{
  if (tmr->status == TMR_START) 
    errexit("Timer already started!\n");

  tmr->time -= seconds();
  tmr->count++;
  tmr->status = TMR_START;
}


/*************************************************************************
* This function stops a timer
**************************************************************************/
void stoptimer(timer *tmr)
{
  if (tmr->status != TMR_START) 
    errexit("Timer has not been started!\n");

  tmr->time += seconds();
  tmr->status = TMR_STOP;

}


/*************************************************************************
* This function prints the contents of a timer 
**************************************************************************/
void printtimer(timer *tmr, char *title)
{
  if (tmr->status == TMR_START) 
    errexit("Timing is in progress!\n");

  printf("\n%s: \t[%6ld, %7.3f]", title, tmr->count, tmr->time);
}

/*************************************************************************
* This function returns the value of a timer 
**************************************************************************/
double gettimer(timer *tmr)
{
  if (tmr->status == TMR_START) 
    errexit("Timing is in progress!\n");

  return tmr->time;
}



/*************************************************************************
* This function returns the seconds
**************************************************************************/
double seconds(void)
{
#ifdef XYZ
  double time;                         /* elapsed time in seconds */
  struct rusage rusage;

  getrusage(RUSAGE_SELF, &rusage);
  time = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
         1.0e-6*(rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));
  return(time);
#endif
  return((double) clock()*1.0e-6);
}
#endif

/*
 * util.c
 *
 * This file contains utility functions
 *
 * Started 8/28/94
 * George
 *
 * $Id: GKlib.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void *errexit(char *f_str,...)
{
  va_list argp;

  va_start(argp, f_str);
  vfprintf(stderr, f_str, argp);
  va_end(argp);

  fprintf(stderr,"\n");
  fflush(stderr);

  exit(0);
}

