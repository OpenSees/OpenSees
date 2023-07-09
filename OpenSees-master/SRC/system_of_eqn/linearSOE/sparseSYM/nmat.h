/*
 * File:  nmat.h
 * =============
 * altered to improve data access
 *
 * Originally written by:  David R. Mackay
 *
 * Modified by:
 *  Jun Peng (junpeng@stanford.edu)
 *  Prof. Kincho H. Law
 *  Stanford University
 * --------------------
 */


#ifndef nmat_h
#define nmat_h

int pfsfct(int neqns, double *diag, double **penv, int nblks, 
	   int *xblk, OFFDBLK **begblk, OFFDBLK *first, int *rowblks);

int pfefct(int neqns, double **penv, double *diag);

void pfsslv(int neqns, double *diag, double **penv, int nblks, 
	    int *xblk, double *rhs, OFFDBLK **begblk);

void pflslv (int neqns, double **penv, double *diag, double *rhs);

void pfuslv(int neqns, double **penv, double *diag, double *rhs);

#endif
