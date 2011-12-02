/*
 * File:  nmat.h
 * =============
 * altered to improve data access
 * Written by:  David Mackay
 *		Kincho Law
 *		Jun Peng
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
