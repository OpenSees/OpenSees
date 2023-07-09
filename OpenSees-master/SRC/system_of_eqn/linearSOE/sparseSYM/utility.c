/*
 * File:  utility.c
 * ================
 * This file contains miscellaneous utility routines used throughout
 * the finite element program.
 *
 * Originally written by:  David R. Mackay
 *
 * Modified by:
 *  Jun Peng (junpeng@stanford.edu)
 *  Prof. Kincho H. Law
 *  Stanford University
 * --------------------
 */
    

#include <stdio.h>
#include <stdlib.h>
#include "utility.h"


void move_real(double *from, double *to, int n)
{
	double *last= to + n;

	while(to < last) *to++ = *from++;
}
/************************* end of function *************************/  

double  dot_real(double *vect_1, double *vect_2, int n)
{
	double *fstop ;
	double  sum    ;

	sum = 0.0  ;
	fstop = vect_1 + n ;
	for( ; vect_1 <fstop ; vect_1++, vect_2++ )
		sum += (*vect_1 * *vect_2) ;

	return(sum) ;
}
/************************* end of function *************************/

int  i_greater(int *p1, int *p2)
{
	return ( (int)(*p1 - *p2) );
}
/************************* end of function *************************/

void  zeroi(int n, int *v)
{
	int  *end ;

	end = v+n;
	for ( ; v < end ; v++)
		*v &= 0;
	return;
}
/************************* end of function *************************/

void  minoni(int n, int *v)
{
	int  *end;
	int  min;

	min = -1;
	end = v+n;
	for ( ; v < end ; v++)
		*v |= min;
	return;
}
/************************* end of function *************************/

void  saxpy(double *v1, double *v2, double alpha, int n)
{
	double  *end;

	end = v1 + n;
	if(n <= 0)
	{
		printf(" n %d\n", n);
		exit(1);
	}
	for( ; v1 < end ; v1++, v2++)
		*v1 += *v2 * alpha;
   
	return;
}
/************************* end of function *************************/

void  copyi(int n, int *from, int *to)
{
	int  i;

	for(i=0 ; i < n ; i++, to++, from++)
		*to = *from;
	return;
}
/************************* end of function *************************/

int  rcomp(double *p1, double *p2)
{
	return( (*p1 < *p2) ? -1 : (*p1 > *p2) );
}
/************************* end of function *************************/


/* UNUSED FUNCTIONS */

void clear_real(double *array, int n)
{
	double *last= array + n;

	while(array < last) *array++ = 0.0;
}
/************************* end of function *************************/

int  icomp(int *p1, int *p2)
{
	return( (*p1 < *p2) ? -1 : (*p1 > *p2) );
}
/************************* end of function *************************/

int fcomp(float *p1, float *p2)
{
	return( (*p1 < *p2) ? -1 : (*p1 > *p2) );
}
