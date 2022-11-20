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

//This file contains the definition for NURBS derivatives
// Written originally by Vinh Phu Nguyen, nvinhphu@gmail.com


#ifndef NurbsDers_h
#define NurbsDers_h

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Subdomain.h>
#include <OPS_Globals.h>

// int      FindSpan(int n, int p, double u, Vector U);
int      FindSpan(int n, int p, double u, Vector& U);
// void     BasisFuns( int i, double u, int p, Vector U, double* N);
void     BasisFuns( int i, double u, int p, Vector& U, Vector& N);
// void     dersBasisFuns(int i, double u, int p, int order, Vector knot, double **ders);
void     dersBasisFuns(int i, double u, int p, int order, Vector& knot, Matrix& ders);
double   OneBasisFun(int p, int m, Vector U, int i, double u);
void     dersOneBasisFuns(int p, int m, Vector U, int i, double u, int n, double* ders);
double** init2DArray(int x, int y);
void     free2Darray(double **array, int x);

#endif