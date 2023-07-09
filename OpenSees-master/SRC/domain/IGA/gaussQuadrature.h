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

// The function quadrature returns a n x 1 column vector W of quadrature
// weights and a n x dim matrix of quadrature points, where n is the
// number of quadrature points.  The function is called as follows:
// 
// [W,Q]=quadrature( nint, type, dim )
// 
// nint is the quadrature order, type is the type of quadrature
// (i.e. gaussian, triangular, etc.. ) and dim is the number of spacial
// dimentions of the problem.  The default for type is GAUSS and the
// default for dim is unity.
// 
// wrQ=quadrature(nint,'TRIANGULAR',2);itten by Jack Chessa
//            j-chessa@northwestern.edu
// Department of Mechanical Engineering
// Northwestern University

// Adapted to OpenSees by Felipe Elgueta and José A. Abell (UANDES, Chile) www.joseabell.com
// Only using gaussian quadrature


#ifndef gaussQuadrature_h
#define gaussQuadrature_h

#include <math.h>      

void gaussQuad(int order, Vector* pt, Vector* wt);
void gaussQuad2dNurbs(int orderU, int orderV, Matrix* quadPoint, Vector* quadWeight);

#endif
