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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-04-04 16:55:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/matrixutil/MatrixUtil.cpp,v $
 
#include <math.h>
                                                                        
#include <stdlib.h>
#include <Vector.h>
#include <MatrixUtil.h>


double invert2by2Matrix(const Matrix &a, Matrix &b)
{
  return a.Invert(b);

  /*
    double det = a(0,0)*a(1,1)-a(0,1)*a(1,0);
    
    if (det != 0.0) {
    b(0,0) =  a(1,1)/det;
    b(1,1) =  a(0,0)/det;
    b(1,0) = -a(0,1)/det;
    b(0,1) = -a(1,0)/det;
    }
    return det;
  */
}



double invert3by3Matrix(const Matrix &a, Matrix &b)
{
  return a.Invert(b);

    /*
   double a00, a01, a02, a11, a12, a22;
   double a01a01, a02a02, a12a12, a00a11, a01a02;
   double det;
               
   a00 = a(0,0);
   a01 = a(0,1);
   a02 = a(0,2);
   a11 = a(1,1);
   a12 = a(1,2);
   a22 = a(2,2);
   a01a01 = a01 * a01;
   a02a02 = a02 * a02;
   a12a12 = a12 * a12;
   a00a11 = a00 * a11;
   a01a02 = a01 * a02;
 
   det =  a00a11*a22 + 2 * a01a02*a12 - (a00*a12a12 + a01a01*a22 + a02a02*a11);

   if (det != 0) {
      b(0,0) = (a11*a22 - a12a12 )/det;
      b(0,1) = (a02*a12 - a01*a22)/det;    
      b(0,2) = (a01*a12 - a02*a11)/det;
      b(1,0) = b(0,1);
      b(1,1) = (a00*a22 - a02a02 )/det;
      b(1,2) = (a01a02  - a00*a12)/det;
      b(2,0) = b(0,2);
      b(2,1) = b(1,2);
      b(2,2) = (a00a11  - a01a01 )/det;
   }
   return det;
    */
}
      

  
      
void invertMatrix(int n, const Matrix &a, Matrix &b)
{
  a.Invert(b);
}



void getCBDIinfluenceMatrix(int nIntegrPts, const Matrix &xi_pt, double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   int i, j, i0, j0;
   double xi;
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nIntegrPts, nIntegrPts);
   Matrix I(nIntegrPts,nIntegrPts);      // an identity matrix for matrix inverse

   for (i = 1; i <= nIntegrPts; i++)
      for (j = 1; j <= nIntegrPts; j++)
      {
         i0 = i - 1;
         j0 = j - 1;
         xi = xi_pt(i0,0);
         G(i0,j0) =  pow(xi,j-1);
         l(i0,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
      }
   
   I.Zero();
   for (i=0; i<nIntegrPts; i++)
     I(i,i) = 1.0;

   //invertMatrix(nIntegrPts, G, Ginv);
   if (G.Solve(I,Ginv) < 0)
     opserr << "LargeDispBeamCol3d::getCBDIinfluenceMatrix() - could not invert G\n";
      
   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}

void getCBDIinfluenceMatrix(int nIntegrPts, double *pts, double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   int i, j, i0, j0;
   double xi;
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nIntegrPts, nIntegrPts);
   Matrix I(nIntegrPts,nIntegrPts);      // an identity matrix for matrix inverse

   for (i = 1; i <= nIntegrPts; i++)
      for (j = 1; j <= nIntegrPts; j++)
      {
         i0 = i - 1;
         j0 = j - 1;
         xi = pts[i0];
         G(i0,j0) =  pow(xi,j-1);
         l(i0,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
      }
   
   I.Zero();
   for (i=0; i<nIntegrPts; i++)
     I(i,i) = 1.0;

   //invertMatrix(nIntegrPts, G, Ginv);
   if (G.Solve(I,Ginv) < 0)
     opserr << "LargeDispBeamCol3d::getCBDIinfluenceMatrix() - could not invert G\n";
      
   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}





