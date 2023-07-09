/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2003-03-04 00:38:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/FindCurvatures.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <FindCurvatures.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>

FindCurvatures::FindCurvatures()
{
}

FindCurvatures::~FindCurvatures()
{
}


int
FindCurvatures::gramSchmidt(const Vector &first, Matrix &R)
{
    // simple GS orthogonalization of input (row) vector first
    int n = first.Size();
    R.resize(n,n);
    R.Zero();
    Matrix R0(R);
    Matrix temp(R);
    
    for (int i = 0; i < n; i++) {
        R0(i,i) = 1;
        R0(n-1,i) = first(i);
        R(n-1,i) = first(i);
    }
    
    // computation of the rows of R
    for (int k = n-2; k >= 0; k--) {
        temp.Zero();
        for (int j = k+1; j < n; j++) {
            double numer = 0;
            double denom = 0;
            for (int i = 0; i < n; i++) {
                numer += R(j,i)*R0(k,i);
                denom += R(j,i)*R(j,i);
            }
            
            for (int i = 0; i < n; i++)
                temp(j,i) = numer/denom * R(j,i);
        }
        
        // sum columns
        for (int j = 0; j < n; j++) {
            double sum_j = 0;
            for (int i = 0; i < n; i++)
                sum_j += temp(i,j);
        
            R(k,j) = R0(k,j) - sum_j;
        }
    }
    
    // normalization of each rows
    for (int k = 0; k < n; k++) {
        double sum = 0;
        for (int j = 0; j < n; j++)
            sum += R(k,j)*R(k,j);
        
        for (int j = 0; j < n; j++)
            R(k,j) = R(k,j)/sqrt(sum);

    }
    
    return 0;
}