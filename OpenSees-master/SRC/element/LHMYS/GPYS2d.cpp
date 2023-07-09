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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-11-08 20:06:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/StringContainer.cpp,v $
//
// Written: Jade Cohen
// Created: 19/05
//
// Description: This file contains the method for GPYS2d.
// GPYS2d is used to build general polynomial yield surface for 2d
// frame elements.
//

#include "GPYS2d.h"
#include <Matrix.h>
#include <Vector.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

void GPYS2d(const Matrix &GPYSC, const Vector &xyref, const Vector &ScVec, double &f, Vector &g, Matrix &h)
{
    if (GPYSC.noCols() != 3) {
        std::cerr<<"GPYSC must have 3 columns\n";
        return;
    }
    
    if ((xyref.Size() != 2) || (ScVec.Size() != 2) || (g.Size() != 2)) {
        std::cerr<<"xy, ScVec and g must be vectors of size 2\n";
        return;
    }
    
    if (h.noCols() != 2 || h.noRows() != 2){
        std::cerr<<"h mut be a 2x2 matrix\n";
        return;
    }
    
    if ((ScVec(0) < 0.0) || (ScVec(1) < 0.0)){
        std::cerr<<"Scaling factors must be positive numbers\n";
        return ;
    }
    
    // Get no. of columns and rows
    int ncol = 2;
    int nrow = GPYSC.noRows();
    
    // Create copy of xy vector
    Vector xy(xyref);
    
    // Scale xy
    for (int i=0; i<ncol; i++){
        xy(i) = xy(i)/ScVec(i);
    }
    
    // Perturbation of Hessian when x = 0 or y = 0 raised to power < 2
    for (int i=0; i<3; i++){
        if ((GPYSC(i,1) < 2) || (GPYSC(i,2) <2)){
            if (xyref(0) == 0.0)
                xy(0) = 1e-6;
            else if (xyref(1) == 0.0)
                xy(1) = 1e-6;
        }
    }
    
    double x = fabs(xy(0));
    double y = fabs(xy(1));
    
    // Get f
    f = 0.0;
    for (int i=0; i<nrow; i++){
        f += GPYSC(i,0)*pow(x,GPYSC(i,1))*pow(y,GPYSC(i,2));
    }
    
    // Get gradient g
    double dfdx = 0.0;
    double dfdy = 0.0;
    for (int i=0; i<nrow; i++){
        if (GPYSC(i,1) != 0.0){
            dfdx += GPYSC(i,0)*GPYSC(i,1)*pow(x,GPYSC(i,1)-1.0)*pow(y,GPYSC(i,2));
        }
        if (GPYSC(i,2) != 0.0){
            dfdy += GPYSC(i,0)*GPYSC(i,2)*pow(x,GPYSC(i,1))*pow(y,GPYSC(i,2)-1.0);
        }
    }
    
    if (xyref(0) >= 0.0)
        g(0) = dfdx/ScVec(0);
    else
        g(0) = -dfdx/ScVec(0);
    
    if (xyref(1) >= 0.0)
        g(1) = dfdy/ScVec(1);
    else
        g(1) = -dfdy/ScVec(1);
    
    // Get hessian h
    Matrix GPYSCx(GPYSC);
    Matrix GPYSCy(GPYSC);
    for (int i=0; i<nrow; i++){
        GPYSCx(i,0) = GPYSC(i,0)*GPYSC(i,1);
        GPYSCx(i,1) = GPYSC(i,1)-1.0;
        GPYSCy(i,0) = GPYSC(i,0)*GPYSC(i,2);
        GPYSCy(i,2) = GPYSC(i,2)-1.0;
    }
    Matrix GPYSCxx(GPYSCx);
    Matrix GPYSCyy(GPYSCy);
    Matrix GPYSCxy(GPYSCx);
    for (int i=0; i<nrow; i++){
        GPYSCxx(i,0) = GPYSCx(i,0)*GPYSCx(i,1);
        GPYSCxx(i,1) = GPYSCx(i,1)-1.0;
        GPYSCyy(i,0) = GPYSCy(i,0)*GPYSCy(i,2);
        GPYSCyy(i,2) = GPYSCy(i,2)-1.0;
        GPYSCxy(i,0) = GPYSCx(i,0)*GPYSCx(i,2);
        GPYSCxy(i,2) = GPYSCx(i,2)-1.0;
    }
    
    double d2fdx2 = 0.0;
    double d2fdy2 = 0.0;
    double d2fdxdy = 0.0;
    for (int i=0; i<nrow; i++){
        if ((GPYSC(i,1) != 0.0) && (GPYSCx(i,1) != 0.0))
            d2fdx2 += GPYSCxx(i,0)*pow(x,GPYSCxx(i,1))*pow(y,GPYSCxx(i,2));
        if ((GPYSC(i,2) != 0.0) && (GPYSCy(i,2) != 0.0))
            d2fdy2 += GPYSCyy(i,0)*pow(x,GPYSCyy(i,1))*pow(y,GPYSCyy(i,2));
        if ((GPYSC(i,1) != 0.0) && (GPYSCx(i,2) != 0.0))
            d2fdxdy += GPYSCxy(i,0)*pow(x,GPYSCxy(i,1))*pow(y,GPYSCxy(i,2));
    }
    
    h(0,0) = d2fdx2/(ScVec(0)*ScVec(0));
    h(1,1) = d2fdy2/(ScVec(1)*ScVec(1));
    if (xyref(0)*xyref(1) >= 0.0)
        h(0,1) = d2fdxdy/(ScVec(0)*ScVec(1));
    else
        h(0,1) = -d2fdxdy/(ScVec(0)*ScVec(1));
    h(1,0)=h(0,1);
    
}
