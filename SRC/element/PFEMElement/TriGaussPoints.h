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

// $Revision $
// $Date$
// $URL$

// Written: Minjie Zhu (zhum@oregonstate.edu)
//
// Description: Gaussian points and weights for triangles
// n - the order of the Gaussian quadrature (n<=12)
// xpts,ypts,wts - x,y coordinates and weights

#ifndef TriGaussPoints_H
#define TriGaussPoints_H

#include <BackgroundDef.h>

class TriGaussPoints
{
public:
    TriGaussPoints();
    ~TriGaussPoints();

    void operator()(int n, VDouble& xpts, VDouble& ypts, VDouble& wts);
};



#endif
