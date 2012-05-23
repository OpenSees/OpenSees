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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-03-13 22:29:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/FindCurvatures.h,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef FindCurvatures_h
#define FindCurvatures_h

#include <Vector.h>
#include <Matrix.h>

class FindCurvatures
{

public:
    FindCurvatures();
    virtual ~FindCurvatures();

    virtual int		computeCurvatures() =0;
    virtual const Vector	&getCurvatures() =0;
    virtual const Vector	&getPrincipalAxes() =0;
    
    int gramSchmidt(const Vector &first, Matrix &R);
  
protected:

private:

};

#endif
