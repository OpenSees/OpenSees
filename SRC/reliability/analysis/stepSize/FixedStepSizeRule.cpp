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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-02-29 19:47:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/FixedStepSizeRule.cpp,v $


//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <FixedStepSizeRule.h>
#include <StepSizeRule.h>
#include <math.h>
#include <Vector.h>


FixedStepSizeRule::FixedStepSizeRule(double passedStepSize)
:StepSizeRule()
{
	stepSize = passedStepSize;
}

FixedStepSizeRule::~FixedStepSizeRule()
{
}



int
FixedStepSizeRule::computeStepSize(const Vector &u, 
				   const Vector &grad_G, 
				   double G, 
				   const Vector &d,
				   int stepNumber, int reschk)
{
	// This method is in fact not necessary for the fixed step size rule. The 
	// user has already given the step size. BUT limit the maximum jump that can occur 

    Vector u_temp(u);
    double udiff = 20.0;
    double stepSizeMod = 4./3.;
    
    // Determine new iteration point (take the step), BUT limit the maximum jump that can occur
    while ( fabs(udiff) > 15 ) {
        u_temp = u;
        stepSizeMod *= 3./4.;
        
        if (stepSizeMod < 1)
            opserr << "FixedStepSizeRule:: reducing stepSize using modification factor of " << stepSizeMod << endln;
        
        //u = u_old + (searchDirection * stepSize);
        u_temp.addVector(1.0, d, stepSize*stepSizeMod);
        
        udiff = u.Norm() - u_temp.Norm();
    }
    
    stepSize = stepSize*stepSizeMod;
    
	return 0;  

}


double
FixedStepSizeRule::getStepSize()
{
	return stepSize;

}


double
FixedStepSizeRule::getInitialStepSize()
{
	return stepSize;

}
