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
                                                                        
// $Revision: 1.6 $
// $Date: 2003-10-27 23:45:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/FixedStepSizeRule.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <FixedStepSizeRule.h>
#include <StepSizeRule.h>
#include <ProbabilityTransformation.h>
#include <math.h>
#include <Vector.h>


FixedStepSizeRule::FixedStepSizeRule(double passedStepSize)
:StepSizeRule()
{
	stepSize = passedStepSize;
	gFunValue = -1;
}

FixedStepSizeRule::~FixedStepSizeRule()
{
}



int
FixedStepSizeRule::computeStepSize(Vector u, 
									Vector grad_G, 
									double G, 
									Vector d,
									int stepNumber)
{
	// This method is in fact not neccesary 
	// for the fixed step size rule. The 
	// user has already given the step size. 

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

double
FixedStepSizeRule::getGFunValue()
{
	return 0.0;
}

