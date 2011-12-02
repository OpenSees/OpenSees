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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-03-04 00:38:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/FirstPrincipalCurvature.cpp,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <FirstPrincipalCurvature.h>
#include <FindCurvatures.h>
#include <LimitStateFunction.h>
#include <Vector.h>


FirstPrincipalCurvature::FirstPrincipalCurvature()
:FindCurvatures(), curvatures(1)
{
}

FirstPrincipalCurvature::~FirstPrincipalCurvature()
{
}


int
FirstPrincipalCurvature::computeCurvatures(ReliabilityDomain *theReliabilityDomain)
{

	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);

	// Get hold of 'u' and 'alpha' at the two last steps
	Vector last_u = theLimitStateFunction->designPoint_u_inStdNormalSpace;
	Vector secondLast_u = theLimitStateFunction->secondLast_u;
	Vector lastAlpha = theLimitStateFunction->normalizedNegativeGradientVectorAlpha;
	Vector secondLastAlpha = theLimitStateFunction->secondLastAlpha;

	// Compute curvature according to Der Kiureghian & De Stefano (1992), Eq.26:

	// Initial computations
	Vector uLastMinus_u = last_u - secondLast_u;
	double signumProduct = secondLastAlpha ^ uLastMinus_u;
	double alphaProduct = secondLastAlpha ^ lastAlpha;
	double sumSquared = 0.0;

	// Compute norm of the difference vector
	for ( int i=0; i<last_u.Size(); i++ ) {
		sumSquared += uLastMinus_u(i)*uLastMinus_u(i);
	}

	double norm_uLastMinus_u = sqrt(sumSquared);

	// Check sign and compute curvature
	if (fabs(signumProduct)==(signumProduct)) {
		curvatures(0) = acos(alphaProduct) / norm_uLastMinus_u;
	}
	else {
		curvatures(0) = -acos(alphaProduct) / norm_uLastMinus_u;
	}

	return 0;
}


Vector
FirstPrincipalCurvature::getCurvatures()
{
	return curvatures;
}


