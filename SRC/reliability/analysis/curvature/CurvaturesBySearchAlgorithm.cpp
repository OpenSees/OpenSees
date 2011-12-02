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
// $Date: 2003-03-04 00:38:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/CurvaturesBySearchAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <CurvaturesBySearchAlgorithm.h>
#include <FindCurvatures.h>
#include <LimitStateFunction.h>
#include <FindDesignPointAlgorithm.h>
#include <RandomVariable.h>
#include <Vector.h>


CurvaturesBySearchAlgorithm::CurvaturesBySearchAlgorithm(int passedNumberOfCurvatures,
									FindDesignPointAlgorithm *passedFindDesignPointAlgorithm)
:FindCurvatures(), curvatures(passedNumberOfCurvatures)
{
	numberOfCurvatures = passedNumberOfCurvatures;
	theFindDesignPointAlgorithm = passedFindDesignPointAlgorithm;
}

CurvaturesBySearchAlgorithm::~CurvaturesBySearchAlgorithm()
{
}


int
CurvaturesBySearchAlgorithm::computeCurvatures(ReliabilityDomain *theReliabilityDomain)
{

	// Initial declarations
	Vector last_u;
	Vector secondLast_u;
	Vector lastAlpha;
	Vector secondLastAlpha;
	Vector lastDirection;
	Vector uLastMinus_u = last_u - secondLast_u;
	double signumProduct;
	double alphaProduct;
	double sumSquared;
	double norm_uLastMinus_u;

	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);

	// The design point in the original space
	Vector x_star = theLimitStateFunction->designPoint_x_inOriginalSpace;

	// Number of random variables
	int nrv = x_star.Size();

	// Declare vector to store all principal axes
	Vector principalAxes( nrv * numberOfCurvatures );

	// Start point of new searches: Perturb x^star by 10% of each standard deviation
	Vector startPoint(nrv);
	RandomVariable *aRandomVariable;
	int numberOfRandomVariables = theReliabilityDomain->getNumberOfRandomVariables();
	int i;
	for (i=0; i<nrv; i++) {
		aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i+1);
		startPoint(i) = aRandomVariable->getStartValue();
	}

	// Compute curvatures by repeated searches
	for (i=0; i<numberOfCurvatures; i++) {

		// Get hold of 'u', 'alpha' and search direction at the two last steps
		last_u = theLimitStateFunction->designPoint_u_inStdNormalSpace;
		secondLast_u = theLimitStateFunction->secondLast_u;
		lastAlpha = theLimitStateFunction->normalizedNegativeGradientVectorAlpha;
		secondLastAlpha = theLimitStateFunction->secondLastAlpha;
		lastDirection = theLimitStateFunction->lastSearchDirection;

		// Compute curvature according to Der Kiureghian & De Stefano (1992), Eq.26:

		// Initial computations
		uLastMinus_u = last_u - secondLast_u;
		signumProduct = secondLastAlpha ^ uLastMinus_u;
		alphaProduct = secondLastAlpha ^ lastAlpha;
		sumSquared = 0.0;

		// Compute norm of the difference vector
		for ( int k=0; k<last_u.Size(); k++ ) {
			sumSquared += uLastMinus_u(k)*uLastMinus_u(k);
		}
		norm_uLastMinus_u = sqrt(sumSquared);

		// Check sign and compute curvature
		if (fabs(signumProduct)==(signumProduct)) {
			curvatures(i) = acos(alphaProduct) / norm_uLastMinus_u;
		}
		else {
			curvatures(i) = -acos(alphaProduct) / norm_uLastMinus_u;
		}

		// Run a new search in the subspace orthogonal to previous principal directions
		if ( i!=(numberOfCurvatures-1) ) {

			// Store all previous principal axes in a vector
			for (int j=0; j<nrv; j++ ) {
				principalAxes((i*nrv)+j) = lastDirection(j);
			}

			// To be completed...
		}
	}

	return 0;
}



Vector
CurvaturesBySearchAlgorithm::getCurvatures()
{
	return curvatures;
}



