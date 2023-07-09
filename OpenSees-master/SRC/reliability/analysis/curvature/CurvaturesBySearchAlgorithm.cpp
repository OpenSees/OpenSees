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
// $Date: 2008-04-10 18:10:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/CurvaturesBySearchAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <CurvaturesBySearchAlgorithm.h>
#include <FindCurvatures.h>
#include <FORMAnalysis.h>
#include <LimitStateFunction.h>
#include <Vector.h>
#include <math.h>


CurvaturesBySearchAlgorithm::CurvaturesBySearchAlgorithm(ReliabilityDomain *passedReliabilityDomain,
                                                         FunctionEvaluator *passedFunctionEvaluator,
                                                         FORMAnalysis *passedFORMAnalysis,
                                                         int passedNumberOfCurvatures)
:FindCurvatures(), curvatures(passedNumberOfCurvatures)
{
    theReliabilityDomain = passedReliabilityDomain;
    theFunctionEvaluator = passedFunctionEvaluator;
    theFORMAnalysis = passedFORMAnalysis;
    numberOfCurvatures = passedNumberOfCurvatures;
}

CurvaturesBySearchAlgorithm::~CurvaturesBySearchAlgorithm()
{
    
}


int
CurvaturesBySearchAlgorithm::computeCurvatures()
{

	// Initial declarations
  /*
	Vector last_u;
	Vector secondLast_u;
	Vector lastAlpha;
	Vector secondLastAlpha;
	Vector lastDirection;
  */
  //Vector uLastMinus_u = last_u - secondLast_u;
	double signumProduct;
	double alphaProduct;
	double sumSquared;
	double norm_uLastMinus_u;

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	// "Download" limit-state function from reliability domain
	int lsfTag = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);

	// The design point in the original space
	// Recorders needed -- MHS 10/7/2011
	//const Vector &x_star = theLimitStateFunction->getFORM_x();
	Vector x_star(nrv);

	// Declare vector to store all principal axes
	Vector principalAxes( nrv * numberOfCurvatures );

	// Start point of new searches: Perturb x^star by 10% of each standard deviation
	Vector startPoint(nrv);
	RandomVariable *aRandomVariable;
	for (int i=0; i<nrv; i++) {
		aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i+1);
		startPoint(i) = aRandomVariable->getStartValue();
	}

	// Compute curvatures by repeated searches
	for (int i=0; i<numberOfCurvatures; i++) {

		// Get hold of 'u', 'alpha' and search direction at the two last steps
	  // Recorders needed -- MHS /10/7/2011
	  /*
	  const Vector &last_u = theLimitStateFunction->getFORM_u();
	  const Vector &secondLast_u = theLimitStateFunction->getSecondLast_u();
	  const Vector &lastAlpha = theLimitStateFunction->getFORM_alpha();
	  const Vector &secondLastAlpha = theLimitStateFunction->getSecondLast_alpha();
	  const Vector &lastDirection = theLimitStateFunction->getLastSearchDirection();
	  */
	  Vector last_u(nrv);
	  Vector secondLast_u(nrv);
	  Vector lastAlpha(nrv);
	  Vector secondLastAlpha(nrv);
	  Vector lastDirection(nrv);

		// Compute curvature according to Der Kiureghian & De Stefano (1992), Eq.26:

		// Initial computations
	  //Vector uLastMinus_u(nrv);
	  //uLastMinus_u = last_u - secondLast_u;
	  //signumProduct = secondLastAlpha ^ uLastMinus_u;
	  signumProduct = secondLastAlpha ^ last_u;
	  signumProduct -= secondLastAlpha ^ secondLast_u;
		alphaProduct = secondLastAlpha ^ lastAlpha;
		sumSquared = 0.0;

		// Compute norm of the difference vector
		for ( int k=0; k<last_u.Size(); k++ ) {
		  //sumSquared += uLastMinus_u(k)*uLastMinus_u(k);
		  double tmp = last_u(k)-secondLast_u(k);
		  sumSquared += tmp*tmp;
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



const Vector &
CurvaturesBySearchAlgorithm::getCurvatures()
{
	return curvatures;
}



const Vector &
CurvaturesBySearchAlgorithm::getPrincipalAxes()
{
	return principalAxes;
}
