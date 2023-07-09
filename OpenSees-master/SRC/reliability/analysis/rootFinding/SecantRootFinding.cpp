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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-11-06 19:32:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/rootFinding/SecantRootFinding.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <SecantRootFinding.h>
#include <RootFinding.h>
#include <FunctionEvaluator.h>
#include <ProbabilityTransformation.h>
#include <ReliabilityDomain.h>
#include <RandomVariable.h>
#include <math.h>
#include <Vector.h>


SecantRootFinding::SecantRootFinding(
						ReliabilityDomain *passedReliabilityDomain,
						ProbabilityTransformation *passedProbabilityTransformation,
						FunctionEvaluator *passedGFunEvaluator,
						int passedMaxIter,
						double ptol,
						double pmaxStepLength)
:RootFinding()
{
	theReliabilityDomain = passedReliabilityDomain;
	theProbabilityTransformation = passedProbabilityTransformation;
	theGFunEvaluator = passedGFunEvaluator;
	maxIter = passedMaxIter;
	tol = ptol;
	maxStepLength = pmaxStepLength;
}

SecantRootFinding::~SecantRootFinding()
{
}


Vector
SecantRootFinding::findLimitStateSurface(int space, double g, Vector pDirection, Vector thePoint)
{
	// Set scale factor for 'g' for convergence check
	double scaleG;
	if (fabs(g)>1.0e-4) { scaleG = g;}
	else { 		scaleG = 1.0;}

	// Normalize the direction vector
	Vector Direction = pDirection/pDirection.Norm();

	// Scale 'maxStepLength' by standard deviation
	// (only if the user has specified to "walk" in original space)
	double perturbation;
	double realMaxStepLength = maxStepLength;
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	if (space == 1) {

		// Go through direction vector and see which element is biggest
		// compared to its standard deviation
		RandomVariable *theRV;
		double stdv, theStdv;
		int theBiggest;
		double maxRatio = 0.0;
		for (int i=0; i<nrv; i++) {
			theRV = theReliabilityDomain->getRandomVariablePtr(i+1);
			stdv = theRV->getStdv();
			if (Direction(i)/stdv > maxRatio) {
				maxRatio = Direction(i)/stdv;
				theStdv = stdv;
				theBiggest = i+1;
			}
		}

		// Now scale so that 'maxStepSize' is related to the real stdv
		perturbation = maxStepLength * theStdv;
		realMaxStepLength = perturbation;
	}
	else {
		perturbation = maxStepLength;
	}

	Vector theTempPoint(nrv);
	double g_old, g_new;
	bool didNotConverge=true;
	double result;
	double tangent;


	int i=0;
	while (i<=maxIter && didNotConverge) {


		// Increment counter right away...
		i++;

		if (i!=1) {

			// Transform the point into x-space if the user has given it in 2-space
			if (space==2) {
			  /*
				result = theProbabilityTransformation->set_u(thePoint);
				if (result < 0) {
					opserr << "GFunVisualizationAnalysis::analyze() - " << endln
						<< " could not set u in the xu-transformation." << endln;
					return -1;
				}

				result = theProbabilityTransformation->transform_u_to_x();
				if (result < 0) {
					opserr << "GFunVisualizationAnalysis::analyze() - " << endln
						<< " could not transform from u to x and compute Jacobian." << endln;
					return -1;
				}
				theTempPoint = theProbabilityTransformation->get_x();
			  */
			  result = theProbabilityTransformation->transform_u_to_x(thePoint, theTempPoint);
			  if (result < 0) {
			    opserr << "GFunVisualizationAnalysis::analyze() - " << endln
				   << " could not transform from u to x and compute Jacobian." << endln;
			    return -1;
			  }
			}
			else {
				theTempPoint = thePoint;
			}


			// Evaluate limit-state function
			result = theGFunEvaluator->runAnalysis();
			if (result < 0) {
				opserr << "GFunVisualizationAnalysis::analyze() - " << endln
					<< " could not run analysis to evaluate limit-state function. " << endln;
				return -1;
			}
			g_new = theGFunEvaluator->evaluateExpression();

		}
		else {
			g_new = g;
		}

		

		// Check convergence
		if (fabs(g_new/scaleG) < tol) {
			didNotConverge = false;
		}
		else {
			if (i==maxIter) {
				opserr << "WARNING: Projection scheme failed to find surface..." << endln;
			}
			else if (i==1) {
				thePoint = thePoint - perturbation * Direction;
				g_old = g_new;
			}
			else {

				// Take a step
				tangent = (g_new-g_old)/perturbation;
				perturbation = -g_new/tangent;
				if (fabs(perturbation) > realMaxStepLength) {
					perturbation = perturbation/fabs(perturbation)*realMaxStepLength;
				}
				thePoint = thePoint - perturbation * Direction;
				g_old = g_new;
			}
		}
	}

	return thePoint;
}


