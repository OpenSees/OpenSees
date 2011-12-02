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
// $Date: 2003-10-27 23:45:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/SearchWithStepSizeAndStepDirection.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef SearchWithStepSizeAndStepDirection_h
#define SearchWithStepSizeAndStepDirection_h

#include <FindDesignPointAlgorithm.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <ProbabilityTransformation.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <HessianApproximation.h>
#include <ReliabilityConvergenceCheck.h>
#include <Matrix.h>
#include <Vector.h>
#include <ReliabilityDomain.h>

#include <fstream>
using std::ofstream;

class SearchWithStepSizeAndStepDirection : public FindDesignPointAlgorithm
{

public:

	// Constructor and destructor
	SearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations, 
					GFunEvaluator *passedGFunEvaluator,
					GradGEvaluator *passedGradGEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					ProbabilityTransformation *passedProbabilityTransformation,
					HessianApproximation *theHessianApproximation,
					ReliabilityConvergenceCheck *theReliabilityConvergenceCheck,
					int printFlag,
					char *fileNamePrint,
					Vector *startPoint);
	~SearchWithStepSizeAndStepDirection();
	
	int findDesignPoint(ReliabilityDomain *theReliabilityDomain);

	Vector get_x();
	Vector get_u();
	Vector get_alpha();
	Vector get_gamma();
	int getNumberOfSteps();
	Vector getSecondLast_u();
	Vector getSecondLast_alpha();
	Vector getLastSearchDirection();
	double getFirstGFunValue();
	double getLastGFunValue();
	Vector getGradientInStandardNormalSpace();
	int    getNumberOfEvaluations();

protected:

private:	

	// The reliability domain and tools for the analysis
	ReliabilityDomain *theReliabilityDomain;
	GFunEvaluator *theGFunEvaluator;
	GradGEvaluator *theGradGEvaluator;
	StepSizeRule *theStepSizeRule;
	SearchDirection *theSearchDirection;
	ProbabilityTransformation *theProbabilityTransformation;
	HessianApproximation *theHessianApproximation;
	ReliabilityConvergenceCheck *theReliabilityConvergenceCheck;

	// Private member functions to do the job
	int doTheActualSearch(bool doRvProjection);
	int doRvProjection(Vector uOld, Vector uNew);

	// Data members set when the object is created
	int maxNumberOfIterations;

	// Data members where the results are to be stored
	Vector x;
	Vector u;
	Vector alpha;
	Vector gradientInStandardNormalSpace;
	Vector gamma;
	Vector uSecondLast;
	Vector alphaSecondLast;
	int i;
	Vector searchDirection;
	double Gfirst;
	double Glast;

	// Data members set through the call when a job is to be done
	Vector *startPoint;
	Vector *designPoint_uStar;

	int printFlag;
	char *fileNamePrint;
	int numberOfEvaluations;
};

#endif
