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
                                                                        
// $Revision: 1.9 $
// $Date: 2008-10-22 16:41:39 $
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
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <ReliabilityConvergenceCheck.h>
#include <Matrix.h>
#include <Vector.h>
#include <Domain.h>
#include <ReliabilityDomain.h>

#include <fstream>
using std::ofstream;

class SearchWithStepSizeAndStepDirection : public FindDesignPointAlgorithm
{

public:

	// Constructor and destructor
	SearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations, 
					ReliabilityDomain *passedReliabilityDomain,
                    Domain *passedOpenSeesDomain,
					FunctionEvaluator *passedGFunEvaluator,
					GradientEvaluator *passedGradGEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					ProbabilityTransformation *passedProbabilityTransformation,
					ReliabilityConvergenceCheck *theReliabilityConvergenceCheck,
					int printFlag,
					char *fileNamePrint);
	~SearchWithStepSizeAndStepDirection();
    
    int gradientStandardNormal(double gFunctionValue);
	int findDesignPoint();

	const Vector &get_x();
	const Vector &get_u();
	const Vector &get_alpha();
	const Vector &get_gamma();
	int getNumberOfSteps();
    double getFirstCurvature();
	const Vector &getLastSearchDirection();
	double getFirstGFunValue();
	double getLastGFunValue();
	const Vector &getGradientInStandardNormalSpace();
	const Vector &getGradientInOriginalSpace();
	int getNumberOfEvaluations();

protected:

private:	

	// The reliability domain and tools for the analysis
	ReliabilityDomain *theReliabilityDomain;
    Domain *theOpenSeesDomain;
	FunctionEvaluator *theFunctionEvaluator;
	GradientEvaluator *theGradientEvaluator;
	StepSizeRule *theStepSizeRule;
	SearchDirection *theSearchDirection;
	ProbabilityTransformation *theProbabilityTransformation;
	ReliabilityConvergenceCheck *theReliabilityConvergenceCheck;

	// Data members set when the object is created
	int maxNumberOfIterations;

	// Data members where the results are to be stored
	Vector *x;
	Vector *u;
	Vector *alpha;
	Vector *gamma;
	Vector *gradientInStandardNormalSpace;
    Vector *gradientInOriginalSpace;
	Vector *uSecondLast;
	Vector *alphaSecondLast;
	Vector *searchDirection;
    Matrix *Jux;
    Matrix *Jxu;
	
	int steps;
	double Gfirst;
	double Glast;
	int printFlag;
	char fileNamePrint[256];
	int numberOfEvaluations;
    double firstCurvature;
};

#endif
