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
// $Date: 2003-10-27 23:45:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/ArmijoStepSizeRule.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef ArmijoStepSizeRule_h
#define ArmijoStepSizeRule_h

#include <StepSizeRule.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <ProbabilityTransformation.h>
#include <MeritFunctionCheck.h>
#include <RootFinding.h>

class ArmijoStepSizeRule : public StepSizeRule
{

public:
	ArmijoStepSizeRule(GFunEvaluator *passedGFunEvaluator, 
		ProbabilityTransformation *theProbabilityTransformation,
		MeritFunctionCheck *theMeritFunctionCheck,
		RootFinding *theRootFindingAlgorithm, 
		double base,
		int maxNumReductions,
		double b0,
		int numberOfShortSteps,
		double radius,
		double surfaceDistance,
		double evolution,
		int printFlag);
		
	~ArmijoStepSizeRule();

	int		computeStepSize(Vector u, Vector grad_G, double G, Vector d, int stepNumber);
	double	getStepSize();
	double	getInitialStepSize();
	double  getGFunValue();

protected:

private:
	double stepSize;
	GFunEvaluator *theGFunEvaluator;
	GradGEvaluator *theGradGEvaluator;
	ProbabilityTransformation *theProbabilityTransformation;
	MeritFunctionCheck *theMeritFunctionCheck;
	RootFinding *theRootFindingAlgorithm;
	double gFunValue;
	double base;
	int maxNumReductions;
	double b0;
	int numberOfShortSteps;
	double radius;
	double surfaceDistance;
	double evolution;
	bool isCloseToSphere;
	int printFlag;

};

#endif
