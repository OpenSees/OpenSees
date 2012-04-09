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
// $Date: 2008-02-29 19:47:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/ArmijoStepSizeRule.h,v $


//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#ifndef ArmijoStepSizeRule_h
#define ArmijoStepSizeRule_h

#include <StepSizeRule.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <ProbabilityTransformation.h>
#include <MeritFunctionCheck.h>
#include <RootFinding.h>

class ArmijoStepSizeRule : public StepSizeRule
{

public:
	ArmijoStepSizeRule(ReliabilityDomain *passedRelDomain, 
                       FunctionEvaluator *passedGFunEvaluator, 
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
    int initialize(void);
    int getNumReductions();
	int computeStepSize(const Vector &u, const Vector &grad_G,
			    double G, const Vector &d, int stepNumber, int reschk = 0);
	double	getStepSize();
	double	getInitialStepSize();

protected:

private:
    ReliabilityDomain *theReliabilityDomain;
	FunctionEvaluator *theGFunEvaluator;
	GradientEvaluator *theGradGEvaluator;
	ProbabilityTransformation *theProbabilityTransformation;
	MeritFunctionCheck *theMeritFunctionCheck;
	RootFinding *theRootFindingAlgorithm;

   	double stepSize;
	double base;
	int maxNumReductions;
	double b0;
	int numberOfShortSteps;
	double radius;
	double surfaceDistance;
	double evolution;
	int printFlag;
	int numReduction;
    
    bool isCloseToSphere;
    bool isOutsideSphere;
	bool isSecondTime;
	bool FEconvergence;

};

#endif
