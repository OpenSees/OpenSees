/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/ReliabilityDirectIntegrationAnalysis.h,v $
                                                                    
#ifndef ReliabilityDirectIntegrationAnalysis_h
#define ReliabilityDirectIntegrationAnalysis_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for 
// DirectIntegrationAnalysis. DirectIntegrationAnalysis is a 
// subclass of TransientAnalysis. It is used to perform a 
// dynamic analysis on the FE\_Model using a direct integration scheme.
//
// What: "@(#) DirectIntegrationAnalysis.h, revA"

#include <TransientAnalysis.h>
// AddingSensitivity:BEGIN //////////////////////////////////
#include <SensitivityAlgorithm.h>
#include <GFunEachStepEvaluator.h>
#include<Integrator.h>
// AddingSensitivity:END ////////////////////////////////////

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class TransientIntegrator;
class LinearSOE;
class EquiSolnAlgo;
class ConvergenceTest;
class Integrator;//Abbas

class ReliabilityDirectIntegrationAnalysis: public TransientAnalysis
{
  public:
    ReliabilityDirectIntegrationAnalysis(Domain &theDomain, 
			      ConstraintHandler &theHandler,
			      DOF_Numberer &theNumberer,
			      AnalysisModel &theModel,
			      EquiSolnAlgo &theSolnAlgo,		   
			      LinearSOE &theSOE,
			      TransientIntegrator &theIntegrator,
			      ConvergenceTest *theTest = 0);

    virtual ~ReliabilityDirectIntegrationAnalysis();

    void clearAll(void);	    
    
    int analyze(int numSteps, double dT);
    int initialize(void);

    int domainChanged(void);

    int setNumberer(DOF_Numberer &theNumberer);    
    int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    int setIntegrator(TransientIntegrator &theIntegrator);
    int setLinearSOE(LinearSOE &theSOE); 
    int setConvergenceTest(ConvergenceTest &theTest);
    
    int checkDomainChange(void);
    EquiSolnAlgo *getAlgorithm(void);
    TransientIntegrator *getIntegrator(void);
    ConvergenceTest *getConvergenceTest(void);

    // AddingSensitivity:BEGIN ///////////////////////////////
    int setSensitivityAlgorithm(/*SensitivityAlgorithm*/ Integrator  *theSensitivityAlgorithm);
	///// added by K Fujimura /////
	bool reliabilityAnalysis(void){ return true;}
    // AddingSensitivity:END /////////////////////////////////
    
	//////////////////////////
	// added by K Fujimura //
	//////////////////////////
	//SensitivityAlgorithm* getSensitivityAlgorithm()
	Integrator* getSensitivityAlgorithm()//Abbas
	{return  theSensitivityAlgorithm;};
	void setGFunEachStepEvaluator(GFunEachStepEvaluator *pGFunEachStepEvaluator)
	{ theGFunEachStepEvaluator=pGFunEachStepEvaluator;}
	void inactivateGFunEachStepEvaluator()
	{ theGFunEachStepEvaluator=0;}
	Matrix* getEachStepResult();
	Matrix* getEachStepConvFlag();

	//////////////////////////
	// added by K Fujimura //
	//////////////////////////
  protected:
    
  private:
    ConstraintHandler 	*theConstraintHandler;    
    DOF_Numberer 	*theDOF_Numberer;
    AnalysisModel 	*theAnalysisModel;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    TransientIntegrator *theIntegrator;
    ConvergenceTest     *theTest;

    int domainStamp;

    // AddingSensitivity:BEGIN ///////////////////////////////
   // SensitivityAlgorithm *theSensitivityAlgorithm;
    Integrator *theSensitivityAlgorithm;

    GFunEachStepEvaluator* theGFunEachStepEvaluator;
    // AddingSensitivity:END ///////////////////////////////

};

#endif

