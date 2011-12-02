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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-05-11 21:32:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/VariableTimeStepDirectIntegrationAnalysis.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/VariableTimeStepDirectIntegrationAnalysis.C
// 
// Written: fmk 
// Created: 10/00
// Revision: A
//
// Description: This file contains the implementation of the
// VariableTimeStepDirectIntegrationAnalysis class.
//
// What: "@(#) VariableTimeStepDirectIntegrationAnalysis.C, revA"

#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include <EquiSolnAlgo.h>
#include <TransientIntegrator.h>
#include <Domain.h>
#include <ConvergenceTest.h>
#include <float.h>
#include <AnalysisModel.h>

// Constructor
VariableTimeStepDirectIntegrationAnalysis::VariableTimeStepDirectIntegrationAnalysis(
			      Domain &the_Domain,
			      ConstraintHandler &theHandler,
			      DOF_Numberer &theNumberer,
			      AnalysisModel &theModel,
			      EquiSolnAlgo &theSolnAlgo,		   
			      LinearSOE &theLinSOE,
			      TransientIntegrator &theTransientIntegrator,
			      ConvergenceTest *theTest)

:DirectIntegrationAnalysis(the_Domain, theHandler, theNumberer, theModel, 
			   theSolnAlgo, theLinSOE, theTransientIntegrator, theTest)
{

}    

VariableTimeStepDirectIntegrationAnalysis::~VariableTimeStepDirectIntegrationAnalysis()
{

}    

int 
VariableTimeStepDirectIntegrationAnalysis::analyze(int numSteps, double dT, double dtMin, double dtMax, int Jd)
{
  // get some pointers
  Domain *theDom = this->getDomainPtr();
  EquiSolnAlgo *theAlgo = this->getAlgorithm();
  TransientIntegrator *theIntegratr = this->getIntegrator();
  ConvergenceTest *theTest = theAlgo->getConvergenceTest();
  AnalysisModel *theModel = this->getModel();


  // set some variables
  int result = 0;  
  double totalTimeIncr = numSteps * dT;
  double currentTimeIncr = 0.0;
  double currentDt = dT;

  // loop until analysis has performed the total time incr requested
  while (currentTimeIncr < totalTimeIncr) {

    if (theModel->analysisStep(currentDt) < 0) {
      opserr << "DirectIntegrationAnalysis::analyze() - the AnalysisModel failed in newStepDomain";
      opserr << " at time " << theDom->getCurrentTime() << endln;
      theDom->revertToLastCommit();
      return -2;
    }

    if (this->checkDomainChange() != 0) {
      opserr << "VariableTimeStepDirectIntegrationAnalysis::analyze() - failed checkDomainChange\n";
      return -1;
    }

    //
    // do newStep(), solveCurrentStep() and commit() as in regular
    // DirectINtegrationAnalysis - difference is we do not return
    // if a failure - we stop the analysis & resize time step if failure
    //

    if (theIntegratr->newStep(currentDt) < 0) {
      result = -2;
    }


    if (result >= 0) {
      result = theAlgo->solveCurrentStep();
      if (result < 0) 
	result = -3;
    }    

    if (result >= 0) {
      result = theIntegratr->commit();
      if (result < 0) 
	result = -4;
    }

    // if the time step was successfull increment delta T for the analysis
    // othewise revert the Domain to last committed state & see if can go on

    if (result >= 0) 
      currentTimeIncr += currentDt;
    else {

      // invoke the revertToLastCommit
      theDom->revertToLastCommit();	    
      theIntegratr->revertToLastStep();

      // if last dT was <= min specified the analysis FAILS - return FAILURE
      if (currentDt <= dtMin) {
	opserr << "VariableTimeStepDirectIntegrationAnalysis::analyze() - ";
	opserr << " failed at time " << theDom->getCurrentTime() << endln;
	return result;
      }

      
      // if still here reset result for next loop
      result = 0;
    }

    // now we determine a new delta T for next loop
    currentDt = this->determineDt(currentDt, dtMin, dtMax, Jd, theTest);
  }


  return 0;
}





double 
VariableTimeStepDirectIntegrationAnalysis::determineDt(double dT, 
						       double dtMin, 
						       double dtMax, 
						       int Jd,
						       ConvergenceTest *theTest)
{
  double newDt = dT;
    
  // get the number of trial steps in the last solveCurrentStep()
  double numLastIter = 1.0;
  if (theTest != 0)
    numLastIter = theTest->getNumTests();
  
  
  // determine new dT based on last dT and Jd and #iter of last step
  double factor = Jd/numLastIter;
  newDt *= factor;
  
  // ensure: dtMin <~~ dT <= dtMax
  if (newDt < dtMin)
    newDt = dtMin - DBL_EPSILON;  // to ensure we get out of the analysis 
                               // loop if can't converge on next step
  else if (newDt > dtMax)
    newDt = dtMax;
    
  return newDt;
}


