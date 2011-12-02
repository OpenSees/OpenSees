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
// $Date: 2000-12-13 04:44:25 $
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

// Constructor
VariableTimeStepDirectIntegrationAnalysis::VariableTimeStepDirectIntegrationAnalysis(
			      Domain &the_Domain,
			      ConstraintHandler &theHandler,
			      DOF_Numberer &theNumberer,
			      AnalysisModel &theModel,
			      EquiSolnAlgo &theSolnAlgo,		   
			      LinearSOE &theLinSOE,
			      TransientIntegrator &theTransientIntegrator)

:DirectIntegrationAnalysis(the_Domain, theHandler, theNumberer, theModel, 
			   theSolnAlgo, theLinSOE, theTransientIntegrator)
{

}    

VariableTimeStepDirectIntegrationAnalysis::~VariableTimeStepDirectIntegrationAnalysis()
{

}    


int 
VariableTimeStepDirectIntegrationAnalysis::analyze(int numSteps, double dT, double dtMin, double dtMax, double perCent)
{
  // get some pointers
  Domain *theDom = this->getDomainPtr();
  EquiSolnAlgo *theAlgo = this->getAlgorithm();
  TransientIntegrator *theIntegratr = this->getIntegrator();
  ConvergenceTest *theTest = theAlgo->getTest();

  // set some variables
  int result = 0;  
  double totalTimeIncr = numSteps * dT;
  double currentTimeIncr = 0.0;
  double currentDt = dT;
  
  // loop until analysis has performed the total time incr requested
  while (currentTimeIncr < totalTimeIncr) {

    // check if domain has undergone change
    if (this->checkDomainChange() < 0) {
      cerr << "VariableTimeStepDirectIntegrationAnalysis::analyze() - failed";
      cerr << " failed at time " << theDom->getCurrentTime() << endl;
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
	cerr << "VariableTimeStepDirectIntegrationAnalysis::analyze() - ";
	cerr << " failed at time " << theDom->getCurrentTime() << endl;
	return result;
      }
      
      // if still here reset result for next loop
      result = 0;
    }

    // now we determine a new delta T for next loop
    currentDt = this->determineDt(currentDt, dtMin, dtMax, perCent, theTest);
  }

  return 0;
}





double 
VariableTimeStepDirectIntegrationAnalysis::determineDt(double dT, 
						       double dtMin, 
						       double dtMax, 
						       double perCent,
						       ConvergenceTest *theTest)
{
  double newDt = dT;
    
  // get the number of trial steps in the last solveCurrentStep()
  double algoPercent = 1.0;
  if (theTest != 0)
    algoPercent = theTest->getRatioNumToMax();
  
  
  // determine new dT based on last dT and Jd and #iter of last step
  double factor = perCent/algoPercent;
  newDt *= factor;
  
  // ensure: dtMin <~~ dT <= dtMax
  if (newDt < dtMin)
    newDt = dtMin - DBL_EPSILON;  // to ensure we get out of the analysis 
                               // loop if can't converge on next step
  else if (newDt > dtMax)
    newDt = dtMax;
    
  return newDt;
}


