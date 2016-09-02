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
                                                                        
// $Revision: 1.10 $
// $Date: 2007-04-02 23:41:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/NewtonRaphson.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/NewtonRaphson.C 
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
// Revision: A 
//

// Description: This file contains the class definition for 
// NewtonRaphson. NewtonRaphson is a class which uses the
// Newton-Raphson solution algorihm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)NewtonRaphson.C, revA"

#include <NewtonRaphson.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>
#include <elementAPI.h>
#include <string>


void* OPS_NewtonRaphsonAlgorithm()
{
    int formTangent = CURRENT_TANGENT;

    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-secant")==0 || strcmp(type,"-Secant")==0) {
	    formTangent = CURRENT_SECANT;
	} else if(strcmp(type,"-initial")==0 || strcmp(type,"-Initial")==0) {
	    formTangent = INITIAL_TANGENT;
	} else if(strcmp(type,"-intialThenCurrent")==0 || strcmp(type,"-intialCurrent")==0) {
	    formTangent = INITIAL_THEN_CURRENT_TANGENT;
	}
    }

    return new NewtonRaphson(formTangent);

}

// Constructor
NewtonRaphson::NewtonRaphson(int theTangentToUse)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 tangent(theTangentToUse)
{

}


NewtonRaphson::NewtonRaphson(ConvergenceTest &theT, int theTangentToUse)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 tangent(theTangentToUse)
{

}

// Destructor
NewtonRaphson::~NewtonRaphson()
{
  

}


int 
NewtonRaphson::solveCurrentStep(void)
{
  
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
    //IncrementalIntegrator *theIntegratorSens=this->getIncrementalIntegratorPtr();//Abbas
    LinearSOE  *theSOE = this->getLinearSOEptr();

    if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
	|| (theTest == 0)){
	opserr << "WARNING NewtonRaphson::solveCurrentStep() - setLinks() has";
	opserr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    if (theIntegrator->formUnbalance() < 0) {
      opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }	    

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "NewtnRaphson::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    int result = -1;
    numIterations = 0;

    do {

      if (tangent == INITIAL_THEN_CURRENT_TANGENT) {
	if (numIterations == 0) {
	  SOLUTION_ALGORITHM_tangentFlag = INITIAL_TANGENT;
	  if (theIntegrator->formTangent(INITIAL_TANGENT) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	  } 
	} else {
	  SOLUTION_ALGORITHM_tangentFlag = CURRENT_TANGENT;
	  if (theIntegrator->formTangent(CURRENT_TANGENT) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	  } 
	}
      }	else {
	
	  SOLUTION_ALGORITHM_tangentFlag = tangent;
	if (theIntegrator->formTangent(tangent) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	}		    
      } 
           if (theSOE->solve() < 0) {
	opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	opserr << "the LinearSysOfEqn failed in solve()\n";	
	return -3;
      }	    

      if (theIntegrator->update(theSOE->getX()) < 0) {
	opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	opserr << "the Integrator failed in update()\n";	
	return -4;
      }	        
      if (theIntegrator->formUnbalance() < 0) {
	opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	opserr << "the Integrator failed in formUnbalance()\n";	
	return -2;
      }	

     
     


      result = theTest->test();
       numIterations++;
      this->record(numIterations);

      ///////////////////*
//ActivateAensitivity() is to account for the reliability effect.
//SensitivityintegratorScheme returns:  true for DC, and false for LC
 if( ((theIntegrator->activateSensitivity())==true) && (theIntegrator->computeSensitivityAtEachIteration())==true)
    {
       
 theIntegrator->computeSensitivities();
 theIntegrator->formUnbalance();

    } 
     
////////////////////

    } while (result == -1);

    if (result == -2) {
      opserr << "NewtnRaphson::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }
//IntegratorSens=theIntegrator;

    if( ((theIntegrator->activateSensitivity())==true) && (theIntegrator->computeSensitivityAtEachIteration())==false)
    {
      theIntegrator->computeSensitivities();//Abbas
    }    
// note - if postive result we are returning what the convergence test returned
    // which should be the number of iterations
    
        return result;
}


int
NewtonRaphson::sendSelf(int cTag, Channel &theChannel)
{
  static ID data(1);
  data(0) = tangent;
  return theChannel.sendID(this->getDbTag(), cTag, data);
}

int
NewtonRaphson::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  static ID data(1);
  theChannel.recvID(this->getDbTag(), cTag, data);
  tangent = data(0);
  return 0;
}


void
NewtonRaphson::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "NewtonRaphson" << endln;
  }
}


int
NewtonRaphson::getNumIterations(void)
{
  return numIterations;
}


