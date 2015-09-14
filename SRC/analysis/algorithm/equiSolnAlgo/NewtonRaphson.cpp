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

void* OPS_NewNewtonRaphsonAlgorithm()
{
    int formTangent = CURRENT_TANGENT;

    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type = OPS_GetString();
	if(type=="-secant" || type=="-Secant") {
	    formTangent = CURRENT_SECANT;
	} else if(type=="-initial" || type=="-Initial") {
	    formTangent = INITIAL_TANGENT;
	} else if(type=="-intialThenCurrent" || type=="-intialCurrent") {
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
	  if (theIntegrator->formTangent(INITIAL_TANGENT) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	  } 
	} else {
	  if (theIntegrator->formTangent(CURRENT_TANGENT) < 0){
	    opserr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	  } 
	}	  
      } else {
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

    } while (result == -1);

    if (result == -2) {
      opserr << "NewtnRaphson::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
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
