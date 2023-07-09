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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-04-02 23:41:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/ModifiedNewton.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/ModifiedNewton.C 
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//
// Description: This file contains the class definition for 
// ModifiedNewton. ModifiedNewton is a class which uses the
// Newton-Raphson solution algorithm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)ModifiedNewton.C, revA"

#include <ModifiedNewton.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
//#include <Timer.h>
#include <elementAPI.h>

void* OPS_ModifiedNewton()
{
    int formTangent = CURRENT_TANGENT;
    double iFactor = 0;
    double cFactor = 1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char* type = OPS_GetString();
      if (strcmp(type,"-secant") == 0) {
	formTangent = CURRENT_SECANT;
      } else if (strcmp(type,"-initial") == 0) {
	formTangent = INITIAL_TANGENT;
      } else if(strcmp(type,"-hall")==0 || strcmp(type,"-Hall")==0) {
	formTangent = HALL_TANGENT;
	iFactor = 0.1;
	cFactor = 0.9;
	if (OPS_GetNumRemainingInputArgs() == 2) {
	  double data[2];
	  int numData = 2;
	  if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
	    opserr << "WARNING invalid data reading 2 hall factors\n";
	    return 0;
	  }
	  iFactor = data[0];
	  cFactor = data[1];
	}
      }
    }
    
    return new ModifiedNewton(formTangent, iFactor, cFactor);

}

// Constructor
ModifiedNewton::ModifiedNewton(int theTangentToUse, double iFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_ModifiedNewton),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact)
{
  
}


ModifiedNewton::ModifiedNewton(ConvergenceTest &theT, int theTangentToUse, double iFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_ModifiedNewton),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact)
{

}

// Destructor
ModifiedNewton::~ModifiedNewton()
{

}


int 
ModifiedNewton::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel       *theAnalysisModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIncIntegratorr = this->getIncrementalIntegratorPtr();
    LinearSOE	        *theSOE = this->getLinearSOEptr();

    if ((theAnalysisModel == 0) || (theIncIntegratorr == 0) || (theSOE == 0)
	|| (theTest == 0)){
	opserr << "WARNING ModifiedNewton::solveCurrentStep() - setLinks() has";
	opserr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    // we form the tangent
    //    Timer timer1;
    // timer1.start();

    if (theIncIntegratorr->formUnbalance() < 0) {
	opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	opserr << "the Integrator failed in formUnbalance()\n";	
	return -2;
    }	

    SOLUTION_ALGORITHM_tangentFlag = tangent;
    if (theIncIntegratorr->formTangent(tangent, iFactor, cFactor) < 0){
	opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	opserr << "the Integrator failed in formTangent()\n";
	return -1;
    }		    

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "ModifiedNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    // repeat until convergence is obtained or reach max num iterations
    int result = -1;
    numIterations = 0;
    do {
      //Timer timer2;
      //timer2.start();
	if (theSOE->solve() < 0) {
	    opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    opserr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
	}	    
	//timer2.pause();
	//opserr << "TIMER::SOLVE()- " << timer2;
	
	if (theIncIntegratorr->update(theSOE->getX()) < 0) {
	    opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    opserr << "the Integrator failed in update()\n";	
	    return -4;
	}	        

	if (theIncIntegratorr->formUnbalance() < 0) {
	    opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    opserr << "the Integrator failed in formUnbalance()\n";	
	    return -2;
	}	

	this->record(numIterations++);
	result = theTest->test();


    } while (result == -1);

    //timer1.pause();
    //opserr << "TIMER::solveCurrentStep - " << timer1;

    if (result == -2) {
      opserr << "ModifiedNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }
    return result;
}

int
ModifiedNewton::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(3);
  data(0) = tangent;
  data(1) = iFactor;
  data(2) = cFactor;
  return theChannel.sendVector(this->getDbTag(), cTag, data);
}

int
ModifiedNewton::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  static Vector data(3);
  theChannel.recvVector(this->getDbTag(), cTag, data);
  tangent = data(0);
  iFactor = data(1);
  cFactor = data(2);
  return 0;
}

void
ModifiedNewton::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) {
	s << "ModifiedNewton";
    }
}

int
ModifiedNewton::getNumIterations(void)
{
  return numIterations;
}
