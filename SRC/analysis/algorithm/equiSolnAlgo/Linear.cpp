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
                                                                        
// $Revision: 1.4 $
// $Date: 2006-09-05 23:02:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/Linear.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/Linear.C 
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
// Revision: A 
//

// Description: This file contains the class definition for 
// Linear. Linear is a class which performs a linear solution algorithm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)Linear.C, revA"

#include <Linear.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <StaticIntegrator.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>

//#include <Timer.h>
#include <elementAPI.h>
#include <string>

void* OPS_LinearAlgorithm()
{
    int formTangent = CURRENT_TANGENT;
    int factorOnce = 0;

    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type = OPS_GetString();
	if(type=="-secant" || type=="-Secant") {
	    formTangent = CURRENT_SECANT;
	} else if(type=="-initial" || type=="-Initial") {
	    formTangent = INITIAL_TANGENT;
	} else if(type=="-factorOnce" || type=="-FactorOnce") {
	    factorOnce = 1;
	}
    }

    return new Linear(formTangent, factorOnce);

}

// Constructor
Linear::Linear(int theTangent, int Fact)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_Linear), incrTangent(theTangent), factorOnce(Fact)
{

}

// Destructor
Linear::~Linear()
{

}

// int run(void)
//    Performs the linear solution algorithm.

int 
Linear::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass

    AnalysisModel *theAnalysisModel = this->getAnalysisModelPtr(); 
    LinearSOE  *theSOE = this->getLinearSOEptr();
    IncrementalIntegrator  *theIncIntegrator = 
	this->getIncrementalIntegratorPtr(); 

    if ((theAnalysisModel == 0) || (theIncIntegrator ==0 ) || (theSOE == 0)){
	opserr << "WARNING Linear::solveCurrentStep() -";
	opserr << "setLinks() has not been called.\n";
	return -5;
    }

	if (factorOnce != 2) {
		if (theIncIntegrator->formTangent(incrTangent) < 0) {
		  opserr << "WARNING Linear::solveCurrentStep() -";
		  opserr << "the Integrator failed in formTangent()\n";
		  return -1;
		}
		if (factorOnce == 1)
			factorOnce = 2;
    }

    
    if (theIncIntegrator->formUnbalance() < 0) {
	opserr << "WARNING Linear::solveCurrentStep() -";
	opserr << "the Integrator failed in formUnbalance()\n";	
	return -2;
    }

    if (theSOE->solve() < 0) {
	opserr << "WARNING Linear::solveCurrentStep() -";
	opserr << "the LinearSOE failed in solve()\n";	
	return -3;
    }

    const Vector &deltaU = theSOE->getX();

    if (theIncIntegrator->update(deltaU) < 0) {
	opserr << "WARNING Linear::solveCurrentStep() -";
	opserr << "the Integrator failed in update()\n";	
	return -4;
    }

    return 0;
}

int
Linear::setConvergenceTest(ConvergenceTest *theNewTest)
{
  return 0;
}

int
Linear::sendSelf(int cTag, Channel &theChannel)
{
  static ID iData(2);
  iData(0) = incrTangent;
  iData(1) = factorOnce;
  return theChannel.sendID(cTag, 0, iData);
}

int
Linear::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID iData(2);
  theChannel.recvID(cTag, 0, iData);
  incrTangent = iData(0);
  factorOnce = iData(1);
  
  return 0;
}


void
Linear::Print(OPS_Stream &s, int flag)
{
    s << "\t Linear algorithm";
}
