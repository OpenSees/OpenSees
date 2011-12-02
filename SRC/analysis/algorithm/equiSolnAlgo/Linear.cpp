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
                                                                        
// $Revision: 1.3 $
// $Date: 2005-11-29 22:42:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/Linear.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/Linear.C 
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
// Revision: A 
//

// Description: This file contains the class definition for 
// Linear. Linear is a class which performs a linear solution algorihm
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

#include <Timer.h>
// Constructor
Linear::Linear()
:EquiSolnAlgo(EquiALGORITHM_TAGS_Linear)
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

    if (theIncIntegrator->formTangent() < 0) {
	opserr << "WARNING Linear::solveCurrentStep() -";
	opserr << "the Integrator failed in formTangent()\n";
	return -1;
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
    return 0;
}

int
Linear::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}


void
Linear::Print(OPS_Stream &s, int flag)
{
    s << "\t Linear algorithm";
}
