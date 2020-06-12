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
                                                                        
// $Revision: 1.0 $
// $Date: 2019-01-28 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/ExpressNewton.cpp,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
//
// Description: This file contains the implementation for the ExpressNewton class.
// ExpressNewton is a class which performs a ExpressNewton solution algorithm
// to solve the equations.
//
// Reference:
// Junjie Xu, Yuli Huang, Zhe Qu,
// An efficient and unconditionally stable numerical algorithm for nonlinear structural dynamics
// International Journal for Numerical Methods in Engineering,
// https://doi.org/10.1002/nme.6456.
// (https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6456)

// What: "@(#)ExpressNewton.h, revA"

#include <ExpressNewton.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <StaticIntegrator.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>

#include <Timer.h>
#include <elementAPI.h>
#include <string>

// Constructor
ExpressNewton::ExpressNewton(int ni, int tg, double km)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_ExpressNewton), tangent(tg), nIter(ni), kMultiplier(km)
{
  factorOnce = 0;
  if (tg == HALL_TANGENT) factorOnce = 1;
}

// Destructor
ExpressNewton::~ExpressNewton()
{

}

int 
ExpressNewton::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass

    AnalysisModel *theAnalysisModel = this->getAnalysisModelPtr(); 
    LinearSOE  *theSOE = this->getLinearSOEptr();
    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();

    if ((theAnalysisModel == 0) || (theIntegrator ==0 ) || (theSOE == 0)){
	opserr << "WARNING ExpressNewton::solveCurrentStep() -";
	opserr << "setLinks() has not been called.\n";
	return -5;
    }

	if (factorOnce != 2) {
		if (theIntegrator->formTangent(tangent, kMultiplier, 0.0) < 0) {
		  opserr << "WARNING ExpressNewton::solveCurrentStep() -";
		  opserr << "the Integrator failed in formTangent()\n";
		  return -1;
		}
		if (factorOnce == 1)
			factorOnce = 2;
    }

    for (int iter = 0; iter <nIter; ++iter)
    {
    if (theIntegrator->formUnbalance() < 0) {
	opserr << "WARNING ExpressNewton::solveCurrentStep() -";
	opserr << "the Integrator failed in formUnbalance()\n";	
	return -2;
    }

    if (theSOE->solve() < 0) {
	opserr << "WARNING ExpressNewton::solveCurrentStep() -";
	opserr << "the LinearSOE failed in solve()\n";	
	return -3;
    }

    if (theIntegrator->update(theSOE->getX()) < 0) {
	opserr << "WARNING ExpressNewton::solveCurrentStep() -";
	opserr << "the Integrator failed in update()\n";	
	return -4;
    }
    }

    return 0;
}

int
ExpressNewton::setConvergenceTest(ConvergenceTest *theNewTest)
{
  return 0;
}

int
ExpressNewton::sendSelf(int cTag, Channel &theChannel)
{
  static ID iData(2);
  iData(0) = tangent;
  iData(1) = factorOnce;
  return theChannel.sendID(cTag, 0, iData);
}

int
ExpressNewton::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID iData(2);
  theChannel.recvID(cTag, 0, iData);
  tangent = iData(0);
  factorOnce = iData(1);
  
  return 0;
}


void
ExpressNewton::Print(OPS_Stream &s, int flag)
{
    s << "\t ExpressNewton algorithm";
}
