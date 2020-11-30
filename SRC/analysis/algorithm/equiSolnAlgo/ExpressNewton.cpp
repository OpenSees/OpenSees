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
ExpressNewton::ExpressNewton(int ni, double km, int tg, int fo)
  :EquiSolnAlgo(EquiALGORITHM_TAGS_ExpressNewton), nIter(ni), factorOnce(fo) 
{
  if (tg == INITIAL_TANGENT) {
    kMultiplier1 = km;
    kMultiplier2 = 0.0;
  }
  else {
    kMultiplier1 = 0.0;
    kMultiplier2 = km;
  }
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
		if (theIntegrator->formTangent(HALL_TANGENT, kMultiplier1, kMultiplier2) < 0) {
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
  static Vector data(4);
  data(0) = nIter;
  data(1) = kMultiplier1;
  data(1) = kMultiplier2;
  data(2) = factorOnce;
  return theChannel.sendVector(this->getDbTag(), cTag, data);


}

int
ExpressNewton::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(4);
  theChannel.recvVector(this->getDbTag(), cTag, data);
  nIter = int(data(0));
  kMultiplier1 = data(1);
  kMultiplier2 = data(2);
  factorOnce = int(data(3));
  return 0;
}


void
ExpressNewton::Print(OPS_Stream &s, int flag)
{
    s << "\t ExpressNewton algorithm";
}
