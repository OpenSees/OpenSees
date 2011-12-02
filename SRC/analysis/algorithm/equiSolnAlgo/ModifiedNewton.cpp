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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/ModifiedNewton.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/ModifiedNewton.C 
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//
// Description: This file contains the class definition for 
// ModifiedNewton. ModifiedNewton is a class which uses the
// Newton-Raphson solution algorihm
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

// Constructor
ModifiedNewton::ModifiedNewton()
:EquiSolnAlgo(EquiALGORITHM_TAGS_ModifiedNewton),
 theTest(0)
{

}


ModifiedNewton::ModifiedNewton(ConvergenceTest &theT)
:EquiSolnAlgo(EquiALGORITHM_TAGS_ModifiedNewton),
 theTest(&theT)
{

}

// Destructor
ModifiedNewton::~ModifiedNewton()
{

}

void 
ModifiedNewton::setTest(ConvergenceTest &newTest)
{
    theTest = &newTest;
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
	cerr << "WARNING ModifiedNewton::solveCurrentStep() - setLinks() has";
	cerr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      cerr << "NewtnRaphson::solveCurrentStep() -";
      cerr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }
    
    // we form the tangent
    
    if (theIncIntegratorr->formTangent() < 0){
	cerr << "WARNING ModifiedNewton::solveCurrentStep() -";
	cerr << "the Integrator failed in formTangent()\n";
	return -1;
    }		    

    if (theIncIntegratorr->formUnbalance() < 0) {
	cerr << "WARNING ModifiedNewton::solveCurrentStep() -";
	cerr << "the Integrator failed in formUnbalance()\n";	
	return -2;
    }	

    // repeat until convergence is obtained or reach max num iterations
    int result = -1;
    do {
	if (theSOE->solve() < 0) {
	    cerr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    cerr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
	}	    

	if (theIncIntegratorr->update(theSOE->getX()) < 0) {
	    cerr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    cerr << "the Integrator failed in update()\n";	
	    return -4;
	}	        

	if (theIncIntegratorr->formUnbalance() < 0) {
	    cerr << "WARNING ModifiedNewton::solveCurrentStep() -";
	    cerr << "the Integrator failed in formUnbalance()\n";	
	    return -2;
	}	

	
	this->record(0);

	result = theTest->test();

    } while (result == -1);

    if (result == -2) {
      cerr << "NewtnRaphson::solveCurrentStep() -";
      cerr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }

    return result;
}

int
ModifiedNewton::sendSelf(int cTag, Channel &theChannel)
{
  int result = 0;
  int dataTag = this->getDbTag();
  ID data(2);
  data(0) = theTest->getClassTag();
  data(1) = theTest->getDbTag();
  result = theChannel.sendID(dataTag, cTag, data);
  if (result != 0) {
    cerr << "ModifiedNewton::sendSelf() - failed to send ID\n";
    return result;
  }

  result = theTest->sendSelf(cTag, theChannel);
  if (result != 0) {
    cerr << "ModifiedNewton::sendSelf() - failed to send CTest object\n";
    return result;
  }
  
  return 0;
}

int
ModifiedNewton::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    ID data(2);
    int result;
    int dataTag = this->getDbTag();

    result = theChannel.recvID(dataTag, cTag, data);    
    if (result != 0) {
      cerr << "ModifiedNewton::recvSelf() - failed to receive ID\n";
      return result;
    }
    int ctType = data(0);
    int ctDb = data(1);
    
    theTest = theBroker.getNewConvergenceTest(ctType);
    theTest->setDbTag(ctDb);
    result = theTest->recvSelf(cTag, theChannel, theBroker);
    if (result != 0) {
      cerr << "ModifiedNewton::recvSelf() - failed to recv CTest object\n";
      return result;
    }
    
    return 0;
}

void
ModifiedNewton::Print(ostream &s, int flag)
{
    if (flag == 0) {
	s << "ModifiedNewton";
    }
}
