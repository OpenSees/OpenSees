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

// Constructor
NewtonRaphson::NewtonRaphson()
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 theTest(0)
{

}


NewtonRaphson::NewtonRaphson(ConvergenceTest &theT)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 theTest(&theT)
{

}

// Destructor
NewtonRaphson::~NewtonRaphson()
{

}

void 
NewtonRaphson::setTest(ConvergenceTest &newTest)
{
    theTest = &newTest;
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
	cerr << "WARNING NewtonRaphson::solveCurrentStep() - setLinks() has";
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

    if (theIntegrator->formUnbalance() < 0) {
      cerr << "WARNING NewtonRaphson::solveCurrentStep() -";
      cerr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }	    

    int result = -1;
    do {
	if (theIntegrator->formTangent() < 0){
	    cerr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    cerr << "the Integrator failed in formTangent()\n";
	    return -1;
	}		    
	
	if (theSOE->solve() < 0) {
	    cerr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    cerr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
	}	    

	if (theIntegrator->update(theSOE->getX()) < 0) {
	    cerr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    cerr << "the Integrator failed in update()\n";	
	    return -4;
	}	        

	if (theIntegrator->formUnbalance() < 0) {
	    cerr << "WARNING NewtonRaphson::solveCurrentStep() -";
	    cerr << "the Integrator failed in formUnbalance()\n";	
	    return -2;
	}	

//	cerr << theSOE->getX();
//	cerr << theSOE->getB();
	

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
NewtonRaphson::sendSelf(int cTag, Channel &theChannel)
{
  int result = 0;
  int dataTag = this->getDbTag();
  ID data(2);
  data(0) = theTest->getClassTag();
  data(1) = theTest->getDbTag();
  result = theChannel.sendID(dataTag, cTag, data);
  if (result != 0) {
    cerr << "NewtonRaphson::sendSelf() - failed to send ID\n";
    return result;
  }

  result = theTest->sendSelf(cTag, theChannel);
  if (result != 0) {
    cerr << "NewtonRaphson::sendSelf() - failed to send CTest object\n";
    return result;
  }
  
  return 0;
}

int
NewtonRaphson::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    ID data(2);
    int result;
    int dataTag = this->getDbTag();

    result = theChannel.recvID(dataTag, cTag, data);    
    if (result != 0) {
      cerr << "NewtonRaphson::recvSelf() - failed to receive ID\n";
      return result;
    }
    int ctType = data(0);
    int ctDb = data(1);
    
    theTest = theBroker.getNewConvergenceTest(ctType);
    theTest->setDbTag(ctDb);
    result = theTest->recvSelf(cTag, theChannel, theBroker);
    if (result != 0) {
      cerr << "NewtonRaphson::recvSelf() - failed to recv CTest object\n";
      return result;
    }
    
    return 0;
}


void
NewtonRaphson::Print(ostream &s, int flag)
{
    if (flag == 0) {
	s << "ModifiedNewton";
    }
}









