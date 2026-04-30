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

// $Revision: 1.6 $
// $Date: 2007-05-04 06:59:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/NewtonLineSearch.cpp,v $

// Written: fmk 
// Created: 11/96 
// Modified: Ed "C++" Love 10/00 to perform the line search
//
// Description: This file contains the implementation for NewtonLineSearch. 
// 
// What: "@(#)NewtonLineSearch.h, revA"

#include <NewtonLineSearch.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>

#include <InitialInterpolatedLineSearch.h>
#include <BisectionLineSearch.h>
#include <SecantLineSearch.h>
#include <RegulaFalsiLineSearch.h>

#include <elementAPI.h>
#include <string.h>

// Shared Tcl/Python OPS_NewtonLineSearch; ConvergenceTest attached by interpreter.
void* OPS_NewtonLineSearch()
{
    double tol        = 0.8;
    int    maxIter    = 10;
    double maxEta     = 10.0;
    double minEta     = 0.1;
    int    pFlag      = 1;
    int    typeSearch = 0;

    int numdata = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();

	if (strcmp(flag, "-tol") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read tol\n";
		return 0;
	    }
	} else if (strcmp(flag, "-maxIter") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    if (OPS_GetIntInput(&numdata, &maxIter) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read maxIter\n";
		return 0;
	    }
	} else if (strcmp(flag, "-pFlag") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    if (OPS_GetIntInput(&numdata, &pFlag) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read pFlag\n";
		return 0;
	    }
	} else if (strcmp(flag, "-minEta") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    if (OPS_GetDoubleInput(&numdata, &minEta) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read minEta\n";
		return 0;
	    }
	} else if (strcmp(flag, "-maxEta") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    if (OPS_GetDoubleInput(&numdata, &maxEta) < 0) {
		opserr << "WARNING NewtonLineSearch failed to read maxEta\n";
		return 0;
	    }
	} else if (strcmp(flag, "-type") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    const char* flag2 = OPS_GetString();
	    if (strcmp(flag2, "Bisection") == 0)
		typeSearch = 1;
	    else if (strcmp(flag2, "Secant") == 0)
		typeSearch = 2;
	    else if (strcmp(flag2, "RegulaFalsi") == 0)
		typeSearch = 3;
	    else if (strcmp(flag2, "LinearInterpolated") == 0)
		typeSearch = 3;
	    else if (strcmp(flag2, "InitialInterpolated") == 0)
		typeSearch = 0;
	}
    }

    LineSearch *theLineSearch = 0;
    if (typeSearch == 0)
	theLineSearch = new InitialInterpolatedLineSearch(tol, maxIter, minEta, maxEta, pFlag);
    else if (typeSearch == 1)
	theLineSearch = new BisectionLineSearch(tol, maxIter, minEta, maxEta, pFlag);
    else if (typeSearch == 2)
	theLineSearch = new SecantLineSearch(tol, maxIter, minEta, maxEta, pFlag);
    else if (typeSearch == 3)
	theLineSearch = new RegulaFalsiLineSearch(tol, maxIter, minEta, maxEta, pFlag);

    return new NewtonLineSearch(theLineSearch);
}


//Null Constructor
NewtonLineSearch::NewtonLineSearch( )
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonLineSearch),
 theTest(0), theOtherTest(0), theLineSearch(0)
{   
}


//Constructor 
NewtonLineSearch::NewtonLineSearch( ConvergenceTest &theT, 
				   LineSearch *theSearch) 
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonLineSearch),
 theTest(&theT), theLineSearch(theSearch)
{
  theOtherTest = theTest->getCopy(10);
  theOtherTest->setEquiSolnAlgo(*this);
}


// Constructor with deferred ConvergenceTest assignment.
// The test is wired in later via setConvergenceTest() (e.g. by the
// interpreter-level setAlgorithm() / StaticAnalysis bookkeeping). Keeps
// argument-parsing factories like OPS_NewtonLineSearch() free from any
// dependency on a specific interpreter's "current test" state.
NewtonLineSearch::NewtonLineSearch(LineSearch *theSearch)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonLineSearch),
 theTest(0), theOtherTest(0), theLineSearch(theSearch)
{
}


// Destructor
NewtonLineSearch::~NewtonLineSearch()
{
  if (theOtherTest != 0)
    delete theOtherTest;
  if (theLineSearch != 0)
    delete theLineSearch;
}

int
NewtonLineSearch::setConvergenceTest(ConvergenceTest *newTest)
{
    theTest = newTest;
    if (theOtherTest != 0)
      delete theOtherTest;
    theOtherTest = theTest->getCopy(10);
    theOtherTest->setEquiSolnAlgo(*this);
    return 0;
}


int 
NewtonLineSearch::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
    LinearSOE  *theSOE = this->getLinearSOEptr();

    if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
	|| (theTest == 0)){
	opserr << "WARNING NewtonLineSearch::solveCurrentStep() - setLinks() has";
	opserr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    theLineSearch->newStep(*theSOE);

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "NewtonLineSearch::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    if (theIntegrator->formUnbalance() < 0) {
      opserr << "WARNING NewtonLineSearch::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }	    

    int result = -1;
    do {

	//residual at this iteration before next solve 
	const Vector &Resid0 = theSOE->getB() ;
	
	//form the tangent
        if (theIntegrator->formTangent() < 0){
	    opserr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	}		    
	
	//solve 
	if (theSOE->solve() < 0) {
	    opserr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    opserr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
	}	    


	//line search direction 
	const Vector &dx0 = theSOE->getX() ;

	//initial value of s
	double s0 = - (dx0 ^ Resid0) ; 

	if (theIntegrator->update(theSOE->getX()) < 0) {
	    opserr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    opserr << "the Integrator failed in update()\n";	
	    return -4;
	}	        

	if (theIntegrator->formUnbalance() < 0) {
	    opserr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    opserr << "the Integrator failed in formUnbalance()\n";	
	    return -2;
	}	

	// do a line search only if convergence criteria not met
	theOtherTest->start();
	result = theOtherTest->test();

	if (result < 1) {
	  //new residual 
	  const Vector &Resid = theSOE->getB() ;
	  
	  //new value of s 
	  double s = - ( dx0 ^ Resid ) ;
	  
	  if (theLineSearch != 0)
	    theLineSearch->search(s0, s, *theSOE, *theIntegrator);
	}

	this->record(0);
	  
	result = theTest->test();

    } while (result == -1);

    if (result == -2) {
      opserr << "NewtonLineSearch::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }

    // note - if positive result we are returning what the convergence test returned
    // which should be the number of iterations
    return result;
}

ConvergenceTest *
NewtonLineSearch::getConvergenceTest(void)
{
  return theTest;
}

int
NewtonLineSearch::sendSelf(int cTag, Channel &theChannel)
{
  static ID data(1);
  data(0) = theLineSearch->getClassTag();
  if (theChannel.sendID(0, cTag, data) < 0) {
    opserr << "NewtonLineSearch::sendSelf(int cTag, Channel &theChannel)   - failed to send date\n";
    return -1;
  }

  if (theLineSearch->sendSelf(cTag, theChannel) < 0) {
    opserr << "NewtonLineSearch::sendSelf(int cTag, Channel &theChannel)   - failed to send line search\n";
    return -1;
  }

  return 0;
}

int
NewtonLineSearch::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  static ID data(1);
  if (theChannel.recvID(0, cTag, data) < 0) {
    opserr << "NewtonLineSearch::recvSelf(int cTag, Channel &theChannel) - failed to recv data\n";
    return -1;
  }

  int lineSearchClassTag = data(0);

  if (theLineSearch == 0 || theLineSearch->getClassTag() != lineSearchClassTag) {
    if (theLineSearch != 0)
      delete theLineSearch;

    theLineSearch = theBroker.getLineSearch(lineSearchClassTag);
    if (theLineSearch == 0) {
      opserr << "NewtonLineSearch::recvSelf(int cTag, Channel &theChannel) - failed to obtain a LineSerach object\n";
      return -1;
    }
  }

  if (theLineSearch->recvSelf(cTag, theChannel, theBroker) < 0) {
      opserr << "NewtonLineSearch::recvSelf(int cTag, Channel &theChannel) - failed to recv the LineSerach object\n";
      return -1;
  }

  return 0;

}


void
NewtonLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) 
    s << "NewtonLineSearch\n";

  if (theLineSearch != 0)
    theLineSearch->Print(s, flag);
}









