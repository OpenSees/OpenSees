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

// $Revision: 1.2 $
// $Date: 2000-12-14 08:39:40 $
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

#define EquiALGORITHM_TAGS_NewtonLineSearch 2011 

//Null Constructor
NewtonLineSearch::NewtonLineSearch( )
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonLineSearch),
 theTest(0), tolerance(0.8), r0(0), x0(0), x(0), xOld(0)
{   
}


//Constructor 
NewtonLineSearch::NewtonLineSearch( ConvergenceTest &theT, 
				    double LineSearchTolerance = 0.8 ) 
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonLineSearch),
 theTest(&theT), tolerance(LineSearchTolerance), r0(0), x0(0), x(0), xOld(0)
{

  tolerance = fabs( LineSearchTolerance ) ;    
  
  if ( tolerance < 0.5 ) 
      tolerance = 0.5 ;

  if ( tolerance > 0.8 ) 
      tolerance = 0.8 ;
}


// Destructor
NewtonLineSearch::~NewtonLineSearch()
{
  if (r0 != 0)
    delete r0 ;

  if (x0 != 0)
    delete x0 ;

  if (x != 0)
    delete x ;

  if (xOld != 0)
    delete xOld ;
}

void 
NewtonLineSearch::setTest(ConvergenceTest &newTest)
{
    theTest = &newTest;
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
	cerr << "WARNING NewtonLineSearch::solveCurrentStep() - setLinks() has";
	cerr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    // check that the vectors we have are of correct size
    int systemSize = ( theSOE->getB() ).Size();
    if ((r0 == 0) || (r0->Size() != systemSize)) {

      if (r0 != 0) {
	delete r0 ;
	delete x ;
	delete x0 ;
	delete xOld ;
      } // end if 

      r0 = new Vector(systemSize) ;
      x = new Vector(systemSize) ;
      x0 = new Vector(systemSize) ;
      xOld = new Vector(systemSize) ;
      
      if (r0 == 0 || x == 0 || x0 == 0 || xOld == 0 ||
	  r0->Size() != systemSize || x->Size() != systemSize ||
	  x0->Size() != systemSize || xOld->Size() != systemSize) {
	cerr << "WARNING NewtonLineSearch::solveCurrentStep() - out of memory";
	cerr << " creating vectors of size " << systemSize << endl;
	if (r0 != 0)
	  delete r0 ;
	if (x0 != 0)
	  delete x0 ;
	if (x != 0)
	  delete x ;
	if (xOld != 0)
	  delete xOld ;	
	return -5;
      } // end if

  } // end if 


    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      cerr << "NewtonLineSearch::solveCurrentStep() -";
      cerr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    if (theIntegrator->formUnbalance() < 0) {
      cerr << "WARNING NewtonLineSearch::solveCurrentStep() -";
      cerr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }	    

    // create some references to the vector pointers
    Vector &Resid0 = *r0 ;
    Vector &dx0 = *x0 ;
    Vector &dx = *x ;
    Vector &dxOld = *xOld ;

    int result = -1;
    do {

	//residual at this iteration before next solve 
	Resid0 = theSOE->getB() ;

	
	//form the tangent
        if (theIntegrator->formTangent() < 0){
	    cerr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    cerr << "the Integrator failed in formTangent()\n";
	    return -1;
	}		    
	
	//solve 
	if (theSOE->solve() < 0) {
	    cerr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    cerr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
	}	    


	//line search direction 
	dx0 = theSOE->getX() ;

	//intial value of s
	double s0 = - (dx0 ^ Resid0) ; //what is this bullshit 
	                               //inner product notation ?
	
	
	if (theIntegrator->update(theSOE->getX()) < 0) {
	    cerr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    cerr << "the Integrator failed in update()\n";	
	    return -4;
	}	        

	if (theIntegrator->formUnbalance() < 0) {
	    cerr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	    cerr << "the Integrator failed in formUnbalance()\n";	
	    return -2;
	}	

	//new residual 
        const Vector &Resid = theSOE->getB() ;

	//new value of s 
	double s = - ( dx0 ^ Resid ) ;

	//intialize r = ratio of residuals 
	double r = 0.0 ;
        double r0 = 0.0 ;
	
        if ( s0 != 0.0 ) 
	    r = fabs( s / s0 ) ;
	
	r0 = r ;

	/******************************************
	if  ( r <= tolerance )
	    cerr << "Line Search Not Required : ";
	    cerr << "Residual Decrease Less Than Tolerance" << endl;
        else
	    cerr << "Line Search, Iteration " << 0 
	         << " : Ratio |s/s0| = " << r 
	         << endl ;
	***************************************/
	
	double eta = 1.0 ; //initial value of line search parameter

	int count = 0 ; //intial value of iteration counter 

	dxOld = dx0 ;

       
	while ( r > tolerance  &&  count < 10 ) {
	
	    count++ ;

	    eta *=  -s0 / (s - s0) ; 

            //-- I don't know if these should be here-----
	    if ( eta > 10.0 )  eta = 10.0 ;
	    if (   r > r0   )  eta =  1.0 ;
	    //--------------------------------------------
	    
	    //dx = ( eta * dx0 ) ; 
	    dx = dx0 ;	
	    dx *= eta ;
	    //Dr. Francis Thomas McKenna wants to save one Vector constructor call
	    
	    if (theIntegrator->update( dx - dxOld ) < 0) {
	      cerr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	      cerr << "the Integrator failed in update()\n";	
	      return -4;
	    }

  	    if (theIntegrator->formUnbalance() < 0) {
	      cerr << "WARNING NewtonLineSearch::solveCurrentStep() -";
	      cerr << "the Integrator failed in formUnbalance()\n";	
	      return -2;
	    }	

	    //new residual
	    const Vector &ResidI = theSOE->getB() ;

	    //new value of s
	    s = - ( dx0 ^ ResidI ) ;

	    //new value of r 
	    r = fabs( s / s0 ) ; 

	    cerr << "Line Search, Iteration " << count 
		 << " : Ratio |s/s0| = " << r << endl ;

	    //swap increments 
	    dxOld = dx ; 

	    if ( count == 10 ) 
		cerr << "Line Search Terminated After 10 Iterations" << endl ;
	    
	} //end while
	
	this->record(0);
	result = theTest->test();

    } while (result == -1);

    if (result == -2) {
      cerr << "NewtonLineSearch::solveCurrentStep() -";
      cerr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }

    // note - if postive result we are returning what the convergence test returned
    // which should be the number of iterations
    return result;
}

ConvergenceTest *
NewtonLineSearch::getTest(void)
{
  return theTest;
}

int
NewtonLineSearch::sendSelf(int cTag, Channel &theChannel)
{
  int result = 0;
  int dataTag = this->getDbTag();
  ID data(2);
  data(0) = theTest->getClassTag();
  data(1) = theTest->getDbTag();
  result = theChannel.sendID(dataTag, cTag, data);
  if (result != 0) {
    cerr << "NewtonLineSearch::sendSelf() - failed to send ID\n";
    return result;
  }

  Vector tol(1);
  tol(0) = tolerance;
  result = theChannel.sendVector(dataTag, cTag, tol);
  if (result != 0) {
    cerr << "NewtonLineSearch::sendSelf() - failed to send Vector\n";
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
NewtonLineSearch::recvSelf(int cTag, 
			Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    ID data(2);
    int result;
    int dataTag = this->getDbTag();

    result = theChannel.recvID(dataTag, cTag, data);    
    if (result != 0) {
      cerr << "NewtonLineSearch::recvSelf() - failed to receive ID\n";
      return result;
    }
    int ctType = data(0);
    int ctDb = data(1);

    
    Vector tol(1);
    result = theChannel.recvVector(dataTag, cTag, tol);
    if (result != 0) {
	cerr << "NewtonLineSearch::sendSelf() - failed to send Vector\n";
	return result;
    }      
    tolerance = tol(0);
    
    theTest = theBroker.getNewConvergenceTest(ctType);
    theTest->setDbTag(ctDb);
    result = theTest->recvSelf(cTag, theChannel, theBroker);
    if (result != 0) {
      cerr << "NewtonLineSearch::recvSelf() - failed to recv CTest object\n";
      return result;
    }
    
    return 0;
}


void
NewtonLineSearch::Print(ostream &s, int flag)
{
    if (flag == 0) 
	s << "NewtonLineSearch :: Line Search Tolerance = " << tolerance << endl ; 
}









