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
// $Date: 2000-12-12 07:58:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestEnergyIncr.cpp,v $
                                                                        
                                                                        
// File: ~/convergenceTest/CTestEnergyIncr.C
//
// Written: fmk 
// Date: 09/98
// Revised:
//
// Purpose: This file contains the class implementation for CTestEnergyIncr.
// A CTestEnergyIncr object tests for convergence using the energy increment,
// which is 0.5 times the absolute value of the product of the rhs and 
// the solution vector of the LinearSOE.

#include <CTestEnergyIncr.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>

CTestEnergyIncr::CTestEnergyIncr()	    	
:ConvergenceTest(CONVERGENCE_TEST_CTestEnergyIncr),
 theSOE(0), tol(0), maxNumIter(0), currentIter(0), printFlag(0),
 norms(1)
{
    
}

CTestEnergyIncr::CTestEnergyIncr(double theTol, int maxIter, int printIt)
:ConvergenceTest(CONVERGENCE_TEST_CTestEnergyIncr),
 theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0),printFlag(printIt),
 norms(maxNumIter)
{

}

CTestEnergyIncr::~CTestEnergyIncr()
{
    
}

void
CTestEnergyIncr::setTolerance(double newTol)
{
    tol = newTol;
}

int
CTestEnergyIncr::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
    if (theSOE == 0) {
	cerr << "WARNING: CTestEnergyIncr::setEquiSolnAlgo - no SOE\n";	
	return -1;
    }
    else
	return 0;
}

    
int
CTestEnergyIncr::test(void)
{
  // check to ensure the SOE has been set - this should not happen if the 
  // return from start() is checked
  if (theSOE == 0)
    return -2;

  // check to ensure the algo does invoke start() - this is needed otherwise
  // may never get convergence later on in analysis!
  if (currentIter == 0) {
    cerr << "WARNING: CTestEnergyIncr::test() - start() was never invoked.\n";	
    return -2;
  }


  // determine the energy & save value in norms vector
  const Vector &b = theSOE->getB();
  const Vector &x = theSOE->getX();    
  double product = x ^ b;
  if (product < 0.0)
    product *= -0.5;
  else
    product *= 0.5;

  if (currentIter <= maxNumIter) 
    norms(currentIter-1) = product;

  // print the data if required
  if (printFlag == 1) {
    cerr << "CTestEnergyIncr::test() - iteration: " << currentIter;
    cerr << " current Product: " << product << " (max permissable: " << tol << ")\n";
  }


  //
  // check if the algorithm converged
  //

  // if converged - print & return ok
  if (product <= tol) {

    // do some printing first
    if (printFlag != 0) {
      if (printFlag == 1) 
	cerr << endl;
      else if (printFlag == 2) {
	cerr << "CTestEnergyIncr::test() - iteration: " << currentIter;
	cerr << " last incr: " << product << " (max permissable: " << tol << ")\n";
      }
    }

    // return the number of times test has been called - SUCCESSFULL
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - but RETURN OK
  else if (printFlag == 5 && currentIter >= maxNumIter) {
    cerr << "WARNING: CTestEnergyIncr::test() - failed to converge but GOING ON\n";
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - return FAILURE -2
  else if (currentIter >= maxNumIter) { // >= in case algorithm does not check
    cerr << "WARNING: CTestEnergyIncr::test() - failed to converge \n";
    cerr << "after: " << currentIter << " iterations\n";	
    currentIter++;    
    return -2;
  } 

  // algo not yet converged - increment counter and return -1
  else {
      currentIter++;    
      return -1;
    }
}

int
CTestEnergyIncr::start(void)
{
    if (theSOE == 0) {
	cerr << "WARNING: CTestEnergyIncr::test() - no SOE returning true\n";
	return -1;
    }

    currentIter = 1;
    norms.Zero();
    return 0;
}


int 
CTestEnergyIncr::getNumTests(void)
{
  return currentIter;
}

int 
CTestEnergyIncr::getMaxNumTests(void)
{
  return maxNumIter;
}

double 
CTestEnergyIncr::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector &
CTestEnergyIncr::getNorms(void) 
{
  return norms;
}


int 
CTestEnergyIncr::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  Vector x(2);
  x(0) = tol;
  x(1) = maxNumIter;

  
  res = theChannel.sendVector(this->getDbTag(), cTag, x);
  if (res < 0) 
    cerr << "CTestEnergyIncr::sendSelf() - failed to send data\n";
    
  return res;
}

int 
CTestEnergyIncr::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  Vector x(2);
  res = theChannel.recvVector(this->getDbTag(), cTag, x);    

  if (res < 0) {
      cerr << "CTestEnergyIncr::sendSelf() - failed to send data\n";
      tol = 1.0e-8;
      maxNumIter = 25;
  }
  else {
      tol = x(0);
      maxNumIter = x(1);
      norms.resize(maxNumIter);
  }
    
  return res;
}






