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
// $Date: 2000-09-15 08:23:17 $
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
 theSOE(0), tol(0), maxNumIter(0), currentIter(0), printFlag(0)
{
    
}

CTestEnergyIncr::CTestEnergyIncr(double theTol, int maxIter, int printIt)
:ConvergenceTest(CONVERGENCE_TEST_CTestEnergyIncr),
 theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0),printFlag(printIt)
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
    if (theSOE == 0) {
	return -2;
    }
    
    const Vector &b = theSOE->getB();
    const Vector &x = theSOE->getX();    
    double product = x ^ b;
    if (product < 0.0)
	product *= -0.5;
    else
	product *= 0.5;

    // print the data if required
    if (printFlag == 1) {
      cerr << "CTestEnergyIncr::test() - iteration: " << currentIter;
      cerr << " current Product: " << product << " (max permissable: " << tol << ")\n";
    }

    if (product <= tol) {
      if (printFlag != 0) {
	if (printFlag == 1) 
	  cerr << endl;
	else if (printFlag == 2) {
	  cerr << "CTestEnergyIncr::test() - iteration: " << currentIter;
	  cerr << " last incr: " << product << " (max permissable: " << tol << ")\n";
	}
      }
      return currentIter;
    }

    else if (currentIter >= maxNumIter) { // >= in case algorithm does not check
	cerr << "WARNING: CTestEnergyIncr::test() - failed to converge \n";
	cerr << "after: " << currentIter << " iterations\n";	
	currentIter++;    
	return -2;
    } 

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
    return 0;
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
  }
    
  return res;
}






