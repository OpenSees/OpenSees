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
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestNormDispIncr.cpp,v $
                                                                        
                                                                        

#include <CTestNormDispIncr.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>

CTestNormDispIncr::CTestNormDispIncr()	    	
:ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
 theSOE(0), tol(0), maxNumIter(0), currentIter(0), printFlag(0), 
 norms(1)
{
    
}

CTestNormDispIncr::CTestNormDispIncr(double theTol, int maxIter, int printIt)
:ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
 theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
 norms(maxIter)
{
    tol = theTol;
}

CTestNormDispIncr::~CTestNormDispIncr()
{
    
}

void
CTestNormDispIncr::setTolerance(double newTol)
{
    tol = newTol;
}

int
CTestNormDispIncr::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
    if (theSOE == 0) {
	cerr << "WARNING: CTestNormDisp::setEquiSolnAlgo() - no SOE\n";	
	return -1;
    }
    else
	return 0;
}

    
int
CTestNormDispIncr::test(void)
{
  // check to ensure the SOE has been set - this should not happen if the 
  // return from start() is checked
  if (theSOE == 0)
    return -2;

  // check to ensure the algo does invoke start() - this is needed otherwise
  // may never get convergence later on in analysis!
  if (currentIter == 0) {
    cerr << "WARNING: CTestNormDisp::test() - start() was never invoked.\n";	
    return -2;
  }

  // get the X vector & determine it's norm & save the value in norms vector
  const Vector &x = theSOE->getX();
  double norm = x.Norm();
  if (currentIter <= maxNumIter) 
    norms(currentIter-1) = norm;


  // print the data if required
  if (printFlag == 1) {
    cerr << "\t CTestNormDispIncr::test() - iteration: " << currentIter;
    cerr << " current Norm: " << norm << " (max permissable: " << tol << ")\n";
  } 

  //
  // check if the algorithm converged
  //

  // if converged - print & return ok
  if (norm <= tol){ 

    // do some printing first
    if (printFlag != 0) {
      if (printFlag == 1) 
	cerr << endl;
      else if (printFlag == 2) {
	cerr << "\t CTestNormDispIncr::test() - iteration: " << currentIter;
	cerr << " current Norm: " << norm << " (max permissable: " << tol << ")\n";
      }
    }

    // return the number of times test has been called
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - but RETURN OK
  else if (printFlag == 5 && currentIter >= maxNumIter) {
    cerr << "WARNING: CTestDispIncr::test() - failed to converge but GOING ON\n";
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - return FAILURE -2
  else if (currentIter >= maxNumIter) { // failes to converge
    cerr << "WARNING: CTestDispIncr::test() - failed to converge \n";
    cerr << "after: " << currentIter << " iterations\n";	
    currentIter++;    
    return -2;
  } 

  // algo not yet converged - increment counter and return -1
  else { // has not yet converged
    currentIter++;    
    return -1;
  }
}


int
CTestNormDispIncr::start(void)
{
    if (theSOE == 0) {
	cerr << "WARNING: CTestNormDispIncr::test() - no SOE returning true\n";
	return -1;
    }

    // set iteration count = 1
    norms.Zero();
    currentIter = 1;
    return 0;
}


int 
CTestNormDispIncr::getNumTests()
{
  return currentIter;
}


int 
CTestNormDispIncr::getMaxNumTests(void)
{
  return maxNumIter;
}

double 
CTestNormDispIncr::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector &
CTestNormDispIncr::getNorms() 
{
  return norms;
}

int 
CTestNormDispIncr::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  Vector x(2);
  x(0) = tol;
  x(1) = maxNumIter;
  res = theChannel.sendVector(this->getDbTag(), cTag, x);
  if (res < 0) 
    cerr << "CTestNormDispIncr::sendSelf() - failed to send data\n";
    
  return res;
}

int 
CTestNormDispIncr::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  Vector x(2);
  res = theChannel.recvVector(this->getDbTag(), cTag, x);    

  if (res < 0) {
      cerr << "CTestNormDispIncr::sendSelf() - failed to send data\n";
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

