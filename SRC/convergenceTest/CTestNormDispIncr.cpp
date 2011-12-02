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
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestNormDispIncr.cpp,v $
                                                                        
                                                                        

#include <CTestNormDispIncr.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>

CTestNormDispIncr::CTestNormDispIncr()	    	
:ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
 theSOE(0), tol(0), maxNumIter(0), currentIter(0), printFlag(0), totalNumIter(1)
{
    
}

CTestNormDispIncr::CTestNormDispIncr(double theTol, int maxIter, int printIt)
:ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
 theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
 totalNumIter(1)
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
	cerr << "WARNING: CTestNormDisp::setEquiSolnAlgo - no SOE\n";	
	return -1;
    }
    else
	return 0;
}

    
int
CTestNormDispIncr::test(void)
{
    // this should not happen if the rturn from start() is checked
    if (theSOE == 0)
       return -2;

    const Vector &x = theSOE->getX();
	
    double norm = x.Norm();
// cerr << "called norm: " << norm << "Vector: " << x << endl;
    // print the data if required
    if (printFlag == 1) {
      cerr << "\t CTestNormDispIncr::test() - iteration: " << currentIter;
      cerr << " current Norm: " << norm << " (max permissable: " << tol << ")\n";
    } 

    if (norm <= tol){ // the algorithm converged
      if (printFlag != 0) {
	if (printFlag == 1) 
	  cerr << endl;
	else if (printFlag == 2) {
	  cerr << "\t CTestNormDispIncr::test() - iteration: " << currentIter;
	  cerr << " current Norm: " << norm << " (max permissable: " << tol << ")\n";
	}
	else if (printFlag == 3) {
	  cerr << "\t CTestNormDispIncr::test() - #iteration: " << currentIter;
	  cerr << "(total #iteration: " << totalNumIter << ") current Norm: ";
	  cerr << norm << " (max permissable: " << tol << ")\n";
	} 

	  
      }
      return currentIter;
    }

	else if (printFlag == 5 && currentIter >= maxNumIter) {
		cerr << "WARNING: CTestDispIncr::test() - failed to converge but GOING ON\n";
		return currentIter;
	}
    else if (currentIter >= maxNumIter) { // failes to converge
	cerr << "WARNING: CTestDispIncr::test() - failed to converge \n";
	cerr << "after: " << currentIter << " iterations\n";	
	currentIter++;    
	return -2;
    } 

    else { // has not yet converged
      currentIter++;    
      totalNumIter++;    
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
    currentIter = 1;
    return 0;
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
  }
  return res;
}

