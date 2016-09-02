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

// $Revision: 1.8 $
// $Date: 2006-01-10 18:13:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestNormDispIncr.cpp,v $


#include <CTestNormDispIncr.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <elementAPI.h>

void* OPS_CTestNormDispIncr()
{
    if(OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"insufficient number of arguments\n";
	return 0;
    }

    // tolerance
    double tol = 1e-6;
    int numData = 1;
    if(OPS_GetDoubleInput(&numData,&tol) < 0) {
	opserr << "WARNING NormUnbalance failed to read tol\n";
	return 0;
    }

    // maxIter
    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 3) numData = 3;
    int data[3] = {0,0,2};
    if(OPS_GetIntInput(&numData,&data[0]) < 0) {
	opserr << "WARNING NormUnbalance failed to read int values\n";
	return 0;
    }

    // maxTol
    double maxTol = OPS_MAXTOL;
    if(OPS_GetNumRemainingInputArgs() > 0) {
	numData = 1;
	if(OPS_GetDoubleInput(&numData,&maxTol) < 0) {
	    opserr << "WARNING NormUnbalance failed to read maxTol\n";
	    return 0;
	}
    }
    
    // create test
    return new CTestNormDispIncr(tol,data[0],data[1],data[2],maxTol);
}

CTestNormDispIncr::CTestNormDispIncr()	    	
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
      theSOE(0), tol(0), maxTol(OPS_MAXTOL), maxNumIter(0), currentIter(0), printFlag(0), 
      norms(25), nType(2)
{
    
}


CTestNormDispIncr::CTestNormDispIncr(double theTol, int maxIter, int printIt, int normType, double max)
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
      theSOE(0), tol(theTol), maxTol(max), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
      nType(normType), norms(maxIter)
{
    
}


CTestNormDispIncr::~CTestNormDispIncr()
{
    
}


ConvergenceTest* CTestNormDispIncr::getCopy(int iterations)
{
    CTestNormDispIncr *theCopy ;
    theCopy = new CTestNormDispIncr(this->tol, iterations, 0, this->nType, this->maxTol) ;
    
    theCopy->theSOE = this->theSOE ;
    
    return theCopy ;
}


void CTestNormDispIncr::setTolerance(double newTol)
{
    tol = newTol;
}


int CTestNormDispIncr::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
   
    return 0;
}


int CTestNormDispIncr::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0) {
		  opserr << "WARNING: CTestNormDispIncr::test() - no SOE set.\n";	
        return -2;
	}
    
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestNormDispIncr::test() - start() was never invoked.\n";	
        return -2;
    }
    
    // get the X vector & determine it's norm & save the value in norms vector
    const Vector &x = theSOE->getX();
    double norm = x.pNorm(nType);
    if (currentIter <= maxNumIter) 
        norms(currentIter-1) = norm;
    
    // print the data if required
    if (printFlag == 1) {
        opserr << "CTestNormDispIncr::test() - iteration: " << currentIter;
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << ")\n";
    } 
    if (printFlag == 4) {
        opserr << "CTestNormDispIncr::test() - iteration: " << currentIter;
        opserr << " current Norm: " << norm << " (max: " << tol << ")\n";
        opserr << "\tNorm deltaX: " << norm << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << endln;
        opserr << "\tdeltaX: " << x << "\tdeltaR: " << theSOE->getB();
    } 
    
    //
    // check if the algorithm converged
    //
    
    // if converged - print & return ok
    if (norm <= tol) { 
        
        // do some printing first
        if (printFlag != 0) {
            if (printFlag == 1 || printFlag == 4) 
                opserr << endln;
            else if (printFlag == 2 || printFlag == 6) {
                opserr << "CTestNormDispIncr::test() - iteration: " << currentIter;
                opserr << " current Norm: " << norm << " (max: " << tol;
                opserr << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << ")\n";
            }
        }
        
        // return the number of times test has been called
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
        opserr << "WARNING: CTestNormDispIncr::test() - failed to converge but going on - ";
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << ")\n";
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter || norm > maxTol) { // failes to converge
        opserr << "WARNING: CTestNormDispIncr::test() - failed to converge \n";
        opserr << "after: " << currentIter << " iterations ";	
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << ")\n";
        currentIter++;    
        return -2;
    } 
    
    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;    
        return -1;
    }
}


int CTestNormDispIncr::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING: CTestNormDispIncr::test() - no SOE returning true\n";
        return -1;
    }
    
    // set iteration count = 1
    norms.Zero();
    currentIter = 1;
    return 0;
}


int CTestNormDispIncr::getNumTests()
{
    return currentIter;
}


int CTestNormDispIncr::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestNormDispIncr::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestNormDispIncr::getNorms() 
{
    return norms;
}


int CTestNormDispIncr::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  Vector x(5);
  x(0) = tol;
  x(1) = maxNumIter;
  x(2) = printFlag;
  x(3) = nType;
  x(4) = maxTol;
  res = theChannel.sendVector(this->getDbTag(), cTag, x);
  if (res < 0) 
    opserr << "CTestNormDispIncr::sendSelf() - failed to send data\n";
    
  return res;
}

int 
CTestNormDispIncr::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  Vector x(5);
  res = theChannel.recvVector(this->getDbTag(), cTag, x);    


  if (res < 0) {
    opserr << "CTestNormDispIncr::sendSelf() - failed to send data\n";
    tol = 1.0e-8;
    maxNumIter = 25;
    printFlag = 0;
    nType = 2;
    norms.resize(maxNumIter);
  } else {
    tol = x(0);
    maxNumIter = (int)x(1);
    printFlag = (int)x(2);
    nType = (int)x(3);
    norms.resize(maxNumIter);
    maxTol = x(4);
  } 
  return res;
}


