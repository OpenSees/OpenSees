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
// $Date: 2005-12-15 00:19:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestNormUnbalance.cpp,v $


#include <CTestNormUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <elementAPI.h>
#include <iostream>
#include <fstream>

void* OPS_CTestNormUnbalance()
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
    if(numData > 4) numData = 4;
    int data[4] = {0,0,2,-1};
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
    return new CTestNormUnbalance(tol,data[0],data[1],data[2],data[3],maxTol);
}

CTestNormUnbalance::CTestNormUnbalance()	    	
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
      theSOE(0), tol(0.0), maxTol(OPS_MAXTOL), maxNumIter(0), currentIter(0), printFlag(0),
      norms(1), nType(2), maxIncr(0), numIncr(0)
{
    
}


CTestNormUnbalance::CTestNormUnbalance(double theTol, int maxIter, int printIt, int normType, int maxincr, double max)
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
      theSOE(0), tol(theTol), maxTol(max), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
      nType(normType), norms(maxNumIter), maxIncr(maxincr), numIncr(0)
{
    if(maxIncr < 0) {
        maxIncr = maxNumIter;
    }
}


CTestNormUnbalance::~CTestNormUnbalance()
{
    
}


ConvergenceTest* CTestNormUnbalance::getCopy(int iterations)
{
    CTestNormUnbalance *theCopy ;
    theCopy = new CTestNormUnbalance(this->tol, iterations, this->printFlag, this->nType, this->maxIncr, this->maxTol) ;
    
    theCopy->theSOE = this->theSOE ;
    
    return theCopy ;
}


void CTestNormUnbalance::setTolerance(double newTol)
{
    tol = newTol;
}


int CTestNormUnbalance::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
    return 0;
}


int CTestNormUnbalance::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0) {
		 opserr << "WARNING: CTestNormUnbalance::test() - no SOE set.\n";
        return -2;
	}
    
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestNormUnbalance::test() - start() was never invoked.\n";	
        return -2;
    }
    
    // get the B vector & determine it's norm & save the value in norms vector
    const Vector &x = theSOE->getB();
    double norm = x.pNorm(nType);
    if (currentIter <= maxNumIter) 
        norms(currentIter-1) = norm;

    if(currentIter > 1) {
        if(norms(currentIter-2) < norm) {
            numIncr++;
        }
    }
    
    // print the data if required
    if (printFlag == 1) {
        opserr << "CTestNormUnbalance::test() - iteration: " << currentIter;
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
    }
    if (printFlag == 4) {
        opserr << "CTestNormUnbalance::test() - iteration: " << currentIter;
        opserr << " current Norm: " << norm << " (max: " << tol << ")\n";
        opserr << "\tNorm deltaX: " << theSOE->getX().pNorm(nType) << ", Norm deltaR: " << norm << endln;
        opserr << "\tdeltaX: " << theSOE->getX() << "\tdeltaR: " << x;
    }

    if (printFlag == 7) {
      std::ofstream outDu;
      std::ofstream outDp;

      if (currentIter == 1) {
	outDu.open ("dX.out",std::ios::out);
	outDp.open ("dP.out", std::ios::out);
      } else {
	outDu.open ("dX.out",std::ios::app);
	outDp.open ("dP.out", std::ios::app);
      }
      const Vector &Du = theSOE->getX();
      const Vector &Dp = theSOE->getB();
      for (int i=0; i<Du.Size(); i++) {
	outDu << Du[i] << " ";
	outDp << Dp[i] << " ";
      }
      outDu << "\n";
      outDp << "\n";
      outDu.close();
      outDp.close();
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
            else if (printFlag == 2 || printFlag == 6 || printFlag == 7) {
                opserr << "CTestNormUnbalance::test() - iteration: " << currentIter;
                opserr << " current Norm: " << norm << " (max: " << tol;
                opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
            }
        }
        
        // return the number of times test has been called
        return currentIter;
    } 
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && (currentIter >= maxNumIter||numIncr>=maxIncr)) {
        opserr << "WARNING: CTestNormUnbalance::test() - failed to converge but going on -";
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter || numIncr >= maxIncr || norm > maxTol) { // the algorithm failed to converge
        opserr << "WARNING: CTestNormUnbalance::test() - failed to converge \n";
        opserr << "after: " << currentIter << " iterations ";	
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
        currentIter++;  // we increment in case analysis does not check for convergence
        return -2;
    } 
    
    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;    
        return -1;
    }
}


int CTestNormUnbalance::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING: CTestNormUnbalance::test() - no SOE returning true\n";
        return -1;
    }
    
    // set iteration count = 1
    norms.Zero();
    currentIter = 1;    
    numIncr = 0;
    return 0;
}


int CTestNormUnbalance::getNumTests()
{
    return currentIter;
}


int CTestNormUnbalance::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestNormUnbalance::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestNormUnbalance::getNorms() 
{
    return norms;
}


int CTestNormUnbalance::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector x(5);
    x(0) = tol;
    x(1) = maxNumIter;
    x(2) = printFlag;
    x(3) = nType;
    x(4) = maxTol;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0) 
        opserr << "CTestNormUnbalance::sendSelf() - failed to send data\n";
    
    return res;
}


int CTestNormUnbalance::recvSelf(int cTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector x(5);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);    
    
    if (res < 0) {
        opserr << "CTestNormUnbalance::sendSelf() - failed to send data\n";
        tol = 1.0e-8;
        maxNumIter = 25;
        printFlag = 0;
        nType = 2;
    }
    else {
        tol = x(0);
        maxNumIter = (int) x(1);
        printFlag = (int) x(2);
        nType = (int) x(3);
        norms.resize(maxNumIter);
	maxTol = x(4);
    }
    return res;
}
