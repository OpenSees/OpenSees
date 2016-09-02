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
#include <elementAPI.h>

void* OPS_CTestEnergyIncr()
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
    return new CTestEnergyIncr(tol,data[0],data[1],data[2],maxTol);
}

CTestEnergyIncr::CTestEnergyIncr()	    	
    : ConvergenceTest(CONVERGENCE_TEST_CTestEnergyIncr),
      theSOE(0), tol(0), maxTol(OPS_MAXTOL), maxNumIter(0), currentIter(0), printFlag(0),
      nType(2), norms(20)
{
    
}


CTestEnergyIncr::CTestEnergyIncr(double theTol, int maxIter, int printIt, int normType, double max)
    : ConvergenceTest(CONVERGENCE_TEST_CTestEnergyIncr),
      theSOE(0), tol(theTol), maxTol(max), maxNumIter(maxIter), currentIter(0),printFlag(printIt),
      nType(normType), norms(maxNumIter)
{
    
}


CTestEnergyIncr::~CTestEnergyIncr()
{
    
}


ConvergenceTest* CTestEnergyIncr::getCopy(int iterations)
{
    CTestEnergyIncr *theCopy ;
    theCopy = new CTestEnergyIncr(this->tol, iterations, this->printFlag, this->nType, this->maxTol);
    
    theCopy->theSOE = this->theSOE ;
    
    return theCopy ;
}


void CTestEnergyIncr::setTolerance(double newTol)
{
    tol = newTol;
}


int CTestEnergyIncr::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
	return 0;
}


int CTestEnergyIncr::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0) {
		opserr << "WARNING: CTestEnergyIncr::test() - no SOE set\n";
        return -2;
	}
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestEnergyIncr::test() - start() was never invoked.\n";	
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
        opserr << "CTestEnergyIncr::test() - iteration: " << currentIter;
        opserr << " current EnergyIncr: " << product << " (max: " << tol << ")\n";
    }
    if (printFlag == 4) {
        opserr << "CTestEnergyIncr::test() - iteration: " << currentIter;
        opserr << " current EnergyIncr: " << product << " (max: " << tol << ")\n";
        opserr << "\tNorm deltaX: " << x.pNorm(nType) << ", Norm deltaR: " << b.pNorm(nType) << endln;
        opserr << "\tdeltaX: " << x << "\tdeltaR: " << b;
    }
    
    //
    // check if the algorithm converged
    //
    
    // if converged - print & return ok
    if (product <= tol) {
        
        // do some printing first
        if (printFlag != 0) {
            if (printFlag == 1 || printFlag == 4) 
                opserr << endln;
            else if (printFlag == 2 || printFlag == 6) {
                opserr << "CTestEnergyIncr::test() - iteration: " << currentIter;
                opserr << " last EnergyIncr: " << product << " (max: " << tol << ")\n";
            }
        }
        
        // return the number of times test has been called - SUCCESSFULL
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
        opserr << "WARNING: CTestEnergyIncr::test() - failed to converge but goin on -";
        opserr << " current EnergyIncr: " << product << " (max: " << tol << ")\n";
        opserr << "\tNorm deltaX: " << x.pNorm(nType) << ", Norm deltaR: " << b.pNorm(nType) << endln;
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter || product > maxTol) { // >= in case algorithm does not check
        opserr << "WARNING: CTestEnergyIncr::test() - failed to converge \n";
        opserr << "after: " << currentIter << " iterations\n";	
        opserr << " current EnergyIncr: " << product << " (max: " << tol << ") ";
        opserr << "\tNorm deltaX: " << x.pNorm(nType) << ", Norm deltaR: " << b.pNorm(nType) << endln;
        currentIter++;    
        return -2;
    } 
    
    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;    
        return -1;
    }
}


int CTestEnergyIncr::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING: CTestEnergyIncr::test() - no SOE returning true\n";
        return -1;
    }
    
    // set iteration count = 1
    currentIter = 1;
    norms.Zero();
    return 0;
}


int CTestEnergyIncr::getNumTests(void)
{
    return currentIter;
}


int CTestEnergyIncr::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestEnergyIncr::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestEnergyIncr::getNorms(void) 
{
    return norms;
}


int CTestEnergyIncr::sendSelf(int cTag, Channel &theChannel)
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
        opserr << "CTestEnergyIncr::sendSelf() - failed to send data\n";
    
    return res;
}


int CTestEnergyIncr::recvSelf(int cTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector x(5);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);    
    
    if (res < 0) {
        opserr << "CTestEnergyIncr::sendSelf() - failed to send data\n";
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
