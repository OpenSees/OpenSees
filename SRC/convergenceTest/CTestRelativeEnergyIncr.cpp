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

// $Revision: 1.3 $
// $Date: 2005-12-15 00:19:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestRelativeEnergyIncr.cpp,v $

// Written: fmk 
// Date: 02/02


#include <CTestRelativeEnergyIncr.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>


CTestRelativeEnergyIncr::CTestRelativeEnergyIncr()	    	
    : ConvergenceTest(CONVERGENCE_TEST_CTestRelativeEnergyIncr),
    theSOE(0), tol(0), maxNumIter(0), currentIter(0), printFlag(0),
    norms(1), norm0(0.0), nType(2)
{
    
}


CTestRelativeEnergyIncr::CTestRelativeEnergyIncr(double theTol, int maxIter, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestRelativeEnergyIncr),
    theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0),printFlag(printIt),
    norms(maxNumIter), norm0(0.0), nType(normType)
{
    
}


CTestRelativeEnergyIncr::~CTestRelativeEnergyIncr()
{
    
}


ConvergenceTest* CTestRelativeEnergyIncr::getCopy(int iterations)
{
    CTestRelativeEnergyIncr *theCopy ;
    theCopy = new CTestRelativeEnergyIncr(this->tol, iterations, this->printFlag, this->nType);
    
    theCopy->theSOE = this->theSOE ;
    
    return theCopy ;
}


void CTestRelativeEnergyIncr::setTolerance(double newTol)
{
    tol = newTol;
}


int CTestRelativeEnergyIncr::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
    if (theSOE == 0) {
        opserr << "WARNING: CTestRelativeEnergyIncr::setEquiSolnAlgo - no SOE\n";	
        return -1;
    }
    else
        return 0;
}


int CTestRelativeEnergyIncr::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0)
        return -2;
    
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestRelativeEnergyIncr::test() - start() was never invoked.\n";	
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
    
    // if first pass through .. set norm0
    if (currentIter == 1) {
        norm0 = product;
    }
    
    // get ratio
    if (norm0 != 0.0)
        product /= norm0;
    
    // print the data if required
    if (printFlag == 1) {
        opserr << "CTestRelativeEnergyIncr::test() - iteration: " << currentIter;
        opserr << " current Ratio (dX*dR/dX1*dR1): " << product << " (max: " << tol << ")\n";
    }
    if (printFlag == 4) {
        opserr << "CTestRelativeEnergyIncr::test() - iteration: " << currentIter;
        opserr << " current Ratio (dX*dR/dX1*dR1): " << product << " (max: " << tol << ")\n";
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
                opserr << "CTestRelativeEnergyIncr::test() - iteration: " << currentIter;
                opserr << " last Ratio (dX*dR/dX1*dR1): " << product << " (max: " << tol << ")\n";
            }
        }
        
        // return the number of times test has been called - SUCCESSFULL
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
        opserr << "WARNING: CTestRelativeEnergyIncr::test() - failed to converge but goin on -";
        opserr << " current Ratio (dX*dR/dX1*dR1): " << product << " (max: " << tol << ")\n";
        opserr << "\tNorm deltaX: " << x.pNorm(nType) << ", Norm deltaR: " << b.pNorm(nType) << endln;
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter) { // >= in case algorithm does not check
        opserr << "WARNING: CTestRelativeEnergyIncr::test() - failed to converge \n";
        opserr << "after: " << currentIter << " iterations\n";	
        currentIter++;    
        return -2;
    } 
    
    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;    
        return -1;
    }
}


int CTestRelativeEnergyIncr::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING: CTestRelativeEnergyIncr::test() - no SOE returning true\n";
        return -1;
    }
    
    currentIter = 1;
    norms.Zero();
    norm0 = 0.0;
    
    return 0;
}


int CTestRelativeEnergyIncr::getNumTests(void)
{
    return currentIter;
}


int CTestRelativeEnergyIncr::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestRelativeEnergyIncr::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestRelativeEnergyIncr::getNorms(void) 
{
    return norms;
}


int CTestRelativeEnergyIncr::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector x(4);
    x(0) = tol;
    x(1) = maxNumIter;
    x(2) = printFlag;
    x(3) = nType;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0) 
        opserr << "CTestRelativeEnergyIncr::sendSelf() - failed to send data\n";
    
    return res;
}


int CTestRelativeEnergyIncr::recvSelf(int cTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector x(4);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);    
    
    if (res < 0) {
        opserr << "CTestRelativeEnergyIncr::sendSelf() - failed to send data\n";
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
    }
    
    return res;
}
