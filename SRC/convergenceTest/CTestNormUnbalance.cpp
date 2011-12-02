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


CTestNormUnbalance::CTestNormUnbalance()	    	
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
    theSOE(0), tol(0.0), maxNumIter(0), currentIter(0), printFlag(0),
    norms(1), nType(2)
{
    
}


CTestNormUnbalance::CTestNormUnbalance(double theTol, int maxIter, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
    theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
    norms(maxNumIter), nType(normType)
{
    
}


CTestNormUnbalance::~CTestNormUnbalance()
{
    
}


ConvergenceTest* CTestNormUnbalance::getCopy(int iterations)
{
    CTestNormUnbalance *theCopy ;
    theCopy = new CTestNormUnbalance(this->tol, iterations, this->printFlag, this->nType) ;
    
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
    if (theSOE == 0) {
        opserr << "WARNING: CTestNormUnbalance::setEquiSolnAlgo - no SOE\n";	
        return -1;
    }
    else
        return 0;
}


int CTestNormUnbalance::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0)
        return -2;
    
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
                opserr << "CTestNormUnbalance::test() - iteration: " << currentIter;
                opserr << " current Norm: " << norm << " (max: " << tol;
                opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
            }
        }
        
        // return the number of times test has been called
        return currentIter;
    } 
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
        opserr << "WARNING: CTestNormUnbalance::test() - failed to converge but going on -";
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter) { // the algorithm failed to converge
        opserr << "WARNING: CTestNormUnbalance::test() - failed to converge \n";
        opserr << "after: " << currentIter << " iterations\n";	
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
    static Vector x(4);
    x(0) = tol;
    x(1) = maxNumIter;
    x(2) = printFlag;
    x(3) = nType;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0) 
        opserr << "CTestNormUnbalance::sendSelf() - failed to send data\n";
    
    return res;
}


int CTestNormUnbalance::recvSelf(int cTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector x(4);
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
    }
    return res;
}
