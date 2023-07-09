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

// $Revision: 1.5 $
// $Date: 2007-05-15 21:55:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestRelativeNormUnbalance.cpp,v $


#include <CTestRelativeNormUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <elementAPI.h>

void* OPS_CTestRelativeNormUnbalance()
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

    // create test
    return new CTestRelativeNormUnbalance(tol,data[0],data[1],data[2]);
}

CTestRelativeNormUnbalance::CTestRelativeNormUnbalance()	    	
    : ConvergenceTest(CONVERGENCE_TEST_CTestRelativeNormUnbalance),
    theSOE(0), tol(0.0), maxNumIter(0), currentIter(0), printFlag(0),
    norms(1), norm0(0.0), nType(2)
{
    
}


CTestRelativeNormUnbalance::CTestRelativeNormUnbalance(double theTol, int maxIter, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestRelativeNormUnbalance),
    theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
    norms(maxNumIter+1), norm0(0.0), nType(normType)
{
    
}


CTestRelativeNormUnbalance::~CTestRelativeNormUnbalance()
{
    
}


ConvergenceTest* CTestRelativeNormUnbalance::getCopy(int iterations)
{
    CTestRelativeNormUnbalance *theCopy ;
    theCopy = new CTestRelativeNormUnbalance(this->tol, iterations, this->printFlag, this->nType) ;
    
    theCopy->theSOE = this->theSOE ;
    
    return theCopy ;
}


void CTestRelativeNormUnbalance::setTolerance(double newTol)
{
    tol = newTol;
}


int CTestRelativeNormUnbalance::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
    
    return 0;
}


int CTestRelativeNormUnbalance::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0)  {
        opserr << "WARNING: CTestRelativeNormUnbalance::test() - no SOE set.\n";
        return -1;
    }
    
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestRelativeNormUnbalance::test() - start() was never invoked.\n";	
        return -2;
    }
    
    // get the B vector & determine it's norm & save the value in norms vector
    const Vector &x = theSOE->getB();
    double norm = x.pNorm(nType);
    if (currentIter <= maxNumIter) 
        norms(currentIter) = norm;
    
    // determine the ratio
    if (norm0 != 0.0)
        norm /= norm0;
    
    // print the data if required
    if (printFlag == 1) {
        opserr << "CTestRelativeNormUnbalance::test() - iteration: " << currentIter;
        opserr << " current Ratio (|dR|/|dR0|): " << norm << " (max: " << tol << ")\n";
    }
    if (printFlag == 4) {
        opserr << "CTestRelativeNormUnbalance::test() - iteration: " << currentIter;
        opserr << " current Ratio (|dR|/|dR0|): " << norm << " (max: " << tol << ")\n";
        opserr << "\tNorm deltaX: " << theSOE->getX().pNorm(nType) << ", Norm deltaR: " << norm << endln;
        opserr << "\tdeltaX: " << theSOE->getX() << "\tdeltaR: " << x;
    }
    
    //
    // check if the algorithm converged
    //
    
    // if converged - print & return ok
    
    if (norm <= tol) { // the algorithm converged
        
        // do some printing first
        if (printFlag != 0) {
            if (printFlag == 1 || printFlag == 4) 
                opserr << endln;
            else if (printFlag == 2 || printFlag == 6) {
                opserr << "CTestRelativeNormUnbalance::test() - iteration: " << currentIter;
                opserr << " current Ratio (|dR|/|dR0|): " << norm << " (max: " << tol << ")\n";
            }
        }
        
        // return the number of times test has been called
        return currentIter;
    } 
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
        opserr << "WARNING: CTestRelativeNormUnbalance::test() - failed to converge but going on -";
        opserr << " current Ratio (dR/dR0): " << norm << " (max: " << tol;
        opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter) { // the algorithm failed to converge
        opserr << "WARNING: CTestRelativeNormUnbalance::test() - failed to converge \n";
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


int CTestRelativeNormUnbalance::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING: CTestRelativeNormUnbalance::test() - no SOE returning true\n";
        return -1;
    }
    
    // set iteration count = 1
    norms.Zero();
    currentIter = 1;    
    norm0 = 0.0;
    
    // determine the initial norm .. the the norm of the initial unbalance
    const Vector &x = theSOE->getB();
    double norm = x.pNorm(nType);
    if (currentIter <= maxNumIter) 
        norms(0) = norm;
    norm0 = norm;
    
    return 0;
}


int CTestRelativeNormUnbalance::getNumTests()
{
    return currentIter;
}


int CTestRelativeNormUnbalance::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestRelativeNormUnbalance::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestRelativeNormUnbalance::getNorms() 
{
    return norms;
}


int CTestRelativeNormUnbalance::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    Vector x(4);
    x(0) = tol;
    x(1) = maxNumIter;  
    x(2) = printFlag;
    x(3) = nType;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0) 
        opserr << "CTestRelativeNormUnbalance::sendSelf() - failed to send data\n";
    
    return res;
}


int CTestRelativeNormUnbalance::recvSelf(int cTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    Vector x(4);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);    
    
    if (res < 0) {
        opserr << "CTestRelativeNormUnbalance::sendSelf() - failed to send data\n";
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
