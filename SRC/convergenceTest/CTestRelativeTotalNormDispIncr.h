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
// $Date: 2005-12-19 22:21:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestRelativeTotalNormDispIncr.h,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 05/05
// Revision: A
//
// Purpose: This file contains the class definition for CTestRelativeTotalNormDispIncr.
// A CTestRelativeTotalNormDispIncr object tests for convergence using the ratio of the 
// current norm to the total norm (the total norm since start was invoked) of the 
// solution vector of the LinearSOE object passed in the constructor
// and a tolerance, set in the constructor.

#ifndef CTestRelativeTotalNormDispIncr_h
#define CTestRelativeTotalNormDispIncr_h

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;


class CTestRelativeTotalNormDispIncr: public ConvergenceTest
{
public:
    // constructors
    CTestRelativeTotalNormDispIncr();
    CTestRelativeTotalNormDispIncr(double tol, int maxNumIter, int printFlag, int normType =2);
    
    // destructor
    ~CTestRelativeTotalNormDispIncr();
    
    ConvergenceTest *getCopy(int iterations);
    
    void setTolerance(double newTol);
    int setEquiSolnAlgo(EquiSolnAlgo &theAlgo);
    
    int test(void);
    int start(void);
    
    int getNumTests(void);
    int getMaxNumTests(void);
    double getRatioNumToMax(void);
    const Vector &getNorms(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
protected:
    
private:
    LinearSOE *theSOE;
    double tol;         // the tol on the norm used to test for convergence
    
    int maxNumIter;     // max number of iterations
    int currentIter;    // number of times test() has been invokes since last start()
    int printFlag;      // a flag indicating if to print on test
    int nType;          // type of norm to use (1-norm, 2-norm, p-norm, max-norm)
    
    Vector norms;       // vector to hold the norms
    double totNorm;     // norm at first iteration of each step
};

#endif
