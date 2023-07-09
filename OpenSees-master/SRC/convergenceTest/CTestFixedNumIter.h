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

// $Revision: 1.1 $
// $Date: 2005-12-15 00:13:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestFixedNumIter.h,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 09/05
// Revision: A

// Purpose: This file contains the class definition for CTestFixedNumIter.
// A CTestFixedNumIter object performs a fixed number of iterations without
// testing for convergence. This test is useful for hybrid simulation where
// the residual error is corrected for.

#ifndef CTestFixedNumIter_h
#define CTestFixedNumIter_h

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;


class CTestFixedNumIter: public ConvergenceTest
{
public:
    // constructors
    CTestFixedNumIter();
    CTestFixedNumIter(int maxNumIter, int printFlag, int normType=2);
    
    // destructor
    ~CTestFixedNumIter();
    
    ConvergenceTest *getCopy(int iterations);
    
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
    
    int maxNumIter;     // max number of iterations
    int currentIter;    // number of times test() has been invokes since last start()
    int printFlag;      // a flag indicating if to print on test
    int nType;          // type of norm to use (1-norm, 2-norm, p-norm, max-norm)
    
    Vector norms;       // vector to hold the norms
};

#endif
