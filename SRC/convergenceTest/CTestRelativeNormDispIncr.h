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
// $Date: 2005-12-15 00:19:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestRelativeNormDispIncr.h,v $

// Written: fmk 
// Date: 02/02
// Modified: 05/05 ahs
//
// Purpose: This file contains the class definition for CTestRelativeNormDispIncr.
// A CTestRelativeNormDispIncr object tests for convergence using the ratio of the 
// current norm to the initial norm (the norm when start is invoked) of the 
// solution vector of the LinearSOE object passed in the constructor
// and a tolerance, set in the constructor.

#ifndef CTestRelativeNormDispIncr_h
#define CTestRelativeNormDispIncr_h

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;


class CTestRelativeNormDispIncr: public ConvergenceTest
{
public:
    // constructors
    CTestRelativeNormDispIncr();	    	
    CTestRelativeNormDispIncr(double tol, int maxNumIter, int printFlag, int normType=2);

    // destructor
    ~CTestRelativeNormDispIncr();
    
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
    double norm0;       // norm at first iteration of each step
};

#endif
