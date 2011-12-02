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
// $Date: 2002-03-02 00:48:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestRelativeEnergyIncr.h,v $
                                                                        
                                                                        
#ifndef CTestRelativeEnergyIncr_h
#define CTestRelativeEnergyIncr_h

// Written: fmk 
// Date: 02/02


// Purpose: This file contains the class definition for CTestRelativeEnergyIncr.
// A CTestRelativeNormEnergyIncr object tests for convergence using the ratio of the 
// current norm to the initial norm (the norm when start is invoked) of the 
// which is 0.5 times the absolute value of the product of the rhs and 
// the solution vector of the LinearSOE.

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;

class CTestRelativeEnergyIncr: public ConvergenceTest
{
  public:
    // constructors and destructor
    CTestRelativeEnergyIncr();	    	
    CTestRelativeEnergyIncr(double tol, int maxNumIter, int printFlag);	
    ~CTestRelativeEnergyIncr();

    ConvergenceTest *getCopy( int iterations ) ;

    void setTolerance(double newTol);
    int setEquiSolnAlgo(EquiSolnAlgo &theAlgo);
    
    int test(void);
    int start(void);
    
    int getNumTests(void);
    int getMaxNumTests(void);        
    double getRatioNumToMax(void);                
    const Vector &getNorms(void);    
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    

  protected:

  private:
    LinearSOE *theSOE;
    double tol;         // the tol on the energy used to test for convergence

    int maxNumIter;     // max number of iterations
    int currentIter;    // number of times test() has been invokes since last start()
    int printFlag;      // a flag indicating if to print on test
    
    Vector norms;       // vector to hold the norms
    double norm0;       // norm at first iteration of each step
};

#endif
