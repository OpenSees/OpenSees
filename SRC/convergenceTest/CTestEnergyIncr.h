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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestEnergyIncr.h,v $
                                                                        
                                                                        
#ifndef CTestEnergyIncr_h
#define CTestEnergyIncr_h

// File: ~/convergenceTest/CTestEnergyIncr.h
//
// Written: fmk 
// Date: 09/98
// Revised:
//
// Purpose: This file contains the class definition for CTestEnergyIncr.
// A CTestEnergyIncr object tests for convergence using the energy increment,
// which is 0.5 times the absolute value of the product of the rhs and 
// the solution vector of the LinearSOE.

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;

class CTestEnergyIncr: public ConvergenceTest
{
  public:
    // constructors and destructor
    CTestEnergyIncr();	    	
    CTestEnergyIncr(double tol, int maxNumIter, int printFlag);	
    ~CTestEnergyIncr();

    void setTolerance(double newTol);
    int setEquiSolnAlgo(EquiSolnAlgo &theAlgo);
    
    int test(void);
    int start(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    

  protected:

  private:
    LinearSOE *theSOE;
    double tol;

    int maxNumIter;
    int currentIter;
    int printFlag;
};

#endif
