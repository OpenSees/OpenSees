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
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestNormDispIncr.h,v $
                                                                        
                                                                        
#ifndef CTestNormDispIncr_h
#define CTestNormDispIncr_h

// File: ~/convergenceTest/CTestNormDispIncr.h
//
// Written: fmk 
// Date: 09/98
// Revised:
//
// Purpose: This file contains the class definition for CTestNormDispIncr.
// A CTestNormDispIncr object tests for convergence using the norm of the 
// solution vector of the LinearSOE object passed in the constructor
// and a tolerance, set in the constructor

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;


class CTestNormDispIncr: public ConvergenceTest
{
  public:
    // constructors and destructor
    CTestNormDispIncr();	    	
    CTestNormDispIncr(double tol, int maxNumIter, int printFlag);	
    ~CTestNormDispIncr();

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
    int totalNumIter;
};

#endif
