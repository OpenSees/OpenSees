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
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestNormUnbalance.h,v $
                                                                        
                                                                        
// File: ~/convergenceTest/CTestNormUnbalance.h
//
// Written: fmk 
// Date: 09/98
// Revised:
//
// Purpose: This file contains the class definition for CTestNormUnbalance.
// A CTestNormUnbalance object tests for convergence using the norm of the 
// right hand side vector of the LinearSOE object passed in the constructor
// and a tolerance, set in the constructor

#ifndef CTestNormUnbalance_h
#define CTestNormUnbalance_h

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;


class CTestNormUnbalance: public ConvergenceTest
{
  public:
    // constructors and destructor
    CTestNormUnbalance();	    	
    CTestNormUnbalance(double tol, int maxNumIter, int printFlag);
    ~CTestNormUnbalance();

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
