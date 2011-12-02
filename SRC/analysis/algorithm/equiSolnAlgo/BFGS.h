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
// $Date: 2003-02-14 23:00:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/BFGS.h,v $
                                                                        
#ifndef BFGS_h
#define BFGS_h

// File: ~/OOP/analysis/algorithm/BFGS.h 
// 
// Written: Ed Love
// Created: 06/01

// Description: This file contains the class definition for
// BFGS.
// 
// What: "@(#)BFGS.h, revA"

#include <EquiSolnAlgo.h>
#include <Vector.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h> 

class BFGS: public EquiSolnAlgo
{
  public:

    BFGS(int tangent = CURRENT_TANGENT, int n = 10);    
    BFGS(ConvergenceTest &theTest, int tangent = CURRENT_TANGENT, int n = 10);
    ~BFGS();

    int solveCurrentStep(void);    

    void setTest(ConvergenceTest &theNewTest);

    ConvergenceTest *getTest(void);     
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:

    ConvergenceTest *theTest;

    ConvergenceTest *localTest;
    int tangent;

    int numberLoops;

    Vector **s;  //displacement increments

    Vector **z;  

    Vector *residOld;  //residuals
    Vector *residNew;

    Vector *du; //displacement increment

    Vector *b;  //current right-hand side

    Vector *temp; //temporary vector 

    double *rdotz;

    double *sdotr;

    void BFGSUpdate(IncrementalIntegrator *theIntegrator,
		    LinearSOE *theSOE,
		    Vector &du, 
		    Vector &b, 
		    int count);
  
};

#endif


