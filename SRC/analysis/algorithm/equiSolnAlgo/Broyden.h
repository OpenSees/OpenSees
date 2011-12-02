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
// $Date: 2005-11-29 22:42:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/Broyden.h,v $
                                                                        
#ifndef Broyden_h
#define Broyden_h

// File: ~/OOP/analysis/algorithm/Broyden.h 
// 
// Written: Ed C++ Love
// Created: 04/01

// Description: This file contains the class definition for
// Broyden.  Broyden is a class which performs the Broyden
// solution algorihm in solving the equations.
// No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)Broyden.h, revA"

#include <EquiSolnAlgo.h>
#include <Vector.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h> 

class Broyden: public EquiSolnAlgo
{
  public:

    Broyden(int tangent = CURRENT_TANGENT, int n = 10 );    

    Broyden(ConvergenceTest &theTest, int tangent = CURRENT_TANGENT, int n = 10 );

    ~Broyden();

    int solveCurrentStep(void);    

    
    int setConvergenceTest(ConvergenceTest *theNewTest);
    ConvergenceTest *getConvergenceTest(void);     
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:

    ConvergenceTest *theTest;

    ConvergenceTest *localTest ;

    int tangent ;

    int numberLoops ; //number of Broyden iterations

    Vector **s ;  //displacement increments

    Vector **z ;  

    Vector *residOld ;  //residuals
    Vector *residNew ;

    Vector *du ; //displacement increment

    Vector *temp ; //temporary vector

    void BroydenUpdate( IncrementalIntegrator *theIntegrator,
                        LinearSOE *theSOE,
			Vector &du, 
			int count ) ;
  
};

#endif


