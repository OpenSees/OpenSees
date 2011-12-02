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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/ModifiedNewton.h,v $
                                                                        
                                                                        
#ifndef ModifiedNewton_h
#define ModifiedNewton_h

// File: ~/OOP/analysis/algorithm/ModifiedNewton.h 
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//
// Description: This file contains the class definition for 
// ModifiedNewton. ModifiedNewton is a class which performs a modified 
// Newton-Raphson  solution algorihm in solving the equations.
// No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)ModifiedNewton.h, revA"

#include <EquiSolnAlgo.h>
#include <Vector.h>

class ConvergenceTest;

class ModifiedNewton: public EquiSolnAlgo
{
  public:
    ModifiedNewton();    
    ModifiedNewton(ConvergenceTest &theTest);
    ~ModifiedNewton();

    int solveCurrentStep(void);    
    void setTest(ConvergenceTest &theNewTest);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(ostream &s, int flag =0);    
    
  protected:
    
  private:
    ConvergenceTest *theTest;
};

#endif


