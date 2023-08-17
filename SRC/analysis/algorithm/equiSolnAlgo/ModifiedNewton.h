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
                                                                        
// $Revision: 1.6 $
// $Date: 2007-04-02 23:41:13 $
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
// Newton-Raphson solution algorithm in solving the equations.
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
  ModifiedNewton(int tangent, double iFactor = 0.0, double cFactor = 1.0, int factOnce=0);
  ModifiedNewton(ConvergenceTest &theTest, int tangent = CURRENT_TANGENT, double iFactor = 0.0, double cFactor = 1.0, int factOnce=0);
  ~ModifiedNewton();

    int solveCurrentStep(void);    
    int getNumIterations(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    int tangent;
    int numIterations;
    int factorOnce;

    double iFactor;
    double cFactor;
};

#endif


