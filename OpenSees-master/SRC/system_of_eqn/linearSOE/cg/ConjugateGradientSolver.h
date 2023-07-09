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
// $Date: 2000-09-15 08:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/cg/ConjugateGradientSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/ConjugateGradientSolver.h
//
// Written: fmk 
// Created: 06/00
// Revision: A
//
// Description: This file contains the class definition for 
// ConjugateGradientSolver. ConjugateGradientSolver is an abstract 
// that implements the method solve and which declares a method
// formAp to be pure virtual.
//
// What: "@(#) ConjugateGradientSolver.h, revA"

#ifndef ConjugateGradientSolver_h
#define ConjugateGradientSolver_h

#include <LinearSOESolver.h>
class LinearSOE;
class Vector;

class ConjugateGradientSolver : public LinearSOESolver
{
  public:
    ConjugateGradientSolver(int classTag, LinearSOE *theLinearSOE, double tol);    
    virtual ~ConjugateGradientSolver();

    virtual int setSize(void);    
    virtual int solve(void);
    virtual int formAp(const Vector &p, Vector &Ap) = 0;    
//    virtual int setLinearSOE(LinearSOE &theSOE) =0;

  protected:
    
  private:
    Vector *r, *p, *Ap, *x;
    LinearSOE *theLinearSOE;
    double tolerance;
};

#endif

