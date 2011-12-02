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
// $Date: 2000-09-15 08:23:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/DomainSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/DomainSolver.h
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for DomainSolver.
// DomainSolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of DomainSolver are used
// to solve a system of equations and to do static condensation operations on the
// linear system of equations.
//
// What: "@(#) DomainSolver.h, revA"

#ifndef DomainSolver_h
#define DomainSolver_h

#include <LinearSOESolver.h>
class Matrix;
class Vector;

class DomainSolver: public LinearSOESolver
{
  public:
    DomainSolver(int classTag);    
    virtual ~DomainSolver();

    virtual int condenseA(int numInt) =0;
    virtual int condenseRHS(int numInt, Vector *u =0) =0;
    virtual int computeCondensedMatVect(int numInt, const Vector &u) =0;
    virtual const Matrix &getCondensedA(void) =0;
    virtual const Vector &getCondensedRHS(void) =0;
    virtual const Vector &getCondensedMatVect(void) =0;

    virtual int setComputedXext(const Vector &) =0;
    virtual int solveXint(void) =0;
    
  protected:
    
  private:
};

#endif

