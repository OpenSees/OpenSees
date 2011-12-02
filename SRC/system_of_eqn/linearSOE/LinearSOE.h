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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/LinearSOE.h,v $
                                                                        
                                                                        
#ifndef LinearSOE_h
#define LinearSOE_h

// File: ~/system_of_eqn/LinearSOE.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for LinearSOE.
// LinearSOE is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. LinearSystemOfEqn is an abstraction 
// of the linear system of equation given by : [A]{X} = {B} - {C},
// where A is a matrix and B,C and X are vectors. To solve the equation means
// given A, B and C to find the unknown X such that the equation is satisfied.
//
// What: "@(#) LinearSOE.h, revA"

#ifndef _bool_h
#include <bool.h>
#endif

#include <SystemOfEqn.h>

class LinearSOESolver;
class Graph;
class Matrix;
class Vector;
class ID;

class LinearSOE : public SystemOfEqn
{
  public:
    LinearSOE(LinearSOESolver &theSolver, int classTag);    
    virtual ~LinearSOE();

    virtual int solve(void);    

    // pure virtual functions
    virtual int setSize(Graph &theGraph) =0;    
    virtual int getNumEqn(void) const =0;
    
    virtual int addA(const Matrix &, const ID &, double fact = 1.0) =0;
    virtual int addB(const Vector &, const ID &, double fact = 1.0) =0;    
    virtual int setB(const Vector &, double fact = 1.0) =0;        

    virtual void zeroA(void) =0;
    virtual void zeroB(void) =0;

    virtual const Vector &getX(void) = 0;
    virtual const Vector &getB(void) = 0;    
    virtual double normRHS(void) = 0;

    virtual void setX(int loc, double value) =0;
    
  protected:
    int setSolver(LinearSOESolver &newSolver);	        
    LinearSOESolver *getSolver(void);
    
  private:
    LinearSOESolver *theSolver;    
};


#endif

