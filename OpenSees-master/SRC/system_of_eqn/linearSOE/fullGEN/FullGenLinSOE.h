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
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/fullGEN/FullGenLinSOE.h,v $
                                                                        
                                                                        
#ifndef FullGenLinSOE_h
#define FullGenLinSOE_h

// File: ~/system_of_eqn/linearSOE/fullGEN/FullGenLinSOE.h
//
// Written: fmk 
// Created: 02/97
// Revision: A
//
// Description: This file contains the class definition for FullGenLinSOE
// FullGenLinSOE is a subclass of LinearSOE. It stores all the components
// of the linear system of equation in 1d arrays.
//
// What: "@(#) FullGenLinSOE.h, revA"

#include <LinearSOE.h>
#include <Vector.h>

class FullGenLinSolver;

class FullGenLinSOE : public LinearSOE
{
  public:
    FullGenLinSOE(FullGenLinSolver &theSolver);        
    FullGenLinSOE(int N, FullGenLinSolver &theSolver);        
    FullGenLinSOE();

    ~FullGenLinSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    int addColA(const Vector &col, int colIndex, double fact = 1.0);
    
    void zeroA(void);
    void zeroB(void);
    
    int formAp(const Vector &p, Vector &Ap);

    const Vector &getX(void);
    const Vector &getB(void);    
    const Matrix *getA(void);

    double normRHS(void);

    void setX(int loc, double value);        
    void setX(const Vector &x);        

    int setFullGenSolver(FullGenLinSolver &newSolver);    

    friend class FullGenLinLapackSolver;    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  protected:
    
  private:
    int size;    
    double *A, *B, *X;
    Vector *vectX;
    Vector *vectB;    
    Matrix *matA;
    int Asize, Bsize;
    bool factored;
};


#endif

