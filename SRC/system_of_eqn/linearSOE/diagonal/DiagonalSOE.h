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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-01-27 22:22:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DiagonalSOE.h,v $
                                                                        
// Written: fmk 
// Created: Jan 2005
// Revision: A
//
// Description: This file contains the class definition for DiagonalSOE
// DiagonalSOE is a subclass of LinearSOE. It stores a diagonal system
// of equation, i.e. just the diagonal

// What: "@(#) DiagonalSOE.h, revA"

#ifndef DiagonalSOE_h
#define DiagonalSOE_h

#include <LinearSOE.h>
#include <Vector.h>
class DiagonalSolver;

class DiagonalSOE : public LinearSOE
{
  public:
    DiagonalSOE(DiagonalSolver &theSolver);
    DiagonalSOE(int N, DiagonalSolver &theSolver);

    ~DiagonalSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);

    int formAp(const Vector &p, Vector &Ap);
    
    void setX(int loc, double value);
    void setX(const Vector &x);

    const Vector &getX(void);
    const Vector &getB(void);
    double normRHS(void);

    int setDiagonalSolver(DiagonalSolver &newSolver);    
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    friend class DiagonalSolver;    
    friend class DiagonalDirectSolver;
    
  protected:
    
  private:
    int size;
    double *A, *B, *X;
    Vector *vectX;
    Vector *vectB;
    bool isAfactored;
};


#endif



