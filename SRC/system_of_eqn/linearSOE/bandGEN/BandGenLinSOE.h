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
// $Date: 2000-09-15 08:23:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandGEN/BandGenLinSOE.h,v $
                                                                        
                                                                        
#ifndef BandGenLinSOE_h
#define BandGenLinSOE_h

// File: ~/system_of_eqn/linearSOE/bandGEN/BandGenLinSOE.h
//
// Written: fmk 
// Created: 02/97
// Revision: A
//
// Description: This file contains the class definition for BandGenLinSOE
// BandGenLinSOE is a subclass of LinearSOE. It uses the LAPACK storage
// scheme to store the components of the A matrix, which is a banded 
// unsymmetric matrix.
//
// What: "@(#) BandGenLinSOE.h, revA"


#include <LinearSOE.h>
#include <Vector.h>

class BandGenLinSolver;

class BandGenLinSOE : public LinearSOE
{
  public:
    BandGenLinSOE(BandGenLinSolver &theSolver);    
    BandGenLinSOE(int N, int numSuperDiagonals, int numSubDiagonal,
		  BandGenLinSolver &theSolver);        
    
    ~BandGenLinSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        

    void zeroA(void);
    void zeroB(void);

    const Vector &getX(void);
    const Vector &getB(void);
    double normRHS(void);

    void setX(int loc, double value);    

    int setBandGenSolver(BandGenLinSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class BandGenLinLapackSolver;

  protected:
    
  private:
    int size, numSuperD, numSubD;    
    double *A, *B, *X;
    Vector *vectX;
    Vector *vectB;
    int Asize, Bsize;
    bool factored;
};


#endif

