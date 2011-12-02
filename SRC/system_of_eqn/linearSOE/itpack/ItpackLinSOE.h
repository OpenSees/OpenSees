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
                                                                        
// $Revision: 1.2 $
// $Date: 2002-06-08 16:17:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/itpack/ItpackLinSOE.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/umfGEN/ItpackGenLinSOE.h
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for 
// ItpackGenLinSOE. It stores the sparse matrix A in a fashion
// required by the ItpackLinSolver object.
//
// What: "@(#) ItpackGenLinSOE.h, revA"


#ifndef ItpackGenLinSOE_h
#define ItpackGenLinSOE_h

#include <LinearSOE.h>
#include <Vector.h>

class ItpackLinSolver;

class ItpackLinSOE : public LinearSOE
{
  public:
    ItpackLinSOE(ItpackLinSolver &theSolver);        
    ItpackLinSOE(int N, int NNZ, int *rowStartA, int *colA,
		    ItpackLinSolver &theSolver);        

    ~ItpackLinSOE();

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
    void setX(const Vector &x);
    int setItpackLinSolver(ItpackLinSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker); 

    friend class ItpackLinSolver;

  protected:
    
  private:
    int size;            // order of A
    int nnz;             // number of non-zeros in A
    double *A, *B, *X;   // 1d arrays containing coefficients of A, B and X
    int *colA, *rowStartA; // int arrays containing info about coeff's in A
    
    Vector *vectX;
    Vector *vectB;    
    int Asize, Bsize;    // size of the 1d array holding A

    bool Aformed;        // false when zeroA() is called
};

#endif
