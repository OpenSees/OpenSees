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
// $Date: 2001-12-07 00:17:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinSOE.h,v $
                                                                        
                                                                        
#ifndef BandSPDLinSOE_h
#define BandSPDLinSOE_h

// File: ~/system_of_eqn/linearSOE/bandSPD/BandSPDLinSOE.h
//
// Written: fmk 
// Created: Febuary 1997
// Revision: A
//
// Description: This file contains the class definition for BandSPDLinSOE
// BandSPDLinSOE is a subclass of LinearSOE. It uses the LAPACK Upper storage
// scheme to store the components of the A matrix.
//
// What: "@(#) BandSPDLinSOE.h, revA"



#include <LinearSOE.h>
#include <Vector.h>

class BandSPDLinSolver;



class BandSPDLinSOE : public LinearSOE
{
  public:
    BandSPDLinSOE(BandSPDLinSolver &theSolver);    
    BandSPDLinSOE(int N, int bandwidth, BandSPDLinSolver &theSolver);        

    ~BandSPDLinSOE();

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
    int setBandSPDSolver(BandSPDLinSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    friend class BandSPDLinSolver;
    friend class BandSPDLinLapackSolver;    
    friend class BandSPDLinThreadSolver;        
    
  protected:
    
  private:
    int size, half_band;    
    double *A, *B, *X;
    Vector *vectX;
    Vector *vectB;    
    int Asize, Bsize;
    int aFactored;
    bool factored;
};


#endif

