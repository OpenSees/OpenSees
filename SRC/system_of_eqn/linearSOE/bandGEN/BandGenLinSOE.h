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
                                                                        
// $Revision: 1.4 $
// $Date: 2009-05-11 20:57:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandGEN/BandGenLinSOE.h,v $
                                                                        
                                                                        
#ifndef BandGenLinSOE_h
#define BandGenLinSOE_h

// Written: fmk 
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
    BandGenLinSOE();
    BandGenLinSOE(int classTag);
    BandGenLinSOE(BandGenLinSolver &theSolver);

    BandGenLinSOE(int N, int numSuperDiagonals, int numSubDiagonal,
		  BandGenLinSolver &theSolver);        
    
    virtual ~BandGenLinSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);
    
    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addColA(const Vector &col, int colIndex, double fact = 1.0);
    virtual int addB(const Vector &, const ID &, double fact = 1.0);    
    virtual int setB(const Vector &, double fact = 1.0);        

    virtual void zeroA(void);
    virtual void zeroB(void);

    virtual const Vector &getX(void);
    virtual const Vector &getB(void);
    virtual double normRHS(void);

    virtual void setX(int loc, double value);    
    virtual void setX(const Vector &x);    

    virtual int setBandGenSolver(BandGenLinSolver &newSolver);    

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class BandGenLinLapackSolver;

  protected:
    int size, numSuperD, numSubD;    
    double *A, *B, *X;
    Vector *vectX;
    Vector *vectB;
    int Asize, Bsize;
    bool factored;
    
  private:
};


#endif

