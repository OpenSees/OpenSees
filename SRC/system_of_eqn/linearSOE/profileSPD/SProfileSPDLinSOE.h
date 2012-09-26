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
// $Date: 2009-05-11 20:58:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.h,v $
                                                                        
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for ProfileSPDLinSOE
// ProfileSPDLinSOE is a subclass of LinearSOE. It uses the LAPACK Upper storage
// scheme to store the components of the A matrix.

// What: "@(#) ProfileSPDLinSOE.h, revA"

#ifndef SProfileSPDLinSOE_h
#define SProfileSPDLinSOE_h

#include <LinearSOE.h>
#include <Vector.h>
class SProfileSPDLinSolver;

class SProfileSPDLinSOE : public LinearSOE
{
  public:
    SProfileSPDLinSOE(SProfileSPDLinSolver &theSolver);
    SProfileSPDLinSOE(SProfileSPDLinSolver &theSolver, int classTag);
    SProfileSPDLinSOE(int classTag);
    SProfileSPDLinSOE(int N, int *iLoc, SProfileSPDLinSolver &theSolver);

    virtual ~SProfileSPDLinSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);
    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addB(const Vector &, const ID &, double fact = 1.0);    
    virtual int setB(const Vector &, double fact = 1.0);
    
    virtual void zeroA(void);
    virtual void zeroB(void);

    virtual void setX(int loc, double value);
    virtual void setX(const Vector &x);
    
    virtual const Vector &getX(void);
    virtual const Vector &getB(void);
    virtual double normRHS(void);

    virtual int setProfileSPDSolver(SProfileSPDLinSolver &newSolver);    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    friend class SProfileSPDLinSolver;    
    
  protected:
    int size, profileSize;    
    float *A, *B, *X;
    double *doubleB, *doubleX;
    Vector *vectX;
    Vector *vectB;
    int *iDiagLoc;
    int Asize, Bsize;
    bool isAfactored, isAcondensed;
    int numInt;
    
  private:
};


#endif



