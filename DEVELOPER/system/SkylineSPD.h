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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for MyProfileSPD_SOE
// SkylineSPD is a subclass of LinearSOE. It stores symmetric positive
// definite matrix in skyline and solves directly using a cholesky LL^t 
// decomposition.

// What: "@(#) SkylineSPD.h, revA"

#ifndef SkylineSPD_h
#define SkylineSPD_h

#include <LinearSOE.h>
#include <Vector.h>

class SkylineSPD : public LinearSOE
{
  public:
    SkylineSPD();
    ~SkylineSPD();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);

    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);
    
    void zeroA(void);
    void zeroB(void);

    void setX(int loc, double value);
    void setX(const Vector &x);
    
    const Vector &getX(void);
    const Vector &getB(void);

    double normRHS(void);

    int solve(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  protected:

  private:
    int size, profileSize;    
    double *A, *B, *X;
    Vector *vectX;
    Vector *vectB;
    int *iDiagLoc;
    int Asize, Bsize;
    bool isAfactored, isAcondensed;
    int numInt;
    double minDiagTol;    
    int *RowTop;
    double **topRowPtr, *invD;
};


#endif



