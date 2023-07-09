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
// $Date: 2000-09-15 08:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectBlockSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectBlockSolver.h
//
// Written: fmk 
// Created: February 1997
// Revision: A
//
// Description: This file contains the class definition for 
// ProfileSPDLinDirectBlockSolver. ProfileSPDLinDirectBlockSolver is a subclass 
// of LinearSOESOlver. It solves a ProfileSPDLinSOE object using
// the LDL^t factorization.

// What: "@(#) ProfileSPDLinDirectBlockSolver.h, revA"

#ifndef ProfileSPDLinDirectBlockSolver_h
#define ProfileSPDLinDirectBlockSolver_h

#include <ProfileSPDLinSolver.h>
class ProfileSPDLinSOE;

class ProfileSPDLinDirectBlockSolver : public ProfileSPDLinSolver
{
  public:
    ProfileSPDLinDirectBlockSolver(double tol=1.0e-12, int blockSize = 4);    
    virtual ~ProfileSPDLinDirectBlockSolver();

    virtual int solve(void);        
    virtual int setSize(void);    

    virtual int setProfileSOE(ProfileSPDLinSOE &theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    
  protected:
    double minDiagTol;
    int blockSize;
    int maxColHeight;
    int size;
    int *RowTop;
    double **topRowPtr, *invD;
    
  private:

};


#endif

