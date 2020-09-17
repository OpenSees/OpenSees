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
// $Date: 2001-02-17 06:32:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.h
//
// Written: fmk 
// Created: February 1997
// Revision: A
//
// Description: This file contains the class definition for 
// ProfileSPDLinDirectSolver. ProfileSPDLinDirectSolver is a subclass 
// of LinearSOESOlver. It solves a ProfileSPDLinSOE object using
// the LDL^t factorization.

// What: "@(#) ProfileSPDLinDirectSolver.h, revA"

#ifndef ProfileSPDLinDirectSolver_h
#define ProfileSPDLinDirectSolver_h

#include "ProfileSPDLinSolver.h"

class ProfileSPDLinSOE;

class ProfileSPDLinDirectSolver : public ProfileSPDLinSolver
{
  public:
    ProfileSPDLinDirectSolver(double tol=1.0e-12);    
    virtual ~ProfileSPDLinDirectSolver();

    virtual int solve(void);        
    virtual int setSize(void);    
    double getDeterminant(void);

    
    virtual int factor(int n);
    virtual int setProfileSOE(ProfileSPDLinSOE &theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
  protected:
    double minDiagTol;
    int size;
    int *RowTop;
    double **topRowPtr, *invD;
    
  private:

};


#endif

