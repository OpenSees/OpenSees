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
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectThreadSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectThreadSolver.h
//
// Written: fmk 
// Created: Febuary 1997
// Revision: A
//
// Description: This file contains the class definition for 
// ProfileSPDLinDirectThreadSolver. ProfileSPDLinDirectThreadSolver is a subclass 
// of LinearSOESOlver. It solves a ProfileSPDLinSOE object using
// the LDL^t factorization.

// What: "@(#) ProfileSPDLinDirectThreadSolver.h, revA"

#ifndef ProfileSPDLinDirectThreadSolver_h
#define ProfileSPDLinDirectThreadSolver_h

#include <ProfileSPDLinSolver.h>
class ProfileSPDLinSOE;

class ProfileSPDLinDirectThreadSolver : public ProfileSPDLinSolver
{
  public:
    ProfileSPDLinDirectThreadSolver();      
    ProfileSPDLinDirectThreadSolver(int numProcessors, int blockSize, double tol);    
    virtual ~ProfileSPDLinDirectThreadSolver();

    virtual int solve(void);        
    virtual int setSize(void);    

    virtual int setProfileSOE(ProfileSPDLinSOE &theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

  protected:
    int NP;
    int running;
    
    double minDiagTol;
    int blockSize;
    int maxColHeight;
    int size;
    int *RowTop;
    double **topRowPtr, *invD;
    
  private:

};


#endif

