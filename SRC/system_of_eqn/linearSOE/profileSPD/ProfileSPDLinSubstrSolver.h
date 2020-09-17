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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSubstrSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSubstrSolver.h
//
// Written: fmk 
// Created: February 1997
// Revision: A
//
// Description: This file contains the class definition for 
// ProfileSPDLinSubstrSolver. ProfileSPDLinSubstrSolver is a 
// subclass of DomainSolver and ProfileSPDlinDirectSolver.
// It perform the operations needed of a substr domain decomp method
// on a ProfileSPDlinSOE object. 

// What: "@(#) ProfileSPDLinSubstrSolver.h, revA"

#ifndef ProfileSPDLinSubstrSolver_h
#define ProfileSPDLinSubstrSolver_h

#include <DomainSolver.h>
#include "ProfileSPDLinDirectSolver.h"

class ProfileSPDLinSOE;

class ProfileSPDLinSubstrSolver : public ProfileSPDLinDirectSolver,
                                  public DomainSolver
{
  public:
    ProfileSPDLinSubstrSolver(double tol=1.0e-12);    
    ~ProfileSPDLinSubstrSolver();

    int solve(void);
    int condenseA(int numInt);
    int condenseRHS(int numInt, Vector *v =0);
    int computeCondensedMatVect(int numInt, const Vector &u);    
    const Matrix &getCondensedA(void);
    const Vector &getCondensedRHS(void);
    const Vector &getCondensedMatVect(void);
    
    int setComputedXext(const Vector &);
    int solveXint(void);

    int setSize(void);
    int getClassTag() const;
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
  protected:
    
  private:
    int dSize;
    double *DU;
    Matrix *Aext;
    Vector *Yext;
    
};


#endif

