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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/PetscSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/Petsc/PetscSolver.h
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// PetscSolver. It solves the FullGenLinSOE object by calling
// Lapack routines.
//
// What: "@(#) PetscSolver.h, revA"

#ifndef PetscSolver_h
#define PetscSolver_h

#include <LinearSOESolver.h>

// extern "C" {
#include <sles.h>
// }

class PetscSOE;

class PetscSolver : public LinearSOESolver
{
  public:
    PetscSolver(KSPType method, PCType preconditioner);    
    PetscSolver(KSPType method, PCType preconditioner, double rTol, double aTol, int maxIts);    
    ~PetscSolver();

    int solve(void);
    int setSize(void);
    virtual int setLinearSOE(PetscSOE &theSOE);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    

    friend class PetscSolver;
    friend class ActorPetscSOE;
    friend class ShadowPetscSOE;
    
  protected:
    PetscSOE *theSOE;

  private:
    SLES sles;
    KSP ksp;
    PC pc;
    int its;
    KSPType method;
    PCType preconditioner;
    double rTol;
    double aTol;
    int maxIts;
};

#endif

