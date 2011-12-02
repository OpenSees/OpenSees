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
// $Date: 2003-02-14 23:02:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/PetscSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/Petsc/PetscSolver.C
//
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the class implementation for 
// FullGenLinLapackSolver. It solves the FullGenLinSOE object by calling
// Lapack routines.
//
// What: "@(#) FullGenLinLapackSolver.h, revA"

#include <PetscSolver.h>
#include <PetscSOE.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>


PetscSolver::PetscSolver(KSPType meth, PCType pre)
  :LinearSOESolver(SOLVER_TAGS_PetscSolver), method(meth), preconditioner(pre),
   rTol(0.0), aTol(0.0), maxIts(0)
{
  int ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRA(ierr);
}

PetscSolver::PetscSolver(KSPType meth, PCType pre, double relTol, double absTol, int maxIterations)
  :LinearSOESolver(SOLVER_TAGS_PetscSolver), method(meth), preconditioner(pre),
   rTol(relTol), aTol(absTol), maxIts(maxIterations)
{
  int ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRA(ierr);
}

PetscSolver::~PetscSolver()
{

}



int
PetscSolver::solve(void)
{
 Mat A = theSOE->A;
 Vec b = theSOE->b;
 Vec x = theSOE->x;
 
 int factored = theSOE->isFactored;

 if (factored == 0) {
   int ierr = SLESSetOperators(sles,A,A,DIFFERENT_NONZERO_PATTERN); 
   CHKERRA(ierr);
   
   ierr = SLESGetKSP(sles,&ksp); CHKERRA(ierr); 
   ierr = SLESGetPC(sles,&pc); CHKERRA(ierr); 

   ierr = KSPSetType(ksp,method); CHKERRA(ierr); 
   ierr = PCSetType(pc,preconditioner); CHKERRA(ierr); 
   if (maxIts != 0) 
     ierr = KSPSetTolerances(ksp,rTol,aTol,PETSC_DEFAULT, maxIts); 
   else
     ierr = KSPSetTolerances(ksp, 1.0e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); 
     
   CHKERRA(ierr); 

   ierr = SLESSetFromOptions(sles); CHKERRA(ierr);

   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr); 
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr); 
 }
 int ierr;
 ierr = VecAssemblyBegin(b); CHKERRA(ierr); 
 ierr = VecAssemblyEnd(b); CHKERRA(ierr); 

 Timer timer;
 timer.start();

 // solve and mark as having been solved
 ierr = SLESSolve(sles,b,x,&its); CHKERRA(ierr); 
 theSOE->isFactored = 1;

 timer.pause();

 opserr << "PetscSolver::solve -real  " << timer.getReal() << "  cpu: " << timer.getCPU() << "  page: " << timer.getNumPageFaults() << endln; 

 opserr << "PetscSolver::solve(void) - number of iterations: " << its << endln;

 return ierr;
}


int
PetscSolver::setSize()
{
    return 0;
}


int 
PetscSolver::setLinearSOE(PetscSOE &theSys)
{
  theSOE = &theSys;
  return 0;
}



PetscSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PetscSolver::recvSelf(int cTag, Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}



