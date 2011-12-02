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
// $Date: 2005-05-18 19:26:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/PetscSparseSeqSolver.cpp,v $
                                                                        
// Written: fmk 
// Created: 04/05
                                                                        
// What: "@(#) PetscSparseSeqSolver.h, revA"

#include <PetscSparseSeqSolver.h>
#include <SparseGenRowLinSOE.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>

PetscSparseSeqSolver::PetscSparseSeqSolver(KSPType meth, PCType pre)
  :SparseGenRowLinSolver(SOLVER_TAGS_PetscSparseSeqSolver), 
   rTol(PETSC_DEFAULT), aTol(PETSC_DEFAULT), dTol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT)
{
  PetscInitialize(0, PETSC_NULL, (char *)0, PETSC_NULL);

  method = meth;
  preconditioner = pre;
}

PetscSparseSeqSolver::PetscSparseSeqSolver(KSPType meth, PCType pre, double relTol, double absTol, double divTol, int maxIterations)
  :SparseGenRowLinSolver(SOLVER_TAGS_PetscSparseSeqSolver), 
   rTol(relTol), aTol(absTol), dTol(divTol), maxIts(maxIterations)
{
  PetscInitialize(0, PETSC_NULL, (char *)0, PETSC_NULL);

  method = meth;
  preconditioner = pre;
}


PetscSparseSeqSolver::~PetscSparseSeqSolver()
{
  MatDestroy(A);
  VecDestroy(x);
  VecDestroy(b);
  KSPDestroy(ksp);
}



int
PetscSparseSeqSolver::solve(void)
{
  PetscErrorCode ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr); 

  return ierr;
}


int
PetscSparseSeqSolver::getNumIterations(void)
{
  PetscErrorCode ierr = KSPGetIterationNumber(ksp, &its);
    
  return its;
}


int
PetscSparseSeqSolver::setSize()
{
  int n = theSOE->size;
  int nnz = theSOE->nnz;
  double *Adata = theSOE->A;
  double *Xdata = theSOE->X;
  double *Bdata = theSOE->B;
  int *colA = theSOE->colA;
  int *rowStartA = theSOE->rowStartA;


  /* 
   * Call Petsc VecCreate & MatCreate; NOTE: using previously allocated storage
   * 
   */

  PetscTruth flg;
  PetscErrorCode ierr = PetscOptionsGetInt(PETSC_NULL,"-n", &n, &flg); CHKERRQ(ierr);
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, n, Xdata, &x); CHKERRQ(ierr); 
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, n, Bdata, &b); CHKERRQ(ierr); 
  ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, n, n, rowStartA, colA,  Adata , &A);

  /* 
   * Create linear solver context
   */

   KSPCreate(PETSC_COMM_WORLD, &ksp);

   /* 
    *  Set operators. NOTE: matrix that defines the linear system
    *  also serves as the preconditioning matrix.
    */

   KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);

   /* 
    *  Set solution scheme & tolerances
    */

   ierr = KSPSetType(ksp, method); CHKERRQ(ierr); 
   ierr = KSPSetTolerances(ksp, rTol, aTol, dTol, maxIts); 

   /* 
    *  Set preconditioning scheme
    */

   KSPGetPC(ksp, &pc);
   ierr = PCSetType(pc,  preconditioner); CHKERRQ(ierr); 

   /* 
    *  Finally mark so that uses last solution as initial guess in each solution
    *    NOTE: maybe change this as another user supplied option
    */

   ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr); 

   return ierr;
}


int 
PetscSparseSeqSolver::setLinearSOE(SparseGenRowLinSOE &theSystem)
{
  theSOE = &theSystem;
  return 0;
}


int
PetscSparseSeqSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PetscSparseSeqSolver::recvSelf(int cTag, Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}



