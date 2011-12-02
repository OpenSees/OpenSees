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
// $Date: 2005-06-24 18:19:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/PetscSolver.cpp,v $
                                                                        
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
#include <ID.h>
#include <Vector.h>

PetscSolver::PetscSolver()
  :LinearSOESolver(SOLVER_TAGS_PetscSolver), 
   rTol(PETSC_DEFAULT), aTol(PETSC_DEFAULT), dTol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT)
{

}

PetscSolver::PetscSolver(KSPType meth, PCType pre)
  :LinearSOESolver(SOLVER_TAGS_PetscSolver), method(meth), preconditioner(pre),
   rTol(PETSC_DEFAULT), aTol(PETSC_DEFAULT), dTol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT)
{

}

PetscSolver::PetscSolver(KSPType meth, PCType pre, double relTol, double absTol, double divTol, int maxIterations)
  :LinearSOESolver(SOLVER_TAGS_PetscSolver), method(meth), preconditioner(pre),
   rTol(relTol), aTol(absTol), dTol(divTol), maxIts(maxIterations)
{

}

PetscSolver::~PetscSolver()
{
  
}



int
PetscSolver::solve(void)
{
  int size = theSOE->size;
  int factored = theSOE->isFactored;
  
  int numProcesses = theSOE->numProcesses;
  int processID = theSOE->processID;
  
  int ierr;
  ierr = MatAssemblyBegin(theSOE->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
  ierr = MatAssemblyEnd(theSOE->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 

  //
  // if parallel, we must zero X & form the total B: each processor has own contributions
  //
  
  static Vector recvVector(1);

  if (numProcesses > 1) {
    Vector *vectX = theSOE->vectX;
    Vector *vectB = theSOE->vectB;

    // zero X
    vectX->Zero();

    //
    // form B on each
    //

    int numChannels = theSOE->numChannels;
    Channel **theChannels = theSOE->theChannels;
    
    if (processID != 0) {
      Channel *theChannel = theChannels[0];
      
      theChannel->sendVector(0, 0, *vectB);
      theChannel->recvVector(0, 0, *vectB);
      
    } else {
      
      if (recvVector.Size() != size) 
	recvVector.resize(size);
      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->recvVector(0, 0, recvVector);
	*vectB += recvVector;
      }
      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->sendVector(0, 0, *vectB);
      }
    }
  }

  //
  // solve and mark as having been solved
  //

  ierr = KSPSolve(ksp, theSOE->b, theSOE->x); CHKERRQ(ierr); 
  theSOE->isFactored = 1;

  //
  // if parallel, we must form the total X: each processor has startRow through endRow-1
  //

  if (numProcesses > 1) {
    Vector *vectX = theSOE->vectX;

    int numChannels = theSOE->numChannels;
    Channel **theChannels = theSOE->theChannels;
    
    if (processID != 0) {
      Channel *theChannel = theChannels[0];
      
      theChannel->sendVector(0, 0, *vectX);
      theChannel->recvVector(0, 0, *vectX);
      
    } else {
      
      if (recvVector.Size() != size) 
	recvVector.resize(size);

      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->recvVector(0, 0, recvVector);
	*vectX += recvVector;
      }
      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->sendVector(0, 0, *vectX);
      }
    }
  }

  return ierr;
}


int
PetscSolver::setSize()
{
  /* 
   * Create linear solver context
   */
  
   KSPCreate(PETSC_COMM_WORLD, &ksp);

   /* 
    *  Set operators. NOTE: matrix that defines the linear system
    *  also serves as the preconditioning matrix.
    */

   KSPSetOperators(ksp, theSOE->A, theSOE->A, DIFFERENT_NONZERO_PATTERN);

   /* 
    *  Set solution scheme & tolerances
    */
   
   int ierr;
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

   //   ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr); 
   return ierr;
}


int 
PetscSolver::setLinearSOE(PetscSOE &theSys)
{
  theSOE = &theSys;
  return 0;
}


int
PetscSolver::sendSelf(int cTag, Channel &theChannel)
{
  static ID idData(3);
  idData(0) = maxIts;
  if (strcmp(method, KSPCG) == 0) 
    idData(1) = 0;
  else if (strcmp(method, KSPBICG) == 0)
    idData(1) = 1;
  else if (strcmp(method, KSPRICHARDSON) == 0)
    idData(1) = 2;
  else if (strcmp(method, KSPCHEBYCHEV) == 0)
    idData(1) = 3;
  else if (strcmp(method, KSPGMRES) ==0) 
    idData(1) = 4;
  else {
    opserr << "PetscSolver::sendSelf() - unknown method set\n";
    return -1;
  }

  if (strcmp(preconditioner, PCJACOBI) == 0)
    idData(2) = 0;
  else if (strcmp(preconditioner, PCILU) == 0)
    idData(2) = 1;
  else if (strcmp(preconditioner, PCICC) == 0)
    idData(2) = 2;
  else if (strcmp(preconditioner, PCBJACOBI) == 0)
    idData(2) = 3;
  else if (strcmp(preconditioner, PCNONE) == 0)
    idData(2) = 4;
  else {
    opserr << "PetscSolver::sendSelf() - unknown preconditioner set\n";
    return -1;
  }

  theChannel.sendID(0, cTag, idData);

  static Vector data(3);
  data(0) = rTol;
  data(1) = aTol;
  data(2) = dTol;

  theChannel.sendVector(0, cTag, data);
  return 0;
}

int
PetscSolver::recvSelf(int cTag, Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  static ID idData(3);
  theChannel.recvID(0, cTag, idData);
  maxIts = idData(0);
  if (idData(1) == 0) 
    method = KSPCG;
  else if (idData(1) == 1)
    method = KSPBICG;
  else if (idData(1) == 2) 
    method = KSPRICHARDSON;
  else if (idData(1) == 3)
    method = KSPCHEBYCHEV;
  else if (idData(1) == 4)
    method = KSPGMRES;
  else {
    opserr << "PetscSolver::recvSelf() - unknown method recvd\n";
    return -1;
  }

  if (idData(2) == 0)
    preconditioner = PCJACOBI;
  else if (idData(2) == 1)
    preconditioner = PCILU;
  else if (idData(2) == 2)
    preconditioner = PCICC;
  else if (idData(2) == 3)
    preconditioner = PCBJACOBI;
  else if (idData(2) == 4)
    preconditioner = PCNONE;
  else {
    opserr << "PetscSolver::sendSelf() - unknown preconditioner set\n";
    return -1;
  }


  static Vector data(3);
  theChannel.recvVector(0, cTag, data);
  rTol = data(0);
  aTol = data(1);
  dTol = data(2);

  return 0;
}



