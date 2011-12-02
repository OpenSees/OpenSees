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
                                                                        
// Description: This file contains the implementation of MumpsSolver

// $Revision: 1.3 $
// $Date: 2006-04-04 22:59:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp,v $

// Written: fmk 
// Created: 02/06

#include <MumpsSolver.h>
#include <MumpsSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

extern "C" {
#include <mpi.h>
}

MumpsSolver::MumpsSolver()
  :LinearSOESolver(SOLVER_TAGS_MumpsSolver),
   theMumpsSOE(0)
{
  init = false;
  id.job=-1; 
  id.par=1; 

  id.comm_fortran=MPI_COMM_WORLD;
}


MumpsSolver::MumpsSolver(int ICNTL7)
  :LinearSOESolver(SOLVER_TAGS_MumpsSolver),
   theMumpsSOE(0)
{
  init = false;
  id.job=-1; 
  id.par=1; 

  id.comm_fortran=MPI_COMM_WORLD;
  id.ICNTL(7)=ICNTL7;;
}

MumpsSolver::~MumpsSolver()
{
  id.job=-2; 
  dmumps_c(&id); /* Terminate instance */
}

int
MumpsSolver::solve(void)
{
  int nnz = theMumpsSOE->nnz;
  int n = theMumpsSOE->size;
  int *rowA = theMumpsSOE->rowA;
  int *colA = theMumpsSOE->colA;

  double *X = theMumpsSOE->X;
  double *B = theMumpsSOE->B;

  // increment row and col A values by 1 for mumps fortran indexing
  for (int i=0; i<nnz; i++) {
    rowA[i]++;
    colA[i]++;
    //    opserr << rowA[i] << " " << colA[i] << " " << theMumpsSOE->A[i] << endln;
  }
  
  for (int i=0; i<n; i++)
    X[i] = B[i];

  int info = 0;
  if (theMumpsSOE->factored == false) {

    // factor the matrix
    id.n   = theMumpsSOE->size; 
    id.nz  = theMumpsSOE->nnz; 
    id.irn = theMumpsSOE->rowA;
    id.jcn = theMumpsSOE->colA;
    id.a   = theMumpsSOE->A; 
    id.rhs = theMumpsSOE->X;

    // No outputs 
    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
    // Call the MUMPS package to factor & solve the system
    id.job = 5;
    dmumps_c(&id);

    theMumpsSOE->factored = true;
  } else {
    // factor the matrix
    id.n   = theMumpsSOE->size; 
    id.nz  = theMumpsSOE->nnz; 
    id.irn = theMumpsSOE->rowA;
    id.jcn = theMumpsSOE->colA;
    id.a   = theMumpsSOE->A; 
    id.rhs = theMumpsSOE->X;

    // No outputs 
    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
    // Call the MUMPS package to factor & solve the system
    id.job = 3;
    dmumps_c(&id);
  }	


  
  info = id.infog[0];
  if (info != 0) {	
    opserr << "WARNING MumpsSolver::solve(void)- ";
    opserr << " Error " << info << " returned in substitution dmumps()\n";
    return info;
  }

  // decrement row and col A values by 1 to return to C++ indexing
  for (int i=0; i<nnz; i++) {
    rowA[i]--;
    colA[i]--;
  }

  return 0;
}


int 
MumpsSOE::setMumpsSolver(MumpsSolver &newSolver)
{
  return this->setSolver(newSolver);
}


int
MumpsSolver::setSize()
{
  if (init == false) {
    id.sym=theMumpsSOE->matType; 
    dmumps_c(&id);
    init = true;
  }

  int nnz = theMumpsSOE->nnz;
  int *rowA = theMumpsSOE->rowA;
  int *colA = theMumpsSOE->colA;

  // increment row and col A values by 1 for mumps fortran indexing
  for (int i=0; i<nnz; i++) {
    rowA[i]++;
    colA[i]++;
  }

  // analyze the matrix
  id.n   = theMumpsSOE->size; 
  id.nz  = theMumpsSOE->nnz; 
  id.irn = theMumpsSOE->rowA;
  id.jcn = theMumpsSOE->colA;
  id.a   = theMumpsSOE->A; 
  id.rhs = theMumpsSOE->X;
  
  // No outputs 
  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
  // Call the MUMPS package to factor & solve the system
  id.job = 1;
  dmumps_c(&id);

  int info = id.infog[0];
  if (info != 0) {	
    opserr << "WARNING MumpsSolver::setSize(void)- ";
    opserr << " Error " << info << " returned in substitution dmumps()\n";
    return info;
  }

  // decrement row and col A values by 1 to return to C++ indexing
  for (int i=0; i<nnz; i++) {
    rowA[i]--;
    colA[i]--;
  }
  return info;
}

int
MumpsSolver::sendSelf(int cTag, Channel &theChannel)
{
  // nothing to do
  return 0;
}

int
MumpsSolver::recvSelf(int ctag,
		      Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  // nothing to do
  return 0;
}

int 
MumpsSolver::setLinearSOE(MumpsSOE &theSOE)
{
  theMumpsSOE = &theSOE;
  return 0;
}













