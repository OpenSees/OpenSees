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

// $Revision: 1.6 $
// $Date: 2008-04-01 00:35:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp,v $

// Written: fmk 
// Created: 02/06

#include <iostream>

#include <MumpsSolver.h>
#include <MumpsSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

#ifdef _PARALLEL_PROCESSING
#include <mpi.h>
#define OPENSEES_INCLUDED_MPI
#endif

#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#define OPENSEES_INCLUDED_MPI
#endif

#ifndef OPENSEES_INCLUDED_MPI
#include <libseq\mpi.h>
#endif

MumpsSolver::MumpsSolver(int ICNTL7, int ICNTL14)
  :LinearSOESolver(SOLVER_TAGS_MumpsSolver),
   theMumpsSOE(0)
{

  std::cerr << "MumpsSOlver - constructor\n";

  memset(&id, 0, sizeof(id));

  init = false;
  id.job=-1; 
  id.par=1; 

#ifdef _OPENMPI
  //  id.comm_fortran=-987654;
  id.comm_fortran=0;
#else
  id.comm_fortran=MPI_COMM_WORLD;
#endif

  id.ICNTL(14) = ICNTL14;
  id.ICNTL(7)=ICNTL7;;

  needsSetSize = false;
}

MumpsSolver::~MumpsSolver()
{
  std::cerr << "MumpsSOlver - destructor\n";
  opserr << "MumpsParallelSOlver::DESTRUCTOR - start\n";
  id.job=-2; 
  dmumps_c(&id); /* Terminate instance */
  opserr << "MumpsParallelSOlver::DESTRUCTOR - end\n";
}

int 
MumpsSolver::initializeMumps(void)
{
  if (needsSetSize == false) {
    return 0;
  }
  else {
    
    if (init == false) {
      std::cerr << "MumpsSOlver - initMumps\n";
      id.sym = theMumpsSOE->matType;
      id.job=-1; 
      dmumps_c(&id);
      init = true;
    }
    
    int nnz = theMumpsSOE->nnz;
    int *rowA = theMumpsSOE->rowA;
    int *colA = theMumpsSOE->colA;
    
    // increment row and col A values by 1 for mumps fortran indexing
    for (int i = 0; i < nnz; i++) {
      rowA[i]++;
      colA[i]++;
    }
    
    // analyze the matrix
    id.n = theMumpsSOE->size;
    id.nz = theMumpsSOE->nnz;
    id.irn = theMumpsSOE->rowA;
    id.jcn = theMumpsSOE->colA;
    id.a = theMumpsSOE->A;
    id.rhs = theMumpsSOE->X;
    
    // No outputs 
    id.ICNTL(1) = -1; id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 0;

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
    for (int i = 0; i < nnz; i++) {
      rowA[i]--;
      colA[i]--;
    }
    
    needsSetSize = false;
    
    return info;
  }
}

int 
MumpsSolver::solveAfterInitialization(void)
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
	switch(info) {
	  case -5:
		opserr << " out of memory allocation error\n";
      case -6:  
		opserr << " cause: Matrix is Singular in Structure: check your model\n";
	  case -7:
		opserr << " out of memory allocation error\n";
	  case -8:
		opserr << "Work array too small; use -ICNTL14 option, the default is -ICNTL 20 make 20 larger\n";
	  case -9:
		opserr << "Work array too small; use -ICNTL14 option, the default is -ICNTL 20 make 20 larger\n";
	  case -10:  
		opserr << " cause: Matrix is Singular Numerically\n";
	  default:
		  ;
	}
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
MumpsSolver::solve(void)
{
	int initializationResult = initializeMumps();

	if (initializationResult == 0)
		return solveAfterInitialization();
	else
		return initializationResult;
}

int 
MumpsSOE::setMumpsSolver(MumpsSolver &newSolver)
{
  return this->setSolver(newSolver);
}


int
MumpsSolver::setSize()
{
	needsSetSize = true;
	return 0;
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













