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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/ActorPetscSOE.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/petsc/ActorPetscSOE.C
//
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the implementation for BandGenLinSOE


#include "ActorPetscSOE.h"
#include "PetscSolver.h"
#include "PetscSOE.h"
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

ActorPetscSOE::ActorPetscSOE(PetscSolver &theSOESolver, int blockSize)
  :myRank(0), theSOE(0), theSolver(&theSOESolver), recvBuffer(0)
{

  MPI_Comm_rank(PETSC_COMM_WORLD, &myRank);
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcessors);
  //  MPI_Comm_dup(PETSC_COMM_WORLD, &theComm);
  if (myRank == 0) {
    opserr << " ActorPetscSOE::ActorPetscSOE - must be rank 0\n";
  }
  recvBuffer = (void *)(&recvData[0]);
  MPI_Barrier(PETSC_COMM_WORLD);

/***************
  // recv data on the PetscSolver & determine KSP and PC type
  MPI_Bcast(recvBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);    

  int method = recvData[0];
  int preconditioner = recvData[1];

  KSPType petscMethod =0;
  PCType petscPre =0;
  if (method == 1)
    petscMethod = KSPCG;
  else if (method == 2)
    petscMethod = KSPGMRES;
  else if (method == 3)
    petscMethod = KSPBCGS;
  else if (method == 4)
    petscMethod = KSPCGS;
  else if (method == 5)
    petscMethod = KSPTFQMR;
  else if (method == 6)
    petscMethod = KSPTCQMR;
  else if (method == 7)
    petscMethod = KSPCR;
  else if (method == 8)
    petscMethod = KSPLSQR;
  else if (method == 9)
    petscMethod = KSPRICHARDSON;
  else if (method == 10)
    petscMethod = KSPCHEBYCHEV;
  else if (method == 11)
    petscMethod = KSPPREONLY;
  else {
    opserr << "ActorPetscSOE::ActorPetscSOE - unknown KSP method\n";
  }

  if (preconditioner == 1)
    petscPre = PCNONE;
  else if (preconditioner == 2)
    petscPre = PCJACOBI;  
  else if (preconditioner == 3)
    petscPre = PCSOR;  
  else if (preconditioner == 4)
    petscPre = PCEISENSTAT;  
  else if (preconditioner == 5)
    petscPre = PCBJACOBI;  
  else if (preconditioner == 6)
    petscPre = PCASM;  
  else if (preconditioner == 7)
    petscPre = PCILU;  
  else if (preconditioner == 8)
    petscPre = PCICC;  
  else if (preconditioner == 9)
    petscPre = PCBGS;  
  else if (preconditioner == 10)
    petscPre = PCMG;  
  else if (preconditioner == 11)
    petscPre = PCSHELL;  
  else if (preconditioner == 12)
    petscPre = PCLU;  
  else {
    opserr << "ActorPetscSOE::ActorPetscSOE - unknown PC method\n";
  }

  // construct the PetscSolver and PetscSOE
  if ((petscMethod != 0) && (petscPre != 0)) {
    theSolver = new PetscSolver(petscMethod, petscPre);
    theSOE = new PetscSOE(*theSolver);
  }
  *********/

  theSOE = new PetscSOE(*theSolver, blockSize);
}


    
ActorPetscSOE::~ActorPetscSOE()
{
  //  if (theSOE != 0)
  //    delete theSOE;
}


int 
ActorPetscSOE::run(void)
{
  int flag = 1;
  int tag, low, high, ierr, n, numRows;
  double *theData;
  int dnz = 0;
  int onz = 0;
  int *dnnz = PETSC_NULL;
  int *onnz = PETSC_NULL;
  MPI_Status status;
  void *buffer = 0;
  
  while (flag != 0) {
    MPI_Bcast(recvBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);  
    flag = recvData[0];
    switch(flag) {

    case 0:
      MPI_Barrier(PETSC_COMM_WORLD);
      break;

    case 1:
      if (theSOE == 0)
	return -1;
      theSOE->isFactored = recvData[1];
      theSOE->solve();
      break;

    case 2:
      if (theSOE == 0)
	return -1;
      tag = 100;
      MPI_Recv(recvBuffer, 3, MPI_INT, 0, tag, PETSC_COMM_WORLD, &status);

      numRows = recvData[1];
      n = recvData[2];
      onnz = new int[numRows];
      dnnz = new int[numRows];
      buffer = (void *)dnnz;
      tag = 101;
      MPI_Recv(buffer, numRows, MPI_INT, 0, tag, PETSC_COMM_WORLD, &status);

      buffer = (void *)onnz;
      tag = 102;
      MPI_Recv(buffer, numRows, MPI_INT, 0, tag, PETSC_COMM_WORLD, &status);

      MPI_Bcast(recvBuffer, 3, MPI_INT, 0, PETSC_COMM_WORLD);  
      theSOE->setSizeParallel(numRows, n, dnz, dnnz, onz, onnz);	
      break;

    case 3:
      if (theSOE == 0)
	return -1;
      theSOE->zeroA();
      break;

    case 4:
      if (theSOE == 0)
	return -1;
      theSOE->zeroB();
      break;

    case 5:
      if (theSOE == 0)
	return -1;
      tag = 99;
      ierr = VecGetOwnershipRange(theSOE->x, &low, &high); CHKERRA(ierr);
      recvData[0] = low; recvData[1] = high;
      MPI_Send(recvBuffer, 2, MPI_INT, 0, tag, PETSC_COMM_WORLD);
      ierr = VecGetArray(theSOE->x, &theData); CHKERRA(ierr);       
      MPI_Send(theData, high-low, MPI_DOUBLE, 0, tag, PETSC_COMM_WORLD); 
      ierr = VecRestoreArray(theSOE->x, &theData); CHKERRA(ierr);             
      break;

    case 6:
      if (theSOE == 0)
	return -1;
      // some work here

    default:
      opserr << "ActorPetscSOE::invalid action " << flag << " received\n";
    }
  }
  return 0;
}




