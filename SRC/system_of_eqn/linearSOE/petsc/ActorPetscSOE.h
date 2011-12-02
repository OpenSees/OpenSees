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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/ActorPetscSOE.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/petsc/ActorPetscSOE.h
//
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the class definition for ActorPetscSOE
// ActorPetscSOE is a subclass of LinearSOE. It uses the LAPACK storage
// scheme to store the components of the A matrix, which is a full matrix.


// What: "@(#) ActorPetscSOE.h, revA"

#ifndef ActorPetscSOE_h
#define ActorPetscSOE_h

#include <LinearSOE.h>
#include <Vector.h>

// extern "C" {
#include <petsc.h>
// }

class PetscSolver;
class PetscSOE;

class ActorPetscSOE
{
  public:
    ActorPetscSOE(PetscSolver &theSolver, int blockSize);     
    
    ~ActorPetscSOE();
    int run(void);
    
  protected:
    
  private:
    MPI_Comm theComm;
    PetscSOE *theSOE;  // the local portion of the SOE
    PetscSolver *theSolver; // created locally via data from process 0
    int myRank;				 
    int recvData[3];				 
    void *recvBuffer;
    int numProcessors;
};


#endif



