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
// $Date: 2007-02-14 20:12:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/ShadowPetscSOE.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/petsc/ShadowPetscSOE.h
//
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the class definition for ShadowPetscSOE
// ShadowPetscSOE is a subclass of LinearSOE. 

// What: "@(#) ShadowPetscSOE.h, revA"

#ifndef ShadowPetscSOE_h
#define ShadowPetscSOE_h

#include <LinearSOE.h>
#include <Vector.h>
#include <PetscSOE.h>

extern "C" {
#include <petsc.h>
}

class PetscSolver;

class ShadowPetscSOE : public LinearSOE
{
  public:
    ShadowPetscSOE(PetscSolver &theSolver, int blockSize);    
    
    ~ShadowPetscSOE();

    int solve(void);    

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        

    void zeroA(void);
    void zeroB(void);

    const Vector &getX(void);
    const Vector &getB(void);
    double normRHS(void);

    void setX(int loc, double value);    

    int setSolver(PetscSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
  protected:
    
  private:
    MPI_Comm theComm; // a comm for communicating to the ActorPetscSOE's
                      // without using PETSC_COMM_WORLD
    PetscSOE theSOE;  // the local portion of the SOE
    PetscSolver *theSolver; // created by the user
    int myRank;
    int numProcessors;
    int sendData[3];
    void *sendBuffer;
    int blockSize;
};


#endif






