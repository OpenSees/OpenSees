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
// $Date: 2007-02-14 20:12:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/PetscSOE.h,v $
                                                                        
                                                                        
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the class definition for PetscSOE
// PetscSOE is a subclass of LinearSOE. It uses the LAPACK storage
// scheme to store the components of the A matrix, which is a full matrix.


// What: "@(#) PetscSOE.h, revA"

#ifndef PetscSOE_h
#define PetscSOE_h

#include <LinearSOE.h>
#include <Vector.h>

//extern "C" {
#include <petscksp.h>
//}

class PetscSolver;

class PetscSOE : public LinearSOE
{
  public:
    PetscSOE(PetscSolver &theSolver, int blockSize=1);    
    
    ~PetscSOE();

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
    void setX(const Vector &x);    

    int setSolver(PetscSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    friend class PetscSolver;
    friend class ActorPetscSOE;
    friend class ShadowPetscSOE;
    
  protected:
    int setChannels(int nChannels, Channel **theChannels);
    
  private:
    int isFactored;
    int size;
    int processID;
    int numProcesses;

    double *B, *X;
    int *indices;
    Vector *vectX;
    Vector *vectB;
    Mat A;
    Vec x, b;
    int blockSize;
    PetscTruth flg;

    int numChannels;
    Channel **theChannels;
    ID **localCol;

    int startRow, endRow;
};


#endif
 
