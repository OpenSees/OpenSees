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


#ifndef PetscSolver_h
#define PetscSolver_h

#include <petscksp.h>
#include <LinearSOESolver.h>


// Uncomment or define -D_DEBUG_PETSCSOLVER to enable debugging
// #define _DEBUG_PETSCSOLVER

#ifdef _DEBUG_PETSCSOLVER
#define PETSCSOLVER_DEBUGOUT cout // or any other ostream
#else
#define PETSCSOLVER_DEBUGOUT 0 && cout
#endif

// Uncomment to enable logger that will creat petsc_log.txt on each
// run, to profile PETSc. 
// #define _ENABLE_PETSC_LOGGER

class PetscSOE;

class PetscSolver : public LinearSOESolver
{
public:
    PetscSolver();
    PetscSolver(KSPType method, PCType preconditioner);
    PetscSolver(KSPType method, PCType preconditioner, double rTol, double aTol, double dTol, int maxIts, MatType mat=MATMPIAIJ);//Guanzhou
    ~PetscSolver();

    int solve(void);
    int setSize(void);
    virtual int setLinearSOE(PetscSOE &theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                    FEM_ObjectBroker &theBroker);

    friend class ActorPetscSOE;
    friend class ShadowPetscSOE;

protected:
    PetscSOE *theSOE;

private:
    KSP ksp;
    PC pc;
    int its;
    KSPType method;
    PCType preconditioner;
    double rTol;
    double aTol;
    double dTol;
    int maxIts;
    MatType matType;

    bool is_KSP_initialized;
};

#endif

