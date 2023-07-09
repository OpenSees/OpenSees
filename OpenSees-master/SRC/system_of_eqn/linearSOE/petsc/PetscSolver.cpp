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


#include <OPS_Globals.h>
#include <PetscSolver.h>
#include <PetscSOE.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>
#include <ID.h>
#include <Vector.h>
#include <string.h>
#include "petscksp.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <petsctime.h>

using namespace std;
//-------------------------------------

PetscSolver::PetscSolver()
    : LinearSOESolver(SOLVER_TAGS_PetscSolver),
      rTol(PETSC_DEFAULT), aTol(PETSC_DEFAULT), dTol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT), matType(MATMPIAIJ),
      is_KSP_initialized(false)
{

}

PetscSolver::PetscSolver(KSPType meth, PCType pre)
    : LinearSOESolver(SOLVER_TAGS_PetscSolver), method(meth), preconditioner(pre),
      rTol(PETSC_DEFAULT), aTol(PETSC_DEFAULT), dTol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT), matType(MATMPIAIJ),
      is_KSP_initialized(false)
{

}

PetscSolver::PetscSolver(KSPType meth, PCType pre, double relTol, double absTol, double divTol, int maxIterations, MatType mat)
    : LinearSOESolver(SOLVER_TAGS_PetscSolver), method(meth), preconditioner(pre),
      rTol(relTol), aTol(absTol), dTol(divTol), maxIts(maxIterations), matType(mat),
      is_KSP_initialized(false)
{

}

PetscSolver::~PetscSolver()
{
    // KSPDestroy(&ksp);
    // CHKERRQ(ierr);
}


int
PetscSolver::solve(void)
{

    int size = theSOE->size;
    // int numProcesses = theSOE->numProcesses;
    int processID = theSOE->processID;
    int ierr;

    if (processID >= 0)
    {

        //MatSetType(theSOE->A, matType);
        PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") MatAssemblyBegin\n";
        ierr = MatAssemblyBegin(theSOE->A, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);

        ierr = MatAssemblyEnd(theSOE->A, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") MatAssemblyEnd\n";


        if (not is_KSP_initialized)
        {

            PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") KSPCreate Begin\n";
            KSPCreate(PETSC_COMM_WORLD, &ksp);
            KSPSetFromOptions(ksp);
            PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") KSPCreate End\n";


            is_KSP_initialized = true;
        }
        
        PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") KSPSetOperators Begin\n";
        KSPSetOperators(ksp, theSOE->A, theSOE->A);
        PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") KSPSetOperators End\n";

    }



    // Everyone sends the assembled vector B to the rank 0 processor. Which sends it back...
    static Vector recvVector(1);


    Vector *vectX = theSOE->vectX;
    Vector *vectB = theSOE->vectB;

    // zero X
    vectX->Zero();

    //
    // form B on each
    //
    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") FormBBegin\n";
    int numChannels = theSOE->numChannels;
    Channel **theChannels = theSOE->theChannels;

    if (processID != 0)
    {
        Channel *theChannel = theChannels[0];

        theChannel->sendVector(0, 0, *vectB);
        theChannel->recvVector(0, 0, *vectB);

        double vectBNorm  = vectB->Norm();
        if (std::isnan(vectBNorm))
        {
            cout << "norm(vectB) = " << vectBNorm << endl;
            return -1;
        }

    }
    else //Master process assembles b vector. This is a bottleneck and can be improved with
        //a collective reduction
    {

        if (recvVector.Size() != size)
        {
            recvVector.resize(size);
        }

        //Better done with ALLREDUCE
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->recvVector(0, 0, recvVector);
            *vectB += recvVector;

            double rvNorm  = recvVector.Norm() ;
            if (std::isnan(rvNorm))
            {
                cout << "norm(recvVector) = " << rvNorm << endl;
                return -1;
            }
        }

        //Better done with a BCast
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->sendVector(0, 0, *vectB);
        }
    }
    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") FormBEnd\n";


    //
    // solve and mark as having been solved
    //

    double t1 = 0;
    double t2 = 0;
    PetscTime(&t1);
    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") SolveBegin\n";
    if (processID >= 0)
    {
        ierr = KSPSolve(ksp, theSOE->b, theSOE->x);
        CHKERRQ(ierr);
        theSOE->isFactored = 1;
    }
    PetscTime(&t2);
    double delta_time = t2 - t1;
    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") SolveEnd dt = " << delta_time << "s \n";





    //
    // if parallel, we must form the total X: each processor has startRow through endRow-1
    //  
    vectX = theSOE->vectX;

    numChannels = theSOE->numChannels;
    theChannels = theSOE->theChannels;
    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") SendXBegin\n";
    if (processID != 0)
    {
        Channel *theChannel = theChannels[0];

        theChannel->sendVector(0, 0, *vectX);
        theChannel->recvVector(0, 0, *vectX);

    }
    else //Again, the master process forms the global X vector which is then sent to all processors
    {

        if (recvVector.Size() != size)
        {
            recvVector.resize(size);
        }

        //Better done with ALLREDUCE
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->recvVector(0, 0, recvVector);
            *vectX += recvVector;
        }

        //Better done with a BCast
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->sendVector(0, 0, *vectX);
        }
    }

    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") SendXEnd\n";
    // Destroy KSP and collect the error at P0
    KSPDestroy(&ksp);
    is_KSP_initialized = false;

    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") Logging Begin\n";

#ifdef _ENABLE_PETSC_LOGGER
        PetscViewer    viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, "petsc_log.txt" , &viewer);
        PetscLogView(viewer);
#endif
        
    PETSCSOLVER_DEBUGOUT << "PetscSolver::solve (" << processID << ") Logging End - Done solve\n";

    return ierr;
}






int
PetscSolver::setSize()
{
    /*
     * Create linear solver context
     */
    if (theSOE->processID >= 0)
    {
        KSPCreate(PETSC_COMM_WORLD, &ksp);
    }


    return 0;
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
   
    return 0;
}

int
PetscSolver::recvSelf(int cTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker)
{
   

    return 0;
}



