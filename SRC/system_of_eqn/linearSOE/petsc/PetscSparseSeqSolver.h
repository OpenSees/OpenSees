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
// $Date: 2005/05/18 19:26:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/PetscSparseSeqSolver.h,v $

// Written: fmk
// Created: 04/05

//
// Description: This file contains the class definition for
// PetscSparseSeqSolver. It solves the SparseGenLinSOE object by calling Petsc routines
//
// What: "@(#) PetscSparseSeqSolver.h, revA"

#ifndef PetscSparseSeqSolver_h
#define PetscSparseSeqSolver_h

#include <petscksp.h>
#include <SparseGenRowLinSolver.h>

#define PetscMinRelTol 1.0e-6

class SparseGenRowLinSOE;

class PetscSparseSeqSolver : public SparseGenRowLinSolver
{
    public:
        PetscSparseSeqSolver(KSPType method, PCType preconditioner);
        PetscSparseSeqSolver(KSPType method, PCType preconditioner, double rTol, double aTol, double dTol, int maxIts);
        ~PetscSparseSeqSolver();

        int solve(void);
        int setSize(void);
        int getNumIterations(void);
        virtual int setLinearSOE(SparseGenRowLinSOE& theSOE);

        int sendSelf(int commitTag, Channel& theChannel);
        int receiveSelf(int commitTag, Channel& theChannel,
                     FEM_ObjectBroker& theBroker);

    protected:
        SparseGenRowLinSOE* theSOE;

    private:
        KSP ksp;                      // linear solver context
        PC  pc;
        int its;
        KSPType method;
        PCType  preconditioner;
        double rTol;
        double aTol;
        double dTol;
        int maxIts;

        Mat A;
        Vec x, b;
};

#endif

