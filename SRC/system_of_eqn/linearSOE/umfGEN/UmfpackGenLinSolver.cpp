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
// $Date: 2009-05-20 17:30:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/98
//
// Description: This file contains the class definition for 
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK5.7.1 routines.
//
// What: "@(#) UmfpackGenLinSolver.C, revA"

#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

void* OPS_UmfpackGenLinSolver()
{
    UmfpackGenLinSolver *theSolver = new UmfpackGenLinSolver();
    return new UmfpackGenLinSOE(*theSolver);  
}

UmfpackGenLinSolver::
UmfpackGenLinSolver()
    :LinearSOESolver(SOLVER_TAGS_UmfpackGenLinSolver), Symbolic(0), theSOE(0)
{
}


UmfpackGenLinSolver::~UmfpackGenLinSolver()
{
    if (Symbolic != 0) {
	umfpack_di_free_symbolic(&Symbolic);
    }
}

int
UmfpackGenLinSolver::solve(void)
{
    int n = theSOE->X.Size();
    int nnz = (int)theSOE->Ai.size();
    if (n == 0 || nnz==0) return 0;
    
    int* Ap = &(theSOE->Ap[0]);
    int* Ai = &(theSOE->Ai[0]);
    double* Ax = &(theSOE->Ax[0]);
    double* X = &(theSOE->X(0));
    double* B = &(theSOE->B(0));

    // check if symbolic is done
    if (Symbolic == 0) {
	opserr<<"WARNING: setSize has not been called -- Umfpackgenlinsolver::solve\n";
	return -1;
    }
    
    // numerical analysis
    void* Numeric = 0;
    int status = umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);

    // check error
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: numeric analysis returns "<<status<<" -- Umfpackgenlinsolver::solve\n";
	return -1;
    }

    // solve
    status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,X,B,Numeric,Control,Info);

    // delete Numeric
    if (Numeric != 0) {
	umfpack_di_free_numeric(&Numeric);
    }
    
    // check error
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: solving returns "<<status<<" -- Umfpackgenlinsolver::solve\n";
	return -1;
    }

    return 0;
}


int
UmfpackGenLinSolver::setSize()
{
    // set default control parameters
    umfpack_di_defaults(Control);
    Control[UMFPACK_PIVOT_TOLERANCE] = 1.0;
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;

    int n = theSOE->X.Size();
    int nnz = (int)theSOE->Ai.size();
    if (n == 0 || nnz==0) return 0;
    
    int* Ap = &(theSOE->Ap[0]);
    int* Ai = &(theSOE->Ai[0]);
    double* Ax = &(theSOE->Ax[0]);

    // symbolic analysis
    if (Symbolic != 0) {
	umfpack_di_free_symbolic(&Symbolic);
    }
    int status = umfpack_di_symbolic(n,n,Ap,Ai,Ax,&Symbolic,Control,Info);

    // check error
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: symbolic analysis returns "<<status<<" -- Umfpackgenlinsolver::setsize\n";
	Symbolic = 0;
	return -1;
    }
    return 0;
}

int
UmfpackGenLinSolver::setLinearSOE(UmfpackGenLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    return 0;
}

int
UmfpackGenLinSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
UmfpackGenLinSolver::recvSelf(int ctag,
			      Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}

