// File: ~/system_of_eqn/linearSOE/symLinSolver/SymSparseLinSolver.C
//
// Written: Jun Peng
// Created: 12/98
// Revision: A
//
// Description: This file contains the class definition for 
// SymSparseinSolver. It solves the SymSparseLinSOE object by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) SymSparseLinSolver.C, revA"


#include "SymSparseLinSOE.h"
#include "SymSparseLinSolver.h"
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

extern "C" {
  #include "FeStructs.h"
  #include "globalVars.h"
  #include "tim.h"
}


SymSparseLinSolver::SymSparseLinSolver()
:LinearSOESolver(SOLVER_TAGS_SymSparseLinSolver),
 theSOE(0), factored(false)
{
    // nothing to do.
}


SymSparseLinSolver::~SymSparseLinSolver()
{ 
    // nothing to do.
}


extern "C" int pfsfct(int neqns, double *diag, double **penv, int nblks, int *xblk,
		      OFFDBLK **begblk, OFFDBLK *first, int *rowblks);

extern "C" void pfsslv(int neqns, double *diag, double **penv, int nblks,
		       int *xblk, double *rhs, OFFDBLK **begblk);


int
SymSparseLinSolver::solve(void)
{ 
    if (theSOE == 0) {
	cerr << "WARNING SymSparseLinSolver::solve(void)- ";
	cerr << " No LinearSOE object has been set\n";
	return -1;
    }

    int      nblks = theSOE->nblks;
    int      *xblk = theSOE->xblk;
    int      *invp = theSOE->invp;
    double   *diag = theSOE->diag;
    double   **penv = theSOE->penv;
    int      *rowblks = theSOE->rowblks;
    OFFDBLK  **begblk = theSOE->begblk;
    OFFDBLK  *first = theSOE->first;
    
    
    int neq = theSOE->size;

    // check for quick return
    if (neq == 0)
	return 0;

    // first copy B into X

    for (int i=0; i<neq; i++) {
        theSOE->X[i] = theSOE->B[i];
    }
    double *Xptr = theSOE->X;

    timer (FACTOR);
    if (factored == false) {

        //factor the matrix
        //call the "C" function to do the numerical factorization.
        int factor;
	factor = pfsfct(neq, diag, penv, nblks, xblk, begblk, first, rowblks);
	if (factor>0) {
	    cerr << "In SymSparseLinSolver: error in factorization.\n";
	    return -1;
	}
	factored = true;
    }

    //do forward and backward substitution.
    //call the "C" function.
    timer (SOLVE);
    pfsslv(neq, diag, penv, nblks, xblk, Xptr, begblk);


    //Since the X we get by solving AX=B is P*X, we need to reordering
    //the Xptr to ge the wanted X.

    double *tempX = new double[neq];
    if (tempX == 0) {
        cerr << "WARNING LawSparseLinSover::SymSparseLinSolver :";
	cerr << " ran out of memory for vectors (tempX) ";
	return -1;
    } 

    for (int m=0; m<neq; m++) {
        tempX[m] = Xptr[invp[m]];
    }
	
    for (int k=0; k<neq; k++) {
        Xptr[k] = tempX[k];
    }
	
    delete [] tempX;

    timer (END);

    return 0;
}


int
SymSparseLinSolver::setSize()
{
    // nothing to do
    return 0;
}


int
SymSparseLinSolver::setLinearSOE(SymSparseLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    return 0;
}


int
SymSparseLinSolver::sendSelf(int cTAg, Channel &theChannel)
{
    // nothing to do
    return 0;
}


int
SymSparseLinSolver::recvSelf(int cTag,
			     Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}




