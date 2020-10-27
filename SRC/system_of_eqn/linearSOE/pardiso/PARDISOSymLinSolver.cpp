// File: ~/system_of_eqn/linearSOE/PARDISOLinSolver.h
//
// Written: M. Salehi
// Created: 02/19
// Revision: A
//
// Description: This file contains the class definition for 
// PARDISOLinSolver. It solves the Sparse General SOE by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) PARDISOSymLinSolver.h, revA"



#include <PARDISOSymLinSolver.h>
#include <PARDISOSymLinSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>

PARDISOSymLinSolver::PARDISOSymLinSolver()
:LinearSOESolver(SOLVER_TAGS_PARDISOSymLinSolver),
 theSOE(0)
{
	
}


PARDISOSymLinSolver::~PARDISOSymLinSolver()
{ 

}



int
PARDISOSymLinSolver::solve(void)
{ 
    if (theSOE == 0) {
	opserr << "WARNING PARDISOLinSolver::solve(void)- ";
	opserr << " No LinearSOE object has been set\n";
	return -1;
    }
	int n = theSOE->size;
	int* ia = theSOE->rowStartA;
	int* ja = theSOE->colA;
	double* a = theSOE->A;

	double* Xptr = theSOE->X;
	double* Bptr = theSOE->B;

	/*opserr << "ia : ";
	for (int i = 0; i < n + 1; i++)
	{
		opserr << ia[i] << " ";
	}

	opserr << "\nja : ";
	for (int i = 0; i < theSOE->nnz; i++)
	{
		opserr << ja[i] << " ";
	}
	opserr << endln;

	opserr << "\na : ";
	for (int i = 0; i < theSOE->nnz; i++)
	{
		opserr << a[i] << " ";
	}
	opserr << endln;*/

	/*1	Real and structurally symmetric	IN
	2	Real and symmetric positive definite
	-2	Real and symmetric indefinite
	3	Complex and structurally symmetric
	4	Complex and Hermitian positive definite
	-4	Complex and Hermitian indefinite
	6	Complex and symmetric matrix
	11	Real and unsymmetric matrix
	13	Complex and unsymmetric matrix*/
	int mtype = 2;
	int nrhs =1;
	void* pt[64];
	int* iparm = new int[64];
	int maxfct, mnum, phase, error, msglvl;

	
	double ddum;			/* Double dummy */
	int idum;			/* Integer dummy. */

	for (int i = 0; i < 64; i++)
	{
		iparm[i] = 0;
		pt[i] = 0;
	}

	iparm[0] = 1;			/* No solver default */
	iparm[1] = 2;			/* Fill-in reordering from METIS */
	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0;			/* No iterative-direct algorithm */
	iparm[4] = 0;			/* No user fill-in reducing permutation */
	iparm[5] = 0;			/* Write solution into x */
	iparm[6] = 0;			/* Not in use */
	iparm[7] = 2;			/* Max numbers of iterative refinement steps */
	iparm[8] = 0;			/* Not in use */
	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;		/* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;		/* Not in use */
	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;		/* Output: Number of perturbed pivots */
	iparm[14] = 0;		/* Not in use */
	iparm[15] = 0;		/* Not in use */
	iparm[16] = 0;		/* Not in use */
	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;		/* Output: Mflops for LU factorization */
	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
	maxfct = 1;			/* Maximum number of numerical factorizations. */
	mnum = 1;			/* Which factorization to use. */
	msglvl = 0;			/* Print statistical information in file */
	error = 0;			/* Initialize error flag */

	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */

//#if _DEBUG
	//opserr << "zise :" << theSOE->size << endln;
	//opserr << "non zero length :" << n << endln;
	//opserr << "A :";
	//for (int i = 0; i < theSOE->nnz; i++)
	//{
	//	opserr << " " << a[i];
	//}
	//opserr << endln << "ja :" << n << endln;
	//for (int i = 0; i < theSOE->nnz; i++)
	//{
	//	opserr << " " << (int)ja[i];
	//}
	//opserr << endln << "ia :" << n << endln;
	//for (int i = 0; i < theSOE->size + 1; i++)
	//{
	//	opserr << " " << (int)ia[i];
	//}
	//opserr << endln;
//#endif

	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

	if (error != 0)
	{
		opserr << "\nERROR during symbolic factorization: " << error;
		return -1;
	}

	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */

	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		opserr << "\nERROR during numerical factorization: " << error;
		return -2;
	}

	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2;			/* Max numbers of iterative refinement steps. */

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);

	if (error != 0)
	{
		opserr << "\nERROR during solution: " << error;
		return -3;
	}

	/* Solve completed ... 
	   The solution of the system is: */

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;			/* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	//if (pt != 0)
	//	delete[] pt;

	//if (iparm != 0)
	//	delete[] iparm;

	/*opserr << endln << "x :" << n << endln;
	for (int i = 0; i < theSOE->size; i++)
	{
		opserr << " " << Xptr[i];
	}
	opserr << endln;*/
    return 0;
}


int
PARDISOSymLinSolver::setSize()
{
    // nothing to do
    return 0;
}


int
PARDISOSymLinSolver::setLinearSOE(PARDISOSymLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    return 0;
}


int
PARDISOSymLinSolver::sendSelf(int cTAg, Channel &theChannel)
{
    // doing nothing
    return 0;
}


int
PARDISOSymLinSolver::recvSelf(int cTag,
			     Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}




