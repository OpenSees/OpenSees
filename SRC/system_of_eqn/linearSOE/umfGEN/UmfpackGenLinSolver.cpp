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
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.C
//
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for 
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK2.2.1 routines.
//
// What: "@(#) UmfpackGenLinSolver.C, revA"

#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#ifdef _WIN32
extern "C" int _stdcall UMD21I(int *keep, double *cntl, int *icntl);
#else
extern "C" int umd21i_(int *keep, double *cntl, int *icntl);
#endif


UmfpackGenLinSolver::
UmfpackGenLinSolver()
:LinearSOESolver(SOLVER_TAGS_UmfpackGenLinSolver),
 copyIndex(0), lIndex(0), work(0), theSOE(0)
{
  // perform the initialisation needed in UMFpack
#ifdef _WIN32
  UMD21I(keep, cntl, icntl);    
#else
  umd21i_(keep, cntl, icntl);
#endif
}


UmfpackGenLinSolver::~UmfpackGenLinSolver()
{
  if (copyIndex != 0)
    delete [] copyIndex;

  if (work != 0)
    delete [] work;
}

#ifdef _WIN32
extern "C" int _stdcall UMD2FA(int *n, int *ne, int *job, logical *transa,
		       int *lvalue, int *lindex, double *value,
		       int *index, int *keep, double *cntl, int *icntl,
		       int *info, double *rinfo);

/*
extern "C" int umd2rf_(int *n, int *ne, int *job, logical *transa,
		       int *lvalue, int *lindex, double *value,
		       int *index, int *keep, double *cntl, int *icntl,
		       int *info, double *rinfo);
*/
extern "C" int _stdcall UMD2SO(int *n, int *job, logical *transa,
		       int *lvalue, int *lindex, double *value,
		       int *index, int *keep, double *b, double *x, 
		       double *w, double *cntl, int *icntl,
		       int *info, double *rinfo);
#else
extern "C" int umd2fa_(int *n, int *ne, int *job, logical *transa,
		       int *lvalue, int *lindex, double *value,
		       int *index, int *keep, double *cntl, int *icntl,
		       int *info, double *rinfo);

extern "C" int umd2so_(int *n, int *job, logical *transa,
		       int *lvalue, int *lindex, double *value,
		       int *index, int *keep, double *b, double *x, 
		       double *w, double *cntl, int *icntl,
		       int *info, double *rinfo);
#endif

int
UmfpackGenLinSolver::solve(void)
{
    if (theSOE == 0) {
	cerr << "WARNING UmfpackGenLinSolver::solve(void)- ";
	cerr << " No LinearSOE object has been set\n";
	return -1;
    }
    
    int n = theSOE->size;
    int ne = theSOE->nnz;
    int lValue = theSOE->lValue;

    // check for quick return
    if (n == 0)
	return 0;

    // first copy B into X
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;
    double *Aptr = theSOE->A;

    int job =0; // set to 1 if wish to do iterative refinment
    logical trans = FALSE_;

    if (theSOE->factored == false) {

      // make a copy of index
      for (int i=0; i<2*ne; i++) {
	copyIndex[i] = theSOE->index[i];
      }

      // factor the matrix
#ifdef _WIN32
      UMD2FA(&n, &ne, &job, &trans, &lValue, &lIndex, Aptr,
	      copyIndex, keep, cntl, icntl, info, rinfo);
#else
      umd2fa_(&n, &ne, &job, &trans, &lValue, &lIndex, Aptr,
	      copyIndex, keep, cntl, icntl, info, rinfo);
#endif      
      
      if (info[0] != 0) {	
	cerr << "WARNING UmfpackGenLinSolver::solve(void)- ";
	cerr << info[0] << " returned in factorization UMD2FA()\n";
	return -info[0];
      }
      theSOE->factored = true;
    }	

    // do forward and backward substitution
#ifdef _WIN32
    UMD2SO(&n, &job, &trans, &lValue, &lIndex, Aptr, copyIndex, 
	   keep, Bptr, Xptr, work, cntl, icntl, info, rinfo);    
#else
    umd2so_(&n, &job, &trans, &lValue, &lIndex, Aptr, copyIndex, 
	    keep, Bptr, Xptr, work, cntl, icntl, info, rinfo);
#endif

    if (info[0] != 0) {	
       cerr << "WARNING UmfpackGenLinSolver::solve(void)- ";
       cerr << info[0] << " returned in substitution dgstrs()\n";
       return -info[0];
    }

    return 0;
}


int
UmfpackGenLinSolver::setSize()
{
    int n = theSOE->size;
    int ne = theSOE->nnz;
    if (n > 0) {
      if (work != 0)
	delete [] work;

      work = new double[4*n];

      lIndex = 37*n + 4*ne + 10;
      if (copyIndex != 0)
	delete [] copyIndex;

      copyIndex = new int[lIndex];
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

