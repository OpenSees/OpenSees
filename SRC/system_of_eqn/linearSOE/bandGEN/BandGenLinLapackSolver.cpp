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
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandGEN/BandGenLinLapackSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/bandGEN/BandGenLinLapackSolver.h
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// BandGenLinLapackSolver. It solves the BandGenLinSOE object by calling
// Lapack routines.
//
// What: "@(#) BandGenLinLapackSolver.h, revA"

#include <BandGenLinLapackSolver.h>
#include <BandGenLinSOE.h>
#include <cmath>

void* OPS_BandGenLinLapack()
{
    BandGenLinSolver *theSolver = new BandGenLinLapackSolver();
    BandGenLinSOE *theSOE = new BandGenLinSOE(*theSolver);
    return theSOE;
}

BandGenLinLapackSolver::BandGenLinLapackSolver()
:BandGenLinSolver(SOLVER_TAGS_BandGenLinLapackSolver),
 iPiv(0), iPivSize(0)
{
    
}

BandGenLinLapackSolver::~BandGenLinLapackSolver()
{
    if (iPiv != 0)
	delete [] iPiv;
}

#ifdef _WIN32

extern "C" int DGBSV(int *N, int *KL, int *KU, int *NRHS, double *A, 
			      int *LDA, int *iPiv, double *B, int *LDB, 
			      int *INFO);

extern "C" int DGBTRS(char *TRANS, 
			       int *N, int *KL, int *KU, int *NRHS,
			       double *A, int *LDA, int *iPiv, 
			       double *B, int *LDB, int *INFO);

#else

extern "C" int dgbsv_(int *N, int *KL, int *KU, int *NRHS, double *A, 
		      int *LDA, int *iPiv, double *B, int *LDB, int *INFO);
		      

extern "C" int dgbtrs_(char *TRANS, int *N, int *KL, int *KU, int *NRHS, 
		       double *A, int *LDA, int *iPiv, double *B, int *LDB, 
		       int *INFO);
#endif
int
BandGenLinLapackSolver::solve(void)
{
    if (theSOE == 0) {
	opserr << "WARNING BandGenLinLapackSolver::solve(void)- ";
	opserr << " No LinearSOE object has been set\n";
	return -1;
    }

    int n = theSOE->size;

    if (n == 0)
	return 0;

    // check iPiv is large enough
    if (iPivSize < n) {
	opserr << "WARNING BandGenLinLapackSolver::solve(void)- ";
	opserr << " iPiv not large enough - has setSize() been called?\n";
	return -1;
    }	    

    int kl = theSOE->numSubD;
    int ku = theSOE->numSuperD;

    // LAPACK requires non-negative bandwidths no larger than n - 1 and
    // enough leading storage for the band matrix. Reject malformed state
    // before it reaches the Fortran routine, which otherwise may crash the
    // host process on some platforms.
    if (kl < 0 || ku < 0 || kl >= n || ku >= n) {
	opserr << "WARNING BandGenLinLapackSolver::solve() - invalid band dimensions" << endln;
	return -1;
    }

    const long long ldAValue = 2LL * kl + ku + 1;
    if (ldAValue <= 0 || ldAValue * n > theSOE->Asize) {
	opserr << "WARNING BandGenLinLapackSolver::solve() - invalid leading dimension" << endln;
	return -1;
    }

    int ldA = static_cast<int>(ldAValue);
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;
    int    *iPIV = iPiv;

    if (Aptr == 0 || Xptr == 0 || Bptr == 0 || iPIV == 0) {
	opserr << "WARNING BandGenLinLapackSolver::solve() - incomplete solver storage" << endln;
	return -1;
    }

    for (long long i = 0; i < static_cast<long long>(ldA) * n; ++i) {
	if (!std::isfinite(Aptr[i])) {
	    opserr << "WARNING BandGenLinLapackSolver::solve() - non-finite value in matrix A" << endln;
	    return -1;
	}
    }
    
    // first copy B into X
    for (int i=0; i<n; i++) {
	*(Xptr++) = *(Bptr++);
	if (!std::isfinite(*(Xptr - 1))) {
	    opserr << "WARNING BandGenLinLapackSolver::solve() - non-finite value in right-hand side" << endln;
	    return -1;
	}
    }
    Xptr = theSOE->X;

    // now solve AX = B

#ifdef _WIN32
    {if (theSOE->factored == false)  
	// factor and solve 
	DGBSV(&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);	
    else  {
	// solve only using factored matrix
	unsigned int sizeC = 1;
	//DGBTRS("N",&sizeC,&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
    char type[] = "N";
	DGBTRS(type,&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
    }}
#else
    {if (theSOE->factored == false)      
	// factor and solve 	
	dgbsv_(&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
    else
	// solve only using factored matrix	
	dgbtrs_("N",&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
    }
#endif
    // check if successful
    if (info != 0) {
      if (info > 0) {
	opserr << "WARNING BandGenLinLapackSolver::solve() -";
	opserr << "factorization failed, matrix singular U(i,i) = 0, i= " << info-1 << endln;
	return -info+1;
      } else {
	opserr << "WARNING BandGenLinLapackSolver::solve() - OpenSees code error\n";
	return info;
      }
    }

    theSOE->factored = true;
    return 0;
}
    


int
BandGenLinLapackSolver::setSize()
{
    if (theSOE->size == 0)
	return 0;

    // if iPiv not big enough, free it and get one large enough
    if (iPivSize < theSOE->size) {
	if (iPiv != 0)
	    delete [] iPiv;
	
	iPiv = new int[theSOE->size];
	if (iPiv == 0) {
	    opserr << "WARNING BandGenLinLapackSolver::setSize() ";
	    opserr << " - ran out of memory for iPiv of size ";
	    opserr << theSOE->size << endln;
	    return -1;
	} else
	    iPivSize = theSOE->size;
    }
	
    return 0;
}

int    
BandGenLinLapackSolver::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
BandGenLinLapackSolver::recvSelf(int commitTag,
				 Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}
