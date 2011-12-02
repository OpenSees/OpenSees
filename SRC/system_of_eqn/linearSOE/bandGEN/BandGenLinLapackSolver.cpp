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
// $Date: 2000-09-15 08:23:28 $
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
#include <f2c.h>
#include <math.h>


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

extern "C" int _stdcall DGBSV(int *N, int *KL, int *KU, int *NRHS, double *A, 
			      int *LDA, int *iPiv, double *B, int *LDB, 
			      int *INFO);

extern "C" int _stdcall DGBTRS(char *TRANS, unsigned int *sizeC, 
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
	cerr << "WARNING BandGenLinLapackSolver::solve(void)- ";
	cerr << " No LinearSOE object has been set\n";
	return -1;
    }

    int n = theSOE->size;    
    // check iPiv is large enough
    if (iPivSize < n) {
	cerr << "WARNING BandGenLinLapackSolver::solve(void)- ";
	cerr << " iPiv not large enough - has setSize() been called?\n";
	return -1;
    }	    

    int kl = theSOE->numSubD;
    int ku = theSOE->numSuperD;
    int ldA = 2*kl + ku +1;
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;
    int    *iPIV = iPiv;
    
    // first copy B into X
    for (int i=0; i<n; i++) {
	*(Xptr++) = *(Bptr++);
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
	DGBTRS("N", &sizeC, &n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
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
    // check if successfull
    if (info != 0) {
	cerr << "WARNING BandGenLinLapackSolver::solve() -";
	cerr << "LAPACK routine returned " << info << endl;
	return -info;
    }

    theSOE->factored = true;
    return 0;
}
    


int
BandGenLinLapackSolver::setSize()
{
    // if iPiv not big enough, free it and get one large enough
    if (iPivSize < theSOE->size) {
	if (iPiv != 0)
	    delete [] iPiv;
	
	iPiv = new int[theSOE->size];
	if (iPiv == 0) {
	    cerr << "WARNING BandGenLinLapackSolver::setSize() ";
	    cerr << " - ran out of memory for iPiv of size ";
	    cerr << theSOE->size << endl;
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
