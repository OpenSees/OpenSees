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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/fullGEN/FullGenLinLapackSolver.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
//
// Description: This file contains the class definition for 
// FullGenLinLapackSolver. It solves the FullGenLinSOE object by calling
// Lapack routines.
//
// What: "@(#) FullGenLinLapackSolver.h, revA"

#include <FullGenLinLapackSolver.h>
#include <FullGenLinSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

void* OPS_FullGenLinLapackSolver()
{
    FullGenLinLapackSolver *theSolver = new FullGenLinLapackSolver();
    return new FullGenLinSOE(*theSolver);
}

FullGenLinLapackSolver::FullGenLinLapackSolver()
:FullGenLinSolver(SOLVER_TAGS_FullGenLinLapackSolver),iPiv(0),sizeIpiv(0)
{
    
}

FullGenLinLapackSolver::~FullGenLinLapackSolver()
{
    if (iPiv != 0)
	delete [] iPiv;
}


#ifdef _WIN32
extern "C" int  DGESV(int *N, int *NRHS, double *A, int *LDA, 
			      int *iPiv, double *B, int *LDB, int *INFO);
			     
extern "C" int  DGETRS(char *TRANS,
			       int *N, int *NRHS, double *A, int *LDA, 
			       int *iPiv, double *B, int *LDB, int *INFO);
#else
extern "C" int dgesv_(int *N, int *NRHS, double *A, int *LDA, int *iPiv, 
		      double *B, int *LDB, int *INFO);

extern "C" int dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, 
		       int *iPiv, double *B, int *LDB, int *INFO);		       
#endif
int
FullGenLinLapackSolver::solve(void)
{
    if (theSOE == 0) {
	opserr << "WARNING FullGenLinLapackSolver::solve(void)- ";
	opserr << " No LinearSOE object has been set\n";
	return -1;
    }
    
    int n = theSOE->size;
    
    // check for quick return
    if (n == 0)
	return 0;
    
    // check iPiv is large enough
    if (sizeIpiv < n) {
	opserr << "WARNING FullGenLinLapackSolver::solve(void)- ";
	opserr << " iPiv not large enough - has setSize() been called?\n";
	return -1;
    }	
	
    int ldA = n;
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;
    int *iPIV = iPiv;
    
    // first copy B into X
    for (int i=0; i<n; i++)
	*(Xptr++) = *(Bptr++);
    Xptr = theSOE->X;

    // now solve AX = Y

#ifdef _WIN32
    {if (theSOE->factored == false)  
	// factor and solve 
	DGESV(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
     else {
	// solve only using factored matrix	 
	 unsigned int sizeC = 1;
	 DGETRS("N", &n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);	 
     }}
#else
    {if (theSOE->factored == false)      
	dgesv_(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
     else
	dgetrs_("N", &n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
    }
#endif
    
    // check if successful
    if (info != 0) {
	opserr << "WARNING FullGenLinLapackSolver::solve()";
	opserr << " - lapack solver failed - " << info << " returned\n";
	return -info;
    }

    
    theSOE->factored = true;
    return 0;
}


int
FullGenLinLapackSolver::setSize()
{
    int n = theSOE->size;
    if (n > 0) {
	if (sizeIpiv < n) {
	    if (iPiv != 0)
		delete [] iPiv;
	    iPiv = new int[n];		
	    if (iPiv == 0) {
		opserr << "WARNING FullGenLinLapackSolver::setSize()";
		opserr << " - ran out of memory\n";
		return -1;
	    }		
	    sizeIpiv = n;
	}
    } else if (n == 0)
	return 0;
    else {
	opserr << "WARNING FullGenLinLapackSolver::setSize()";
	opserr << " - ran out of memory\n";
	return -1;	
    }
	
    return 0;
}

int
FullGenLinLapackSolver::sendSelf(int commitTag,
				 Channel &theChannel)

{
    // nothing to do
    return 0;
}

int
FullGenLinLapackSolver::recvSelf(int commitTag,
				 Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}



