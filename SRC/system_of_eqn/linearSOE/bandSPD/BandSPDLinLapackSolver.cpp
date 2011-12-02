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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-04-02 22:02:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinLapackSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/bandSPD/BandSPDLinLapackSolver.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of BandSPDLinLapackSolver.
//
// What: "@(#) BandSPDLinLapackSolver.h, revA"

#include <BandSPDLinLapackSolver.h>
#include <BandSPDLinSOE.h>
#include <f2c.h>
#include <math.h>


BandSPDLinLapackSolver::BandSPDLinLapackSolver()
:BandSPDLinSolver(SOLVER_TAGS_BandSPDLinLapackSolver)
{
    
}

BandSPDLinLapackSolver::~BandSPDLinLapackSolver()
{
    
}


#ifdef _WIN32
extern "C" int  DPBSV(char *UPLO,
			      int *N, int *KD, int *NRHS, 
			      double *A, int *LDA, double *B, int *LDB, 
			      int *INFO);

extern "C" int  DPBTRS(char *UPLO,
			       int *N, int *KD, int *NRHS, 
			       double *A, int *LDA, double *B, int *LDB, 
			       int *INFO);
#else

extern "C" int dpbsv_(char *UPLO, int *N, int *KD, int *NRHS, 
		      double *A, int *LDA, double *B, int *LDB, 
		      int *INFO);

extern "C" int dpbtrs_(char *UPLO, int *N, int *KD, int *NRHS, 
		       double *A, int *LDA, double *B, int *LDB, 
		       int *INFO);

#endif
		       

int
BandSPDLinLapackSolver::solve(void)
{
    if (theSOE == 0) {
	opserr << "WARNING BandSPDLinLapackSolver::solve(void)- ";
	opserr << " No LinearSOE object has been set\n";
	return -1;
    }

    int n = theSOE->size;
    int kd = theSOE->half_band -1;
    int ldA = kd +1;
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;

    // first copy B into X
    for (int i=0; i<n; i++)
	*(Xptr++) = *(Bptr++);
    Xptr = theSOE->X;

    // now solve AX = Y

	

#ifdef _WIN32
    if (theSOE->factored == false) {
	// factor and solve 	
	unsigned int sizeC = 1;
	DPBSV("U", &n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);	
    }
      else {
	// solve only using factored matrix	  
	unsigned int sizeC = 1;	
	//DPBTRS("U", sizeC, &n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);	
	DPBTRS("U", &n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);
    }
#else	
    { if (theSOE->factored == false)          
	dpbsv_("U",&n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);
      else
	dpbtrs_("U",&n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);
    }
#endif    

    // check if successfull
    if (info != 0) {
	opserr << "WARNING BandSPDLinLapackSolver::solve() - the LAPACK";
	opserr << " routines returned " << info << endln;
	return -info;
    }

    theSOE->factored = true;
    return 0;
}
    


int
BandSPDLinLapackSolver::setSize()
{
    // nothing to do    
    return 0;
}

int
BandSPDLinLapackSolver::sendSelf(int cTag,
				 Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
BandSPDLinLapackSolver::recvSelf(int tag,
				 Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}



