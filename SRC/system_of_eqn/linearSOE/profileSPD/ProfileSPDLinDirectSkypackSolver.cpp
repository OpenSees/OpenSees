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
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSkypackSolver.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 03/98

//
// Description: This file contains the class definition for 
// ProfileSPDLinDirectSkypackSolver. ProfileSPDLinDirectSkypackSolver 
// is a subclass of LinearSOESOlver. It solves a ProfileSPDLinSOE object using
// the Skypack library developed by Osni Marques, software available at
//   http://www.nersc.gov/research/SCG/Osni/marques_software.html

// What: "@(#) ProfileSPDLinDirectSkypackSolver.C, revA"

#include <ProfileSPDLinDirectSkypackSolver.h>
#include <ProfileSPDLinSOE.h>
#include <math.h>

#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Timer.h>

ProfileSPDLinDirectSkypackSolver::ProfileSPDLinDirectSkypackSolver()
:ProfileSPDLinSolver(SOLVER_TAGS_ProfileSPDLinDirectSkypackSolver),
 mCols(0), mRows(0),rw(0),tw(0), index(0),
 size(0), invD(0)
{

}

ProfileSPDLinDirectSkypackSolver::ProfileSPDLinDirectSkypackSolver(int Mcols, int Mrows)
:ProfileSPDLinSolver(SOLVER_TAGS_ProfileSPDLinDirectSkypackSolver),
 mCols(Mcols), mRows(Mrows),rw(0),tw(0), index(0),
 size(0), invD(0)
{
    if (mCols != 0 && mRows != 0) {
	rw = new double[mRows*mCols];
	tw = new double[mRows*mRows];
	if (mCols > mRows)
	    index = new int[mCols];
	else
	    index = new int[mRows];
    } else { // make sure mCols and mRows == 0
	mCols = 0;
	mRows = 0;
    }

    // check we got the space requested
    if (rw == 0 || tw == 0 || index == 0) {
	opserr << "WARNING ProfileSPDLinDirectSkypackSolver::ProfileSPDLinDirectSkypack";
	opserr << "Solver() - ran out of memory for work areas, setting mCols and mRows = 0\n";
	mCols = 0; mRows = 0;
    }
}

    
ProfileSPDLinDirectSkypackSolver::~ProfileSPDLinDirectSkypackSolver()
{
 if (invD != 0) delete [] invD;
 if (rw != 0) delete [] rw;
 if (tw != 0) delete [] tw;
 if (index != 0) delete [] index;
}

int
ProfileSPDLinDirectSkypackSolver::setSize(void)
{
    int result = 0;

    if (theSOE == 0) {
	opserr << "WARNING ProfileSPDLinDirectSkypackSolver::setSize()";
	opserr << " No system has been set\n";
	return -1;
    }

    // check for quick return 
    if (theSOE->size == 0)
	return 0;
    
    size = theSOE->size;
    block[0] = 1;
    block[1] = size;
    block[2] = 1;

    // delete the old work areas
    if (invD != 0) 
	delete [] invD;

    // create the new work areas
    invD = new double[size];
    if (invD == 0) {
      opserr << "Warning :ProfileSPDLinDirectSkypackSolver::setSize():";
      opserr << " ran out of memory for work area invD\n";
      result = -2;;
    }    

    return result;
}

extern "C" int skysf2_(double *A, double *D, int *DGPNT, int *JMIN, int *JMAX); 

extern "C" int skypf2_(int *index, int *JMAX, int *JMIN, int *MCOLS, int *MROWS, int *DGPNT, 
		       double *A, double *RW, double *TW);

extern "C" int skyss_(int *LDX, int *N, int *NRHS, 
		      double *A, double *D, double *X, int *DGPNT,
		      int *BLOCK, int *NBLOCK,
		      char *FNAME, int *FUNIT, int *INFO);

int 
ProfileSPDLinDirectSkypackSolver::solve(void)
{
    // check for quick returns
    if (theSOE == 0) {
	opserr << "ProfileSPDLinDirectSkypackSolver::solve(void): ";
	opserr << " - No ProfileSPDSOE has been assigned\n";
	return -1;
    }
    
    if (theSOE->size == 0)
	return 0;    
    

    // check that work area invD has been created
    if (invD == 0) {
	opserr << "ProfileSPDLinDirectSkypackSolver::solve(void): ";
	opserr << " - no space for invD - has setSize() been called?\n";
	return -1;
    }	

    // set some pointers
    double *A = theSOE->A;
    double *B = theSOE->B;
    double *X = theSOE->X;
    int *iDiagLoc = theSOE->iDiagLoc;
    int theSize = theSOE->size;
    // copy B into X
    for (int ii=0; ii<theSize; ii++)
	X[ii] = B[ii];

    char *FILE = "INCORE";
    
    if (theSOE->isAfactored == false)  {
      
	// FACTOR 
	if (mRows == 0 || mCols == 0) { // factor using skysf2_
	    int JMIN =1;
	    int JMAX = size;
	    skysf2_(A, invD, iDiagLoc, &JMIN, &JMAX); 
	}
	else { // factor using skypf2_
	    int JMIN =1;
	    int JMAX = size;
	    int MCOLS = mCols;
	    int MROWS = mRows;
	    double *RW = rw;
	    double *TW = tw;
	    int *INDEX = index;
	    skypf2_(INDEX, &JMAX, &JMIN, &MCOLS, &MROWS, iDiagLoc, A, RW, TW);

	    // set up invD
	    for (int i=0; i<size; i++)
		invD[i] = 1.0/A[iDiagLoc[i]-1]; // iDiagLoc has fortran array indexing
	}
      
	// mark the system as having been factored
	theSOE->isAfactored = true;
	theSOE->numInt = 0;      
    }

    /*
     * now do the forward and back substitution
     */

    int LDX = theSize;
    int NRHS = 1;
    int numBlock = 1;
    int fileFD = 0;
    int INFO;
    A = theSOE->A;
    X = theSOE->X;
    int *BLOCK = &block[0];
    iDiagLoc = theSOE->iDiagLoc;
    
    skyss_(&LDX, &theSize, &NRHS, A, invD, X, iDiagLoc, BLOCK, &numBlock, 
	   FILE,  &fileFD, &INFO);
      
    // return
    if (INFO < 0) {
	opserr << "WARNING ProfileSPDLinDirectSkypackSolver::solve()";
	opserr << " error value returned from skyss()\n";
    }    

    return INFO;

}


int 
ProfileSPDLinDirectSkypackSolver::setProfileSOE(ProfileSPDLinSOE &theNewSOE)
{
    theSOE = &theNewSOE;
    return 0;
}
	
int
ProfileSPDLinDirectSkypackSolver::sendSelf(int cTag,
					   Channel &theChannel)
{
    if (size != 0)
	opserr << "ProfileSPDLinDirectSkypackSolver::sendSelf - does not send itself YET\n"; 
    return 0;
}


int 
ProfileSPDLinDirectSkypackSolver::recvSelf(int cTag,
					   Channel &theChannel, 
					   FEM_ObjectBroker &theBroker)
{
    return 0;
}


