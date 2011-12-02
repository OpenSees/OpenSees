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
// $Date: 2000-09-15 08:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectBlockSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/ProfileSPD/ProfileSPDLinDirectBlockSolver.C
//
// Written: fmk 
// Created: Mar 1998
// Revision: A
//
// Description: This file contains the class definition for 
// ProfileSPDLinDirectBlockSolver. ProfileSPDLinDirectBlockSolver is a subclass 
// of LinearSOESOlver. It solves a ProfileSPDLinSOE object using
// the LDL^t factorization and a block approach

// What: "@(#) ProfileSPDLinDirectBlockSolver.C, revA"

#include <ProfileSPDLinDirectBlockSolver.h>
#include <ProfileSPDLinSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Timer.h>

ProfileSPDLinDirectBlockSolver::ProfileSPDLinDirectBlockSolver(double tol, int blckSize)
:ProfileSPDLinSolver(SOLVER_TAGS_ProfileSPDLinDirectBlockSolver),
 minDiagTol(tol), blockSize(blckSize), maxColHeight(0), 
 size(0), RowTop(0), topRowPtr(0), invD(0)
{

}

    
ProfileSPDLinDirectBlockSolver::~ProfileSPDLinDirectBlockSolver()
{
    if (RowTop != 0) delete [] RowTop;
    if (topRowPtr != 0) free((void *)topRowPtr);
    if (invD != 0) delete [] invD;
}

int
ProfileSPDLinDirectBlockSolver::setSize(void)
{

    if (theSOE == 0) {
	cerr << "ProfileSPDLinDirectBlockSolver::setSize()";
	cerr << " No system has been set\n";
	return -1;
    }

    // check for quick return 
    if (theSOE->size == 0)
	return 0;
    if (size != theSOE->size) {    
      size = theSOE->size;
    
      if (RowTop != 0) delete [] RowTop;
      if (topRowPtr != 0) delete [] topRowPtr;
      if (invD != 0) delete [] invD;

      RowTop = new int[size];

      // we cannot use topRowPtr = new (double *)[size] with the cxx compiler
      topRowPtr = (double **)malloc(size *sizeof(double *));

      invD = new double[size]; 
	
      if (RowTop == 0 || topRowPtr == 0 || invD == 0) {
	cerr << "Warning :ProfileSPDLinDirectBlockSolver::ProfileSPDLinDirectBlockSolver :";
	cerr << " ran out of memory for work areas \n";
	return -1;
      }
    }


    // set some pointers
    double *A = theSOE->A;
    int *iDiagLoc = theSOE->iDiagLoc;

    // set RowTop and topRowPtr info

    maxColHeight = 0;
    RowTop[0] = 0;
    topRowPtr[0] = A;
    for (int j=1; j<size; j++) {
	int icolsz = iDiagLoc[j] - iDiagLoc[j-1];
        if (icolsz > maxColHeight) maxColHeight = icolsz;
	RowTop[j] = j - icolsz +  1;
	topRowPtr[j] = &A[iDiagLoc[j-1]]; // FORTRAN array indexing in iDiagLoc
    }

    size = theSOE->size;
    return 0;
}


int 
ProfileSPDLinDirectBlockSolver::solve(void)
{
    // check for quick returns
    if (theSOE == 0) {
	cerr << "ProfileSPDLinDirectBlockSolver::solve(void): ";
	cerr << " - No ProfileSPDSOE has been assigned\n";
	return -1;
    }
    
    if (theSOE->size == 0)
	return 0;

    // set some pointers
    double *B = theSOE->B;
    double *X = theSOE->X;
    int n = theSOE->size;

    // copy B into X
    for (int ii=0; ii<n; ii++)
	X[ii] = B[ii];
    
    if (theSOE->isAfactored == false)  {

	// FACTOR 
	invD[0] = 1.0/theSOE->A[0];	
        int startRow = 0;
	int lastRow = startRow+blockSize-1;
        int lastColEffected = lastRow+maxColHeight -1;
        int nBlck = n/blockSize;
	if ((n % blockSize) != 0)
	  nBlck++;

	// for every block across      
	for (int i=0; i<nBlck; i++) {

	  // first factor the diagonal block int Ui,i and Di
          int j;
	  for (j=0; j<blockSize; j++) {
            int currentRow = startRow + j;

	    if (currentRow < n) { // this is for case when size%blockSize != 0

	      int rowjTop = RowTop[currentRow];
	      double *akjPtr = topRowPtr[currentRow];
	      int maxRowijTop;
	      if (rowjTop < startRow) {
		akjPtr += startRow-rowjTop; // pointer to start of block row
		maxRowijTop = startRow;
	      } else
		maxRowijTop = rowjTop;

	      int k;
	      for (k=maxRowijTop; k<currentRow; k++) {
		double tmp = *akjPtr;
		int rowkTop = RowTop[k];
		int maxRowkjTop;
		double *alkPtr, *aljPtr;
		if (rowkTop < rowjTop) {
		  alkPtr = topRowPtr[k] + (rowjTop - rowkTop);
		  aljPtr = topRowPtr[currentRow];
		  maxRowkjTop = rowjTop;
		} else {
		  alkPtr = topRowPtr[k];
		  aljPtr = topRowPtr[currentRow] + (rowkTop - rowjTop);
		  maxRowkjTop = rowkTop;
		}

		for (int l = maxRowkjTop; l<k; l++) 
		  tmp -= *alkPtr++ * *aljPtr++;
		
		*akjPtr++ = tmp;
	      }

	      double ajj = *akjPtr;
	      akjPtr = topRowPtr[currentRow];
	      double *bjPtr  = &X[rowjTop];  
	      double tmp = 0;	    

	      for (k=rowjTop; k<currentRow; k++){
		double akj = *akjPtr;
		double lkj = akj * invD[k];
		tmp -= lkj * *bjPtr++; 		
		*akjPtr++ = lkj;
		ajj = ajj -lkj * akj;
	      }

	      X[currentRow] += tmp;

	      // check that the diag > the tolerance specified
	      if (ajj <= 0.0) {
		cerr << "ProfileSPDLinDirectBlockSolver::solve() - ";
		cerr << " aii < 0 (i, aii): (" << currentRow << ", " << ajj << ")\n"; 
		return(-2);
	      }
	      if (ajj <= minDiagTol) {
		cerr << "ProfileSPDLinDirectBlockSolver::solve() - ";
		cerr << " aii < minDiagTol (i, aii): (" << currentRow;
		cerr << ", " << ajj << ")\n"; 
		return(-2);
	      }		
	      invD[currentRow] = 1.0/ajj; 

	    } else 
	      j = blockSize;

	  }

	  // now do rest of i'th block row doing a block of columns at a time
          // forming Ui,j*Di
	  int currentCol = startRow + blockSize;
	  for (j=i+1; j<nBlck; j++) 
	    for (int k=0; k<blockSize; k++) {

	      if (currentCol < n) { // this is for case when size%blockSize != 0

		int rowkTop = RowTop[currentCol];
		double *alkPtr = topRowPtr[currentCol];
		int maxRowikTop;
		if (rowkTop < startRow) {
		  alkPtr += startRow-rowkTop; // pointer to start of block row
		  maxRowikTop = startRow;
		} else
		  maxRowikTop = rowkTop;

		for (int l=maxRowikTop; l<=lastRow; l++) {
		  double tmp = *alkPtr;
		  int rowlTop = RowTop[l];
		  int maxRowklTop;
		  double *amlPtr, *amkPtr;
		  if (rowlTop < rowkTop) {
		    amlPtr = topRowPtr[l] + (rowkTop - rowlTop);
		    amkPtr = topRowPtr[currentCol];
		    maxRowklTop = rowkTop;
		  } else {
		    amlPtr = topRowPtr[l];
		    amkPtr = topRowPtr[currentCol] + (rowlTop - rowkTop);
		    maxRowklTop = rowlTop;
		  }
		  
		  for (int m = maxRowklTop; m<l; m++) 
		    tmp -= *amkPtr++ * *amlPtr++;
		  
		  *alkPtr++ = tmp;
		}
		currentCol++;
		if (currentCol > lastColEffected) {
		  k = blockSize;
		  j = nBlck;
		}

	      } else
		k = blockSize;

	    }

	  startRow += blockSize;
	  lastRow = startRow + blockSize -1;
          lastColEffected = lastRow + maxColHeight -1;
	}

	theSOE->isAfactored = true;
	theSOE->numInt = 0;
	
	// divide by diag term 
	double *bjPtr = X; 
	double *aiiPtr = invD;
	for (int j=0; j<n; j++) 
	  *bjPtr++ = *aiiPtr++ * X[j];

    
	// now do the back substitution storing result in X
	for (int k=(n-1); k>0; k--) {
      
	  int rowktop = RowTop[k];
	  double bk = X[k];
	  double *ajiPtr = topRowPtr[k]; 		

	  for (int j=rowktop; j<k; j++) 
	    X[j] -= *ajiPtr++ * bk;
	}   	 
    } else { // just do forward and back substitution
      
      // do forward substitution 
      for (int i=1; i<n; i++) {
	    
	int rowitop = RowTop[i];	    
	double *ajiPtr = topRowPtr[i];
	double *bjPtr  = &X[rowitop];  
	double tmp = 0;	    
	    
	for (int j=rowitop; j<i; j++) 
	  tmp -= *ajiPtr++ * *bjPtr++; 
	    
	X[i] += tmp;
      }

      // divide by diag term 
      double *bjPtr = X; 
      double *aiiPtr = invD;
      for (int j=0; j<n; j++) 
	*bjPtr++ = *aiiPtr++ * X[j];

    
      // now do the back substitution storing result in X
      for (int k=(n-1); k>0; k--) {
      
	int rowktop = RowTop[k];
	double bk = X[k];
	double *ajiPtr = topRowPtr[k]; 		

	for (int j=rowktop; j<k; j++) 
	  X[j] -= *ajiPtr++ * bk;
      }   	 
    }    
    return 0;
}

int 
ProfileSPDLinDirectBlockSolver::setProfileSOE(ProfileSPDLinSOE &theNewSOE)
{
    if (theSOE != 0) {
	cerr << "ProfileSPDLinDirectBlockSolver::setProfileSOE() - ";
	cerr << " has already been called \n";	
	return -1;
    }
    
    theSOE = &theNewSOE;
    return 0;
}
	
int
ProfileSPDLinDirectBlockSolver::sendSelf(int cTag, Channel &theChannel)
{
    if (size != 0)
	cerr << "ProfileSPDLinDirectBlockSolver::sendSelf - does not send itself YET\n"; 
    return 0;
}


int 
ProfileSPDLinDirectBlockSolver::recvSelf(int cTag,
					 Channel &theChannel, 
					 FEM_ObjectBroker &theBroker)
{
    return 0;
}


