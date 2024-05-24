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
                                                                        
// Written: fmk 
// Created: 2012
// Revision: A
//
// Description: This file contains the implementation for ProfileSPDLinSOESolver
// Description: This file contains the class definition for 
// ProfileSPDLinDirectSolver. ProfileSPDLinDirectSolver is a subclass 
// of LinearSOESOlver. It solves a ProfileSPDLinSOE object using
// the LDL^t factorization.

#include <SProfileSPDLinSolver.h>
#include <SProfileSPDLinSOE.h>
#include <math.h>
#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>

void* OPS_SProfileSPDLinSolver()
{
    SProfileSPDLinSolver *theSolver = new SProfileSPDLinSolver(); 	
    return new SProfileSPDLinSOE(*theSolver); 
}

SProfileSPDLinSolver::SProfileSPDLinSolver(double tol)
:LinearSOESolver(SOLVER_TAGS_SProfileSPDLinSolver),
 theSOE(0), minDiagTol(tol), size(0), RowTop(0), topRowPtr(0), invD(0)
{

}

    
SProfileSPDLinSolver::~SProfileSPDLinSolver()
{
    if (RowTop != 0) delete [] RowTop;
    if (topRowPtr != 0) free((void *)topRowPtr);
    if (invD != 0) delete [] invD;
}

int
SProfileSPDLinSolver::setSize(void)
{

    if (theSOE == 0) {
	opserr << "SProfileSPDLinSolver::setSize()";
	opserr << " No system has been set\n";
	return -1;
    }

    // check for quick return 
    if (theSOE->size == 0)
	return 0;
    
    size = theSOE->size;
    
    if (RowTop != 0) delete [] RowTop;
    if (topRowPtr != 0) free((void *)topRowPtr);
    if (invD != 0) delete [] invD;

    RowTop = new int[size];

    // we cannot use topRowPtr = new (double *)[size] with the cxx compiler
    topRowPtr = (float **)malloc(size *sizeof(double *));

    invD = new float[size]; 
	
    if (RowTop == 0 || topRowPtr == 0 || invD == 0) {
	opserr << "Warning :SProfileSPDLinSolver::SProfileSPDLinSolver :";
	opserr << " ran out of memory for work areas \n";
	return -1;
    }


    // set some pointers
    float *A = theSOE->A;
    int *iDiagLoc = theSOE->iDiagLoc;

    // set RowTop and topRowPtr info

    RowTop[0] = 0;
    topRowPtr[0] = A;
    for (int j=1; j<size; j++) {
	int icolsz = iDiagLoc[j] - iDiagLoc[j-1];
	RowTop[j] = j - icolsz +  1;
	topRowPtr[j] = &A[iDiagLoc[j-1]]; // FORTRAN array indexing in iDiagLoc
    }

    size = theSOE->size;
    return 0;
}


int 
SProfileSPDLinSolver::solve(void)
{
  Timer timer;
  timer.start();
    // check for quick returns
    if (theSOE == 0) {
	opserr << "SProfileSPDLinSolver::solve(void): ";
	opserr << " - No ProfileSPDSOE has been assigned\n";
	return -1;
    }
    
    if (theSOE->size == 0)
	return 0;

    // set some pointers
    float *B = theSOE->B;
    float *X = theSOE->X;
    double *doubleB = theSOE->doubleB;
    double *doubleX = theSOE->doubleX;
    int theSize = theSOE->size;
    // copy B into X
    for (int ii=0; ii<theSize; ii++)  {
	X[ii] = B[ii];
	doubleB[ii] = B[ii];
    }

    /*
      for (int iii=0; iii<theSize; iii++) {
      int rowiiitop = RowTop[iii];
      float *ajiptr = topRowPtr[iii];
      opserr << "\n COLUMN " << iii << " TopRow " << rowiiitop << " -> ";
      for (int jjj = rowiiitop; jjj <=iii; jjj++)
      opserr << *ajiptr++ << " ";
      }
      opserr << endln;
      
      for (int iii=0; iii<theSOE->size; iii++) {
      opserr << "COLUMN " << iii << " Biii -> " << X[iii] << endln;
      }
      opserr << endln;
    */
    
    
    if (theSOE->isAfactored == false)  {
      
      // FACTOR & SOLVE
      float *ajiPtr, *akjPtr, *akiPtr, *bjPtr;    
      
      // if the matrix has not been factored already factor it into U^t D U
      // storing D^-1 in invD as we go
      
      float a00 = theSOE->A[0];
      if (a00 <= 0.0) {
	opserr << "SProfileSPDLinSolver::solve() - ";
	opserr << " aii < 0 (i, aii): (0,0)\n"; 
	return(-2);
      }    
      
      invD[0] = 1.0/theSOE->A[0];	
      
      // for every col across 
      for (int i=1; i<theSize; i++) {
	
	int rowitop = RowTop[i];
	ajiPtr = topRowPtr[i];
	
	for (int j=rowitop; j<i; j++) {
	  float tmp = *ajiPtr;
	  int rowjtop = RowTop[j];
	  
	  if (rowitop > rowjtop) {
	    
	    akjPtr = topRowPtr[j] + (rowitop-rowjtop);
	    akiPtr = topRowPtr[i];
	    
	    for (int k=rowitop; k<j; k++) 
	      tmp -= *akjPtr++ * *akiPtr++ ;
	    
	    *ajiPtr++ = tmp;
	  }
	  else {
	    akjPtr = topRowPtr[j];
	    akiPtr = topRowPtr[i] + (rowjtop-rowitop);
	    
	    for (int k=rowjtop; k<j; k++) 
	      tmp -= *akjPtr++ * *akiPtr++ ;
	    
	    *ajiPtr++ = tmp;
	  }
	}
	
	/* now form i'th col of [U] and determine [dii] */
	
	float aii = theSOE->A[theSOE->iDiagLoc[i] -1]; // FORTRAN ARRAY INDEXING
	ajiPtr = topRowPtr[i];
	float *bjPtr  = &X[rowitop];  
	float tmp = 0;	    
	
	for (int jj=rowitop; jj<i; jj++) {
	  float aji = *ajiPtr;
	  float lij = aji * invD[jj];
	  tmp -= lij * *bjPtr++; 		
	  *ajiPtr++ = lij;
	  aii = aii - lij*aji;
	}
	
	// check that the diag > the tolerance specified
	if (aii == 0.0) {
	  opserr << "SProfileSPDLinSolver::solve() - ";
	  opserr << " aii < 0 (i, aii): (" << i << ", " << aii << ")\n"; 
	  return(-2);
	}
	if (fabs(aii) <= minDiagTol) {
	  opserr << "SProfileSPDLinSolver::solve() - ";
	  opserr << " aii < minDiagTol (i, aii): (" << i;
	  opserr << ", " << aii << ")\n"; 
	  return(-2);
	}		
	invD[i] = 1.0/aii; 
	X[i] += tmp;	    
      }
      
      theSOE->isAfactored = true;
      theSOE->numInt = 0;
      
      
      // divide by diag term 
      bjPtr = X; 
      float *aiiPtr = invD;
      for (int j=0; j<theSize; j++) 
	*bjPtr++ = *aiiPtr++ * X[j];
      
      
      // now do the back substitution storing result in X
      for (int k=(theSize-1); k>0; k--) {
	
	int rowktop = RowTop[k];
	float bk = X[k];
	float *ajiPtr = topRowPtr[k]; 		
	
	for (int j=rowktop; j<k; j++) 
	  X[j] -= *ajiPtr++ * bk;
      }   	 	
    }
    
    else {
      
      // JUST DO SOLVE
      
      // do forward substitution 
      for (int i=1; i<theSize; i++) {
	
	int rowitop = RowTop[i];	    
	float *ajiPtr = topRowPtr[i];
	float *bjPtr  = &X[rowitop];  
	float tmp = 0;	    
	
	for (int j=rowitop; j<i; j++) 
	  tmp -= *ajiPtr++ * *bjPtr++; 
	
	X[i] += tmp;
      }
      
      // divide by diag term 
      float *bjPtr = X; 
      float *aiiPtr = invD;
      for (int j=0; j<theSize; j++) 
	*bjPtr++ = *aiiPtr++ * X[j];
      
      
      // now do the back substitution storing result in X
      for (int k=(theSize-1); k>0; k--) {
	
	int rowktop = RowTop[k];
	float bk = X[k];
	float *ajiPtr = topRowPtr[k]; 		
	
	for (int j=rowktop; j<k; j++) 
	  X[j] -= *ajiPtr++ * bk;
      }   	 
    }  
  
    for (int i=0; i<size; i++) {
      doubleX[i] = X[i];
    }

    /*
      opserr << "BBBB " << theSOE->getB();
      opserr << "XXXX " << theSOE->getX();
    */

  timer.pause();
  timer.Print(opserr);
    return 0;
}

double
SProfileSPDLinSolver::getDeterminant(void) 
{
   int theSize = theSOE->size;
   float determinant = 1.0;
   for (int i=0; i<theSize; i++)
     determinant *= invD[i];
   determinant = 1.0/determinant;
   return determinant;
}

int 
SProfileSPDLinSolver::setLinearSOE(SProfileSPDLinSOE &theNewSOE)
{
    if (theSOE != 0) {
	opserr << "SProfileSPDLinSolver::setProfileSOE() - ";
	opserr << " has already been called \n";	
	return -1;
    }
    
    theSOE = &theNewSOE;
    return 0;
}
	

int 
SProfileSPDLinSolver::factor(int n)
{

    // check for quick returns
    if (theSOE == 0) {
	opserr << "SProfileSPDLinSolver::factor: ";
	opserr << " - No ProfileSPDSOE has been assigned\n";
	return -1;
    }

    int theSize = theSOE->size;    
    if (n > theSize) {
	opserr << "SProfileSPDLinSolver::factor: ";
	opserr << " - n " << n << " greater than size of system" << theSize << endln;
	return -1;
    }

    if (theSize == 0 || n == 0)
	return 0;


    // set some pointers
    if (theSOE->isAfactored == false)  {

	// FACTOR & SOLVE
	float *ajiPtr, *akjPtr, *akiPtr;    
	
	// if the matrix has not been factored already factor it into U^t D U
	// storing D^-1 in invD as we go
    
	invD[0] = 1.0/theSOE->A[0];	
	
	// for every col across 
	for (int i=1; i<n; i++) {

	    int rowitop = RowTop[i];
	    ajiPtr = topRowPtr[i];

	    for (int j=rowitop; j<i; j++) {
		float tmp = *ajiPtr;
		int rowjtop = RowTop[j];

		if (rowitop > rowjtop) {

		    akjPtr = topRowPtr[j] + (rowitop-rowjtop);
		    akiPtr = topRowPtr[i];

		    for (int k=rowitop; k<j; k++) 
			tmp -= *akjPtr++ * *akiPtr++ ;

		    *ajiPtr++ = tmp;
		}
		else {
		    akjPtr = topRowPtr[j];
		    akiPtr = topRowPtr[i] + (rowjtop-rowitop);

		    for (int k=rowjtop; k<j; k++) 
			tmp -= *akjPtr++ * *akiPtr++ ;

		    *ajiPtr++ = tmp;
		}
	    }

	    /* now form i'th col of [U] and determine [dii] */

	    float aii = theSOE->A[theSOE->iDiagLoc[i] -1]; // FORTRAN ARRAY INDEXING
	    ajiPtr = topRowPtr[i];
	    
	    for (int jj=rowitop; jj<i; jj++) {
		float aji = *ajiPtr;
		float lij = aji * invD[jj];
		*ajiPtr++ = lij;
		aii = aii - lij*aji;
	    }
	    
	    // check that the diag > the tolerance specified
	    if (aii <= 0.0) {
		opserr << "SProfileSPDLinSolver::solve() - ";
		opserr << " aii < 0 (i, aii): (" << i << ", " << aii << ")\n"; 
		return(-2);
	    }
	    if (aii <= minDiagTol) {
		opserr << "SProfileSPDLinSolver::solve() - ";
		opserr << " aii < minDiagTol (i, aii): (" << i;
		opserr << ", " << aii << ")\n"; 
		return(-2);
	    }		
	    

	    invD[i] = 1.0/aii; 
	}

	theSOE->isAfactored = true;
	theSOE->numInt = n;
	
    }	
    return 0;
}


/*
    float *ajiPtr, *akjPtr, *akiPtr, *bjPtr;    
	
    // if the matrix has not been factored already factor it into U^t D U
    // storing D^-1 in invD as we go
    
    if (theSOE->isAfactored == false) {
	
	invD[0] = 1.0/theSOE->A[0];	
	
	// for every col across 
	for (int i=1; i<n; i++) {

	    int rowitop = RowTop[i];
	    ajiPtr = topRowPtr[i];

	    for (int j=rowitop; j<i; j++) {
		float tmp = *ajiPtr;
		int rowjtop = RowTop[j];

		if (rowitop > rowjtop) {

		    akjPtr = topRowPtr[j] + (rowitop-rowjtop);
		    akiPtr = topRowPtr[i];

		    for (int k=rowitop; k<j; k++) 
			tmp -= *akjPtr++ * *akiPtr++ ;

		    *ajiPtr++ = tmp;
		}
		else {
		    akjPtr = topRowPtr[j];
		    akiPtr = topRowPtr[i] + (rowjtop-rowitop);

		    for (int k=rowjtop; k<j; k++) 
			tmp -= *akjPtr++ * *akiPtr++ ;

		    *ajiPtr++ = tmp;
		}
	    }

	    // now form i'th col of [U] and determine [dii] 

	    float aii = theSOE->A[theSOE->iDiagLoc[i] -1]; // FORTRAN ARRAY INDEXING
	    ajiPtr = topRowPtr[i];
	    
	    for (int jj=rowitop; jj<i; jj++) {
		float aji = *ajiPtr;
		float lij = aji * invD[jj];
		*ajiPtr++ = lij;
		aii = aii - lij*aji;
	    }
	    
	    // check that the diag > the tolerance specified
	    if (aii <= 0.0)     return(-2);
	    if (aii <= minDiagTol) return(-1);
	    invD[i] = 1.0/aii; 
	}
    }
    
    theSOE->isAfactored = true;
    if (n == theSize)
	theSOE->numInt = 0;
    else
	theSOE->numInt = n;
    
    return 0;
}
*/

int
SProfileSPDLinSolver::sendSelf(int cTag,
				    Channel &theChannel)
{
    return 0;
}


int 
SProfileSPDLinSolver::recvSelf(int cTag,
				    Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
    return 0;
}


