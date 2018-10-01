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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation for ProfileSPDLinSOE

#include "SkylineSPD.h"
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>


#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

// init function called once on first loading of package

OPS_Export void
localInit() 
{
  OPS_Error("SkylineSPD - Written by fmk UC Berkeley Copyright 2009\n", 1);
}

OPS_Export void *
OPS_SkylineSPD()
{
  //
  // parse args
  //

  // none to parse for this soe

  //
  // create the SOE
  //

  LinearSOE *theSOE = new SkylineSPD();

  // return it

  return theSOE;

}


//

SkylineSPD::SkylineSPD()
  :LinearSOE(0),
   size(0), profileSize(0), A(0), B(0), X(0), vectX(0), vectB(0),
   iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), isAcondensed(false),
   numInt(0) , minDiagTol(1e-18), RowTop(0), topRowPtr(0), invD(0)
{

}

SkylineSPD::~SkylineSPD()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (iDiagLoc != 0) delete [] iDiagLoc;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;    
    if (RowTop != 0) delete [] RowTop;
    if (topRowPtr != 0) free((void*)topRowPtr);
    if (invD != 0) delete [] invD;
}


int 
SkylineSPD::getNumEqn(void) const
{
    return size;
}

int 
SkylineSPD::setSize(Graph &theGraph)
{
    int oldSize = size;
    int result = 0;
    size = theGraph.getNumVertex();

    // check we have enough space in iDiagLoc and iLastCol
    // if not delete old and create new
    if (size > Bsize) { 
      if (iDiagLoc != 0) delete [] iDiagLoc;
      if (RowTop != 0) delete [] RowTop;
      if (topRowPtr != 0) delete [] topRowPtr;
      if (invD != 0) delete [] invD;
      // if (topRowPtr != 0) free((void *)topRowPtr);
      
      iDiagLoc = new int[size];
      RowTop = new int[size];
      invD = new double[size]; 
      //      topRowPtr = new double *[size] ;
      topRowPtr = (double **)malloc(size *sizeof(double *));


      if (iDiagLoc == 0 || RowTop == 0 || topRowPtr == 0 || invD == 0) {
	opserr << "WARNING SkylineSPD::setSize() : ";
	opserr << " - ran out of memory for iDiagLoc\n";
	size = 0; Asize = 0;
	result = -1;
      }
    }
    
    // zero out iDiagLoc 
    for (int i=0; i<size; i++) {
      iDiagLoc[i] = 0;
    }
    
    // now we go through the vertices to find the height of each col and
    // width of each row from the connectivity information.
    
    Vertex *vertexPtr;
    VertexIter &theVertices = theGraph.getVertices();
    
    while ((vertexPtr = theVertices()) != 0) {
      int vertexNum = vertexPtr->getTag();
      const ID &theAdjacency = vertexPtr->getAdjacency();
      int iiDiagLoc = iDiagLoc[vertexNum];
      int *iiDiagLocPtr = &(iDiagLoc[vertexNum]);
      
      for (int i=0; i<theAdjacency.Size(); i++) {
	int otherNum = theAdjacency(i);
	int diff = vertexNum-otherNum;
	if (diff > 0) {
	  if (iiDiagLoc < diff) {
	    iiDiagLoc = diff;
	    *iiDiagLocPtr = diff;
	  }
	} 
      }
    }
    
    
    // now go through iDiagLoc, adding 1 for the diagonal element
    // and then adding previous entry to give current location.
    if (iDiagLoc != 0)
      iDiagLoc[0] = 1; // NOTE FORTRAN ARRAY LOCATION
    
    for (int j=1; j<size; j++)
      iDiagLoc[j] = iDiagLoc[j] + 1 + iDiagLoc[j-1];
    
    if (iDiagLoc != 0)       
      profileSize = iDiagLoc[size-1];
    
    // check if we need more space to hold A
    // if so then go get it
    if (profileSize > Asize) { 
      
      // delete old space
      if (A != 0)
	delete [] A;
      
      // get new space
      A = new double[profileSize];
      
      if (A == 0) {
	opserr << "SkylineSPD::SkylineSPD :";
	opserr << " ran out of memory for A (size,Profile) (";
	opserr << size <<", " << profileSize << ") \n";
	size = 0;  Asize = 0;  profileSize = 0;
	result = -1;
      }
      else 
	Asize = profileSize;
    }
    
    // zero the matrix
    for (int k=0; k<profileSize; k++)
      A[k] = 0;
    
    isAfactored = false;
    isAcondensed = false;    
    
    if (size > Bsize) { // we have to get another space for A
      
      // delete the old
      if (B != 0) delete [] B;
      if (X != 0) delete [] X;
      
      // create the new	
      B = new double[size];
      X = new double[size];
      
      if (B == 0 || X == 0 ) {
	opserr << "SkylineSPD::SkylineSPD :";
	opserr << " ran out of memory for vectors (size) (";
	opserr << size << ") \n";
	size = 0; Bsize = 0;
	result = -1;
      }
    }
    
    // zero the vectors
    for (int l=0; l<size; l++) {
      B[l] = 0;
      X[l] = 0;
    }
    
    if (size != oldSize) {
      
      if (vectX != 0)
	delete vectX;
      if (vectB != 0)
	delete vectB;
      
      vectX = new Vector(X,size);
      vectB = new Vector(B,size);
      
      if (size > Bsize)
	Bsize = size;
    }
    
    RowTop[0] = 0;
    topRowPtr[0] = A;
    for (int j=1; j<size; j++) {
      int icolsz = iDiagLoc[j] - iDiagLoc[j-1];
      RowTop[j] = j - icolsz +  1;
      topRowPtr[j] = &A[iDiagLoc[j-1]]; // FORTRAN array indexing in iDiagLoc
    }
    
    return result;
}

int 
SkylineSPD::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;
    
    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "SkylineSPD::addA()	- Matrix and ID not of similar sizes\n";
	return -1;
    }

    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
		int minColRow;
		if (col == 0)
		    minColRow = 0;
		else
		    minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0 && 
			row <= col && row >= minColRow) { 

			// we only add upper and inside profile
			double *APtr = coliiPtr + (row-col);
			 *APtr += m(j,i);
		     }
		}  // for j
	    } 
	}  // for i
    } else {
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 		
		int minColRow;
		if (col == 0)
		    minColRow = 0;
		else
		    minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0 && 
			row <= col && row >= minColRow) { 

			// we only add upper and inside profile
			double *APtr = coliiPtr + (row-col);
			 *APtr += m(j,i) * fact;
		     }
		}  // for j
	    } 
	}  // for i
    
    }
    return 0;
}

    
int 
SkylineSPD::addB(const Vector &v, const ID &id, double fact)
{
    
    // check for a quick return 
    if (fact == 0.0)  return 0;

    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "SkylineSPD::addB() -";
	opserr << " Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	
    return 0;
}


int
SkylineSPD::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<size; i++) {
	    B[i] = v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<size; i++) {
	    B[i] = -v(i);
	}
    } else {
	for (int i=0; i<size; i++) {
	    B[i] = v(i) * fact;
	}
    }	
    return 0;
}

void 
SkylineSPD::zeroA(void)
{
    double *Aptr = A;
    for (int i=0; i<Asize; i++)
	*Aptr++ = 0;
    
    isAfactored = false;
}
	
void 
SkylineSPD::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}


void 
SkylineSPD::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void 
SkylineSPD::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
SkylineSPD::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL SkylineSPD::getX - vectX == 0";
	exit(-1);
    }    
    return *vectX;
}


const Vector &
SkylineSPD::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL SkylineSPD::getB - vectB == 0";
	exit(-1);
    }        
    return *vectB;
}

double 
SkylineSPD::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    


int 
SkylineSPD::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
SkylineSPD::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}


int 
SkylineSPD::solve(void)
{
  if (size == 0)
    return 0;
  
  // copy B into X
  for (int ii=0; ii<size; ii++)
    X[ii] = B[ii];
  
  if (isAfactored == false)  {

    // FACTOR & SOLVE
    double *ajiPtr, *akjPtr, *akiPtr, *bjPtr;    
    
    // if the matrix has not been factored already factor it into U^t D U
    // storing D^-1 in invD as we go
    
    double a00 = A[0];
    if (a00 <= 0.0) {
      opserr << "ProfileSPDLinDirectSolver::solve() - ";
      opserr << " aii < 0 (i, aii): (0,0)\n"; 
      return(-2);
    }    
    
    invD[0] = 1.0/A[0];	
    
    // for every col across 
    for (int i=1; i<size; i++) {
      
      int rowitop = RowTop[i];
      ajiPtr = topRowPtr[i];
      
      for (int j=rowitop; j<i; j++) {
	double tmp = *ajiPtr;
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
      
      double aii = A[iDiagLoc[i] -1]; // FORTRAN ARRAY INDEXING
      ajiPtr = topRowPtr[i];
      double *bjPtr  = &X[rowitop];  
      double tmp = 0;	    
      
      for (int jj=rowitop; jj<i; jj++) {
	double aji = *ajiPtr;
	double lij = aji * invD[jj];
	tmp -= lij * *bjPtr++; 		
	*ajiPtr++ = lij;
	aii = aii - lij*aji;
      }
      
      // check that the diag > the tolerance specified
      if (aii == 0.0) {
	opserr << "ProfileSPDLinDirectSolver::solve() - ";
	opserr << " aii < 0 (i, aii): (" << i << ", " << aii << ")\n"; 
	return(-2);
      }
      if (fabs(aii) <= minDiagTol) {
	opserr << "ProfileSPDLinDirectSolver::solve() - ";
	opserr << " aii < minDiagTol (i, aii): (" << i;
	opserr << ", " << aii << ")\n"; 
	return(-2);
      }		
      invD[i] = 1.0/aii; 
      X[i] += tmp;	    
    }
    
    isAfactored = true;
    numInt = 0;
    
    // divide by diag term 
    bjPtr = X; 
    double *aiiPtr = invD;
    for (int j=0; j<size; j++) 
      *bjPtr++ = *aiiPtr++ * X[j];
    
    
    // now do the back substitution storing result in X
    for (int k=(size-1); k>0; k--) {
      
      int rowktop = RowTop[k];
      double bk = X[k];
      double *ajiPtr = topRowPtr[k]; 		
      
      for (int j=rowktop; j<k; j++) 
	X[j] -= *ajiPtr++ * bk;
    }   	 	
  }
  
  else {
    
    // JUST DO SOLVE
    
    // do forward substitution 
    for (int i=1; i<size; i++) {
      
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
    for (int j=0; j<size; j++) 
      *bjPtr++ = *aiiPtr++ * X[j];
    
    
    // now do the back substitution storing result in X
    for (int k=(size-1); k>0; k--) {
      
      int rowktop = RowTop[k];
      double bk = X[k];
      double *ajiPtr = topRowPtr[k]; 		
      
      for (int j=rowktop; j<k; j++) 
	X[j] -= *ajiPtr++ * bk;
    }   	 
  }    
    
  return 0;
}
