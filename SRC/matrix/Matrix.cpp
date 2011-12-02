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
// $Date: 2000-12-12 07:14:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/Matrix.cpp,v $
                                                                        
                                                                        
// File: ~/matrix/Matrix.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for Matrix.
//
// What: "@(#) Matrix.h, revA"

#include "Matrix.h"
#include "Vector.h"
#include "ID.h"

#include <stdlib.h>

#define MATRIX_WORK_AREA 400

double Matrix::MATRIX_NOT_VALID_ENTRY =0.0;
double *Matrix::matrixWork = new double [MATRIX_WORK_AREA];
//double *Matrix::matrixWork = (double *)malloc(400*sizeof(double));

//
// CONSTRUCTORS
//

Matrix::Matrix()
:numRows(0), numCols(0),dataSize(0),data(0),fromFree(0)
{
    // does nothing
}


Matrix::Matrix(int nRows,int nCols)
:numRows(nRows),numCols(nCols),fromFree(0)
{
#ifdef _G3DEBUG
    if (nRows < 0) {
      cerr << "WARNING: Matrix::Matrix(int,int): tried to init matrix ";
      cerr << "with num rows: " << nRows << " <0\n";
      numRows = 0; numCols =0; dataSize =0; data = 0;
    }
    if (nCols < 0) {
      cerr << "WARNING: Matrix::Matrix(int,int): tried to init matrix";
      cerr << "with num cols: " << nCols << " <0\n";
      numRows = 0; numCols =0; dataSize =0; data = 0;
    }
#endif
    dataSize = numRows * numCols;
    data = 0;

    if (dataSize > 0) {
      data = new double[dataSize];
      //data = (double *)malloc(dataSize*sizeof(double));
      if (data == 0) {
	cerr << "WARNING:Matrix::Matrix(int,int): Ran out of memory on init ";
	cerr << "of size " << dataSize << "\n";
	numRows = 0; numCols =0; dataSize =0;
      } else {
	// zero the data
	double *dataPtr = data;
	for (int i=0; i<dataSize; i++)
	  *dataPtr++ = 0.0;
      }
    }
}

Matrix::Matrix(double *theData, int row, int col) 
:numRows(row),numCols(col),dataSize(row*col),data(theData),fromFree(1)
{
#ifdef _G3DEBUG
    if (row < 0) {
      cerr << "WARNING: Matrix::Matrix(int,int): tried to init matrix with numRows: ";
      cerr << row << " <0\n";
      numRows = 0; numCols =0; dataSize =0; data = 0;
    }
    if (col < 0) {
      cerr << "WARNING: Matrix::Matrix(int,int): tried to init matrix with numCols: ";
      cerr << col << " <0\n";
      numRows = 0; numCols =0; dataSize =0; data = 0;
    }    
#endif

    // does nothing
}

Matrix::Matrix(const Matrix &other)
:fromFree(0)
{
    numRows = other.numRows;
    numCols = other.numCols;
    dataSize = other.dataSize;
    data = 0;
    
    if (dataSize != 0) {
      data = new double[dataSize];
      // data = (double *)malloc(dataSize*sizeof(double));
      if (data == 0) {
	cerr << "WARNING:Matrix::Matrix(Matrix &): ";
	cerr << "Ran out of memory on init of size " << dataSize << "\n"; 
	numRows = 0; numCols =0; dataSize = 0;
      } else {
	// copy the data
	double *dataPtr = data;
	double *otherDataPtr = other.data;
	for (int i=0; i<dataSize; i++)
	  *dataPtr++ = *otherDataPtr++;
      }
    }
}


//
// DESTRUCTOR
//

Matrix::~Matrix()
{
  if (data != 0) 
      if (fromFree == 0)
	  delete [] data; 
  //  if (data != 0) free((void *) data);
}
    

//
// METHODS - Zero, Assemble, Solve
//

void
Matrix::Zero(void)
{
  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ = 0;
}


int
Matrix::resize(int rows, int cols) {

  int newSize = rows*cols;

  if (newSize <= 0) {
    g3ErrorHandler->warning("Matrix::resize) - rows %d or cols %d specified <= 0\n",
			    rows, cols); 
    return -1;
  }

  else if (newSize > dataSize) {

    // free the old space
    if (data != 0) 
      if (fromFree == 0)
	delete [] data; 
    //  if (data != 0) free((void *) data);

    fromFree = 0;
    // create new space
    data = new double[newSize];
    // data = (double *)malloc(dataSize*sizeof(double));
    if (data == 0) {
      g3ErrorHandler->warning("Matrix::resize(%d %d) - out of memory\n",
			      rows, cols); 
      numRows = 0; numCols =0; dataSize = 0;
      return -2;
    }
    dataSize = newSize;
    numRows = rows;
    numCols = cols;
  }

  // just reset the cols and rows - save two memory calls at expense of holding 
  // onto extra memory
  else {
    numRows = rows;
    numCols = cols;
  }

  return 0;
}





int
Matrix::Assemble(const Matrix &V, const ID &rows, const ID &cols, double fact) 
{
  int pos_Rows, pos_Cols;
  int res = 0;

  for (int i=0; i<cols.Size(); i++) {
    pos_Cols = cols(i);
    for (int j=0; j<rows.Size(); j++) {
      pos_Rows = rows(j);
      
      if ((pos_Cols >= 0) && (pos_Rows >= 0) && (pos_Rows < numRows) &&
	  (pos_Cols < numCols) && (i < V.numCols) && (j < V.numRows))
	(*this)(pos_Rows,pos_Cols) += V(j,i)*fact;
      else {
	cerr << "WARNING: Matrix::Assemble(const Matrix &V, const ID &l): ";
	cerr << " - position (" << pos_Rows << "," << pos_Cols << ") outside bounds \n";
	res = -1;
      }
    }
  }

  return res;
}



int
Matrix::Solve(const Vector &b, Vector &x) const
{

    int n = numRows;

#ifdef _G3DEBUG    
    if (numRows != numCols) {
      g3ErrorHandler->warning("Matrix::Solve(b,x) - %s [%d %d] is not square\n",
			      "the matrix of dimensions",
			      numRows, numCols);
      return -1;
    }

    if (n != x.Size()) {
      g3ErrorHandler->warning("Matrix::Solve(b,x) - dimension of x, %d, %s %d\n",
			      numRows, "is not same as matrix", x.Size());
      return -2;
    }

    if (n != b.Size()) {
      g3ErrorHandler->warning("Matrix::Solve(b,x) - dimension of b, %d, %s %d\n",
			      numRows, "is not same as matrix", b.Size());
      return -2;
    }
#endif
    
    // check work area can hold all the data
    if (dataSize > MATRIX_WORK_AREA) {
       g3ErrorHandler->warning("Matrix::Solve(b,x) - matrix dimension [%d %d] %s %d\n",
			       numRows, numCols, "larger than work area",
			       MATRIX_WORK_AREA);
       return -3;
    }

    // copy the data
    int i;
    for (i=0; i<dataSize; i++)
      matrixWork[i] = data[i];

    // set x equal to b
    x = b;

    double aii,lji,xi;
    for (i=0; i<n-1; i++) {
      int j;
      double *aiiPtr = &matrixWork[i+i*numRows];
      aii = *aiiPtr;
      if (aii == 0.0) {
	g3ErrorHandler->warning("Matrix::Solve(b,x) - 0 %s %d\n",
				"diagonal in factorization at column", i);
				
	return -4;
      }

      xi = x(i);	

      // compute column i of L, store in A
      // do forward substution
      double *ajiPtr = aiiPtr+1;
      for (j= i+1; j<n; j++) {
	lji = *ajiPtr / aii;
	x(j) -= lji * xi;
	*ajiPtr++ = lji;
      }
	
      // update rest of A
      //    A(j,k) -= A(j,i)*A(i,k);

      for (int k=i+1; k<n; k++) {
	ajiPtr = aiiPtr +1;
	double *ajkPtr = &matrixWork[i+1+k*numRows];
	double aik = matrixWork[i + k*numRows];
	for (j=i+1; j<n; j++) {
	  *ajkPtr++ -= *ajiPtr++ * aik;
	}
      }
	  

    }
    
    // now do backward substitution
    for (i=n-1; i>=0; i--) {
      double *aiiPtr = &matrixWork[i+i*numRows];
      xi = x(i) /  *aiiPtr--;
      x(i) = xi;
      for (int j = i-1; j>=0; j--)
	x(j) -= *aiiPtr-- * xi;
    }

    return 0;
}



int
Matrix::Solve(const Matrix &b, Matrix &x) const
{

    int n = numRows;
    int nrhs = x.numCols;

#ifdef _G3DEBUG    
    if (numRows != numCols) {
      g3ErrorHandler->warning("Matrix::Solve(B,X) - the matrix of dimensions [%d %d] is not square\n",
			      numRows, numCols);
      return -1;
    }

    if (n != x.numRows) {
      g3ErrorHandler->warning("Matrix::Solve(B,X) - #rows of X, %d, is not same as matrix %d\n",
			      numRows, x.numRows);
      return -2;
    }

    if (n != b.numRows) {
      g3ErrorHandler->warning("Matrix::Solve(B,X) - #rows of B, %d, is not same as matrix %d\n",
			      numRows, b.numRows);
      return -2;
    }

    if (x.numCols != b.numCols) {
      g3ErrorHandler->warning("Matrix::Solve(B,X) - #cols of B, %d, is not same as that of X, %d\n",
			      b.numCols, x.numCols);
      return -3;
    }
#endif
    
    // check work area can hold all the data
    if (dataSize > MATRIX_WORK_AREA) {
      g3ErrorHandler->warning("Matrix::Solve(b,x) - matrix dimension [%d %d] larger than work area %d\n",
			      numRows, numCols, MATRIX_WORK_AREA);
      return -3;
    }

    // copy the data
    int i;
    for (i=0; i<dataSize; i++)
      matrixWork[i] = data[i];

    // set x equal to b
    x = b;

    double aii,lji,xi;

    for (i=0; i<n-1; i++) {
      int j;
      double *aiiPtr = &matrixWork[i+i*numRows];
      aii = *aiiPtr;
      if (aii == 0.0) {
	g3ErrorHandler->warning("Matrix::Solve(b,x) - 0 diagonal in factorization at column %d\n",
				i);
	return -4;
      }

      // compute column i of L, store in A
      // do forward substution
      double *ajiPtr = aiiPtr+1;
      for (j= i+1; j<n; j++) {
	lji = *ajiPtr / aii;
	for (int k=0; k<nrhs; k++)
	  x(j,k) -= lji * x(i,k);
	*ajiPtr++ = lji;
      }
	
      // update rest of A
      //    A(j,k) -= A(j,i)*A(i,k);

      for (int k=i+1; k<n; k++) {
	ajiPtr = aiiPtr +1;
	double *ajkPtr = &matrixWork[i+1+k*numRows];
	double aik = matrixWork[i + k*numRows];
	for (j=i+1; j<n; j++) {
	  *ajkPtr++ -= *ajiPtr++ * aik;
	}
      }
    }
    
    // now do backward substitution
    for (i=n-1; i>=0; i--) {
      double *aiiPtrA = &matrixWork[i+i*numRows];
      for (int k=0; k<nrhs; k++) {
	double *aiiPtr = aiiPtrA;
	xi = x(i,k) /  *aiiPtr--;
	x(i,k) = xi;
	for (int j = i-1; j>=0; j--)
	  x(j,k) -= *aiiPtr-- * xi;
      }
    }

    return 0;
}
    
		    

int
Matrix::addMatrix(double factThis, const Matrix &other, double factOther)
{
    if (factThis == 1.0 && factOther == 0.0)
      return 0;

#ifdef _G3DEBUG
    if ((other.numRows != numRows) || (other.numCols != numCols)) {
      g3ErrorHandler->warning("Matrix::addMatrix(): incompatable matrices, this[%d %d] other[%d %d]\n",
			      numRows, numCols, other.numRows, other.numCols);
      return -1;
    }
#endif

    if (factThis == 1.0) {

      // want: this += other * factOther
      if (factOther == 1.0) {
	double *dataPtr = data;
	double *otherDataPtr = other.data;		    
	for (int i=0; i<dataSize; i++)
	  *dataPtr++ += *otherDataPtr++;
      } else {
	double *dataPtr = data;
	double *otherDataPtr = other.data;		    
	for (int i=0; i<dataSize; i++)
	  *dataPtr++ += *otherDataPtr++ * factOther;
      }
    } 

    else if (factThis == 0.0) {

      // want: this = other * factOther
      if (factOther == 1.0) {
	double *dataPtr = data;
	double *otherDataPtr = other.data;		    
	for (int i=0; i<dataSize; i++)
	  *dataPtr++ = *otherDataPtr++;
      } else {
	double *dataPtr = data;
	double *otherDataPtr = other.data;		    
	for (int i=0; i<dataSize; i++)
	  *dataPtr++ = *otherDataPtr++ * factOther;
      }
    } 

    else {

      // want: this = this * thisFact + other * factOther
      if (factOther == 1.0) {
	double *dataPtr = data;
	double *otherDataPtr = other.data;		    
	for (int i=0; i<dataSize; i++) {
	  double value = *dataPtr * factThis + *otherDataPtr++;
	  *dataPtr++ = value;
	}
      } else {
	double *dataPtr = data;
	double *otherDataPtr = other.data;		    
	for (int i=0; i<dataSize; i++) {
	  double value = *dataPtr * factThis + *otherDataPtr++ * factOther;
	  *dataPtr++ = value;
	}
      }
    } 

    // successfull
    return 0;
}



int
Matrix::addMatrixProduct(double thisFact, 
			 const Matrix &B, 
			 const Matrix &C, 
			 double otherFact)
{
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;
#ifdef _G3DEBUG
    if ((B.numRows != numRows) || (C.numCols != numCols) || (B.numCols != C.numRows)) {
      g3ErrorHandler->warning("Matrix::addMatrixProduct():%s[%d %d] A[%d &d] B{%d %d]\n", 
			      "incompatable matrices, this",
			      numRows, numCols, B.numRows, B.numCols, C.numRows, C.numCols);
      return -1;
    }
#endif
    // NOTE: looping as per blas3 dgemm_: j,k,i
    if (thisFact == 1.0) {

      // want: this += B * C  otherFact
      int numColB = B.numCols;
      double *ckjPtr  = &(C.data)[0];
      for (int j=0; j<numCols; j++) {
	double *aijPtrA = &data[j*numRows];
	for (int k=0; k<numColB; k++) {
	  double tmp = *ckjPtr++ * otherFact;
	  double *aijPtr = aijPtrA;
	  double *bikPtr = &(B.data)[k*numRows];
	  for (int i=0; i<numRows; i++)
	    *aijPtr++ += *bikPtr++ * tmp;
	}
      }
    }

    else if (thisFact == 0.0) {

      // want: this = B * C  otherFact
      double *dataPtr = data;
      for (int i=0; i<dataSize; i++)
	  *dataPtr++ = 0.0;
      int numColB = B.numCols;
      double *ckjPtr  = &(C.data)[0];
      for (int j=0; j<numCols; j++) {
	double *aijPtrA = &data[j*numRows];
	for (int k=0; k<numColB; k++) {
	  double tmp = *ckjPtr++ * otherFact;
	  double *aijPtr = aijPtrA;
	  double *bikPtr = &(B.data)[k*numRows];
	  for (int i=0; i<numRows; i++)
	    *aijPtr++ += *bikPtr++ * tmp;
	}
      }
    } 

    else {
      // want: this = B * C  otherFact
      double *dataPtr = data;
      for (int i=0; i<dataSize; i++)
	  *dataPtr++ *= thisFact;
      int numColB = B.numCols;
      double *ckjPtr  = &(C.data)[0];
      for (int j=0; j<numCols; j++) {
	double *aijPtrA = &data[j*numRows];
	for (int k=0; k<numColB; k++) {
	  double tmp = *ckjPtr++ * otherFact;
	  double *aijPtr = aijPtrA;
	  double *bikPtr = &(B.data)[k*numRows];
	  for (int i=0; i<numRows; i++)
	    *aijPtr++ += *bikPtr++ * tmp;
	}
      }
    } 

    return 0;
}


// to perform this += T' * B * T
int
Matrix::addMatrixTripleProduct(double thisFact, 
			       const Matrix &T, 
			       const Matrix &B, 
			       double otherFact)
{
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;
#ifdef _G3DEBUG
    if ((numCols != numRows) || (B.numCols != B.numRows) || (T.numCols != numRows) ||
	(T.numRows != B.numCols)) {
      g3ErrorHandler->warning("Matrix::addMatrixTripleProduct():%s[%d %d] T[%d %d] B[%d %d]\n", 
			      "incompatable matrices, this",
			      numRows, numCols, T.numRows, T.numCols, B.numRows, B.numCols);
      return -1;
    }
#endif

    // cheack work area can hold the temporary matrix
    int dimB = B.numCols;
    int sizeWork = dimB * numCols;

    if (sizeWork > MATRIX_WORK_AREA) {
      this->addMatrix(thisFact, T^B*T, otherFact);
      return 0;
    }

    // zero out the work area
    double *matrixWorkPtr = matrixWork;
    for (int l=0; l<sizeWork; l++)
      *matrixWorkPtr++ = 0.0;

    // now form B * T * fact store in matrixWork == A area
    // NOTE: looping as per blas3 dgemm_: j,k,i

    double *tkjPtr  = &(T.data)[0];
    for (int j=0; j<numCols; j++) {
      double *aijPtrA = &matrixWork[j*dimB];
      for (int k=0; k<dimB; k++) {
	double tmp = *tkjPtr++ * otherFact;
	double *aijPtr = aijPtrA;
	double *bikPtr = &(B.data)[k*dimB];
	for (int i=0; i<dimB; i++) 
	  *aijPtr++ += *bikPtr++ * tmp;
      }
    }

    // now form T' * matrixWork
    // NOTE: looping as per blas3 dgemm_: j,i,k
    if (thisFact == 1.0) {
      double *dataPtr = &data[0];
      for (int j=0; j< numCols; j++) {
	double *workkjPtrA = &matrixWork[j*dimB];
	for (int i=0; i<numRows; i++) {
	  double *ckiPtr = &(T.data)[i*dimB];
	  double *workkjPtr = workkjPtrA;
	  double aij = 0.0;
	  for (int k=0; k< dimB; k++)
	    aij += *ckiPtr++ * *workkjPtr++;
	  *dataPtr++ += aij;
	}
      }
    } else if (thisFact == 0.0) {
      double *dataPtr = &data[0];
      for (int j=0; j< numCols; j++) {
	double *workkjPtrA = &matrixWork[j*dimB];
	for (int i=0; i<numRows; i++) {
	  double *ckiPtr = &(T.data)[i*dimB];
	  double *workkjPtr = workkjPtrA;
	  double aij = 0.0;
	  for (int k=0; k< dimB; k++)
	    aij += *ckiPtr++ * *workkjPtr++;
	  *dataPtr++ = aij;
	}
      }

    } else {
      double *dataPtr = &data[0];
      for (int j=0; j< numCols; j++) {
	double *workkjPtrA = &matrixWork[j*dimB];
	for (int i=0; i<numRows; i++) {
	  double *ckiPtr = &(T.data)[i*dimB];
	  double *workkjPtr = workkjPtrA;
	  double aij = 0.0;
	  for (int k=0; k< dimB; k++)
	    aij += *ckiPtr++ * *workkjPtr++;
	  double value = *dataPtr * thisFact + aij;
	  *dataPtr++ = value;
	}
      }
    }

    return 0;
}



//
// OVERLOADED OPERATOR () to CONSTRUCT A NEW MATRIX
//

Matrix
Matrix::operator()(const ID &rows, const ID & cols) const
{
    int nRows, nCols;
    nRows = rows.Size();
    nCols = cols.Size();
    Matrix result(nRows,nCols);
    double *dataPtr = result.data;
    for (int i=0; i<nCols; i++)
	for (int j=0; j<nRows; j++)
	    *dataPtr++ = (*this)(rows(j),cols(i));

    return result;
}
		
// Matrix &operator=(const Matrix  &V):
//      the assignment operator, This is assigned to be a copy of V. if sizes
//      are not compatable this.data [] is deleted. The data pointers will not
//      point to the same area in mem after the assignment.
//



Matrix &
Matrix::operator=(const Matrix &other)
{
  // first check we are not trying other = other
  if (this == &other) 
    return *this;

/*
#ifdef _G3DEBUG    
  if ((numCols != other.numCols) || (numRows != other.numRows)) {
    g3ErrorHandler->warning("Matrix::operator=() - matrix dimensions do not match: [%d %d] != [%d %d]\n",
			    numRows, numCols, other.numRows, other.numCols);
    return *this;
  }
#endif
*/
  if ((numCols != other.numCols) || (numRows != other.numRows)) {
#ifdef _G3DEBUG    
      g3ErrorHandler->warning("Matrix::operator=() - matrix dimensions do not match: [%d %d] != [%d %d]\n",
			      numRows, numCols, other.numRows, other.numCols);
#endif
      if (this->data != 0)
	  delete [] this->data;
      
      int theSize = other.numCols*other.numRows;
      
      data = new double[theSize];
      
      this->dataSize = theSize;
      this->numCols = other.numCols;
      this->numRows = other.numRows;
  }


  // now copy the data
  double *dataPtr = data;
  double *otherDataPtr = other.data;		    
  for (int i=0; i<dataSize; i++)
      *dataPtr++ = *otherDataPtr++;
  
  return *this;
}





Matrix &
Matrix::operator=(const Tensor &V)
{
  int rank = V.rank();
  if (rank != 4) {
      g3ErrorHandler->warning("Matrix::operator=() - tensor must be of rank 4\n");
      return *this;
  }
  int dim = V.dim(1);
  if (dim != V.dim(2) != V.dim(3) != V.dim(4)) {
      g3ErrorHandler->warning("Matrix::operator=() - tensor must have square dimensions\n");
      return *this;
  }

  if (dim != 2 || dim != 3 || dim != 1) {
      g3ErrorHandler->warning("Matrix::operator=() - tensor must be of dimension 2 or 3\n");
      return *this;
  }      

  if (dim == 1) {
      if ((numCols != 1) || (numRows != 1)) {      
	  g3ErrorHandler->warning("Matrix::operator=() - matrix must be %s\n",
				  "1x1 for tensor of dimension 3\n");
	  return *this;
      }      	  
      (*this)(0,0) = V.cval(1,1,1,1);
      
  } else if (dim == 2) {
      if ((numCols != 3) || (numRows != 3)) {      
	  g3ErrorHandler->warning("Matrix::operator=() - matrix must be %s\n",
				  "1x1 for tensor of dimension 3\n");      
	  return *this;
      }
      (*this)(0,0) = V.cval(1,1,1,1);
      (*this)(0,1) = V.cval(1,1,2,2);
      (*this)(0,2) = V.cval(1,1,1,2);      

      (*this)(1,0) = V.cval(2,2,1,1);
      (*this)(1,1) = V.cval(2,2,2,2);
      (*this)(1,2) = V.cval(2,2,1,2);      

      (*this)(2,0) = V.cval(1,2,1,1);
      (*this)(2,1) = V.cval(1,2,2,2);
      (*this)(2,2) = V.cval(1,2,1,2);            

  } else {
      if ((numCols != 6) || (numRows != 6)) {      
	  g3ErrorHandler->warning("Matrix::operator=() - matrix must be %s\n",
				  "1x1 for tensor of dimension 3\n");      
	  return *this;
      }      
      (*this)(0,0) = V.cval(1,1,1,1);
      (*this)(0,1) = V.cval(1,1,2,2);
      (*this)(0,2) = V.cval(1,1,3,3);      
      (*this)(0,3) = V.cval(1,1,1,2);
      (*this)(0,4) = V.cval(1,1,1,3);
      (*this)(0,5) = V.cval(1,1,2,3);      
      
      (*this)(1,0) = V.cval(2,2,1,1);
      (*this)(1,1) = V.cval(2,2,2,2);
      (*this)(1,2) = V.cval(2,2,3,3);      
      (*this)(1,3) = V.cval(2,2,1,2);
      (*this)(1,4) = V.cval(2,2,1,3);
      (*this)(1,5) = V.cval(2,2,2,3);            
      
      (*this)(2,0) = V.cval(3,3,1,1);
      (*this)(2,1) = V.cval(3,3,2,2);
      (*this)(2,2) = V.cval(3,3,3,3);      
      (*this)(2,3) = V.cval(3,3,1,2);
      (*this)(2,4) = V.cval(3,3,1,3);
      (*this)(2,5) = V.cval(3,3,2,3);                  
      
      (*this)(3,0) = V.cval(1,2,1,1);
      (*this)(3,1) = V.cval(1,2,2,2);
      (*this)(3,2) = V.cval(1,2,3,3);      
      (*this)(3,3) = V.cval(1,2,1,2);
      (*this)(3,4) = V.cval(1,2,1,3);
      (*this)(3,5) = V.cval(1,2,2,3);                        
      
      (*this)(4,0) = V.cval(1,3,1,1);
      (*this)(4,1) = V.cval(1,3,2,2);
      (*this)(4,2) = V.cval(1,3,3,3);      
      (*this)(4,3) = V.cval(1,3,1,2);
      (*this)(4,4) = V.cval(1,3,1,3);
      (*this)(4,5) = V.cval(1,3,2,3);                              
      
      (*this)(5,0) = V.cval(2,3,1,1);
      (*this)(5,1) = V.cval(2,3,2,2);
      (*this)(5,2) = V.cval(2,3,3,3);      
      (*this)(5,3) = V.cval(2,3,1,2);
      (*this)(5,4) = V.cval(2,3,1,3);
      (*this)(5,5) = V.cval(2,3,2,3);                                    
  }
  return *this;
}



// virtual Matrix &operator+=(double fact);
// virtual Matrix &operator-=(double fact);
// virtual Matrix &operator*=(double fact);
// virtual Matrix &operator/=(double fact); 
//	The above methods all modify the current matrix. If in
//	derived matrices data kept in data and of sizeData no redef necessary.

Matrix &
Matrix::operator+=(double fact)
{
  // check if quick return
  if (fact == 0.0)
    return *this;

  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ += fact;
  
  return *this;
}




Matrix &
Matrix::operator-=(double fact)
{
  // check if quick return
  if (fact == 0.0)
    return *this;
  
  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ -= fact;

  return *this;
}


Matrix &
Matrix::operator*=(double fact)
{
  // check if quick return
  if (fact == 1.0)
    return *this;
  
  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ *= fact;
  
  return *this;
}

Matrix &
Matrix::operator/=(double fact)
{
    // check if quick return
    if (fact == 1.0)
	return *this;

    if (fact != 0.0) {
      double val = 1.0/fact;

      double *dataPtr = data;
      for (int i=0; i<dataSize; i++)
	*dataPtr++ *= val;

      return *this;
    } else {
      // print out the warining message
      cerr << "WARNING:Matrix::operator/= - 0 factor specified all values in Matrix set to ";
      cerr << MATRIX_VERY_LARGE_VALUE << endl;

      double *dataPtr = data;
      for (int i=0; i<dataSize; i++)
	*dataPtr++ = MATRIX_VERY_LARGE_VALUE;
      
      return *this;
    }
}


//    virtual Matrix operator+(double fact);
//    virtual Matrix operator-(double fact);
//    virtual Matrix operator*(double fact);
//    virtual Matrix operator/(double fact);
//	The above methods all return a new full general matrix.

Matrix
Matrix::operator+(double fact) const
{
    Matrix result(*this);
    result += fact;
    return result;
}

Matrix
Matrix::operator-(double fact) const
{
    Matrix result(*this);
    result -= fact;
    return result;
}

Matrix
Matrix::operator*(double fact) const
{
    Matrix result(*this);
    result *= fact;
    return result;
}

Matrix
Matrix::operator/(double fact) const
{
    if (fact == 0.0) {
	cerr << "Matrix::operator/(const double &fact): ERROR divide-by-zero\n";
	exit(0);
    }
    Matrix result(*this);
    result /= fact;
    return result;
}


//
// MATRIX_VECTOR OPERATIONS
//

Vector
Matrix::operator*(const Vector &V) const
{
    Vector result(numRows);
    
    if (V.Size() != numCols) {
	cerr << "Matrix::operator*(Vector): incompatable sizes\n";
	return result;
    } 
    
    double *dataPtr = data;
    for (int i=0; i<numCols; i++)
      for (int j=0; j<numRows; j++)
	result(j) += *dataPtr++ * V(i);

    /*
    cerr << "HELLO: " << V;
    for (int i=0; i<numRows; i++) {
	double sum = 0.0;
	for (int j=0; j<numCols; j++) {
	    sum += (*this)(i,j) * V(j);
	    if (i == 9) cerr << "sum: " << sum << " " << (*this)(i,j)*V(j) << " " << V(j) << endl;
	}
	result(i) += sum;
    }
    cerr << *this;
    cerr << "HELLO result: " << result;    
    */

    return result;
}

Vector
Matrix::operator^(const Vector &V) const
{
    Vector result(numCols);
    
    if (V.Size() != numRows) {
      cerr << "Matrix::operator*(Vector): incompatable sizes\n";
      return result;
    } 

    double *dataPtr = data;
    for (int i=0; i<numCols; i++)
      for (int j=0; j<numRows; j++)
	result(i) += *dataPtr++ * V(j);

    return result;
}


//
// MATRIX - MATRIX OPERATIONS
//
	    

Matrix
Matrix::operator+(const Matrix &M) const
{
    Matrix result(*this);
    result.addMatrix(1.0,M,1.0);    
    return result;
}
	    
Matrix
Matrix::operator-(const Matrix &M) const
{
    Matrix result(*this);
    result.addMatrix(1.0,M,-1.0);    
    return result;
}
	    
    
Matrix
Matrix::operator*(const Matrix &M) const
{
    Matrix result(numRows,M.numCols);
    
    if (numCols != M.numRows || result.numRows != numRows) {
	cerr << "Matrix::operator*(Matrix): incompatable sizes\n";
	return result;
    } 

    double *resDataPtr = result.data;	    

    int innerDim = numCols;
    int nCols = result.numCols;
    for (int i=0; i<nCols; i++) {
      double *aStartRowDataPtr = data;
      double *bStartColDataPtr = &(M.data[i*innerDim]);
      for (int j=0; j<numRows; j++) {
	double *bDataPtr = bStartColDataPtr;
	double *aDataPtr = aStartRowDataPtr +j;	    
	double sum = 0.0;
	for (int k=0; k<innerDim; k++) {
	  sum += *aDataPtr * *bDataPtr++;
	  aDataPtr += numRows;
	}
	*resDataPtr++ = sum;
      }
    }
    return result;
}



// Matrix operator^(const Matrix &M) const
//	We overload the * operator to perform matrix^t-matrix multiplication.
//	reults = (*this)transposed * M.

Matrix
Matrix::operator^(const Matrix &M) const
{
  Matrix result(numCols,M.numCols);
  
  if (numRows != M.numRows || result.numRows != numCols) {
    cerr << "Matrix::operator*(Matrix): incompatable sizes\n";
    return result;
  } 

    double *resDataPtr = result.data;	    

    int innerDim = numRows;
    int nCols = result.numCols;
    for (int i=0; i<nCols; i++) {
      double *aDataPtr = data;
      double *bStartColDataPtr = &(M.data[i*innerDim]);
      for (int j=0; j<numCols; j++) {
	double *bDataPtr = bStartColDataPtr;
	double sum = 0.0;
	for (int k=0; k<innerDim; k++) {
	  sum += *aDataPtr++ * *bDataPtr++;
	}
	*resDataPtr++ = sum;
      }
    }

    return result;
}
    



Matrix &
Matrix::operator+=(const Matrix &M)
{
#ifdef _G3DEBUG
  if (numRows != M.numRows || numCols != M.numCols) {
    g3ErrorHandler->warning("Matrix::operator+=(const Matrix &M) - matrices incompatable [%d %d] [%d %d]\n",
			    numRows, numCols, M.numRows, M.numCols);
    return *this;
  }
#endif

  double *dataPtr = data;
  double *otherData = M.data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ += *otherData++;
  
  return *this;
}

Matrix &
Matrix::operator-=(const Matrix &M)
{
#ifdef _G3DEBUG
  if (numRows != M.numRows || numCols != M.numCols) {
    g3ErrorHandler->warning("Matrix::operator-=(const Matrix &M) - matrices incompatable [%d %d] [%d %d]\n",
			    numRows, numCols, M.numRows, M.numCols);
    return *this;
  }
#endif

  double *dataPtr = data;
  double *otherData = M.data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ -= *otherData++;
  
  return *this;
}


//
// Input/Output Methods
//

void 
Matrix::Output(ostream &s) const
{
    for (int i=0; i<noRows(); i++) {
	for (int j=0; j<noCols(); j++)
	    s << (*this)(i,j) << " ";
	s << "\n";
    }
}


void 
Matrix::Input(istream &s)
{
    for (int i=0; i<noRows(); i++)
	for (int j=0; j<noCols(); j++)
	    s >> (*this)(i,j);
}	


//
// friend stream functions for input and output
//

ostream &operator<<(ostream &s, const Matrix &V)
{
    s << "\n";
    V.Output(s);
    s << "\n";        
    return s;
}

	
	
istream &operator>>(istream &s, Matrix &V)
{
    V.Input(s);
    return s;
}






int
Matrix::Assemble(const Matrix &V, int init_row, int init_col, double fact) 
{
  int pos_Rows, pos_Cols;
  int res = 0;
  
  int VnumRows = V.numRows;
  int VnumCols = V.numCols;
  
  int final_row = init_row + VnumRows - 1;
  int final_col = init_col + VnumCols - 1;
  
  if ((init_row >= 0) && (final_row < numRows) && (init_col >= 0) && (final_col < numCols))
  {
     for (int i=0; i<VnumCols; i++) 
     {
        pos_Cols = init_col + i;
        for (int j=0; j<VnumRows; j++) 
        {
           pos_Rows = init_row + j;
      
	   (*this)(pos_Rows,pos_Cols) += V(j,i)*fact;
        }
     }
  }  
  else 
  {
     cerr << "WARNING: Matrix::Assemble(const Matrix &V, int init_row, int init_col, double fact): ";
     cerr << "position outside bounds \n";
     res = -1;
  }

  return res;
}




int
Matrix::AssembleTranspose(const Matrix &V, int init_row, int init_col, double fact) 
{
  int pos_Rows, pos_Cols;
  int res = 0;
  
  int VnumRows = V.numRows;
  int VnumCols = V.numCols;
  
  int final_row = init_row + VnumCols - 1;
  int final_col = init_col + VnumRows - 1;
  
  if ((init_row >= 0) && (final_row < numRows) && (init_col >= 0) && (final_col < numCols))
  {
     for (int i=0; i<VnumRows; i++) 
     {
        pos_Cols = init_col + i;
        for (int j=0; j<VnumCols; j++) 
        {
           pos_Rows = init_row + j;
      
	   (*this)(pos_Rows,pos_Cols) += V(i,j)*fact;
        }
     }
  }  
  else 
  {
     cerr << "WARNING: Matrix::AssembleTranspose(const Matrix &V, int init_row, int init_col, double fact): ";
     cerr << "position outside bounds \n";
     res = -1;
  }

  return res;
}




int
Matrix::Extract(const Matrix &V, int init_row, int init_col, double fact) 
{
  int pos_Rows, pos_Cols;
  int res = 0;
  
  int VnumRows = V.numRows;
  int VnumCols = V.numCols;
  
  int final_row = init_row + numRows - 1;
  int final_col = init_col + numCols - 1;
  
  if ((init_row >= 0) && (final_row < VnumRows) && (init_col >= 0) && (final_col < VnumCols))
  {
     for (int i=0; i<numCols; i++) 
     {
        pos_Cols = init_col + i;
        for (int j=0; j<numRows; j++) 
        {
           pos_Rows = init_row + j;
      
	   (*this)(j,i) = V(pos_Rows,pos_Cols)*fact;
        }
     }
  }  
  else 
  {
     cerr << "WARNING: Matrix::Extract(const Matrix &V, int init_row, int init_col, double fact): ";
     cerr << "position outside bounds \n";
     res = -1;
  }

  return res;
}


Matrix operator*(double a, const Matrix &V)
{
  return V * a;
}




