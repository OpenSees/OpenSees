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
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/Vector.cpp,v $
                                                                        
                                                                        
// File: ~/matrix/Vector.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for Vector.
//
// What: "@(#) Vector.C, revA"

#include "Vector.h"
#include "Matrix.h"
#include "ID.h"

#include <stdlib.h>
#include <math.h>

double Vector::VECTOR_NOT_VALID_ENTRY =0.0;

// Vector():
//	Standard constructor, sets size = 0;

Vector::Vector()
: sz(0), theData(0), fromFree(0)
{

}

// Vector(int size):









//	Constructor used to allocate a vector of size size.

Vector::Vector(int size)
: sz(size), theData(0), fromFree(0)
{
#ifdef _G3DEBUG
  if (sz <= 0) {
    g3ErrorHandler->warning("Vector::Vector(int) - size %d specified <= 0\n",size);
    sz = 1;
  }
#endif

  // get some space for the vector
  theData = (double *)malloc(size*sizeof(double));

  if (theData == 0) {
    g3ErrorHandler->fatal("Vector::Vector(int) - out of memory creating vector of size %d\n",size);
    
    sz = 0; // set this should fatal error handler not kill process!!
  }
    
  // zero the components
  for (int i=0; i<sz; i++)
    theData[i] = 0.0;
}


// Vector::Vector(double *data, int size)

Vector::Vector(double *data, int size)
: sz(size),theData(data),fromFree(1)
{
#ifdef _G3DEBUG
  if (sz <= 0) {
    g3ErrorHandler->warning("Vector::Vector(double *, size) - size %d specified <= 0\n",size);
    sz = 0;
  }
#endif
}
 


// Vector(const Vector&):
//	Constructor to init a vector from another.

Vector::Vector(const Vector &other)
: sz(other.sz),theData(0),fromFree(0)
{
#ifdef _G3DEBUG
  if (sz < 0) {
    g3ErrorHandler->warning("Vector::Vector(int) - size %d specified <= 0\n",sz);
    sz = 1;
  }
#endif

  theData = (double *)malloc(other.sz*sizeof(double));    
  
  if (theData == 0) {
    g3ErrorHandler->fatal("Vector::Vector(int) - out of memory creating vector of size %d\n",sz);
    sz = 0;
  }

  // copy the component data
  for (int i=0; i<sz; i++)
    theData[i] = other.theData[i];
}	






// ~Vector():
// 	destructor, deletes the [] data

Vector::~Vector()
{
  if (sz != 0 && fromFree == 0) 
    free((void *)theData);
}



int 
Vector::setData(double *newData, int size){
  if (sz != 0 && fromFree == 0) 
    free((void *)theData);
  sz = size;
  theData = newData;
  fromFree = 1;

  if (sz <= 0) {
    g3ErrorHandler->warning("Vector::Vector(double *, size) - size %d specified <= 0\n",size);
    sz = 0;
  }

  return 0;
}


// Assemble(Vector &x, ID &y, double fact ):
//	Method to assemble into object the Vector V using the ID l.
//	If ID(x) does not exist program writes error message if
//	VECTOR_CHECK defined, otherwise ignores it and goes on.

int 
Vector::Assemble(const Vector &V, const ID &l, double fact )
{
  int result = 0;
  int pos;
  for (int i=0; i<l.Size(); i++) {
    pos = l(i);
    
    if (pos < 0)
      ;
    else if ((pos < sz) && (i < V.Size()))
      // assemble into vector
      theData[pos] += V.theData[i] *fact;
    else {
      result = -1;
      if (pos < sz)
	g3ErrorHandler->warning("Vector::Assemble() %d out of range [1 %d]\n",
				pos, sz-1);
      else
	g3ErrorHandler->warning("Vector::Assemble() %d out of range [1 %d]\n",
				pos, V.Size()-1);
    }
  }
  return result;
}
    


int
Vector::addVector(double thisFact, const Vector &other, double otherFact )
{
  // check if quick return
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0; 

  // if sizes are compatable add
#ifdef _G3DEBUG
  if (sz != other.sz) {
    // else sizes are incompatable, do nothing but warning
    g3ErrorHandler->warning( "WARNING Vector::addVector() - incompatable Vector sizes\n");
    return -1;
  }
#endif


  if (thisFact == 1.0) {

    // want: this += other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ += *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ -= *otherDataPtr++;
    } else 
      for (int i=0; i<sz; i++) 
	*dataPtr++ += *otherDataPtr++ * otherFact;
  } 

  else if (thisFact == 0.0) {

    // want: this = other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ = *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ = *otherDataPtr++;
    } else 
      for (int i=0; i<sz; i++) 
	*dataPtr++ = *otherDataPtr++ * otherFact;
  }

  else {

    // want: this = this * thisFact + other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) {
	double value = *dataPtr * thisFact + *otherDataPtr++;
	*dataPtr++ = value;
      }
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) {
	double value = *dataPtr * thisFact - *otherDataPtr++;
	*dataPtr++ = value;
      }
    } else 
      for (int i=0; i<sz; i++) {
	double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
	*dataPtr++ = value;
      }
  } 

  // successfull
  return 0;
}
	    
	
int
Vector::addMatrixVector(double thisFact, const Matrix &m, const Vector &v, double otherFact )
{
  // see if quick return
  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

  // check the sizes are compatable
#ifdef _G3DEBUG
  // check the sizes are compatable
  if ((sz != m.noRows()) && (m.noCols() != v.sz)) {
    // otherwise incompatable sizes
    g3ErrorHandler->warning("Vector::addMatrixVector() %s [%d] = [%d %d] [%d]\n",
			    "- incompatable sizes",
			    sz, m.noRows(), m.noCols(), v.sz);
    return -1;    
  }
#endif

  if (thisFact == 1.0) {

    // want: this += m * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    } 
    else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] -= *matrixDataPtr++ * otherData;
      }
    } 
    else { // have to do the multiplication
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++ * otherFact;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else if (thisFact == 0.0) {
    
    // want: this = m * v * otherFact
    for (int i=0; i<sz; i++)
      theData[i] = 0.0;

    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    } 
    else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] -= *matrixDataPtr++ * otherData;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++ * otherFact;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else {

    // want: this = this * thisFact + m * v * otherFact
    for (int i=0; i<sz; i++)
      theData[i] *= thisFact;

    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] -= *matrixDataPtr++ * otherData;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++ * otherFact;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }
  
  // successfull
  return 0;
}



int
Vector::addMatrixTransposeVector(double thisFact, 
				 const Matrix &m, 
				 const Vector &v, 
				 double otherFact )
{
  // see if quick return
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0;

#ifdef _G3DEBUG
  // check the sizes are compatable
  if ((sz != m.noRows()) && (m.noRows() != v.sz)) {
    // otherwise incompatable sizes
    g3ErrorHandler->warning("Vector::addMatrixTransposeVector() %s [%d] = [%d %d] [%d]\n",
			    "- incompatable sizes",
			    sz, m.noCols(), m.noRows(), v.sz);
    return -1;    
  }
#endif

  if (thisFact == 1.0) {

    // want: this += m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] += sum;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] -= sum;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] += sum * otherFact;
      }
    }
  }

  else if (thisFact == 0.0) {

    // want: this = m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] = sum;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] = -sum;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] = sum * otherFact;
      }
    }
  } 

  else {

    // want: this = this * thisFact + m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	double value = theData[i] * thisFact + sum;
	theData[i] = value;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	double value = theData[i] * thisFact - sum;
	theData[i] = value;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	double value = theData[i] * thisFact + sum * otherFact;
	theData[i] = value;
      }
    }
}

  return 0;
}
	
	



// double Norm();
//	Method to return the norm of  vector. (non-const as may save norm for later)

double
Vector::Norm(void) const
{
  double value = 0;
  for (int i=0; i<sz; i++) {
    double data = theData[i];
    value += data*data;
  }
  return sqrt(value);
}


double &
Vector::operator[](int x) 
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
      g3ErrorHandler->warning("Vector::() - x %d outside range [0, %d]\n",x,sz-1);
      return VECTOR_NOT_VALID_ENTRY;
  }
#endif

  return theData[x];
}

double Vector::operator[](int x) const
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
    g3ErrorHandler->warning("Vector::() - x %d outside range [0, %d]\n",x,sz-1);
    return VECTOR_NOT_VALID_ENTRY;
  }
#endif

  return theData[x];
}


// operator()(const ID &rows) const
//	Method to return a vector whose components are the components of the
//	current vector located in positions given by the ID rows.


Vector 
Vector::operator()(const ID &rows) const
{
  // create a new Vector to be returned
  Vector result(rows.Size());

  // check if obtained VEctor of correct size
  if (result.Size() != rows.Size()) {
    g3ErrorHandler->warning("Vector::()(ID) - new Vector could not be constructed\n");
    return result;
  }

  // copy the appropraite contents from current to result     
  int pos;
  for (int i=0; i<rows.Size(); i++) {
    pos = rows(i);
    if (pos <0 || pos >= sz) {
      g3ErrorHandler->warning("Vector::()(ID) - invalid location %d outside range [0,%d]\n",
			      pos,sz-1);
    } else
      result(i) = (*this)(pos);
  }
  return result;
}


// Vector &operator=(const Vector  &V):
//	the assignment operator, This is assigned to be a copy of V. if sizes
//	are not compatable this.theData [] is deleted. The data pointers will not
//	point to the same area in mem after the assignment.
//

Vector &
Vector::operator=(const Vector &V) 
{
  // first check we are not trying v = v
  if (this != &V) {

	  /*
#ifdef _G3DEBUG
    // check size compatability, if different warning
    if (sz != V.sz) 
      g3ErrorHandler->warning("Vector::operator=() - vectors of differing sizes\n");
#endif
	  */

      if (sz != V.sz)  {

#ifdef _G3DEBUG
	  g3ErrorHandler->warning("Vector::operator=() - vectors of differing sizes\n");
#endif

	  delete [] this->theData;
	  this->sz = V.sz;
	  theData = new double[sz];
      }


      //	 copy the data
      for (int i=0; i<sz; i++)
	  theData[i] = V.theData[i];
  }

  return *this;
}




Vector &
Vector::operator=(const Tensor &V) 
{
  int rank = V.rank();
  if (rank != 2) {
      g3ErrorHandler->warning("Vector::operator=() - tensor must be of rank 2\n");
      return *this;
  }
  int dim = V.dim(1);
  if (dim != V.dim(2)) {
      g3ErrorHandler->warning("Vector::operator=() - tensor must have square dimensions\n");
      return *this;
  }

  if (dim != 2 || dim != 3 || dim != 1) {
      g3ErrorHandler->warning("Vector::operator=() - tensor must be of dimension 2 or 3\n");
      return *this;
  }      
  
  if (dim == 1) {
      if (sz != 1) {
	  g3ErrorHandler->warning("Vector::operator=() - Vector size must be 1\n"); 
	  return *this;
      }
      theData[0] = V.cval(1,1);
  } else if (dim == 2) {
      if (sz != 3) {
	  g3ErrorHandler->warning("Vector::operator=() - Vector size must be 3\n"); 
	  return *this;
      }
      theData[0] = V.cval(1,1);
      theData[1] = V.cval(2,2);
      theData[2] = V.cval(1,2);
  } else {
      if (sz != 6) {
	  g3ErrorHandler->warning("Vector::operator=() - Vector size must be 6\n"); 
	  return *this;
      }      
      theData[0] = V.cval(1,1);
      theData[1] = V.cval(2,2);
      theData[2] = V.cval(3,3);
      theData[3] = V.cval(1,2);
      theData[4] = V.cval(1,3);
      theData[5] = V.cval(2,3);
  }
  return *this;
}




// Vector &operator+=(double fact):
//	The += operator adds fact to each element of the vector, data[i] = data[i]+fact.

Vector &Vector::operator+=(double fact)
{
  if (fact != 0.0)
    for (int i=0; i<sz; i++)
      theData[i] += fact;
  return *this;    
}



// Vector &operator-=(double fact)
//	The -= operator subtracts fact from each element of the vector, data[i] = data[i]-fact.

Vector &Vector::operator-=(double fact)
{
  if (fact != 0.0)
    for (int i=0; i<sz; i++)
      theData[i] -= fact;
  return *this;    
}



// Vector &operator*=(double fact):
//	The -= operator subtracts fact from each element of the vector, data[i] = data[i]-fact.

Vector &Vector::operator*=(double fact)
{
  for (int i=0; i<sz; i++)
    theData[i] *= fact;
  return *this;
}



// Vector &operator/=(double fact):
//	The /= operator divides each element of the vector by fact, theData[i] = theData[i]/fact.
//	Program exits if divide-by-zero would occur with warning message.

Vector &Vector::operator/=(double fact)
{
  if (fact == 0.0) { // instead of divide-by-zero error set to VECTOR_VERY_LARGE_VALUE
    for (int i=0; i<sz; i++)
      theData[i] = VECTOR_VERY_LARGE_VALUE;
  } else {
    for (int i=0; i<sz; i++)
      theData[i] /= fact;
  }
  
  return *this;
}




// Vector operator+(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]+fact;

Vector 
Vector::operator+(double fact) const
{
  Vector result(*this);
  if (result.Size() != sz) 
    g3ErrorHandler->warning("Vector::operator+(double) - ran out of memory for new Vector\n");

  result += fact;
  return result;
}



// Vector operator-(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]-fact;

Vector 
Vector::operator-(double fact) const
{
    Vector result(*this);
    if (result.Size() != sz) 
      g3ErrorHandler->warning("Vector::operator-(double) - ran out of memory for new Vector\n");

    result -= fact;
    return result;
}



// Vector operator*(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]*fact;

Vector 
Vector::operator*(double fact) const
{
    Vector result(*this);
    if (result.Size() != sz) 
      g3ErrorHandler->warning("Vector::operator-(double) - ran out of memory for new Vector\n");

    result *= fact;
    return result;
}


// Vector operator/(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]/fact; Exits if divide-by-zero error.

Vector 
Vector::operator/(double fact) const
{
    if (fact == 0.0) 
      g3ErrorHandler->warning("Vector::operator/(double fact) - divide-by-zero error coming\n");

    Vector result(*this);
    if (result.Size() != sz) 
      g3ErrorHandler->warning("Vector::operator-(double) - ran out of memory for new Vector\n");

    result /= fact;
    return result;
}



// Vector &operator+=(const Vector &V):
//	The += operator adds V's data to data, data[i]+=V(i). A check to see if
//	vectors are of same size is performed if VECTOR_CHECK is defined.

Vector &
Vector::operator+=(const Vector &other)
{
#ifdef _G3DEBUG
  if (sz != other.sz) {
    g3ErrorHandler->warning("WARNING Vector::operator+=(Vector):Vectors not of same sizes: %d %d\n",
			     sz, other.sz);
    return *this;
  }    
#endif

  for (int i=0; i<sz; i++)
    theData[i] += other.theData[i];
  return *this;	    
}



// Vector &operator-=(const Vector &V):
//	The -= operator subtracts V's data from  data, data[i]+=V(i). A check 
//   	to see if vectors are of same size is performed if VECTOR_CHECK is defined.

Vector &
Vector::operator-=(const Vector &other)
{
#ifdef _G3DEBUG
  if (sz != other.sz) {
    g3ErrorHandler->warning("WARNING Vector::operator-=(Vector): Vectors not of same sizes: %d %d\n",
			    sz, other.sz);
    return *this;
  }
#endif
  
  for (int i=0; i<sz; i++)
    theData[i] -= other.theData[i];
  return *this;    
}



// Vector operator+(const Vector &V):
//	The + operator checks the two vectors are of the same size if VECTOR_CHECK is defined.
// 	Then returns a Vector whose components are the vector sum of current and V's data.

Vector 
Vector::operator+(const Vector &b) const
{
#ifdef _G3DEBUG
  if (sz != b.sz) {
    g3ErrorHandler->warning( "Vector::operator+(Vector): Vectors not of same size, sizes: %d %d",
			     sz, b.sz);
			     
    return *this;
  }
#endif

    Vector result(*this);

    // check new Vector of correct size
  if (result.Size() != sz) {
    g3ErrorHandler->warning("Vector::operator-(Vector): new Vector not of correct size ");
    return result;
  }
  result += b;
  return result;
}


// Vector operator-(const Vector &V):
//	The - operator checks the two vectors are of the same size and then returns a Vector
//	whose components are the vector difference of current and V's data.

Vector 
Vector::operator-(const Vector &b) const
{
#ifdef _G3DEBUG
  if (sz != b.sz) {
    g3ErrorHandler->warning("WARNING Vector::operator-(Vector): Vectors not of same sizes: %d %d\n",
			    sz, b.sz);
    return *this;
  }
#endif

  Vector result(*this);

  // check new Vector of correct size
  if (result.Size() != sz) {
    g3ErrorHandler->warning("Vector::operator-(Vector): new Vector not of correct size ");
    return result;
  }

  result -= b;
  return result;
}



// double operator^(const Vector &V) const;
//	Method to perform (Vector)transposed * vector.
double
Vector::operator^(const Vector &V) const
{
#ifdef _G3DEBUG
  if (sz != V.sz) {
    g3ErrorHandler->warning( "WARNING Vector::operator-(Vector): Vectors not of same sizes: %d %d\n",
			    sz, V.sz);
    return 0.0;
  }
#endif

  double result = 0.0;
  double *dataThis = theData;
  double *dataV = V.theData;
    for (int i=0; i<sz; i++)
      result += *dataThis++ * *dataV++;

    return result;
}


// Vector operator/(const Matrix &M) const;    
//	Method to return inv(M)*this

Vector
Vector::operator/(const Matrix &M) const
{
  Vector res(M.noRows());
    
  if (M.noRows() != M.noCols()) { // if not square do least squares solution
    Matrix A(M^M);
    A.Solve(*this, res);    
  }
  else {
    M.Solve(*this, res);
  }
  return res;
}
    
	

	



// friend ostream &operator<<(ostream &s, const Vector &V)
//	A function is defined to allow user to print the vectors using ostreams.

ostream &operator<<(ostream &s, const Vector &V)
{
  for (int i=0; i<V.Size(); i++) 
      s << V(i) << " ";

  return s << "\n";
}

// friend istream &operator>>(istream &s, Vector &V)
//	A function is defined to allow user to input the data into a Vector which has already
//	been constructed with data, i.e. Vector(int) or Vector(const Vector &) constructors.

istream &operator>>(istream &s, Vector &V)
{
  for (int i=0; i<V.Size(); i++) 
    s >> V(i);
  
    return s;
}



int
Vector::Assemble(const Vector &V, int init_pos, double fact) 
{
  int res = 0;
  int cur_pos   = init_pos;  
  int final_pos = init_pos + V.sz - 1;
  
  if ((init_pos >= 0) && (final_pos < sz))
  {
     for (int j=0; j<V.sz; j++) 
        (*this)(cur_pos++) += V(j)*fact;
  }
  else 
  {
     cerr << "WARNING: Vector::Assemble(const Vector &V, int init_pos, double fact): ";
     cerr << "position outside bounds \n";
     res = -1;
  }

  return res;
}




int
Vector::Extract(const Vector &V, int init_pos, double fact) 
{
  int res = 0;
  int cur_pos   = init_pos;  
  int final_pos = init_pos + sz - 1;
  
  if ((init_pos >= 0) && (final_pos < V.sz))
  {
     for (int j=0; j<sz; j++) 
        (*this)(j) = V(cur_pos++)*fact;
  }
  else 
  {
     cerr << "WARNING: Vector::Assemble(const Vector &V, int init_pos, double fact): ";
     cerr << "position outside bounds \n";
     res = -1;
  }

  return res;
}
