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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-14 23:01:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Information.cpp,v $
                                                                        
                                                                        
// File: ~/element/Information.C
//
// Written: fmk 10/99
// Revised:
//
// Purpose: This file contains the class implementation for Information.

#include <Information.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Tensor.h>

Information::Information() 
  :theType(UnknownType),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
    // does nothing
}

Information::Information(int val) 
  :theType(IntType), theInt(val),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
    // does nothing
}

Information::Information(double val) 
  :theType(DoubleType), theDouble(val),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
  // does nothing
}

Information::Information(const ID &val) 
  :theType(IdType),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
  // Make a copy
  theID = new ID(val);

  if (theID == 0)
    opserr << "Information::Information -- failed to allocate\n";
}

Information::Information(const Vector &val) 
  :theType(VectorType),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
  // Make a copy
  theVector = new Vector(val);
  
  if (theVector == 0)
    opserr << "Information::Information -- failed to allocate Vector\n";
}

Information::Information(const Matrix &val) 
  :theType(MatrixType),
   theID(0), theVector(0), theMatrix(0), theTensor(0)
{
  // Make a copy
  theMatrix = new Matrix(val);
  
  if (theMatrix == 0)
    opserr << "Information::Information -- failed to allocate Matrix\n";
}

Information::Information(const Tensor &val) 
  :theType(TensorType),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
  // Make a copy
  theTensor = new Tensor(val);
  
  if (theTensor == 0)
    opserr << "Information::Iformation -- failed to allocate Tensor\n";
}

Information::~Information() 
{
  if (theID != 0)
    delete theID;
  
  if (theVector != 0)
    delete theVector;
  
  if (theMatrix != 0)
    delete theMatrix;
  
  if (theTensor != 0)
    delete theTensor;
}

int 
Information::setInt(int newInt)
{
  theInt = newInt;
  
  return 0;
}

int 
Information::setDouble(double newDouble)
{
  theDouble = newDouble;
  
  return 0;
}

int 
Information::setID(const ID &newID)
{
  if (theID != 0) {
    *theID = newID;
    return 0;
  } else
    return -1;
}

int 
Information::setVector(const Vector &newVector)
{
  if (theVector != 0) {
    *theVector = newVector;
    return 0;
  }
  else
    return -1;
}

int 
Information::setMatrix(const Matrix &newMatrix)
{
  if (theMatrix != 0) {
    *theMatrix = newMatrix;
    return 0;
  }
  else
    return -1;
}

int 
Information::setTensor(const Tensor &newTensor)
{
  if (theTensor != 0) {
    *theTensor = newTensor;
    return 0;
  }
  else
    return -1;
}

void 
Information::Print(OPS_Stream &s, int flag)
{
  if (theType == IntType)
    s << theInt << " ";
  else if (theType == DoubleType)
    s << theDouble << " ";
  else if (theType == IdType && theID != 0)
    for (int i=0; i<theID->Size(); i++)
      s << (*theID)(i) << " ";
  else if (theType == VectorType && theVector != 0)
    for (int i=0; i<theVector->Size(); i++)
      s << (*theVector)(i) << " ";
  else if (theType == MatrixType && theMatrix != 0) {
    for (int i=0; i<theMatrix->noRows(); i++) {
      for (int j=0; j<theMatrix->noCols(); j++)
	s <<  (*theMatrix)(i,j) << " ";
	s << endln;
    }
  } else if (theType == TensorType && theTensor != 0)
    // No overloaded << for Tensors yet!
    //s << *theTensor;
    s << "No Tensor output";
  else
    return;
}


void 
Information::Print(ofstream &s, int flag)
{
  if (theType == IntType)
    s << theInt << " ";
  else if (theType == DoubleType)
    s << theDouble << " ";
  else if (theType == IdType && theID != 0)
    for (int i=0; i<theID->Size(); i++)
      s << (*theID)(i) << " ";
  else if (theType == VectorType && theVector != 0)
    for (int i=0; i<theVector->Size(); i++)
      s << (*theVector)(i) << " ";
  else if (theType == MatrixType && theMatrix != 0) {
    for (int i=0; i<theMatrix->noRows(); i++) {
      for (int j=0; j<theMatrix->noCols(); j++)
	s <<  (*theMatrix)(i,j) << " ";
	s << endln;
    }
  }
  else if (theType == TensorType && theTensor != 0)
    // No overloaded << for Tensors yet!
    //s << *theTensor;
    s << "No Tensor output";
  else
    return;
}

const Vector &
Information::getData(void) 
{
  if (theType == IntType) {
    if (theVector == 0) 
      theVector = new Vector(1);
    (*theVector)(0) = theInt;
  } else if (theType == DoubleType) {
    if (theVector == 0) 
      theVector = new Vector(1);
    (*theVector)(0) = theDouble;
  } else if (theType == IdType && theID != 0) {
    if (theVector == 0) 
      theVector = new Vector(theID->Size());
    for (int i=0; i<theID->Size(); i++)
      (*theVector)(i) =  (*theID)(i);
  } else if (theType == MatrixType && theMatrix != 0) {
    int noRows = theMatrix->noRows();
    int noCols = theMatrix->noCols();
    if (theVector == 0) 
      theVector = new Vector(noRows * noCols);
    for (int i=0; i<noRows; i++)
      for (int j=0; j<noCols; j++)
	(*theVector)(i) = (*theMatrix)(i,j);
  }
  
  return *theVector;
}
