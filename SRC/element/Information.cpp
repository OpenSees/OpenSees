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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-18 10:29:13 $
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
		g3ErrorHandler->warning("%s -- failed to allocate ID",
			"Information::Information");
}

Information::Information(const Vector &val) 
  :theType(VectorType),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
	// Make a copy
    theVector = new Vector(val);

	if (theVector == 0)
		g3ErrorHandler->warning("%s -- failed to allocate Vector",
			"Information::Information");
}

Information::Information(const Matrix &val) 
  :theType(MatrixType),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
	// Make a copy
    theMatrix = new Matrix(val);

	if (theMatrix == 0)
		g3ErrorHandler->warning("%s -- failed to allocate Matrix",
			"Information::Information");
}

Information::Information(const Tensor &val) 
  :theType(TensorType),
  theID(0), theVector(0), theMatrix(0), theTensor(0)
{
	// Make a copy
    theTensor = new Tensor(val);

	if (theTensor == 0)
		g3ErrorHandler->warning("%s -- failed to allocate Tensor",
			"Information::Information");
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
	}
	else
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
Information::Print(ostream &s, int flag)
{
    if (theType == IntType)
		s << theInt << endl;
	else if (theType == DoubleType)
		s << theDouble << endl;
	else if (theType == IdType && theID != 0)
		s << *theID;
	else if (theType == VectorType && theVector != 0)
		s << *theVector;
	else if (theType == MatrixType && theMatrix != 0)
		s << *theMatrix;
	else if (theType == TensorType && theTensor != 0)
		// No overloaded << for Tensors yet!
		//s << *theTensor;
		s << "No Tensor output";
	else
		return;

}
