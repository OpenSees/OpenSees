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
// $Source: /usr/local/cvs/OpenSees/SRC/element/Information.h,v $
                                                                        
                                                                        
// File: ~/element/Information.h
//
// Written: fmk 
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for Information.
// Information is a class in which all data members are public, i.e. basically
// a struct.
//
// What: "@(#) Information.h, revA"

#ifndef Information_h
#define Information_h

class Matrix;
class Vector;
class ID;
#include <Tensor.h>

class ostream;

enum InfoType {UnknownType, IntType, DoubleType, 
	       IdType, VectorType, MatrixType, TensorType};
		   
class Information
{
  public:
    Information();
	Information(int val);
	Information(double val);
	Information(const ID &val);
	Information(const Vector &val);
	Information(const Matrix &val);
	Information(const Tensor &val);

    virtual ~Information();

	virtual int setInt(int newInt);
	virtual int setDouble(double newDouble);
    virtual int setID(const ID &newID);
	virtual int setVector(const Vector &newVector);
	virtual int setMatrix(const Matrix &newMatrix);
	virtual int setTensor(const Tensor &newTensor);

	virtual void Print(ostream &s, int flag = 0);

    // data that is stored in the information object
    InfoType	theType;   // information about data type
    int		theInt;    // an integer value
    double	theDouble; // a double value
    ID		*theID;    // pointer to an ID object, created elsewhere
    Vector 	*theVector;// pointer to a Vector object, created elsewhere
    Matrix	*theMatrix;// pointer to a Matrix object, created elsewhere
	Tensor  *theTensor;// pointer to a Tensor object, created elsewhere

  protected:
    
  private:    
};

#endif

