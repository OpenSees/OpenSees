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
                                                                        
// $Revision: 1.8 $
// $Date: 2010-01-21 21:43:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Information.h,v $
                                                                        
                                                                        
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


#include <OPS_Globals.h>
class ID;
class Vector;


#include <fstream>
using std::ofstream;

enum InfoType {UnknownType, IntType, DoubleType, 
	       IdType, VectorType, MatrixType};
		   
class Information
{
  public:
    Information();
    Information(int val);
    Information(double val);
    Information(const ID &val);
    Information(const Vector &val);
    Information(const Matrix &val);
    Information(const ID &val1, const Vector &val2);
    
    virtual ~Information();
    
    virtual int setInt(int newInt);
    virtual int setDouble(double newDouble);
    virtual int setID(const ID &newID);
    virtual int setVector(const Vector &newVector);
    virtual int setMatrix(const Matrix &newMatrix);
    virtual int setString(const char *theString);
    
    virtual void Print(OPS_Stream &s, int flag = 0);
    virtual void Print(ofstream &s, int flag = 0);
    virtual const Vector &getData(void);

    // data that is stored in the information object
    InfoType	theType;   // information about data type
    int		theInt;    // an integer value
    double	theDouble; // a double value
    ID		*theID;    // pointer to an ID object, created elsewhere
    Vector 	*theVector;// pointer to a Vector object, created elsewhere
    Matrix	*theMatrix;// pointer to a Matrix object, created elsewhere
    char        *theString;// pointer to string

  protected:
    
  private:        

};

#endif

