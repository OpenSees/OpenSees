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
// What: "@(#) InformationProxy.h, revA"

#ifndef InformationProxy_h
#define InformationProxy_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>



class InformationProxy
{
public:
	InformationProxy(Element* elemenet);
	double getDouble(int id);
	Vector* getVector(int id);
	Matrix* getMatrix(int id);
	ID* getID(int id);
	int getInt(int id);
	char* getString(int id);
	~InformationProxy();

	//// data that is stored in the information object
	//int		theInt;    // an integer value
	//double	theDouble; // a double value
	//ID		*theID;    // pointer to an ID object, created elsewhere
	//Vector 	*theVector;// pointer to a Vector object, created elsewhere
	//Matrix	*theMatrix;// pointer to a Matrix object, created elsewhere
	//char        *theString;// pointer to string

protected:

private:
	Element* theEle;
};

#endif

