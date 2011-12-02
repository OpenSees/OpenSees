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
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/SectionForceDeformation.cpp,v $
                                                                        
                                                                        
// File: ~/material/SectionForceDeformation.C
//
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for SectionForceDeformation.
//
// What: "@(#) SectionForceDeformation.C, revA"

#include <SectionForceDeformation.h>
#include <Information.h>
#include <Matrix.h>
#include <Vector.h>

#include <string.h>

double invert2by2Matrix(const Matrix &a, Matrix &b);
double invert3by3Matrix(const Matrix &a, Matrix &b);
void invertMatrix(int n, const Matrix &a, Matrix &b);

SectionForceDeformation::SectionForceDeformation(int tag, int classTag)
:Material(tag,classTag), fDefault(0)
{

}

SectionForceDeformation::~SectionForceDeformation()
{
	if (fDefault != 0)
		delete fDefault;
}

const Matrix&
SectionForceDeformation::getSectionFlexibility ()
{
	int order = this->getOrder();

	if (fDefault == 0) {		
		fDefault = new Matrix(order,order);
		if (fDefault == 0)
			g3ErrorHandler->fatal("%s -- failed to allocate flexibility matrix",
				"SectionForceDeformation::getSectionFlexibility");
	}

	const Matrix &k = this->getSectionTangent();

	switch(order) {
	case 1:
		if (k(0,0) != 0.0)
			(*fDefault)(0,0) = 1.0/k(0,0);
		break;
	case 2:
		invert2by2Matrix(k,*fDefault);
		break;
	case 3:
		invert3by3Matrix(k,*fDefault);
		break;
	default:
		invertMatrix(order,k,*fDefault);
		break;
	}

	return *fDefault;
}

int 
SectionForceDeformation::setResponse(char **argv, int argc, Information &sectInfo)
{
    // deformations
    if ((strcmp(argv[0],"deformations") ==0) || 
	(strcmp(argv[0],"deformation") ==0)) {

	Vector *theVector = new Vector(this->getOrder());
	if (theVector == 0) {
	    cerr << "WARNING SectionForceDeformation::setResponse() - out of memory\n";
	    return -1;
	} 
	sectInfo.theVector = theVector;
	sectInfo.theType = VectorType;	
	return 1;
    } 

    // stress resultants
    else if ((strcmp(argv[0],"forces") ==0) ||
	     (strcmp(argv[0],"force") ==0)) {

	Vector *theVector = new Vector(this->getOrder());
	if (theVector == 0) {
	    cerr << "WARNING SectionForceDeformation::setResponse() - out of memory\n";
	    return -1;
	} 
	sectInfo.theVector = theVector;
	sectInfo.theType = VectorType;	
	return 2;
    } 

	// tangent stiffness
	else if (strcmp(argv[0],"stiff") == 0 ||
		strcmp(argv[0],"stiffness") == 0) {
		int order = this->getOrder();
		Matrix *newMatrix = new Matrix(order,order);
		if (newMatrix == 0) {
			cerr << "WARNING SectionForceDeformation::setResponse() - out of memory\n";
			return -1;
		} 
		sectInfo.theMatrix = newMatrix;
		sectInfo.theType = MatrixType;	
		return 3;
	}

    // otherwise response quantity is unknown for the Section class
    else
	return -1;    
}

int 
SectionForceDeformation::getResponse(int responseID, Information &sectInfo)
{
  switch (responseID) {
    case -1:
      return -1;
      
    case 1:
      if (sectInfo.theVector != 0)
	  *(sectInfo.theVector) = this->getSectionDeformation();
     return 0;
      
    case 2:
      if (sectInfo.theVector != 0)
	  *(sectInfo.theVector) = this->getStressResultant();
      return 0;      

	case 3:
		if (sectInfo.theMatrix != 0)
			*(sectInfo.theMatrix) = this->getSectionTangent();
		return 0;

    default:
      return -1;
  }
}
