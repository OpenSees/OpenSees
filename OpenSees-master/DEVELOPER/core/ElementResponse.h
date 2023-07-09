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
// $Date: 2009-12-17 23:50:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/ElementResponse.h,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the ElementResponse class interface

#ifndef ElementResponse_h
#define ElementResponse_h

#include <Response.h>
#include <Information.h>

class Element;

class ID;
class Vector;
class Matrix;

class ElementResponse : public Response
{
public:
	ElementResponse(Element *ele, int id);
	ElementResponse(Element *ele, int id, int val);
	ElementResponse(Element *ele, int id, double val);
	ElementResponse(Element *ele, int id, const ID &val);
	ElementResponse(Element *ele, int id, const Vector &val);
	ElementResponse(Element *ele, int id, const Matrix &val);
	ElementResponse(Element *ele, int id, const Vector &val1, const ID &val2);

	~ElementResponse();

	int getResponse(void);
	int getResponseSensitivity(int gradNumber);

private:
	Element *theElement;
	int responseID;
};

#endif
