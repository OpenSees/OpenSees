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
// $Date: 2000-12-18 11:35:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/ElementResponse.cpp,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the ElementResponse class implementation

#include <ElementResponse.h>
#include <Element.h>

ElementResponse::ElementResponse(Element *ele, int id):
Response(), theElement(ele), responseID(id), eleInfo()
{

}

ElementResponse::ElementResponse(Element *ele, int id, int val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, double val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const ID &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Vector &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Matrix &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Tensor &val):
Response(), theElement(ele), responseID(id), eleInfo(val)
{

}

ElementResponse::~ElementResponse()
{

}

int
ElementResponse::getResponse(void)
{
	return theElement->getResponse(responseID, eleInfo);
}

void
ElementResponse::Print(ostream &s, int flag)
{
	eleInfo.Print(s, flag);
}
