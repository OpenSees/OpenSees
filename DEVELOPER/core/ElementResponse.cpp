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
                                                                        
// $Revision: 1.6 $
// $Date: 2009-12-17 23:50:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/ElementResponse.cpp,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the ElementResponse class implementation

#include <ElementResponse.h>
#include <Element.h>

ElementResponse::ElementResponse(Element *ele, int id):
Response(), theElement(ele), responseID(id)
{

}

ElementResponse::ElementResponse(Element *ele, int id, int val):
Response(val), theElement(ele), responseID(id)
{

}

ElementResponse::ElementResponse(Element *ele, int id, double val):
Response(val), theElement(ele), responseID(id)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const ID &val):
Response(val), theElement(ele), responseID(id)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Vector &val):
Response(val), theElement(ele), responseID(id)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Matrix &val):
Response(val), theElement(ele), responseID(id)
{

}

ElementResponse::ElementResponse(Element *ele, int id, const Vector &val1, const ID &val2)
 :Response(val1, val2), theElement(ele), responseID(id)
{

}

ElementResponse::~ElementResponse()
{

}

int
ElementResponse::getResponse(void)
{
  return theElement->getResponse(responseID, myInfo);
}

int
ElementResponse::getResponseSensitivity(int gradNumber)
{
  return theElement->getResponseSensitivity(responseID, gradNumber, myInfo);
}
