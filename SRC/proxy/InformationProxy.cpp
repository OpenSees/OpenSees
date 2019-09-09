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
// $Source: /usr/local/cvs/OpenSees/SRC/element/Information.cpp,v $
                                                                        
                                                                        
// Written: fmk 10/99
// Revised:
//
// Purpose: This file contains the class implementation for Information.

#include <InformationProxy.h>
#include <Information.h>
#include <Element.h>

InformationProxy::InformationProxy(Element* elemenet) : theEle(elemenet)
{
    // does nothing
}

double
InformationProxy::getDouble(int id)
{
	Information* info = new Information();
	theEle->getResponse(id, *info);
	return info->theDouble;
}

Vector*
InformationProxy::getVector(int id)
{
	Information* info = new Information();
	theEle->getResponse(id, *info);
	return info->theVector;
}

Matrix*
InformationProxy::getMatrix(int id)
{
	Information* info = new Information();
	theEle->getResponse(id, *info);
	return info->theMatrix;
}

ID*
InformationProxy::getID(int id)
{
	Information* info = new Information();
	theEle->getResponse(id, *info);
	return info->theID;
}

int
InformationProxy::getInt(int id)
{
	Information* info = new Information();
	theEle->getResponse(id, *info);
	return info->theInt;
}

char*
InformationProxy::getString(int id)
{
	Information* info = new Information();
	theEle->getResponse(id, *info);
	return info->theString;
}

InformationProxy::~InformationProxy()
{
    
}



