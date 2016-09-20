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

// $Revision: 4952 $
// $Date: 2012-08-08 22:56:05 -0700 (Wed, 08 Aug 2012) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/frictionBearing/frictionModel/FrictionResponse.cpp $

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the FrictionResponse class implementation

#include <FrictionResponse.h>
#include <FrictionModel.h>


FrictionResponse::FrictionResponse(FrictionModel *frn, int id)
    : Response(), theFriction(frn), responseID(id)
{

}


FrictionResponse::FrictionResponse(FrictionModel *frn, int id, int val)
    : Response(val), theFriction(frn), responseID(id)
{

}


FrictionResponse::FrictionResponse(FrictionModel *frn, int id, double val)
    : Response(val), theFriction(frn), responseID(id)
{

}


FrictionResponse::FrictionResponse(FrictionModel *frn, int id, const ID &val)
    : Response(val), theFriction(frn), responseID(id)
{

}


FrictionResponse::FrictionResponse(FrictionModel *frn, int id, const Vector &val)
    : Response(val), theFriction(frn), responseID(id)
{

}


FrictionResponse::FrictionResponse(FrictionModel *frn, int id, const Matrix &val)
    : Response(val), theFriction(frn), responseID(id)
{

}


FrictionResponse::~FrictionResponse()
{

}


int FrictionResponse::getResponse()
{
    return theFriction->getResponse(responseID, myInfo);
}
