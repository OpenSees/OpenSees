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

// $Revision: 1.1 $
// $Date: 2009-04-17 23:02:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/FrictionResponse.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the FrictionResponse class implementation
//
// What: "@(#) FrictionResponse.cpp, revA"

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


FrictionResponse::FrictionResponse(FrictionModel *frn, int id, const Tensor &val)
    : Response(val), theFriction(frn), responseID(id)
{

}


FrictionResponse::~FrictionResponse()
{

}


int FrictionResponse::getResponse(void)
{
    return theFriction->getResponse(responseID, myInfo);
}
