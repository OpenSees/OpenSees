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
// $Date: 2004-09-01 03:54:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/DamageResponse.cpp,v $
                                                                        
// Written: Arash Altoontash, Gregory Deierlein
// Created: Aug 2002
//
// Description: This file contains the DamageResponse class implementation

#include <DamageResponse.h>
#include <DamageModel.h>

DamageResponse::DamageResponse(DamageModel *dmg, int id):
Response(), theDamage(dmg), responseID(id)
{

}

DamageResponse::DamageResponse(DamageModel *dmg, int id, int val):
Response(val), theDamage(dmg), responseID(id)
{

}

DamageResponse::DamageResponse(DamageModel *dmg, int id, double val):
Response(val), theDamage(dmg), responseID(id)
{

}

DamageResponse::DamageResponse(DamageModel *dmg, int id, const ID &val):
Response(val), theDamage(dmg), responseID(id)
{

}

DamageResponse::DamageResponse(DamageModel *dmg, int id, const Vector &val):
Response(val), theDamage(dmg), responseID(id)
{

}

DamageResponse::DamageResponse(DamageModel *dmg, int id, const Matrix &val):
Response(val), theDamage(dmg), responseID(id)
{

}

DamageResponse::~DamageResponse()
{

}

int
DamageResponse::getResponse(void)
{
  return theDamage->getResponse(responseID, myInfo);
}
