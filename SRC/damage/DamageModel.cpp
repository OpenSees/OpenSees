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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-04-14 22:38:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/DamageModel.cpp,v $                                                                        
                                                                        
// Written: Arash Altoontash, Gregory Deierlein
// Created: 08/02
// Revision: AA
//
// Description: This file contains the class implementation for DamageModel

#include <DamageModel.h>
#include <Response.h>
#include <DamageResponse.h>
#include <string.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>
#include <api/runtimeAPI.h>

static MapOfTaggedObjects theDamageModelObjects;

bool OPS_addDamageModel(DamageModel *newComponent) {
  return theDamageModelObjects.addComponent(newComponent);
}

DamageModel *OPS_getDamageModel(int tag) {

  TaggedObject *theResult = theDamageModelObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "DamageModel *getDamageModel(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  DamageModel *theMat = (DamageModel *)theResult;

  return theMat;
}

void
OPS_ADD_RUNTIME_VXV(OPS_clearAllDamageModel)
{
  theDamageModelObjects.clearAll();
}

DamageModel::DamageModel(int tag, int clasTag)
:TaggedObject(tag), MovableObject(clasTag)
{

}


DamageModel::~DamageModel()
{
  // does nothing


}

Response*
DamageModel::setResponse(const char **argv, int argc, OPS_Stream &stream)
{
  if ( strcmp(argv[0],"damage") == 0 || strcmp(argv[0],"damageindex") == 0 )
    return new DamageResponse( this , 1 , 0.0 );
  
  else 
    return 0;
  
}

int 
DamageModel::getResponse(int responseID, Information &info)
{
  switch (responseID) {
  case -1:
    return -1;
    
  case 1:
    return info.setDouble( this->getDamage() );
    
  default:
    return -1;
  }
}

