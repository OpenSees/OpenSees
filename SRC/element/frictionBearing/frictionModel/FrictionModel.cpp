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
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/FrictionModel.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class implementation for FrictionModel.
//
// What: "@(#) FrictionModel.cpp, revA"

#include <FrictionModel.h>
#include <FrictionResponse.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theFrictionModelObjects;

bool OPS_addFrictionModel(FrictionModel *newComponent) {
  return theFrictionModelObjects.addComponent(newComponent);
}

FrictionModel *OPS_getFrictionModel(int tag) {

  TaggedObject *theResult = theFrictionModelObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "FrictionModel *getFrictionModel(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  FrictionModel *theMat = (FrictionModel *)theResult;

  return theMat;
}

void OPS_clearAllFrictionModel(void) {
  theFrictionModelObjects.clearAll();
}


FrictionModel::FrictionModel(int tag, int classTag)
    : TaggedObject(tag), MovableObject(classTag),
    trialN(0.0), trialVel(0.0)
{
    
}


FrictionModel::~FrictionModel()
{
    // does nothing
}


double FrictionModel::getNormalForce(void)
{
    return trialN;
}


double FrictionModel::getVelocity(void)
{
    return trialVel;
}


Response* FrictionModel::setResponse(char **argv, int argc, Information &info)
{
    if (argc == 0) 
        return 0;
    
    // normal force
    else if (strcmp(argv[0],"normalForce") == 0)
        return new FrictionResponse(this, 2, this->getNormalForce());
    
    // velocity
    if (strcmp(argv[0],"velocity") == 0)
        return new FrictionResponse(this, 1, this->getVelocity());
    
    // friction force
    else if (strcmp(argv[0],"frictionForce") == 0)
        return new FrictionResponse(this, 3, this->getFrictionForce());
    
    // friction coeff
    else if (strcmp(argv[0],"frictionCoeff") == 0)
        return new FrictionResponse(this, 4, this->getFrictionCoeff());

    // otherwise unknown
    else
        return 0;
}


int FrictionModel::getResponse(int responseID, Information &info)
{
    // each subclass must implement its own stuff    
    switch (responseID)  {
    case 1:
        info.setDouble(trialN);
        return 0;
        
    case 2:
        info.setDouble(trialVel);
        return 0;      
        
    case 3:
        info.setDouble(this->getFrictionForce());
        return 0;      
                
    case 4:
        info.setDouble(this->getFrictionCoeff());
        return 0;      

    default:      
        return -1;
    }
}
