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
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/frictionBearing/frictionModel/FrictionModel.cpp $

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class implementation for FrictionModel.

#include <FrictionModel.h>
#include <FrictionResponse.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theFrictionModelObjects;


bool OPS_addFrictionModel(FrictionModel *newComponent)
{
    return theFrictionModelObjects.addComponent(newComponent);
}


FrictionModel *OPS_getFrictionModel(int tag)
{
    TaggedObject *theResult = theFrictionModelObjects.getComponentPtr(tag);
    if (theResult == 0) {
        opserr << "FrictionModel *getFrictionModel(int tag) - none found with tag: " << tag << endln;
        return 0;
    }
    FrictionModel *theFrnMdl = (FrictionModel *)theResult;
    
    return theFrnMdl;
}


void OPS_clearAllFrictionModel()
{
    theFrictionModelObjects.clearAll();
}


FrictionModel::FrictionModel(int tag, int classTag)
    : TaggedObject(tag), MovableObject(classTag),
    trialN(0.0), trialVel(0.0)
{
    // does nothing
}


FrictionModel::~FrictionModel()
{
    // does nothing
}


double FrictionModel::getNormalForce()
{
    return trialN;
}


double FrictionModel::getVelocity()
{
    return trialVel;
}


Response* FrictionModel::setResponse(const char **argv, int argc,
    OPS_Stream &theOutputStream)
{
    Response *theResponse = 0;
    
    theOutputStream.tag("FrictionModelOutput");
    theOutputStream.attr("frnMdlType", this->getClassType());
    theOutputStream.attr("frnMdlTag", this->getTag());
    
    // normal force
    if (strcmp(argv[0],"normalForce") == 0 || strcmp(argv[0],"N") == 0 ||
        strcmp(argv[0],"normalFrc") == 0)  {
        theOutputStream.tag("ResponseType", "N");
        return new FrictionResponse(this, 1, this->getNormalForce());
    }
    // velocity
    else if (strcmp(argv[0],"velocity") == 0 || strcmp(argv[0],"vel") == 0)  {
        theOutputStream.tag("ResponseType", "vel");
        return new FrictionResponse(this, 2, this->getVelocity());
    }
    // friction force
    else if (strcmp(argv[0],"frictionForce") == 0 || strcmp(argv[0],"Ff") == 0 ||
             strcmp(argv[0],"frnForce") == 0 || strcmp(argv[0],"frnFrc") == 0)  {
        theOutputStream.tag("ResponseType", "frnFrc");
        return new FrictionResponse(this, 3, this->getFrictionForce());
    }
    // coefficient of friction
    else if (strcmp(argv[0],"frictionCoeff") == 0 || strcmp(argv[0],"mu") == 0 ||
             strcmp(argv[0],"frnCoeff") == 0 || strcmp(argv[0],"COF") == 0)  {
        theOutputStream.tag("ResponseType", "COF");
        return new FrictionResponse(this, 4, this->getFrictionCoeff());
    }
    
    theOutputStream.endTag();
    return theResponse;
}


int FrictionModel::getResponse(int responseID, Information &info)
{
    // each subclass must implement its own stuff
    switch (responseID)  {
    case 1:
        info.setDouble(this->getNormalForce());
        return 0;
        
    case 2:
        info.setDouble(this->getVelocity());
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
