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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: Massimo Petracca - ASDEA Software, Italy
// Created: 03/2024

#include <stdlib.h>
#include <math.h>

#include <InitStrainNDMaterial.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <Parameter.h>
#include <OPS_Globals.h>

#include <elementAPI.h>

namespace {

    // a base index to try not to shadow the adpted material's parameters...
    constexpr int param_base_index = 111000;

}

void*
OPS_InitStrainNDMaterial(void)
{
    // Pointer to an ND material that will be returned
    NDMaterial* theMaterial = 0;
    NDMaterial* theOtherMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs != 3 && numArgs != 8) {
        opserr << "Error while creating InitStrain material:" << endln
            << "Only 3 or 8 arguments are required, not " << numArgs << endln
            << "nDMaterial InitStrain tag? otherTag? eps0_11? <eps0_22 eps0_33 eps0_12 eps0_23 eps0_13>" << endln;
        return theMaterial;
    }

    int iData[2];
    int numData = 2;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "Error while creating InitStrain material:" << endln
            << "Cannot read $tag and $otherTag" << numArgs << endln;
        return theMaterial;
    }

    theOtherMaterial = OPS_getNDMaterial(iData[1]);
    if (theOtherMaterial == 0) {
        opserr << "Error while creating InitStrain material:" << endln 
            << "Could not find material with tag: " << iData[1] << endln;
        return theMaterial;
    }

    double eps0[6] = { 0,0,0, 0,0,0 };
    numData = numArgs - 2;
    if (OPS_GetDoubleInput(&numData, eps0) != 0) {
        opserr << "Error while creating InitStrain material:" << endln
            << "Could not read initial strain" << endln;
        return theMaterial;
    }

    // Parsing was successful, allocate the material
    if (numArgs == 3) {
        theMaterial = new InitStrainNDMaterial(iData[0], *theOtherMaterial, eps0[0]);
    }
    else {
        Vector E0(6);
        for (int i = 0; i < 6; ++i) E0(i) = eps0[i];
        theMaterial = new InitStrainNDMaterial(iData[0], *theOtherMaterial, E0);
    }

    if (theMaterial == 0) {
        opserr << "WARNING could not create NDMaterial of type InitStrain\n";
        return theMaterial;
    }

    return theMaterial;
}

InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial& material, const Vector& eps0)
    : NDMaterial(tag, ND_TAG_InitStrainNDMaterial)
    , theMaterial(nullptr)
    , epsInit(eps0)
{
    // get copy of the main material
    theMaterial = material.getCopy("ThreeDimensional");
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- failed to get copy of material (a 3D material is required)\n";
        exit(-1);
    }

    // make sure the input strain vector is of correct size
    if (epsInit.Size() != 6) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- input eps0 vector of incorrect size\n";
        exit(-1);
    }
}

InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial& material, double eps0)
    : NDMaterial(tag, ND_TAG_InitStrainNDMaterial)
    , theMaterial(0)
{
    // get copy of the main material
    theMaterial = material.getCopy("ThreeDimensional");
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- failed to get copy of material (a 3D material is required)\n";
        exit(-1);
    }

    // initialize epsInit
    epsInit.resize(6);
    epsInit.Zero();
    for (int i = 0; i < 3; ++i)
        epsInit(i) = eps0;
}

InitStrainNDMaterial::InitStrainNDMaterial()
    : NDMaterial(0, ND_TAG_InitStrainNDMaterial)
    , theMaterial(nullptr)
    , epsInit(6)
{

}

InitStrainNDMaterial::~InitStrainNDMaterial()
{
    if (theMaterial)
        delete theMaterial;
}

int
InitStrainNDMaterial::setTrialStrain(const Vector& strain)
{
    static Vector total_strain(6);
    total_strain = strain;
    total_strain.addVector(1.0, epsInit, 1.0);
    return theMaterial->setTrialStrain(total_strain);
}

int
InitStrainNDMaterial::setTrialStrain(const Vector& strain, const Vector& /*strainRate*/)
{
    return setTrialStrain(strain);
}

int
InitStrainNDMaterial::setTrialStrainIncr(const Vector& strain)
{
    static Vector strain_from_ele(6);
    strain_from_ele = theMaterial->getStrain();
    strain_from_ele.addVector(1.0, epsInit, -1.0);
    strain_from_ele.addVector(1.0, strain, 1.0);
    return setTrialStrain(strain_from_ele);
}

int
InitStrainNDMaterial::setTrialStrainIncr(const Vector& strain, const Vector& /*strainRate*/)
{
    return setTrialStrainIncr(strain);
}

const Vector&
InitStrainNDMaterial::getStress(void)
{
    return theMaterial->getStress();
}

const Matrix&
InitStrainNDMaterial::getTangent(void)
{
    return theMaterial->getTangent();
}

const Matrix&
InitStrainNDMaterial::getInitialTangent(void)
{
    return theMaterial->getInitialTangent();
}

const Vector&
InitStrainNDMaterial::getStrain(void)
{
    return theMaterial->getStrain();
}

int
InitStrainNDMaterial::commitState(void)
{
    return theMaterial->commitState();
}

int
InitStrainNDMaterial::revertToLastCommit(void)
{
    return theMaterial->revertToLastCommit();
}

int
InitStrainNDMaterial::revertToStart(void)
{
    return theMaterial->revertToStart();
}

double
InitStrainNDMaterial::getRho(void)
{
    return theMaterial->getRho();
}

NDMaterial*
InitStrainNDMaterial::getCopy(void)
{
    InitStrainNDMaterial* theCopy = new InitStrainNDMaterial(getTag(), *theMaterial, epsInit);
    return theCopy;
}

NDMaterial *
InitStrainNDMaterial::getCopy(const char *type)
{
    if (strcmp(type, "ThreeDimensional") == 0)
        return getCopy();
    return NDMaterial::getCopy(type);
}

const char*
InitStrainNDMaterial::getType(void) const
{
    return theMaterial->getType();
}

int InitStrainNDMaterial::getOrder(void) const
{
    return 6;
}

int
InitStrainNDMaterial::sendSelf(int cTag, Channel& theChannel)
{
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::sendSelf() - theMaterial is null, nothing to send" << endln;
        return -1;
    }

    int dbTag = this->getDbTag();

    static ID dataID(3);
    dataID(0) = this->getTag();
    dataID(1) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        theMaterial->setDbTag(matDbTag);
    }
    dataID(2) = matDbTag;
    if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send the ID\n";
        return -1;
    }

    if (theChannel.sendVector(dbTag, cTag, epsInit) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send epsInit\n";
        return -2;
    }

    if (theMaterial->sendSelf(cTag, theChannel) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send the Material\n";
        return -3;
    }

    return 0;
}

int
InitStrainNDMaterial::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int dbTag = this->getDbTag();

    static ID dataID(3);
    if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the ID\n";
        return -1;
    }
    setTag(dataID(0));

    // as no way to change material, don't have to check classTag of the material 
    if (theMaterial == 0) {
        int matClassTag = dataID(1);
        theMaterial = theBroker.getNewNDMaterial(matClassTag);
        if (theMaterial == 0) {
            opserr << "InitStrainNDMaterial::recvSelf() - failed to create Material with classTag "
                << matClassTag << endln;
            return -2;
        }
    }
    theMaterial->setDbTag(dataID(2));

    epsInit.resize(6);
    if (theChannel.recvVector(dbTag, cTag, epsInit) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the epsInit vector\n";
        return -3;
    }

    if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the Material\n";
        return -4;
    }

    return 0;
}

void
InitStrainNDMaterial::Print(OPS_Stream& s, int flag)
{
    s << "InitStrainNDMaterial tag: " << this->getTag() << endln;
    s << "\tMaterial: " << theMaterial->getTag() << endln;
    s << "\tinitital strain: " << epsInit << endln;
}

int
InitStrainNDMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    if (argc > 0) {
        if (strcmp(argv[0], "eps0") == 0) {
            // eps0 is assumed to impose with a single scalar a volumetric stress = I*eps0
            param.setValue(epsInit(0));
            return param.addObject(param_base_index, this);
        }
        else if (strcmp(argv[0], "eps0_11") == 0) {
            param.setValue(epsInit(0));
            return param.addObject(param_base_index + 1, this);
        }
        else if (strcmp(argv[0], "eps0_22") == 0) {
            param.setValue(epsInit(1));
            return param.addObject(param_base_index + 2, this);
        }
        else if (strcmp(argv[0], "eps0_33") == 0) {
            param.setValue(epsInit(2));
            return param.addObject(param_base_index + 3, this);
        }
        else if (strcmp(argv[0], "eps0_12") == 0) {
            param.setValue(epsInit(3));
            return param.addObject(param_base_index + 4, this);
        }
        else if (strcmp(argv[0], "eps0_23") == 0) {
            param.setValue(epsInit(4));
            return param.addObject(param_base_index + 5, this);
        }
        else if (strcmp(argv[0], "eps0_13") == 0) {
            param.setValue(epsInit(5));
            return param.addObject(param_base_index + 6, this);
        }
    }
    return theMaterial->setParameter(argv, argc, param);
}

int InitStrainNDMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {
    case param_base_index:
        epsInit(0) = epsInit(1) = epsInit(2) = info.theDouble;
        return 0;
    case param_base_index + 1:
        epsInit(0) = info.theDouble;
        return 0;
    case param_base_index + 2:
        epsInit(1) = info.theDouble;
        return 0;
    case param_base_index + 3:
        epsInit(2) = info.theDouble;
        return 0;
    case param_base_index + 4:
        epsInit(3) = info.theDouble;
        return 0;
    case param_base_index + 5:
        epsInit(4) = info.theDouble;
        return 0;
    case param_base_index + 6:
        epsInit(5) = info.theDouble;
        return 0;
    default:
        return -1;
    }
}

Response* InitStrainNDMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    return theMaterial->setResponse(argv, argc, output);
}

const Vector&
InitStrainNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
    return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int
InitStrainNDMaterial::commitSensitivity(const Vector& depsdh,
    int gradIndex, int numGrads)
{
    return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
