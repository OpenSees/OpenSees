/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

// Written: Alex Baker
// Created: 04/21
//
// Description: This file contains the material testing commands

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <elementAPI.h>
#include <Matrix.h>

static int testType = 0;
static UniaxialMaterial* theTestingUniaxialMaterial = 0;
static NDMaterial* theTestingNDMaterial = 0;
static SectionForceDeformation* theTestingSection = 0;

int OPS_testUniaxialMaterial()
{
    if (OPS_GetNumRemainingInputArgs() != 1) {
        opserr << "testUniaxialMaterial - You must provide a material tag.\n";
        return -1;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
        opserr << "invalid int value\n";
        return -1;
    }

    UniaxialMaterial* mat = OPS_getUniaxialMaterial(tag);

    if (mat == 0) {
        opserr << "testUniaxialMaterial - Material Not Found.\n";
        return -1;
    }

    theTestingUniaxialMaterial = mat;
    testType = 1;

    return 0;
}

int OPS_testNDMaterial()
{
    if (OPS_GetNumRemainingInputArgs() != 1) {
        opserr << "testNDMaterial - You must provide a material tag.\n";
        return -1;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
        opserr << "invalid int value\n";
        return -1;
    }

    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == 0) {
        opserr << "testNDMaterial - Material Not Found.\n";
        return -1;
    }

    theTestingNDMaterial = mat;
    testType = 2;

    return 0;
}

int OPS_testSection()
{
    if (OPS_GetNumRemainingInputArgs() != 1) {
        opserr << "testSection - You must provide a section tag.\n";
        return -1;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
        opserr << "invalid int value\n";
        return -1;
    }

    SectionForceDeformation* mat = OPS_getNDMaterial(tag);

    if (mat == 0) {
        opserr << "testSection - Section Not Found.\n";
        return -1;
    }

    theTestingSection = mat;
    testType = 3;

    return 0;
}

int OPS_setStrain()
{
    if (testType == 0) {
        opserr << "setStrain WARNING no active test\n";
        return -1;
    }
    OPS_setTrialStrain()
    OPS_commitStrain()
}

int OPS_setTrialStrain()
{
    if (testType == 0) {
        opserr << "setTrialStrain WARNING no active test\n";
        return -1;
    }
    else if (testType == 1) {
        if (OPS_GetNumRemainingInputArgs() != 1) {
            opserr << "testUniaxialMaterial - You must provide a strain value.\n";
            return -1;
        }
        double strain;
        int numData = 1;
        if (OPS_GetDoubleInput(&numData, &strain) < 0) {
            opserr << "invalid double value\n";
            return -1;
        }
        theTestingUniaxialMaterial->setTrialStrain(strain);
    }
    else if (testType == 2) {
        opserr << "nDMaterial test - not implemented yet\n";
    }
    else if (testType == 3) {
        opserr << "section test - not implemented yet\n";
    }

    return 0;
}

int OPS_commitStrain()
{
    if (testType == 0) {
        opserr << "commitStrain WARNING no active test\n";
        return -1;
    }
    else if (testType == 1) {
        theTestingUniaxialMaterial->commitState();
    }
    else if (testType == 2) {
        theTestingNDMaterial->commitState();
    }
    else if (testType == 3) {
        theTestingSection->commitState();
    }
}

int OPS_getStrain()
{
    if (testType == 0) {
        opserr << "getStrain WARNING no active test\n";
        return -1;
    }
    else if (testType == 1) {
        double strain = theTestingUniaxialMaterial->getStrain();
    }
    else if (testType == 2) {
        double strain = theTestingNDMaterial->getStrain();
    }
    else if (testType == 3) {
        double strain = theTestingSection->getStrain();
    }

    double strain = material->getStrain();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &strain, true) < 0) {
        opserr << "failed to set strain\n";
        return -1;
    }

    return 0;
}

int OPS_getStress()
{
    UniaxialMaterial* material = theTestingUniaxialMaterial;
    if (material == 0) {
        opserr << "getStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
        return -1;
    }

    double stress = material->getStress();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &stress, true) < 0) {
        opserr << "failed to set stress\n";
        return -1;
    }

    return 0;
}

int OPS_getTangent()
{
    UniaxialMaterial* material = theTestingUniaxialMaterial;
    if (material == 0) {
        opserr << "getStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
        return -1;
    }

    double tangent = material->getTangent();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &tangent, true) < 0) {
        opserr << "failed to set tangent\n";
        return -1;
    }

    return 0;
}

int OPS_getDampTangent()
{
    UniaxialMaterial* material = theTestingUniaxialMaterial;
    if (material == 0) {
        opserr << "getStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command.\n";
        return -1;
    }

    double tangent = material->getDampTangent();

    int numData = 1;

    if (OPS_SetDoubleOutput(&numData, &tangent, true) < 0) {
        opserr << "failed to set damp tangent\n";
        return -1;
    }

    return 0;
}