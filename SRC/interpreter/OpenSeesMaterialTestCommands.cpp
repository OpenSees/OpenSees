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
#include <Vector.h>
#include <DummyStream.h>
#include <Response.h>

static int testType = 0;
static UniaxialMaterial* theTestingUniaxialMaterial = 0;
static NDMaterial* theTestingNDMaterial = 0;
static SectionForceDeformation* theTestingSection = 0;

void OPS_clearAllTestMaterials(void) 
{
    testType = 0;
    theTestingUniaxialMaterial = 0;
    theTestingNDMaterial = 0;
    theTestingSection = 0;
}

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

    SectionForceDeformation* mat = OPS_getSectionForceDeformation(tag);

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
    // check if test is active
    if (testType == 0) {
        opserr << "setStrain WARNING no active test\n";
        return -1;
    }
    // check if any input values
    if (OPS_GetNumRemainingInputArgs() == 0) {
        opserr << "setStrain WARNING must provide strain values.\n";
        return -1;
    }
    // get strain values from input
    double strainData;
    int numData = OPS_GetNumRemainingInputArgs();
    if (OPS_GetDoubleInput(&numData, &strainData) < 0) {
        opserr << "invalid double value\n";
        return -1;
    }
    // switch for test type
    if (testType == 1) {
        if (numData != 1) {
            opserr << "setTrialStrain WARNING wrong number of arguments.\n";
            return -1;
        }
        theTestingUniaxialMaterial->setTrialStrain(strainData);
        theTestingUniaxialMaterial->commitState();
    }
    else if (testType == 2) {
        Vector strains(&strainData, numData);
        theTestingNDMaterial->setTrialStrain(strains);
        theTestingNDMaterial->commitState();
    }
    else if (testType == 3) {
        Vector strains(&strainData, numData);
        theTestingSection->setTrialSectionDeformation(strains);
        theTestingSection->commitState();
    }

    return 0;
}

int OPS_setTrialStrain()
{
    // check if test is active
    if (testType == 0) {
        opserr << "setTrialStrain WARNING no active test\n";
        return -1;
    }
    // check if any input values
    if (OPS_GetNumRemainingInputArgs() == 0) {
        opserr << "setTrialStrain WARNING must provide strain values.\n";
        return -1;
    }
    // get strain values from input
    double strainData;
    int numData = OPS_GetNumRemainingInputArgs();
    if (OPS_GetDoubleInput(&numData, &strainData) < 0) {
        opserr << "invalid double value\n";
        return -1;
    }
    // switch for test type
    if (testType == 1) {
        if (numData != 1) {
            opserr << "setTrialStrain WARNING wrong number of arguments.\n";
            return -1;
        }
        theTestingUniaxialMaterial->setTrialStrain(strainData);
    }
    else if (testType == 2) {
        Vector strains(&strainData, numData);
        theTestingNDMaterial->setTrialStrain(strains);
    }
    else if (testType == 3) {
        Vector strains(&strainData, numData);
        theTestingSection->setTrialSectionDeformation(strains);
    }

    return 0;
}

int OPS_commitStrain()
{
    // switch for test type
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

    return 0;
}

int OPS_getStrain()
{
    // switch for test type
    if (testType == 0) {
        opserr << "getStrain WARNING no active test\n";
        return -1;
    }
    else if (testType == 1) {
        double strain = theTestingUniaxialMaterial->getStrain();
        int numData = 1;
        // set output
        if (OPS_SetDoubleOutput(&numData, &strain, true) < 0) {
            opserr << "failed to get strain\n";
            return -1;
        }
    }
    else if (testType == 2) {
        const Vector& strain = theTestingNDMaterial->getStrain();
        int numData = strain.Size();
        // convert to double
        double* data = new double[numData];
        for (int i = 0; i < numData; i++)
            data[i] = strain(i);
        // set output
        if (OPS_SetDoubleOutput(&numData, data, true) < 0) {
            opserr << "failed to get strain\n";
            return -1;
        }
        delete[] data;
    }
    else if (testType == 3) {
        const Vector& strain = theTestingSection->getSectionDeformation();
        int numData = strain.Size();
        // convert to double
        double* data = new double[numData];
        for (int i = 0; i < numData; i++)
            data[i] = strain(i);
        // set output
        if (OPS_SetDoubleOutput(&numData, data, true) < 0) {
            opserr << "failed to get strain\n";
            return -1;
        }
        delete[] data;
    }

    return 0;
}

int OPS_getStress()
{
    // switch for test type
    if (testType == 0) {
        opserr << "getStress WARNING no active test\n";
        return -1;
    }
    else if (testType == 1) {
        double stress = theTestingUniaxialMaterial->getStress();
        int numData = 1;
        // set output
        if (OPS_SetDoubleOutput(&numData, &stress, true) < 0) {
            opserr << "failed to get stress\n";
            return -1;
        }
    }
    else if (testType == 2) {
        const Vector& stress = theTestingNDMaterial->getStress();
        int numData = stress.Size();
        // convert to double
        double* data = new double[numData];
        for (int i = 0; i < numData; i++)
            data[i] = stress(i);
        // set output
        if (OPS_SetDoubleOutput(&numData, data, true) < 0) {
            opserr << "failed to get stress\n";
            return -1;
        }
        delete[] data;
    }
    else if (testType == 3) {
        const Vector& stress = theTestingSection->getStressResultant();
        int numData = stress.Size();
        // convert to double
        double* data = new double[numData];
        for (int i = 0; i < numData; i++)
            data[i] = stress(i);
        // set output
        if (OPS_SetDoubleOutput(&numData, data, true) < 0) {
            opserr << "failed to get stress\n";
            return -1;
        }
        delete[] data;
    }

    return 0;
}

int OPS_getTangent()
{
    // switch for test type
    if (testType == 0) {
        opserr << "getTangent WARNING no active test\n";
        return -1;
    }
    else if (testType == 1) {
        double tangent = theTestingUniaxialMaterial->getTangent();
        int numData = 1;
        // set output
        if (OPS_SetDoubleOutput(&numData, &tangent, true) < 0) {
            opserr << "failed to get tangent\n";
            return -1;
        }
    }
    else if (testType == 2) {
        const Matrix& tangent = theTestingNDMaterial->getTangent();
        int numData = tangent.noCols() * tangent.noRows();
        // convert to double
        double* data = new double[numData];
        int k = 0;
        for (int i = 0; i < tangent.noRows(); i++)
            for (int j = 0; j < tangent.noCols(); j++) {
                data[k] = tangent(i, j);
                k++;
            }
        // set output
        if (OPS_SetDoubleOutput(&numData, data, true) < 0) {
            opserr << "failed to get tangent\n";
            return -1;
        }
        delete[] data;
    }
    else if (testType == 3) {
        const Matrix& tangent = theTestingSection->getSectionTangent();
        int numData = tangent.noCols() * tangent.noRows();
        // convert to double
        double* data = new double[numData];
        int k = 0;
        for (int i = 0; i < tangent.noRows(); i++)
            for (int j = 0; j < tangent.noCols(); j++) {
                data[k] = tangent(i, j);
                k++;
            }
        // set output
        if (OPS_SetDoubleOutput(&numData, data, true) < 0) {
            opserr << "failed to get tangent\n";
            return -1;
        }
        delete[] data;
    }

    return 0;
}

int OPS_getDampTangent()
{
    // switch for test type
    if (testType == 0) {
        opserr << "getDampTangent WARNING no active test\n";
        return -1;
    }
    else if (testType == 1) {
        double dampTangent = theTestingUniaxialMaterial->getDampTangent();
        int numData = 1;
        if (OPS_SetDoubleOutput(&numData, &dampTangent, true) < 0) {
            opserr << "failed to get damping tangent\n";
            return -1;
        }
    }
    else if (testType == 2) {
        opserr << "damping tangent not implemented for nDMaterials\n";
        return -1;
    }
    else if (testType == 3) {
        opserr << "damping tangent not implemented for sections\n";
        return -1;
    }

    return 0;
}

static Vector responseData(0);

int OPS_getResponse()
{
    // switch for test type
    if (testType == 0) {
        opserr << "getResponse WARNING no active test\n";
        return -1;
    }
    // check number of input args
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc > 0) {
        char** argv = new char*[argc];
        char buffer[128];
        for (int i = 0; i < argc; i++) {
            argv[i] = new char[128];
            // Turn everything in to a string for setResponse
            strcpy(argv[i], OPS_GetStringFromAll(buffer, 128));
        }
        // set the response
        DummyStream dummy;
        Response* theResponse{};
        if (testType == 1) {
            theResponse = theTestingUniaxialMaterial->setResponse((const char**)argv, argc, dummy);
        }
        else if (testType == 2) {
            theResponse = theTestingNDMaterial->setResponse((const char**)argv, argc, dummy);
        }
        else if (testType == 3) {
            theResponse = theTestingSection->setResponse((const char**)argv, argc, dummy);
        }
        // get clean up and check for null cases
        for (int i = 0; i < argc; i++)
            delete [] argv[i];
        delete [] argv;
        if (theResponse == 0) {
            return 0;
        }
        if (theResponse->getResponse() < 0) {
            delete theResponse;
            return 0;
        }
        // get data from response
        Information &responseInfo = theResponse->getInformation();
        responseData = responseInfo.getData();
        delete theResponse;
        int size = responseData.Size();
        // if data not empty, set as output
        if (size != 0) {
            double* newdata = new double[size];
            for (int i = 0; i < size; i++) {
                newdata[i] = responseData(i);
            }
            if (OPS_SetDoubleOutput(&size, newdata, false) < 0) {
                opserr << "WARNING failed to get test response\n";
                delete [] newdata;
                return -1;
            }
            delete [] newdata;
        }
        else {
            double* newdata = 0;
            if (OPS_SetDoubleOutput(&size, newdata, false) < 0) {
                opserr << "WARNING failed to get test response\n";
                return -1;
            }
        }
    }
    else {
        int size = 0;
        double* newdata = 0;
        if (OPS_SetDoubleOutput(&size, newdata, false) < 0) {
            opserr << "WARNING failed to get test response\n";
            return -1;
        }
    }

    return 0;
}