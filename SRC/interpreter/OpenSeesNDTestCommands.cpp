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
                                                                        
// Written: Minjie

// Description: command to create timeseries

#include "TimeSeries.h"
#include "map"
#include "string.h"
#include "elementAPI.h"
#include "vector"
#include "PathTimeSeries.h"
#include "PathSeries.h"
#include "Vector.h"
#include "Matrix.h"
#include "NDMaterial.h"
#include "Information.h"


void* OPS_NDSetStrain()
{	
	int tag = 0;
	double strains[6];

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
	
	numdata = 6;
	if (OPS_GetDoubleInput(&numdata, strains) < 0) return 0;


    NDMaterial* mat = OPS_getNDMaterial(tag);

	Vector new_strain(6);

	new_strain(0) = strains[0];
	new_strain(1) = strains[1];
	new_strain(2) = strains[2];
	new_strain(3) = strains[3];
	new_strain(4) = strains[4];
	new_strain(5) = strains[5];

	mat->setTrialStrain(new_strain);

	return 0;
}

void* OPS_NDPrintStrain()
{
	int tag = 0;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    const Vector &strain = mat->getStrain();

	opserr << "strain = " << strain;

	return 0;
}

void* OPS_NDPrintStress()
{
	int tag = 0;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    const Vector &stress = mat->getStress();

	opserr << "stress = " << stress;

	return 0;
}


void* OPS_NDGetStrain()
{
	int tag = 0;
	int size = 6;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    const Vector &strain = mat->getStrain();

	std::vector<double> values(size);
	for (int i=0; i<6; i++) {
	    values[i] = strain(i);
	}
	if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	    opserr<<"WARNING OPS_NDGetStress - failed to set double inputs\n";
	    return 0;
	}
    

	return 0;
}



void* OPS_NDGetStress()
{
	int tag = 0;
	int size = 6;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    const Vector &stress = mat->getStress();

	std::vector<double> values(size);
	for (int i=0; i<6; i++) {
	    values[i] = stress(i);
	}
	if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	    opserr<<"WARNING OPS_NDGetStress - failed to set double inputs\n";
	    return 0;
	}
    

	return 0;
}


void* OPS_NDGetTangentStiffness()
{
	int tag = 0;
	int size = 36;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    const Matrix &stiffness = mat->getTangent();

	std::vector<double> values(size);
	for (int i=0; i<6; i++) {
		for (int j = 0; j < 6; j++)
		{
			double one_value = stiffness(i,j);
	    	values[6*i+j] = one_value;
		}
	}
	if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	    opserr<<"WARNING OPS_NDGetStress - failed to set double inputs\n";
	    return 0;
	}
    

	return 0;
}



void* OPS_NDCommitState()
{	
	int tag = 0;
	double strains[6];
	int size = 6;
	double stressdata[6] = {1 , 2 , 3, 4, 5, 6};

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
	

    NDMaterial* mat = OPS_getNDMaterial(tag);

	mat->commitState();

	return 0;
}


void* OPS_NDUpdateIntegerParameter()
{	
	int tag = 0;
	int responseID = 0;
	int theNewIntegerParameterValue = 0;
	
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
	if (OPS_GetIntInput(&numdata, &responseID) < 0) return 0;
	if (OPS_GetIntInput(&numdata, &theNewIntegerParameterValue) < 0) return 0;
	

    NDMaterial* mat = OPS_getNDMaterial(tag);

    Information info;

    info.theInt = theNewIntegerParameterValue;

	mat->updateParameter(responseID, info);

	return 0;
}


void* OPS_NDUpdateDoubleParameter()
{	
	int tag = 0;
	int responseID = 0;
	double theNewDoubleParameterValue = 0;
	
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
	if (OPS_GetIntInput(&numdata, &responseID) < 0) return 0;
	if (OPS_GetDoubleInput(&numdata, &theNewDoubleParameterValue) < 0) return 0;
	

    NDMaterial* mat = OPS_getNDMaterial(tag);

    Information info;

    info.theDouble = theNewDoubleParameterValue;

	mat->updateParameter(responseID, info);

	return 0;
}


namespace {
    
    struct char_cmp { 
	bool operator () (const char *a,const char *b) const 
	    {
		return strcmp(a,b)<0;
	    } 
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;


    static OPS_ParsingFunctionMap functionMap;

    static int setUpFunctions(void)
    {
	functionMap.insert(std::make_pair("SetStrain", &OPS_NDSetStrain));
	functionMap.insert(std::make_pair("CommitState", &OPS_NDCommitState));
	functionMap.insert(std::make_pair("PrintStress", &OPS_NDPrintStress));
	functionMap.insert(std::make_pair("PrintStrain", &OPS_NDPrintStrain));
	functionMap.insert(std::make_pair("GetStrain", &OPS_NDGetStrain));
	functionMap.insert(std::make_pair("GetStress", &OPS_NDGetStress));
	functionMap.insert(std::make_pair("GetTangentStiffness", &OPS_NDGetTangentStiffness));
	functionMap.insert(std::make_pair("UpdateIntegerParameter", &OPS_NDUpdateIntegerParameter));
	functionMap.insert(std::make_pair("UpdateDoubleParameter", &OPS_NDUpdateDoubleParameter));
      
	return 0;
    }
}

int
OPS_NDTest()
{
    static bool initDone = false;
    if (initDone == false) {
	setUpFunctions();
	initDone = true;
    }

    // Identify what specific command of Patch we're calling
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"WARNING too few arguments: NDTest cmd? \n";
	opserr<<" available commands: SetStrain|CommitState|GetStrain|GetStress \n";
	return -1;
    }

    const char* type = OPS_GetString();
    
    OPS_ParsingFunctionMap::const_iterator iter = functionMap.find(type);
    if (iter == functionMap.end()) {
	opserr<<"WARNING NDTest type " << type << " is unknown\n";
	return -1;
    }

    // Call the function
    (*iter->second)();
    
    return 0;

}
