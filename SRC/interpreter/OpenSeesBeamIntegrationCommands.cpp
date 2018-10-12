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

// Description: command to create beam integration

#include <string.h>
#include <elementAPI.h>
#include <map>
#include <ID.h>
#include <BeamIntegration.h>

void* OPS_BeamIntegration(int&,ID&);
void* OPS_LobattoBeamIntegration(int&,ID&);
void* OPS_LegendreBeamIntegration(int&,ID&);
void* OPS_NewtonCotesBeamIntegration(int&,ID&);
void* OPS_RadauBeamIntegration(int&,ID&);
void* OPS_TrapezoidalBeamIntegration(int&,ID&);
void* OPS_CompositeSimpsonBeamIntegration(int&,ID&);
void* OPS_UserDefinedBeamIntegration(int&,ID&);
void* OPS_FixedLocationBeamIntegration(int&,ID&);
void* OPS_LowOrderBeamIntegration(int&,ID&);
void* OPS_MidDistanceBeamIntegration(int&,ID&);
void* OPS_UserHingeBeamIntegration(int&,ID&);
void* OPS_HingeMidpointBeamIntegration(int&,ID&);
void* OPS_HingeRadauBeamIntegration(int&,ID&);
void* OPS_HingeRadauTwoBeamIntegration(int&,ID&);
void* OPS_HingeEndpointBeamIntegration(int&,ID&);

namespace {
    struct char_cmp { 
	bool operator () (const char *a,const char *b) const 
	    {
		return strcmp(a,b)<0;
	    } 
    };

    typedef std::map<const char *, void *(*)(int&,ID&), char_cmp> OPS_ParsingFunctionMap;


    static OPS_ParsingFunctionMap functionMap;

    int setUpFunctions(void) {
	functionMap.insert(std::make_pair("Lobatto", &OPS_LobattoBeamIntegration));
	functionMap.insert(std::make_pair("Legendre", &OPS_LegendreBeamIntegration));
	functionMap.insert(std::make_pair("NewtonCotes", &OPS_NewtonCotesBeamIntegration));
	functionMap.insert(std::make_pair("Radau", &OPS_RadauBeamIntegration));
	functionMap.insert(std::make_pair("Trapezoidal", &OPS_TrapezoidalBeamIntegration));
	functionMap.insert(std::make_pair("CompositeSimpson", &OPS_CompositeSimpsonBeamIntegration));
	functionMap.insert(std::make_pair("UserDefined", &OPS_UserDefinedBeamIntegration));
	functionMap.insert(std::make_pair("FixedLocation", &OPS_FixedLocationBeamIntegration));
	functionMap.insert(std::make_pair("LowOrder", &OPS_LowOrderBeamIntegration));
	functionMap.insert(std::make_pair("MidDistance", &OPS_MidDistanceBeamIntegration));
	functionMap.insert(std::make_pair("UserHinge", &OPS_UserHingeBeamIntegration));
	functionMap.insert(std::make_pair("HingeMidpoint", &OPS_HingeMidpointBeamIntegration));
	functionMap.insert(std::make_pair("HingeRadau", &OPS_HingeRadauBeamIntegration));
	functionMap.insert(std::make_pair("HingeRadauTwo", &OPS_HingeRadauTwoBeamIntegration));
	functionMap.insert(std::make_pair("HingeEndpoint", &OPS_HingeEndpointBeamIntegration));
	return 0;
    }
    
}

int OPS_BeamIntegration()
{
    static bool initDone = false;
    if (initDone == false) {
	setUpFunctions();
	initDone = true;
    }
    
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"WARNING too few arguments: beamIntegration type? tag? ...\n";
	return -1;
    }

    const char* type = OPS_GetString();

    OPS_ParsingFunctionMap::const_iterator iter = functionMap.find(type);
    if (iter == functionMap.end()) {
	opserr<<"WARNING beam integration type " << type << " is unknown\n";
	return -1;
    }

    int iTag;
    ID secTags;
    BeamIntegration* bi = (BeamIntegration*)(*iter->second)(iTag,secTags);
    if (bi == 0) {
	return -1;
    }
    BeamIntegrationRule* rule = new BeamIntegrationRule(iTag,bi,secTags);
    if (rule == 0) {
	return -1;
    }

    // add it
    if (OPS_addBeamIntegrationRule(rule) == false) {
	opserr << "WARNING failed to add BeamIntegration\n";
	delete rule;
	return -1;
    }

    return 0;
}
