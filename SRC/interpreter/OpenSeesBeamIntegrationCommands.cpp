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
	functionMap.insert(std::make_pair("NewtoCotes", &OPS_NewtonCotesBeamIntegration));
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
