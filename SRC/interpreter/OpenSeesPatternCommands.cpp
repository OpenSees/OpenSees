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

// Description: command to create pattern

#include <elementAPI.h>
#include <Domain.h>
#include <NodalLoad.h>
#include <Beam2dPartialUniformLoad.h>
#include <Beam2dUniformLoad.h>
#include <Beam3dUniformLoad.h>
#include <Beam2dPointLoad.h>
#include <Beam3dPointLoad.h>
#include <BrickSelfWeight.h>
#include <SurfaceLoader.h>
#include <SelfWeight.h>
#include <Beam2dThermalAction.h>
#include <Beam2dTempLoad.h>
#include <SP_Constraint.h>
#include <LoadPattern.h>
#include <MultiSupportPattern.h>
#include <UniformExcitation.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <TimeSeriesIntegrator.h>
#include <TimeSeries.h>
#include <GroundMotion.h>
#include <vector>
#include <InterpolatedGroundMotion.h>

void* OPS_LoadPattern();
void* OPS_UniformExcitationPattern();
void* OPS_MultiSupportPattern();
void* OPS_TimeSeriesIntegrator();

namespace {
    static LoadPattern* theActiveLoadPattern = 0;
    static UniformExcitation* theActiveUniformPattern = 0;
    static MultiSupportPattern* theActiveMultiSupportPattern = 0;
    static int eleLoadTag = 0;
}


int OPS_Pattern()
{
    theActiveMultiSupportPattern = 0;
    theActiveUniformPattern = 0;
    theActiveLoadPattern = 0;

    // num args
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"WARNING insufficient args: pattern type ...\n";
	return -1;
    }

    // create pattern
    const char* type = OPS_GetString();
    LoadPattern* pattern = 0;
    if (strcmp(type, "Plain") == 0) {

	theActiveLoadPattern = (LoadPattern*)OPS_LoadPattern();
	pattern = theActiveLoadPattern;

    } else if (strcmp(type, "UniformExcitation") == 0) {

	theActiveUniformPattern = (UniformExcitation*)OPS_UniformExcitationPattern();
	pattern = theActiveUniformPattern;

    } else if (strcmp(type, "MultipleSupport") == 0) {

	theActiveMultiSupportPattern = (MultiSupportPattern*)OPS_MultiSupportPattern();
	pattern = theActiveMultiSupportPattern;

    } else {
	opserr<<"WARNING unknown pattern type"<<type<<"\n";
	return -1;
    }

    if (pattern == 0) {
	opserr<<"WARNING failed to create pattern\n";
	return -1;
    }

    // add to domain
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) {
	opserr<<"WARNING no domain is created\n";
	return -1;
    }
    if (theDomain->addLoadPattern(pattern) == false) {
	opserr<<"WARNING failed to add pattern to domain\n";
	delete pattern;
	pattern = 0;
	return -1;
    }

    return 0;
}

int OPS_NodalLoad()
{
    Domain* theDomain = OPS_GetDomain();
    if(theDomain == 0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }

    int ndf = OPS_GetNumRemainingInputArgs()-1;
    if(ndf < 1) {
	opserr<<"insufficient number of args\n";
	return -1;
    }

    // get node tag
    int ndtag;
    int numData = 1;
    if(OPS_GetIntInput(&numData, &ndtag) < 0) {
	opserr << "WARNING invalid node tag\n";
	return -1;
    }

    // get load vector
    Vector forces(ndf);
    if(OPS_GetDoubleInput(&ndf, &forces(0)) < 0) {
	opserr << "WARNING invalid load vector\n";
	return -1;
    }

    // get options
    bool isLoadConst = false;
    bool userPattern = false;
    int loadPatternTag = 0;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-const") == 0) {
	    isLoadConst = true;
	} else if(strcmp(type,"-pattern") == 0) {
	    int numData = 1;
	    if(OPS_GetIntInput(&numData, &loadPatternTag) < 0) {
		return -1;
	    }
	    userPattern = true;
	}
    }

    // get the current pattern tag
    LoadPattern* currPattern = theActiveLoadPattern;
    if(userPattern == false) {
	if(currPattern==0) {
	    opserr<<"WARNING: no current load pattern is set\n";
	    return -1;
	}
	loadPatternTag = currPattern->getTag();
    }

    // create the load
    static int nodeLoadTag = 0;
    NodalLoad* theLoad = new NodalLoad(nodeLoadTag++, ndtag, forces, isLoadConst);
    if(theLoad == 0) return -1;

    // add load to domain
    if(theDomain->addNodalLoad(theLoad,loadPatternTag) == false) {
	opserr<<"WARNING: failed to add nodal load to domain\n";
	delete theLoad;
	return -1;
    }

    return 0;

}

int OPS_ElementalLoad()
{

    if (theActiveLoadPattern == 0) {
	opserr << "WARNING no active load pattern - eleLoad\n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int ndm = OPS_GetNDM();
    ElementalLoad* theLoad = 0;

    // look for -type, -ele, -range
    int locType = -1;
    int locEle = -1;
    int locRange = -1;
    int currentLoc = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* opt = OPS_GetString();
	if (strcmp(opt, "-type") == 0) {
	    locType = currentLoc;
	} else if (strcmp(opt, "-ele") == 0) {
	    locEle = currentLoc;
	} else if (strcmp(opt, "-range") == 0) {
	    locRange = currentLoc;
	}
	currentLoc++;
    }
    if (locType < 0) {
	opserr << "WARNING eleLoad - no -type option\n";
	return -1;
    }

    // we first create an ID containing the ele tags of all elements
    // for which the load applies.
    ID theEleTags(0,16);
    if (locEle > 0) {
	OPS_ResetCurrentInputArg(locEle+1);
	while(OPS_GetNumRemainingInputArgs() > 0) {
	    int tag;
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &tag) < 0) {
		break;
	    }
	    theEleTags.insert(tag);
	}
    }
    if (locRange > 0) {
	OPS_ResetCurrentInputArg(locRange+1);
	if (OPS_GetNumRemainingInputArgs() > 1) {
	    int tags[2];
	    int numdata = 2;
	    if (OPS_GetIntInput(&numdata, tags) < 0) {
		opserr<<"WARNING failed to read tag range\n";
		return -1;
	    }
	    for (int tag=tags[0]; tag<=tags[1]; tag++) {
		theEleTags.insert(tag);
	    }
	}
    }

    // we then create the load
    OPS_ResetCurrentInputArg(locType+1);
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"WARNING no load type is given\n";
	return -1;
    }
    const char* type = OPS_GetString();
    if (strcmp(type,"-beamUniform") == 0 ||
	strcmp(type,"beamUniform") == 0) {

	if (ndm == 2) {
	    // wta, waa, aL, bL, wtb, wab
        double data[6] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
	    int numdata = OPS_GetNumRemainingInputArgs();
	    if (numdata < 1) {
		opserr<<"WARNING eleLoad - beamUniform want Wya <Wxa> <aL> <bL> <Wyb> <Wxb>\n";
		return -1;
	    }
	    if (numdata > 6) numdata = 6;
	    if (OPS_GetDoubleInput(&numdata, data) < 0) {
		opserr<<"WARNING eleLoad - invalid value for beamUniform\n";
		return -1;
	    }
	    for (int i=0; i<theEleTags.Size(); i++) {
		if (numdata == 3 || numdata == 4) {
		  data[4] = data[0];
		  data[5] = data[1];
		}
		if (data[2] > 0.0 || data[3] < 1.0 || numdata > 4)
		    theLoad = new Beam2dPartialUniformLoad(eleLoadTag, data[0], data[4], data[1], data[5], data[2], data[3], theEleTags(i));
		else
		    theLoad = new Beam2dUniformLoad(eleLoadTag, data[0], data[1], theEleTags(i));

		if (theLoad == 0) {
		    opserr << "WARNING eleLoad - out of memory creating load of type " << type;
		    return -1;
		}

		// get the current pattern tag if no tag given in i/p
		int loadPatternTag = theActiveLoadPattern->getTag();

		// add the load to the domain
		if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
		    opserr << "WARNING eleLoad - could not add following load to domain:\n ";
		    opserr << theLoad;
		    delete theLoad;
		    return -1;
		}
		eleLoadTag++;
	    }

	    return 0;
	}

	else if (ndm == 3) {
	    // wy, wz, wx
	    double data[3] = {0.0, 0.0, 0.0};
	    int numdata = OPS_GetNumRemainingInputArgs();
	    if (numdata < 2) {
		opserr<<"WARNING eleLoad - beamUniform want Wy Wz <Wx>\n";
		return -1;
	    }
	    if (numdata > 3) numdata = 3;
	    if (OPS_GetDoubleInput(&numdata, data) < 0) {
		opserr<<"WARNING eleLoad - invalid value for beamUniform\n";
		return -1;
	    }
	    for (int i=0; i<theEleTags.Size(); i++) {
		theLoad = new Beam3dUniformLoad(eleLoadTag, data[0], data[1], data[2], theEleTags(i));

		if (theLoad == 0) {
		    opserr << "WARNING eleLoad - out of memory creating load of type " << type;
		    return -1;
		}

		// get the current pattern tag if no tag given in i/p
		int loadPatternTag = theActiveLoadPattern->getTag();

		// add the load to the domain
		if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
		    opserr << "WARNING eleLoad - could not add following load to domain:\n ";
		    opserr << theLoad;
		    delete theLoad;
		    return -1;
		}
		eleLoadTag++;
	    }

	    return 0;

	}
	else {
	    opserr << "WARNING eleLoad beamUniform currently only valid only for ndm=2 or 3\n";
	    return -1;
	}

    } else if (strcmp(type,"-beamPoint") == 0 ||
	       strcmp(type,"beamPoint") == 0 ) {

	if (ndm == 2) {
	    // P, x, N
	    double data[3] = {0.0, 0.0, 0.0};
	    int numdata = OPS_GetNumRemainingInputArgs();
	    if (numdata < 2) {
		opserr<<"WARNING eleLoad - beamPoint want Py xL <Px>\n";
		return -1;
	    }
	    if (numdata > 3) numdata = 3;
	    if (OPS_GetDoubleInput(&numdata, data) < 0) {
		opserr<<"WARNING eleLoad - invalid value for beamPoint\n";
		return -1;
	    }

	    if (data[1] < 0.0 || data[1] > 1.0) {
		opserr << "WARNING eleLoad - invalid xDivL of " << data[1];
		opserr << " for beamPoint (valid range [0.0, 1.0]\n";
		return -1;
	    }


	    for (int i=0; i<theEleTags.Size(); i++) {
		theLoad = new Beam2dPointLoad(eleLoadTag, data[0], data[1],
					      theEleTags(i), data[2]);

		if (theLoad == 0) {
		    opserr << "WARNING eleLoad - out of memory creating load of type " << type;
		    return -1;
		}

		// get the current pattern tag if no tag given in i/p
		int loadPatternTag = theActiveLoadPattern->getTag();

		// add the load to the domain
		if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
		    opserr << "WARNING eleLoad - could not add following load to domain:\n ";
		    opserr << theLoad;
		    delete theLoad;
		    return -1;
		}
		eleLoadTag++;
	    }

	    return 0;

	}
	else if (ndm == 3) {
	    // Py, Pz, x, N
	    double data[4] = {0.0, 0.0, 0.0, 0.0};
	    int numdata = OPS_GetNumRemainingInputArgs();
	    if (numdata < 3) {
		opserr<<"WARNING eleLoad - beamPoint want Py Pz xL <Px>\n";
		return -1;
	    }
	    if (numdata > 4) numdata = 4;
	    if (OPS_GetDoubleInput(&numdata, data) < 0) {
		opserr<<"WARNING eleLoad - invalid value for beamPoint\n";
		return -1;
	    }
	    if (data[2] < 0.0 || data[2] > 1.0) {
		opserr << "WARNING eleLoad - invalid xDivL of " << data[2];
		opserr << " for beamPoint (valid range [0.0, 1.0]\n";
		return -1;
	    }

	    for (int i=0; i<theEleTags.Size(); i++) {
		theLoad = new Beam3dPointLoad(eleLoadTag, data[0], data[1], data[2], theEleTags(i), data[3]);

		if (theLoad == 0) {
		    opserr << "WARNING eleLoad - out of memory creating load of type " << type;
		    return -1;
		}

		// get the current pattern tag if no tag given in i/p
		int loadPatternTag = theActiveLoadPattern->getTag();

		// add the load to the domain
		if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
		    opserr << "WARNING eleLoad - could not add following load to domain:\n ";
		    opserr << theLoad;
		    delete theLoad;
		    return -1;
		}
		eleLoadTag++;
	    }
	    return 0;
	}
	else {
	    opserr << "WARNING eleLoad beamPoint type currently only valid only for ndm=2 or 3\n";
	    return -1;
	}
    }
    // Added Joey Yang UC Davis
    else if (strcmp(type,"-BrickW") == 0) {

	for (int i=0; i<theEleTags.Size(); i++) {
	    theLoad = new BrickSelfWeight(eleLoadTag, theEleTags(i));

	    if (theLoad == 0) {
		opserr << "WARNING eleLoad - out of memory creating load of type " << type;
		return -1;
	    }

	    // get the current pattern tag if no tag given in i/p
	    int loadPatternTag = theActiveLoadPattern->getTag();

	    // add the load to the domain
	    if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
		opserr << "WARNING eleLoad - could not add following load to domain:\n ";
		opserr << theLoad;
		delete theLoad;
		return -1;
	    }
	    eleLoadTag++;
	}
	return 0;
    }
    // Added: C.McGann, U.Washington
    else if (strcmp(type,"-surfaceLoad") == 0 || strcmp(type,"-SurfaceLoad") == 0) {

	for (int i=0; i<theEleTags.Size(); i++) {
	    theLoad = new SurfaceLoader(eleLoadTag, theEleTags(i));

	    if (theLoad == 0) {
		opserr << "WARNING eleLoad - out of memory creating load of type " << type;
		return -1;
	    }

	    // get the current pattern tag if no tag given in i/p
	    int loadPatternTag = theActiveLoadPattern->getTag();

	    // add the load to the domain
	    if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
		opserr << "WARNING eleLoad - could not add following load to domain:\n ";
		opserr << theLoad;
		delete theLoad;
		return -1;
	    }
	    eleLoadTag++;
	}
	return 0;
    }
// Added: C.McGann, U.Washington
    else if ((strcmp(type,"-selfWeight") == 0) || (strcmp(type,"-SelfWeight") == 0)) {
	// xf, yf, zf
	double data[3] = {0.0, 0.0, 0.0};
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 2) {
	    opserr<<"WARNING eleLoad - selfWeight want xf, yf, <zf>\n";
	    return -1;
	}
	if (numdata > 3) numdata = 3;
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
	    opserr<<"WARNING eleLoad - invalid value for SelfWeight\n";
	    return -1;
	}
	for (int i=0; i<theEleTags.Size(); i++) {
	    theLoad = new SelfWeight(eleLoadTag, data[0], data[1], data[2], theEleTags(i));

	    if (theLoad == 0) {
		opserr << "WARNING eleLoad - out of memory creating load of type " << type;
		return -1;
	    }

	    // get the current pattern tag if no tag given in i/p
	    int loadPatternTag = theActiveLoadPattern->getTag();

	    // add the load to the domain
	    if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
		opserr << "WARNING eleLoad - could not add following load to domain:\n ";
		opserr << theLoad;
		delete theLoad;
		return -1;
	    }
	    eleLoadTag++;
	}
	return 0;
    }

    //--Adding identifier for Beam2dThermalAction:[BEGIN] by UoE OpenSees Group--//
    else if (strcmp(type,"-beamThermal") == 0) {

	// get the current pattern tag if no tag given in i/p
	int loadPatternTag = theActiveLoadPattern->getTag();

	if (ndm == 2) {
	    //so far three kinds of temperature distribution
	    //(1) 9 temperature points, i.e. 8 layers
	    //(2) 5 temperature points, i.e. 4 layers
	    //(3) 2 temperature points, i.e. 1 layers: linear or uniform

	    //double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5,
	    // t6, locY6, t7, locY7, t8, locY8, t9, locY9;
	    // 9 temperature points are given,i.e. 8 layers are defined; Also the 9 corresponding vertical coordinate is given.
	    // the temperature at each fiber is obtained by interpolating of temperatures at the nearby temperature points.
	    int numdata = OPS_GetNumRemainingInputArgs();
	    double data[18];
        double Temp[9]; double Loc[9];
	    if (numdata == 18) {
		if (OPS_GetDoubleInput(&numdata, data) < 0) {
		    opserr << "WARNING eleLoad - invalid input\n";
		    return -1;
		}
        for (int i = 0; i < 9; i++) {
            Temp[i] = data[2 * i];
            Loc[i] = data[2 * i + 1];
        }
	    }

	    // 5 temperatures are given, i.e. 4 layers are defined.
	    else if (numdata == 10){
		if (OPS_GetDoubleInput(&numdata, data) < 0) {
		    opserr << "WARNING eleLoad - invalid input\n";
		    return -1;
		}

        Temp[0] = data[0]; Temp[2] = data[2]; Temp[4] = data[4]; Temp[6] = data[6]; Temp[8] = data[8];
        Loc[0]  = data[1]; Loc[2]  = data[3]; Loc[4]  = data[5]; Loc[6]  = data[7]; Loc[8]  = data[9];
        for (int i = 1; i < 5; i++) {
            Temp[2 * i - 1] = (Temp[2 * i - 2] + Temp[2 * i]) / 2;
            Loc[2 * i - 1] = (Loc[2 * i - 2] + Loc[2 * i]) / 2;
        }
	    }

	    // two temperature is given,
	    //if the two temperatures are equal,i.e. uniform Temperature change in element
	    //if the two temperatures are different,i.e. linear Temperature change in element
	    else if (numdata == 4){
		if (OPS_GetDoubleInput(&numdata, data) < 0) {
		    opserr << "WARNING eleLoad - invalid input\n";
		    return -1;
		}

        Temp[0] = data[0]; Temp[8] = data[2];
        Loc[0]  = data[1]; Loc[8]  = data[3];
        for (int i = 1; i < 8; i++) {
            Temp[i] = Temp[0] - i*(Temp[0] - Temp[8]) / 8;
            Loc[i] = Loc[0] - i*(Loc[0] - Loc[8]) / 8;
        }
	    }

        //finish the temperature arguments
        else {
            opserr << "WARNING eleLoad -beamThermalAction invalid number of temperature aguments,/n looking for 0, 2, 5 or 9 arguments.\n";
        }

        for (int i = 0; i<theEleTags.Size(); i++) {
            theLoad = new Beam2dThermalAction(eleLoadTag,
                Temp[0], Loc[0], Temp[1], Loc[1],
                Temp[2], Loc[2], Temp[3], Loc[3],
                Temp[4], Loc[4], Temp[5], Loc[5],
                Temp[6], Loc[6], Temp[7], Loc[7],
                Temp[8], Loc[8], theEleTags(i));

            if (theLoad == 0) {
                opserr << "WARNING eleLoad - out of memory creating load of type " << type;
                return -1;
            }

            // add the load to the domain
            if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
                opserr << "WARNING eleLoad - could not add following load to domain:\n ";
                opserr << theLoad;
                delete theLoad;
                return -1;
            }
            eleLoadTag++;
        }
	} // for the if (ndm==2)
	else {//if (ndm=3)
	    opserr << "WARNING eleLoad -beamThermalAction type currently only valid only for ndm=2\n";
	    return -1;
	}
    }
    //--Adding identifier for Beam2dThermalAction:[END] by UoE OpenSees Group--//


    // Added by Scott R. Hamilton   - Stanford
    else if (strcmp(type,"-beamTemp") == 0) {

	if (ndm == 2) {
	    int numdata = OPS_GetNumRemainingInputArgs();
	    double data[4];
	    // double temp1, temp2, temp3, temp4;

	    // Four temps given, Temp change at top node 1, bottom node 1, top node 2, bottom node 2.
	    if (numdata == 4){
		if (OPS_GetDoubleInput(&numdata, data) < 0) {
		    opserr << "WARNING eleLoad - invalid input\n";
		    return -1;
		}

		for (int i=0; i<theEleTags.Size(); i++) {
		    theLoad = new Beam2dTempLoad(eleLoadTag, data[0],
						 data[1], data[2],
						 data[3], theEleTags(i));

		    if (theLoad == 0) {
			opserr << "WARNING eleLoad - out of memory creating load of type " << type;
			return -1;
		    }

		    // get the current pattern tag if no tag given in i/p
		    int loadPatternTag = theActiveLoadPattern->getTag();

		    // add the load to the domain
		    if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
			opserr << "WARNING eleLoad - could not add following load to domain:\n ";
			opserr << theLoad;
			delete theLoad;
			return -1;
		    }
		    eleLoadTag++;
		}

		return 0;

	    }
	    // Two temps given, temp change at top, temp at bottom of element
	    else if (numdata == 2) {
		if (OPS_GetDoubleInput(&numdata, data) < 0) {
		    opserr << "WARNING eleLoad - invalid input\n";
		    return -1;
		}
		for (int i=0; i<theEleTags.Size(); i++) {
		    theLoad = new Beam2dTempLoad(eleLoadTag, data[0], data[1], theEleTags(i));

		    if (theLoad == 0) {
			opserr << "WARNING eleLoad - out of memory creating load of type " << type;
			return -1;
		    }

		    // get the current pattern tag if no tag given in i/p
		    int loadPatternTag = theActiveLoadPattern->getTag();

		    // add the load to the domain
		    if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
			opserr << "WARNING eleLoad - could not add following load to domain:\n ";
			opserr << theLoad;
			delete theLoad;
			return -1;
		    }
		    eleLoadTag++;
		}
	    }
	    // One twmp change give, uniform temp change in element
	    else if (numdata == 1) {
		if (OPS_GetDoubleInput(&numdata, data) < 0) {
		    opserr << "WARNING eleLoad - invalid input\n";
		    return -1;
		}


		theLoad=0;

		for (int i=0; i<theEleTags.Size(); i++) {
		    theLoad = new Beam2dTempLoad(eleLoadTag, data[0], theEleTags(i));

		    if (theLoad == 0) {
			opserr << "WARNING eleLoad - out of memory creating load of type " << type;
			return -1;
		    }

		    // get the current pattern tag if no tag given in i/p
		    int loadPatternTag = theActiveLoadPattern->getTag();

		    // add the load to the domain
		    if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
			opserr << "WARNING eleLoad - could not add following load to domain:\n ";
			opserr << theLoad;
			delete theLoad;
			return -1;
		    }
		    eleLoadTag++;
		}

		return 0;

	    }

	    else {
		opserr << "WARNING eleLoad -beamTempLoad invalid number of temperature aguments,/n looking for 0, 1, 2 or 4 arguments.\n";
		return -1;
	    }

	} else {
	    opserr << "WARNING eleLoad -beamTempLoad type currently only valid only for ndm=2\n";
	    return -1;
	}
    }

    // if get here we have sucessfully created the load and added it to the domain

    return 0;
}


int OPS_SP()
{
    Domain* theDomain = OPS_GetDomain();
    if(theDomain == 0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }

    if(OPS_GetNumRemainingInputArgs() < 3) {
	opserr<<"insufficient number of args\n";
	return -1;
    }

    SP_Constraint *theSP = 0;

    // get tags
    int tags[2];
    int numData = 2;
    if(OPS_GetIntInput(&numData, &tags[0]) < 0) {
	opserr << "WARNING invalid int tags\n";
	return -1;
    }

    // get node
    Node* theNode = theDomain->getNode(tags[0]);
    if(theNode == 0) {
	opserr<<"ERROR node "<<tags[0]<<"does not exsit\n";
	return -1;
    }
    int ndf = theNode->getNumberDOF();
    if(tags[1] > ndf || tags[1] < 0) {
	opserr<<"WARNING invalid dof\n";
	return -1;
    }

    // get value
    double value;
    numData = 1;
    if(OPS_GetDoubleInput(&numData,&value) < 0) {
	opserr << "WARNING invalid double value\n";
	return -1;
    }

    // get sp const
    bool isSpConst = false;
    bool userPattern = false;
    int loadPatternTag = 0;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type, "-const") == 0) {
	    isSpConst = true;
	} else if(strcmp(type, "-pattern") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		int numData = 1;
		if(OPS_GetIntInput(&numData, &loadPatternTag) < 0) {
		    opserr << "WARNING invalid pattern tag\n";
		    return -1;
		}
		userPattern = true;
	    }
	}
    }

    // get the current pattern tag
    LoadPattern* currPattern = theActiveLoadPattern;
    if(userPattern == false) {
	if(currPattern == 0) {
	    opserr<<"WARNING: no current pattern is set\n";
	    return -1;
	}
	loadPatternTag = currPattern->getTag();
    }

    // create pattern
    theSP = new SP_Constraint(tags[0], tags[1]-1, value, isSpConst);
    if(theSP == 0) return -1;

    // add load to domain
    if(theDomain->addSP_Constraint(theSP,loadPatternTag) == false) {
	opserr<<"WARNING: failed to add SP_Constraint to domain\n";
	delete theSP;
	return -1;
    }

    return 0;
}

int OPS_ImposedMotionSP()
{
    // check number of arguments
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING bad command - want: imposedMotion nodeId dofID gMotionID\n";
	return -1;
    }

    // get the nodeID, dofId and value of the constraint
    int nodeId, dofId, gMotionID;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &nodeId) < 0) {
	opserr << "WARNING invalid nodeId: ";
	opserr << " - imposedMotion nodeId dofID gMotionID\n";
	return -1;
    }

    if (OPS_GetIntInput(&numdata, &dofId) < 0) {
	opserr << "WARNING invalid dofId: imposedMotion ";
	opserr << nodeId << " dofID gMotionID\n";
	return -1;
    }
    dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

    if (OPS_GetIntInput(&numdata, &gMotionID) < 0) {
	opserr << "WARNING invalid gMotionID:  -  imposedMotion ";
	opserr << nodeId << " dofID gMotionID\n";
	return -1;
    }

    bool alt = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	if (strcmp(flag,"-other") == 0) {
	    alt = true;
	}
    }

    //
    // check valid node & dof
    //
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(nodeId);
    if (theNode == 0) {
	opserr << "WARNING invalid node " << nodeId << " node not found\n ";
	return -1;
    }
    int nDof = theNode->getNumberDOF();
    if (dofId < 0 || dofId >= nDof) {
	opserr << "WARNING invalid dofId: " << dofId << " dof specified cannot be <= 0 or greater than num dof at nod\n ";
	return -2;
    }


    MultiSupportPattern *thePattern = theActiveMultiSupportPattern;
    if (thePattern == 0) {
	opserr << "WARNING no active multi support pattern - imposedMotion\n";
	return -1;
    }
    int loadPatternTag = thePattern->getTag();

    // create a new ImposedMotionSP
    SP_Constraint *theSP;
    if (alt == true) {
	theSP = new ImposedMotionSP1(nodeId, dofId, loadPatternTag, gMotionID);
    }
    else {
	theSP = new ImposedMotionSP(nodeId, dofId, loadPatternTag, gMotionID);
    }

    if (theSP == 0) {
	opserr << "WARNING ran out of memory for ImposedMotionSP ";
	opserr << " -  imposedMotion ";
	opserr << nodeId << " " << dofId++ << " " << gMotionID << endln;
	return -1;
    }
    if (thePattern->addSP_Constraint(theSP) == false) {
	opserr << "WARNING could not add SP_Constraint to pattern ";
	delete theSP;
	return -1;
    }

    // if get here we have sucessfully created the node and added it to the domain
    return 0;
}

int OPS_groundMotion()
{
    GroundMotion *theMotion = 0;
    int gMotionTag;

    MultiSupportPattern *thePattern = theActiveMultiSupportPattern;
    if (thePattern == 0) {
	opserr << "WARNING no active multi support pattern - groundMotion\n";
	return -1;
    }

    // make sure at least one other argument to contain integrator
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING invalid command - want: groundMotion tag type <args>\n";
	opserr << "           valid types: AccelRecord and Interpolated \n";
	return -1;
    }

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &gMotionTag) < 0) {
	opserr << "WARNING invalid tag: groundMotion tag  type <args>\n";
	return -1;
    }

    const char* type = OPS_GetString();
    if ((strcmp(type,"Series") == 0) ||
	(strcmp(type,"Plain") == 0)) {

	TimeSeries *accelSeries = 0;
	TimeSeries *velSeries = 0;
	TimeSeries *dispSeries = 0;
	TimeSeriesIntegrator *seriesIntegrator = 0;

	double dtInt = 0.01;
	double fact = 1.0;
	while (OPS_GetNumRemainingInputArgs() > 1) {
	    const char* flag = OPS_GetString();
	    if ((strcmp(flag,"-accel") == 0) ||
		(strcmp(flag,"-acceleration") == 0)) {

		int tsTag;
		if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
		    opserr << "WARNING failed to get accel time series tag\n";
		    return -1;
		}

		accelSeries = OPS_getTimeSeries(tsTag);

		if (accelSeries == 0) {
		    opserr << "WARNING invalid accel series: " << tsTag;
		    opserr << " groundMotion tag Series -accel {series}\n";
		    return -1;
		}


	    } else if ((strcmp(flag,"-vel") == 0) ||
		       (strcmp(flag,"-velocity") == 0)) {

		int tsTag;
		if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
		    opserr << "WARNING failed to get accel time series tag\n";
		    return -1;
		}

		velSeries = OPS_getTimeSeries(tsTag);

		if (velSeries == 0) {
		    opserr << "WARNING invalid vel series: " << tsTag;
		    opserr << " groundMotion tag Series -vel {series}\n";
		    return -1;
		}


	    } else if ((strcmp(flag,"-disp") == 0) ||
		       (strcmp(flag,"-displacement") == 0)) {

		int tsTag;
		if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
		    opserr << "WARNING failed to get accel time series tag\n";
		    return -1;
		}

		dispSeries = OPS_getTimeSeries(tsTag);

		if (dispSeries == 0) {
		    opserr << "WARNING invalid disp series: " << tsTag;
		    opserr << " groundMotion tag Series -disp {series}\n";
		    return -1;
		}

	    } else if ((strcmp(flag,"-int") == 0) ||
		       (strcmp(flag,"-integrator") == 0)) {

		seriesIntegrator = (TimeSeriesIntegrator*) OPS_TimeSeriesIntegrator();

		if (seriesIntegrator == 0) {
		    opserr << "WARNING invalid series integrator: ";
		    opserr << " - groundMotion tag Series -int {Series Integrator}\n";
		    return -1;
		}

	    } else if ((strcmp(flag,"-dtInt") == 0) ||
		       (strcmp(flag,"-dtIntegrator") == 0) ||
		       (strcmp(flag,"-deltaT") == 0)) {

		if (OPS_GetDoubleInput(&numdata, &dtInt) < 0) {
		    opserr << "WARNING invalid dtInt: ";
		    opserr << " - groundMotion tag Series -dtInt dt\n";
		    return -1;
		}


	    } else if ((strcmp(flag,"-fact") == 0) ||
		       (strcmp(flag,"-factor") == 0)) {

		if (OPS_GetDoubleInput(&numdata, &fact) < 0) {
		    opserr << "WARNING invalid factor: ";
		    opserr << " - groundMotion tag Series -fact factor\n";
		    return -1;
		}
	    }

	}

	theMotion = new GroundMotion(dispSeries, velSeries,
				     accelSeries, seriesIntegrator, dtInt, fact);

	if (theMotion == 0) {
	    opserr << "WARNING ran out of memory creating ground motion - pattern UniformExcitation ";
	    opserr << gMotionTag << endln;

	    return -1;
	}
    }

    else if (strcmp(type,"Interpolated") == 0) {

	int numMotions = (OPS_GetNumRemainingInputArgs()-1) / 2;
	if (numMotions == 0) {
	    opserr << "WARNING no gMotionTags want :";
	    opserr << " pattern MultiSupport gMotion1? gMotion? .. ";
	    opserr << "-fact fact1? fact2? .. \n";
	    return -1;
	}

	std::vector<GroundMotion*> theMotions(numMotions);

	for (int i=0; i<numMotions; i++) {
	    int motionID;
	    if (OPS_GetIntInput(&numdata, &motionID) < 0) {
		opserr << "WARNING invalid motion id\n";
		return -1;
	    }

	    theMotions[i] = thePattern->getMotion(motionID);
	    if (theMotions[i] == 0) {
		opserr << "WARNING no groundMotion with tag " << motionID <<" :";
		opserr << " pattern MultiSupport gMotion1? gMotion? .. ";
		opserr << "-fact fact1? fact2? .. \n";
		return -1;
	    }
	}

	const char* flag = OPS_GetString();
	if (strcmp(flag, "-fact") != 0) {
	    opserr << "WARNING want -fact flag here\n";
	    return -1;
	}

	Vector facts(numMotions);
	for (int i=0; i<numMotions; i++) {
	    double fact;
	    if (OPS_GetDoubleInput(&numdata, &fact) < 0) {
		opserr << "WARNING invalid fact\n";
		return -1;
	    }
	    facts[i] = fact;
	}

	theMotion = new InterpolatedGroundMotion(&theMotions[0], facts, false);

    }

    else {
	opserr << "WARNING unknown pattern type " << type;
	opserr << " - want: pattern patternType " << gMotionTag ;
	opserr << " \t valid types: Plain, UniformExcitation \n";
	return -1;
    }


    // now add the load pattern to the modelBuilder
    if (theMotion != 0) {
	if (thePattern->addMotion(*theMotion, gMotionTag) < 0) {
	    opserr << "WARNING could not add ground motion with tag " << gMotionTag;
	    opserr << " to pattern\n ";
	    delete theMotion; // free up the memory, pattern destroys the time series
	    return -1;
	}
    }

    return 0;
}
