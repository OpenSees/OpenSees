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

#include <TimeSeries.h>
#include <map>
#include <string.h>
#include <elementAPI.h>
#include <vector>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <Vector.h>


void* OPS_ConstantSeries();
void* OPS_LinearSeries();
void* OPS_TriangleSeries();
void* OPS_TrigSeries();
void* OPS_RampSeries();
void* OPS_RectangularSeries();
void* OPS_PulseSeries();
void* OPS_MPAccSeries();   //Tang.S
void* OPS_DiscretizedRandomProcessSeries();
void* OPS_SimulatedRandomProcessSeries();
void* OPS_PathSeries();

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
	functionMap.insert(std::make_pair("Constant", &OPS_ConstantSeries));
	functionMap.insert(std::make_pair("ConstantSeries", &OPS_ConstantSeries));
	functionMap.insert(std::make_pair("RampSeries", &OPS_RampSeries));
	functionMap.insert(std::make_pair("Ramp", &OPS_RampSeries));
	functionMap.insert(std::make_pair("Trig", &OPS_TrigSeries));
	functionMap.insert(std::make_pair("TrigSeries", &OPS_TrigSeries));
	functionMap.insert(std::make_pair("Sine", &OPS_TrigSeries));
	functionMap.insert(std::make_pair("SineSeries", &OPS_TrigSeries));
	functionMap.insert(std::make_pair("Linear", &OPS_LinearSeries));
	functionMap.insert(std::make_pair("LinearSeries", &OPS_LinearSeries));
	functionMap.insert(std::make_pair("Rectangular", &OPS_RectangularSeries));
	functionMap.insert(std::make_pair("Pulse", &OPS_PulseSeries));
	functionMap.insert(std::make_pair("PulseSeries", &OPS_PulseSeries));
	functionMap.insert(std::make_pair("Triangle", &OPS_TriangleSeries));
	functionMap.insert(std::make_pair("TriangleSeries", &OPS_TriangleSeries));
	functionMap.insert(std::make_pair("Path", &OPS_PathSeries));
	functionMap.insert(std::make_pair("Series", &OPS_PathSeries));
	functionMap.insert(std::make_pair("MPAcc", &OPS_MPAccSeries));  //Tang.S
	functionMap.insert(std::make_pair("MPAccSeries", &OPS_MPAccSeries));
	functionMap.insert(std::make_pair("DiscretizedRandomProcess", &OPS_DiscretizedRandomProcessSeries));
	functionMap.insert(std::make_pair("SimulatedRandomProcess", &OPS_SimulatedRandomProcessSeries));
      
	return 0;
    }
}

int
OPS_TimeSeries()
{
    static bool initDone = false;
    if (initDone == false) {
	setUpFunctions();
	initDone = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"WARNING too few arguments: timeSeries type? tag? ...\n";
	return -1;
    }

    const char* type = OPS_GetString();
    
    OPS_ParsingFunctionMap::const_iterator iter = functionMap.find(type);
    if (iter == functionMap.end()) {
	opserr<<"WARNING timeSeries type " << type << " is unknown\n";
	return -1;
    }

    TimeSeries* ts = (TimeSeries*) (*iter->second)();
    if (ts == 0) {
	return -1;
    }

    // Now add the timeseries to the domain
    if (OPS_addTimeSeries(ts) == false) {
	opserr<<"ERROR could not add timeseries to domain.\n";
	delete ts;
	return -1;
    }

    return 0;

}
