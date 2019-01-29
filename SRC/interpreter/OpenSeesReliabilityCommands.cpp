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

// Description: all reliability APIs are defined or declared here
//

#include "OpenSeesReliabilityCommands.h"
#include <elementAPI.h>

#include <RandomVariable.h>
#include <NormalRV.h>
#include <LognormalRV.h>
#include <GammaRV.h>
#include <ShiftedExponentialRV.h>
#include <ShiftedRayleighRV.h>
#include <ExponentialRV.h>
#include <RayleighRV.h>
#include <UniformRV.h>
#include <BetaRV.h>
#include <Type1LargestValueRV.h>
#include <Type1SmallestValueRV.h>
#include <Type2LargestValueRV.h>
#include <Type3SmallestValueRV.h>
#include <ChiSquareRV.h>
#include <GumbelRV.h>
#include <WeibullRV.h>
#include <UserDefinedRV.h>
#include <LaplaceRV.h>
#include <ParetoRV.h>



// active object
static OpenSeesReliabilityCommands* cmds = 0;

OpenSeesReliabilityCommands::OpenSeesReliabilityCommands(Domain* structuralDomain)
    :theDomain(0)
{
    if (structuralDomain != 0) {
	theDomain = new ReliabilityDomain(structuralDomain);	
    }

    cmds = this;
}

OpenSeesReliabilityCommands::~OpenSeesReliabilityCommands()
{
    if (theDomain != 0) delete theDomain;
    cmds = 0;
}

ReliabilityDomain*
OpenSeesReliabilityCommands::getDomain()
{
    return theDomain;
}

int OPS_randomVariable()
{
    RandomVariable *theRandomVariable = 0;
    int tag;
    double mean = 0;
    double stdv = 0;
    double startPt = 0;
    int use_start_pt = 0;
	
    double param = 0;
    Vector parameters;

    // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
    if (OPS_GetNumRemainingInputArgs() < 4) {
	opserr << "ERROR: invalid number of arguments to randomVariable command : randomVariable tag dist -mean mean -stdv stdv -startPoint startPoint -parameters pram1 pram2 ...\n";
	return -1;
    }

    // GET TAG NUMBER
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "ERROR: invalid input: tag \n";
	return -1;

    }

    // GET DISTRIBUTION
    const char* dist = OPS_GetString();

    // read options
    while (OPS_GetNumRemainingInputArgs() > 0) {

	const char* arg = OPS_GetString();
	
	// user specified mean directly
	if (strcmp(arg,"-mean") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING not enough args, need -mean mean??\n";
		opserr << "for random variable: " << tag << "\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numdata, &mean) < 0) {
		opserr << "WARNING invalid mean\n";
		opserr << " for random variable: " << tag << "\n";
		return -1;
	    }
	}
		
	// user specified standard deviation directly
	else if (strcmp(arg,"-stdv") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING not enough args, need -stdv stdv??\n";
		opserr << " for random variable: " << tag << "\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numdata, &stdv) < 0) {
		opserr << "WARNING invalid standard deviation\n";
		opserr << " for random variable: " << tag << "\n";
		return -1;
	    }
	}
		
	// user specified starting point directly
	else if (strcmp(arg,"-startPoint") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING not enough args, need -startPoint startPt??\n";
		opserr << " for random variable: " << tag << "\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numdata, &startPt) < 0) {
		opserr << "WARNING invalid starting point\n";
		opserr << " for random variable: " << tag << "\n";
		return -1;
	    }
            use_start_pt = 1;
	}
		
	// user input distribution specific parameters directly
	else if (strcmp(arg,"-parameters") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING not enough args, need -parameters param1 ...??\n";
		opserr << " for random variable: " << tag << "\n";
		return -1;
	    }

	    while (true) {
		if (OPS_GetDoubleInput(&numdata, &param) < 0) {
		    OPS_ResetCurrentInputArg(-1);
		    break;
		}
		parameters[parameters.Size()] = param;
	    }
	}
    }
    
    // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
    int param_indx = parameters.Size();
    if (strcmp(dist,"normal") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new NormalRV(tag, parameters);
	else
	    theRandomVariable = new NormalRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"lognormal") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new LognormalRV(tag, parameters);
	else
	    theRandomVariable = new LognormalRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"gamma") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new GammaRV(tag, parameters);
	else
	    theRandomVariable = new GammaRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"shiftedExponential") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new ShiftedExponentialRV(tag, parameters);
	else
	    theRandomVariable = new ShiftedExponentialRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"shiftedRayleigh") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new ShiftedRayleighRV(tag, parameters);
	else
	    theRandomVariable = new ShiftedRayleighRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"exponential") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new ExponentialRV(tag, parameters);
	else
	    theRandomVariable = new ExponentialRV(tag, mean, stdv);
    }

    else if (strcmp(dist,"rayleigh") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new RayleighRV(tag, parameters);
	else {
	    opserr << "Rayleigh random variable with tag " << tag << " cannot be created with only mean/stdv." << "\n";
	    return TCL_ERROR;
	}
    }
	
    else if (strcmp(dist,"uniform") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new UniformRV(tag, parameters);
	else
	    theRandomVariable = new UniformRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"beta") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new BetaRV(tag, parameters);
	else {
	    opserr << "Beta random variable with tag " << tag << " cannot be created with only mean/stdv." << "\n";
	    return TCL_ERROR;
	}
    }
	
    else if (strcmp(dist,"type1LargestValue") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new Type1LargestValueRV(tag, parameters);
	else
	    theRandomVariable = new Type1LargestValueRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"type1SmallestValue") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new Type1SmallestValueRV(tag, parameters);
	else
	    theRandomVariable = new Type1SmallestValueRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"type2LargestValue") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new Type2LargestValueRV(tag, parameters);
	else
	    theRandomVariable = new Type2LargestValueRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"type3SmallestValue") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new Type3SmallestValueRV(tag, parameters);
	else {
	    opserr << "T3S random variable with tag " << tag << " cannot be created with only mean/stdv." << "\n";
	    return TCL_ERROR;
	}
    }
		
    else if (strcmp(dist,"chiSquare") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new ChiSquareRV(tag, parameters);
	else
	    theRandomVariable = new ChiSquareRV(tag, mean, stdv);
    }
		
    else if (strcmp(dist,"gumbel") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new GumbelRV(tag, parameters);
	else
	    theRandomVariable = new GumbelRV(tag, mean, stdv);
    }

    else if (strcmp(dist,"weibull") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new WeibullRV(tag, parameters);
	else
	    theRandomVariable = new WeibullRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"laplace") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new LaplaceRV(tag, parameters);
	else
	    theRandomVariable = new LaplaceRV(tag, mean, stdv);
    }
	
    else if (strcmp(dist,"pareto") == 0) {
	if (param_indx > 0)
	    theRandomVariable = new ParetoRV(tag, parameters);
	else {
	    opserr << "Pareto random variable with tag " << tag << " cannot be created with only mean/stdv." << "\n";
	    return TCL_ERROR;
	}
    }

    else if (strcmp(dist,"userdefined") == 0) {
	// note userdefined is a special case and will not have any input read from the command line yet
	// unless user defined mean and standard deviation for some reason, which will break this input
	// because we assume argi starts at 3 here
        
        // KRM 4/22/2012 userdefined currently not implemented.
	// Vector xPoints;
	// Vector PDFpoints;
	// int numPoints = 0;
		
	// if (strcmp(argv[3],"-list") == 0) {
			
	//     numPoints = (argc-4) % 2;
	//     Vector temp_xPoints(numPoints);
	//     Vector temp_PDFpoints(numPoints);
			
	//     double x = 0.0;
	//     double pdf = 0.0;
	//     double x_old = 0.0;
			
	//     // Read the points
	//     for (int i=0; i < numPoints; i++) {
	// 	if (Tcl_GetDouble(interp, argv[4+2*i], &x) != TCL_OK) {
	// 	    opserr << "ERROR: Invalid x point to user-defined random variable." << "\n";
	// 	    return TCL_ERROR;
	// 	}
	// 	if (Tcl_GetDouble(interp, argv[5+2*i], &pdf) != TCL_OK) {
	// 	    opserr << "ERROR: Invalid PDF value point to user-defined random variable." << "\n";
	// 	    return TCL_ERROR;
	// 	}
	// 	if (i>0 && x<=x_old) {
	// 	    opserr << "ERROR: x-points to user-defined random variable must be consequtive!" << "\n";
	// 	    return TCL_ERROR;
	// 	}
	// 	temp_xPoints(i) = x;
	// 	temp_PDFpoints(i) = pdf;
	// 	x_old = x;
	//     }
			
	//     xPoints = temp_xPoints;
	//     PDFpoints = temp_PDFpoints;
			
	// }
	// else if (strcmp(argv[3],"-file") == 0) {
			
	//     // Open file where the vectors are given
	//     ifstream inputFile( argv[4], ios::in );
	//     if (inputFile.fail()) {
	// 	opserr << "File " << argv[4] << " could not be opened. " << "\n";
	// 	return TCL_ERROR;
	//     }
			
	//     // Loop through file to see how many entries there are
	//     double dummy;
	//     numPoints = 0;
	//     while (inputFile >> dummy) {
	// 	inputFile >> dummy;
	// 	numPoints++;
	//     }
	//     if (numPoints == 0) {
	// 	opserr << "ERROR: No entries in the direction file read by " << "\n"
	// 	       << "user-defined random variable, number " << tag << "\n";
	// 	return TCL_ERROR;
	//     }
			
	//     // rewind
	//     inputFile.clear(); 
	//     inputFile.seekg(0); 
			
	//     // Allocate vectors of correct size
	//     Vector temp_xPoints(numPoints);
	//     Vector temp_PDFpoints(numPoints);
						
	//     // Store the vector
	//     for (int i=0; i<numPoints; i++) {
	// 	inputFile >> temp_xPoints(i);
	// 	inputFile >> temp_PDFpoints(i);
	//     }
	//     inputFile.close();
			
	//     xPoints = temp_xPoints;
	//     PDFpoints = temp_PDFpoints;
	// }
	// else {
	//     opserr << "ERROR: Invalid argument to user-defined random variable, number " << tag << "\n";
	//     return TCL_ERROR;
	// }
		
	//theRandomVariable = new UserDefinedRV(tag, xPoints, PDFpoints);
		
    }
	
    else {
	opserr << "ERROR: unknown random variable type: " << dist << " provided. Must be one of " << "\n";
	return -1;
    }

    if (theRandomVariable == 0) {
	opserr << "ERROR: could not create random variable number " << tag << "\n";
	return -1;
    }
	
    // set start point on object if user provided
    if (use_start_pt == 1) {
	theRandomVariable->setStartValue(startPt);
        theRandomVariable->setCurrentValue(startPt);
    }
    else {
	theRandomVariable->setStartValue(theRandomVariable->getMean());
        theRandomVariable->setCurrentValue(theRandomVariable->getMean());
    }

    // Add the random variable to the domain
    ReliabilityDomain* theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain->addRandomVariable(theRandomVariable) == false) {
	opserr << "ERROR: failed to add random variable to the domain (wrong number of arguments?)\n";
	opserr << "random variable: " << tag << "\n";
	delete theRandomVariable; // otherwise memory leak
	return -1;
    }

    //RVParameter *theRVParam = new RVParameter(tag, theRandomVariable);
    //theStructuralDomain->addParameter(theRVParam);

    return 0;

}
