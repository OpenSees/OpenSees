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

#include <vector>

#include <RandomVariable.h>
#include <RandomVariableIter.h>
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
#include <PythonRV.h>

#include <AllIndependentTransformation.h>
#include <NatafProbabilityTransformation.h>


// active object
static OpenSeesReliabilityCommands* cmds = 0;

OpenSeesReliabilityCommands::OpenSeesReliabilityCommands(Domain* structuralDomain)
  :theDomain(0), theProbabilityTransformation(0)
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
    double mean = 0.0;
    double stdv = 1.0;
    double startPt = 0.0;
    int use_start_pt = 0;
	
    double param = 0;
    Vector parameters;

    // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
    if (OPS_GetNumRemainingInputArgs() < 2) {
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

    char *filename = 0;
    char *functionname = 0;
    
    // read options
    while (OPS_GetNumRemainingInputArgs() > 0) {

	const char* arg = OPS_GetString();

	if (strcmp(arg,"-file") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING not enough args, need -file filename??\n";
		opserr << "for random variable: " << tag << "\n";
		return -1;
	    }
	    filename = (char*)OPS_GetString();
	}
	if (strcmp(arg,"-function") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING not enough args, need -function functionname??\n";
		opserr << "for random variable: " << tag << "\n";
		return -1;
	    }
	    functionname = (char*)OPS_GetString();
	}	
	
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

    else if (strcmp(dist,"python") == 0) {
      if (filename == 0 || functionname == 0) {
	opserr << "ERROR: PythonRV filename or functionname not specified" << endln;
	return -1;
      }
      if (param_indx > 0)
	theRandomVariable = new PythonRV(tag, parameters, filename, functionname);
      else
	theRandomVariable = new PythonRV(tag, mean, stdv, filename, functionname);
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

int OPS_getRVTags()
{
  ReliabilityDomain* theReliabilityDomain = cmds->getDomain();
  if (theReliabilityDomain == 0)
    return -1;

  std::vector<int> rvTags;
  RandomVariable *theRV;
  RandomVariableIter &rvIter = theReliabilityDomain->getRandomVariables();
  while ((theRV = rvIter()) != 0)
    rvTags.push_back(theRV->getTag());

  int size = (int)rvTags.size();
  int *data = &rvTags[0];

  if (OPS_SetIntOutput(&size,data) < 0) {
    opserr << "ERROR: failed to set outputs in getRVTags" << endln;
    return -1;
  }

  return 0;
}

int OPS_getRVMean()
{
  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "ERROR: invalid number of arguments to getMean command : getMean rvTag\n";
    return -1;
  }

  // GET TAG NUMBER
  int rvTag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &rvTag) < 0) {
    opserr << "ERROR: invalid input to getMean: tag \n";
    return -1;
  }
    
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  RandomVariable *rv = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (rv == 0) {
    opserr << "ERROR: getMean - random variable with tag " << rvTag << " not found" << endln;
    return -1;
  }

  double mean = rv->getMean();
  if (OPS_SetDoubleOutput(&numData, &mean) < 0) {
    opserr << "ERROR: getMean - failed to set double output\n";
    return -1;
  }
  
  return 0;
}

int OPS_getRVStdv()
{
  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "ERROR: invalid number of arguments to getStdv command : getStdv rvTag\n";
    return -1;
  }

  // GET TAG NUMBER
  int rvTag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &rvTag) < 0) {
    opserr << "ERROR: invalid input to getStdv: tag \n";
    return -1;
  }
    
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  RandomVariable *rv = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (rv == 0) {
    opserr << "ERROR: getStdv - random variable with tag " << rvTag << " not found" << endln;
    return -1;
  }

  double stdv = rv->getStdv();
  if (OPS_SetDoubleOutput(&numData, &stdv) < 0) {
    opserr << "ERROR: getStdv - failed to set double output\n";
    return -1;
  }
  
  return 0;
}

int OPS_getRVPDF()
{
  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "ERROR: invalid number of arguments to getPDF command : getPDF rvTag X\n";
    return -1;
  }

  // GET TAG NUMBER
  int rvTag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &rvTag) < 0) {
    opserr << "ERROR: invalid input to getPDF: tag \n";
    return -1;
  }

  double x;
  if (OPS_GetDoubleInput(&numData, &x) < 0) {
    opserr << "ERROR: invalid input to getPDF: x \n";
    return -1;
  }
  
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  RandomVariable *rv = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (rv == 0) {
    opserr << "ERROR: getPDF - random variable with tag " << rvTag << " not found" << endln;
    return -1;
  }

  double pdf = rv->getPDFvalue(x);
  if (OPS_SetDoubleOutput(&numData, &pdf) < 0) {
    opserr << "ERROR: getPDF - failed to set double output\n";
    return -1;
  }
  
  return 0;
}

int OPS_getRVCDF()
{
  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "ERROR: invalid number of arguments to getCDF command : getCDF rvTag X\n";
    return -1;
  }

  // GET TAG NUMBER
  int rvTag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &rvTag) < 0) {
    opserr << "ERROR: invalid input to getCDF: tag \n";
    return -1;
  }

  double x;
  if (OPS_GetDoubleInput(&numData, &x) < 0) {
    opserr << "ERROR: invalid input to getCDF: x \n";
    return -1;
  }
  
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  RandomVariable *rv = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (rv == 0) {
    opserr << "ERROR: getCDF - random variable with tag " << rvTag << " not found" << endln;
    return -1;
  }

  double cdf = rv->getCDFvalue(x);
  if (OPS_SetDoubleOutput(&numData, &cdf) < 0) {
    opserr << "ERROR: getCDF - failed to set double output\n";
    return -1;
  }
  
  return 0;
}

int OPS_getRVInverseCDF()
{
  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "ERROR: invalid number of arguments to getInverseCDF command : getInverseCDF rvTag p\n";
    return -1;
  }

  // GET TAG NUMBER
  int rvTag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &rvTag) < 0) {
    opserr << "ERROR: invalid input to getInverseCDF: tag \n";
    return -1;
  }

  double p;
  if (OPS_GetDoubleInput(&numData, &p) < 0) {
    opserr << "ERROR: invalid input to getInverseCDF: p \n";
    return -1;
  }
  
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  RandomVariable *rv = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (rv == 0) {
    opserr << "ERROR: getInverseCDF - random variable with tag " << rvTag << " not found" << endln;
    return -1;
  }

  double invcdf = rv->getInverseCDFvalue(p);
  if (OPS_SetDoubleOutput(&numData, &invcdf) < 0) {
    opserr << "ERROR: getInverseCDF - failed to set double output\n";
    return -1;
  }
  
  return 0;
}

int OPS_addCorrelate()
{
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "ERROR: Wrong number of arguments to correlate command" << endln;
    return -1;
  }

  int rvTag[2];
  double correlationValue;

  int numData = 2;
  if (OPS_GetIntInput(&numData, rvTag) < 0) {
    opserr << "ERROR: invalid input to correlate: tag" << endln;;
    return -1;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &correlationValue) < 0) {
    opserr << "ERROR: invalid input to correlate: value" << endln;;
    return -1;    
  }

  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  int tag = theReliabilityDomain->getNumberOfCorrelationCoefficients();
  CorrelationCoefficient *theCorrelationCoefficient = new CorrelationCoefficient(tag+1, rvTag[0], rvTag[1], correlationValue);
  if (theCorrelationCoefficient == 0) {
    opserr << "ERROR: failed to add correlation coefficient to domain" << endln;
    return -1;
  }

  if (theReliabilityDomain->addCorrelationCoefficient(theCorrelationCoefficient) == false) {
    opserr << "ERROR: failed to add correlation coefficient to domain\n";
    opserr << "tag, rv1, rv2: " << tag << ' ' << rvTag[0] << ' ' << rvTag[1] << endln;
    return -1;
  }

  return 0;
}

void
OpenSeesReliabilityCommands::setProbabilityTransformation(ProbabilityTransformation *transform)
{
  if (theProbabilityTransformation != 0) {
    delete theProbabilityTransformation;
    theProbabilityTransformation = 0;
  }

  theProbabilityTransformation = transform;
  if (transform == 0)
    return;
}

int OPS_probabilityTransformation()
{
  //
  //
  // Check for replacement
  //
  //

  if (OPS_GetNumRemainingInputArgs() != 1 && OPS_GetNumRemainingInputArgs() != 3) {
    opserr << "ERROR: wrong number of arguments to probabilityTransformation" << endln;
    return -1;
  }

  // Get transformation type
  const char* arg = OPS_GetString();

  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  ProbabilityTransformation *theTransf = 0;
  if (strcmp(arg,"Nataf") == 0) {
    theTransf = new NatafProbabilityTransformation(theReliabilityDomain,0);
  }

  else if (strcmp(arg,"AllIndependent") == 0) {
    theTransf = new AllIndependentTransformation(theReliabilityDomain,0);
  }

  else {
    opserr << "ERROR: unrecognized type of probabilityTransformation" << endln;
    return -1;
  }

  if (theTransf == 0) {
    opserr << "ERROR: could not create probabilityTransformation" << endln;
    return -1;
  } else {
    if (cmds != 0)
      cmds->setProbabilityTransformation(theTransf);
  }

  
  return 0;
}

int OPS_transformUtoX()
{
  ProbabilityTransformation *theTransf = cmds->getProbabilityTransformation();
  if (theTransf == 0) {
    opserr << "ERROR: probability transformation has not been set" << endln;
    return -1;
  }

  ReliabilityDomain* theReliabilityDomain = cmds->getDomain();
  int nrv = theReliabilityDomain->getNumberOfRandomVariables();

  int numData = 1;
  double val;
  Vector u(nrv);
  int loc = 0;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    if (OPS_GetDoubleInput(&numData,&val) < 0) {
      OPS_ResetCurrentInputArg(-1);
      break;
    }
    u(loc) = val;
    loc++;
  }
  
  Vector x(nrv);
  theTransf->transform_u_to_x(u, x);

  if (OPS_SetDoubleOutput(&nrv, &x[0]) < 0) {
    opserr << "ERROR: failed to set output in transformUtoX" << endln;
    return -1;
  }
  
  return 0;
}
