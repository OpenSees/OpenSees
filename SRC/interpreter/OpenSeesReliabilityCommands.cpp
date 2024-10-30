/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California
(Regents). All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are
those of the authors and should not be interpreted as representing official
policies, either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

***************************************************************************
*/

// Written: Minjie

// Description: all reliability APIs are defined or declared here
//

#include "OpenSeesReliabilityCommands.h"

#include <BetaRV.h>
#include <ChiSquareRV.h>
#include <ExponentialRV.h>
#include <GammaRV.h>
#include <GumbelRV.h>
#include <ID.h>
#include <LaplaceRV.h>
#include <LognormalRV.h>
#include <NormalRV.h>
#include <ParetoRV.h>
#include <RandomVariable.h>
#include <RandomVariableIter.h>
#include <RayleighRV.h>
#include <ShiftedExponentialRV.h>
#include <ShiftedRayleighRV.h>
#include <Type1LargestValueRV.h>
#include <Type1SmallestValueRV.h>
#include <Type2LargestValueRV.h>
#include <Type3SmallestValueRV.h>
#include <UniformRV.h>
#include <UserDefinedRV.h>
#include <WeibullRV.h>
#include <elementAPI.h>

#include <fstream>
#include <vector>
#ifdef _PYTHON3
#include <PythonEvaluator.h>
#include <PythonRV.h>
#endif

#include <StandardLinearOscillatorDisplacementFilter.h>
#include <StandardLinearOscillatorVelocityFilter.h>
#include <StandardLinearOscillatorAccelerationFilter.h>
#include <KooFilter.h>
#include <GammaModulatingFunction.h>
#include <ConstantModulatingFunction.h>
#include <TrapezoidalModulatingFunction.h>
#include <KooModulatingFunction.h>
#include <JonswapSpectrum.h>
#include <NarrowBandSpectrum.h>
#include <PointsSpectrum.h>

#include <AdkZhangMeritFunctionCheck.h>
#include <AllIndependentTransformation.h>
#include <ArmijoStepSizeRule.h>
#include <CStdLibRandGenerator.h>
#include <FiniteDifferenceGradient.h>
#include <FixedStepSizeRule.h>
#include <GradientProjectionSearchDirection.h>
#include <HLRFSearchDirection.h>
#include <ImplicitGradient.h>
#include <LimitStateFunctionIter.h>
#include <NatafProbabilityTransformation.h>
#include <OptimalityConditionReliabilityConvergenceCheck.h>
#include <PolakHeSearchDirectionAndMeritFunction.h>
#include <SQPsearchDirectionMeritFunctionAndHessian.h>
#include <SearchWithStepSizeAndStepDirection.h>
#include <SecantRootFinding.h>
#include <StandardReliabilityConvergenceCheck.h>
#include <FirstPrincipalCurvature.h>
#include <CurvaturesBySearchAlgorithm.h>

// active object
static OpenSeesReliabilityCommands *cmds = 0;

// other static objects
PolakHeSearchDirectionAndMeritFunction *thePolakHeDualPurpose;
SQPsearchDirectionMeritFunctionAndHessian *theSQPtriplePurpose;

OpenSeesReliabilityCommands::OpenSeesReliabilityCommands(
    Domain *structuralDomain)
    : theDomain(0),
      theStructuralDomain(structuralDomain),
      theProbabilityTransformation(0),
      theRandomNumberGenerator(0),
      theReliabilityConvergenceCheck(0),
      theSearchDirection(0),
      theMeritFunctionCheck(0),
      theStepSizeRule(0),
      theRootFinding(0),
      theFindDesignPointAlgorithm(0),
      theFindCurvatures(0),
      theFunctionEvaluator(0),
      theGradientEvaluator(0),
      thePolakHeDualPurpose(0),
      theSQPtriplePurpose(0),
      theFOSMAnalysis(0),
      theFORMAnalysis(0),
      theSORMAnalysis(0),
      theImportanceSamplingAnalysis(0),
      theMonteCarloAnalysis(0),      
      theSensAlgo(0) {
    if (structuralDomain != 0) {
        theDomain = new ReliabilityDomain(structuralDomain);
    }

    cmds = this;
}

OpenSeesReliabilityCommands::~OpenSeesReliabilityCommands() {
  this->wipe();
  if (theDomain != 0)
    delete theDomain;
  cmds = 0;
}

ReliabilityDomain *OpenSeesReliabilityCommands::getDomain() {
    return theDomain;
}

Domain *OpenSeesReliabilityCommands::getStructuralDomain() {
    return theStructuralDomain;
}

void OpenSeesReliabilityCommands::wipe() {
    // wipe reliability domain
    if (theDomain != 0) {
        theDomain->clearAll();
    }

    if (theProbabilityTransformation != 0) {
        delete theProbabilityTransformation;
        theProbabilityTransformation = 0;
    }
    if (theRandomNumberGenerator != 0) {
        delete theRandomNumberGenerator;
        theRandomNumberGenerator = 0;
    }
    if (theReliabilityConvergenceCheck != 0) {
        delete theReliabilityConvergenceCheck;
        theReliabilityConvergenceCheck = 0;
    }
    if (theSearchDirection != 0) {
        delete theSearchDirection;
        theSearchDirection = 0;
    }
    if (theMeritFunctionCheck != 0) {
        delete theMeritFunctionCheck;
        theMeritFunctionCheck = 0;
    }
    if (theStepSizeRule != 0) {
        delete theStepSizeRule;
        theStepSizeRule = 0;
    }
    if (theRootFinding != 0) {
        delete theRootFinding;
        theRootFinding = 0;
    }
    if (theFindDesignPointAlgorithm != 0) {
        delete theFindDesignPointAlgorithm;
        theFindDesignPointAlgorithm = 0;
    }
    if (theFindCurvatures != 0) {
        delete theFindCurvatures;
        theFindCurvatures = 0;
    }
    if (theFunctionEvaluator != 0) {
        delete theFunctionEvaluator;
        theFunctionEvaluator = 0;
    }
    if (theGradientEvaluator != 0) {
        delete theGradientEvaluator;
        theGradientEvaluator = 0;
    }
    if (theFOSMAnalysis != 0) {
        delete theFOSMAnalysis;
        theFOSMAnalysis = 0;
    }
    if (theFORMAnalysis != 0) {
        delete theFORMAnalysis;
        theFORMAnalysis = 0;
    }
    if (theSORMAnalysis != 0) {
        delete theSORMAnalysis;
        theSORMAnalysis = 0;
    }    
    if (theImportanceSamplingAnalysis != 0) {
        delete theImportanceSamplingAnalysis;
        theImportanceSamplingAnalysis = 0;
    }
    if (theMonteCarloAnalysis != 0) {
        delete theMonteCarloAnalysis;
        theMonteCarloAnalysis = 0;
    }    
}

int OPS_wipeReliability() {
    if (cmds != 0) {
        cmds->wipe();
    }
    return 0;
}

int OPS_performanceFunction() {
    LimitStateFunction *theLSF = 0;
    int tag;

    // Check enough arguments
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr
            << "ERROR: invalid number of arguments to performanceFunction "
               "command: performanceFunction tag ...\n";
        return -1;
    }

    // Get tag
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "ERROR: invalid tag for performanceFunction: tag \n";
        return -1;
    }

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();

    // Evaluate performance function
    FunctionEvaluator *theEvaluator = cmds->getFunctionEvaluator();
    if (OPS_GetNumRemainingInputArgs() < 1) {
        if (theEvaluator == 0) {
            opserr << "Function evaluator must be defined in order to "
                      "evaluate "
                      "performance function"
                   << endln;
            return -1;
        }
        theReliabilityDomain->setTagOfActiveLimitStateFunction(tag);
        double g = theEvaluator->evaluateExpression();
        if (OPS_SetDoubleOutput(&numdata, &g, true) < 0) {
            opserr << "ERROR: performanceFunction - failed to set double "
                      "output\n";
            return -1;
        }
        return 0;
    }

    if (theEvaluator != 0) {
        opserr
            << "ERROR: A limit-state function should not be created after "
               "the GFunEvaluator has been instantiated.\n";
        return -1;
    }

    // Create new performance function
    const char *lsf = OPS_GetString();
    theLSF = new LimitStateFunction(tag, lsf);

    // Add the performance function to the domain
    if (theReliabilityDomain->addLimitStateFunction(theLSF) == false) {
        opserr << "ERROR: failed to add performance function to the "
                  "reliability domain\n";
        opserr << "performanceFunction: " << tag << "\n";
        delete theLSF;  // otherwise memory leak
        return -1;
    }

    return 0;
}

int OPS_gradPerformanceFunction() {
    LimitStateFunction *theLimitStateFunction = 0;

    // Check enough arguments
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "ERROR: invalid number of arguments -- "
                  "command: gradPerformanceFunction lsfTag rvTag expr\n";
        return -1;
    }

    // Get tags
    int numdata = 2;
    int tags[2];
    if (OPS_GetIntInput(&numdata, &tags[0]) < 0) {
        opserr << "ERROR: invalid tag for gradPerformanceFunction: lsfTag "
                  "rvTag \n";
        return -1;
    }

    // get domain
    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "ERROR: reliability domain is invalid\n";
        return -1;
    }

    // Evaluate performance function
    // FunctionEvaluator *theFunctionEvaluator =
    // cmds->getFunctionEvaluator(); if (theFunctionEvaluator != 0) {
    // opserr << "ERROR: A limit-state function should not be created
    // after the FunctionEvaluator has been instantiated." << endln;
    // return TCL_ERROR;
    // }

    // GET LSF pointer
    theLimitStateFunction =
        theReliabilityDomain->getLimitStateFunctionPtr(tags[0]);
    if (theLimitStateFunction == 0) {
        opserr << "ERROR: limit state function with tag " << tags[0]
               << " does not exist\n";
        return -1;
    }

    // ADD THE OBJECT TO THE LSF
    const char *expr = OPS_GetString();
    int ok = theLimitStateFunction->addGradientExpression(expr, tags[1]);
    if (ok < 0) {
        opserr << "ERROR: could not add gradient of LSF " << tags[0]
               << " for random variable " << tags[1] << endln;
        return -1;
    }

    return 0;
}

int OPS_filter() {
  Filter *theFilter = 0;
  int tag;
  double period_Tn, damping;//, dtpulse;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "ERROR: Invalid number of arguments to filter "
      "command: filter tag type arg1 arg2 ..." << endln;
    return -1;
  }

  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "ERROR: invalid tag for filter command" << endln;
    return -1;
  }

  const char *type = OPS_GetString();
  if (strcmp(type,"standard") == 0 || strcmp(type,"Koo") == 0 ||
      strcmp(type,"standardDisplacement") == 0 ||
      strcmp(type,"standardVelocity") == 0 ||
      strcmp(type,"standardAcceleration") == 0) {

    if (OPS_GetDoubleInput(&numdata, &period_Tn) < 0) {
      opserr << "ERROR: invalid period for filter" << endln;
      return -1;
    }
    if (OPS_GetDoubleInput(&numdata, &damping) < 0) {
      opserr << "ERROR: invalid damping for filter" << endln;
      return -1;
    }

    if (strcmp(type,"standard") == 0 || strcmp(type,"standardDisplacement") == 0)
      theFilter = new StandardLinearOscillatorDisplacementFilter(tag, period_Tn, damping);
    if (strcmp(type,"standardVelocity") == 0)
      theFilter = new StandardLinearOscillatorVelocityFilter(tag, period_Tn, damping);
    if (strcmp(type,"standardAcceleration") == 0)
      theFilter = new StandardLinearOscillatorAccelerationFilter(tag, period_Tn, damping);
    if (strcmp(type,"Koo") == 0)
      theFilter = new KooFilter(tag, period_Tn, damping);    
  }
  else {
    opserr << "Unknown filter type: " << type << endln;
    return -1;
  }

  if (theFilter == 0) {
    opserr << "ERROR: ran out of memory creating filter \n";
    opserr << type << ' ' << tag << endln;
    return -1;
  }
  
  // ADD THE OBJECT TO THE DOMAIN
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  if (theReliabilityDomain->addFilter(theFilter) == false) {
    opserr << "ERROR: failed to add filter to the domain\n";
    opserr << type << ' ' << tag << endln;
    delete theFilter; // otherwise memory leak
    return -1;
  }
  
  return 0;
}

int OPS_spectrum() {
  Spectrum *theSpectrum = 0;
  int tag;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "ERROR: Invalid number of arguments to spectrum "
      "command: spectrum tag type arg1 arg2 ..." << endln;
    return -1;
  }

  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "ERROR: invalid tag for spectrum command" << endln;
    return -1;
  }

  const char *type = OPS_GetString();
  if (strcmp(type,"jonswap") == 0) {
    if (OPS_GetNumRemainingInputArgs() < 5) {
      opserr << "ERROR: insufficient arguments for jonswap spectrum" << endln;
      return -1;
    }
    double data[5];
    numdata = 5;
    if (OPS_GetDoubleInput(&numdata,data) < 0) {
      opserr << "ERROR: invalid double data for spectrum with tag " << tag << endln;
      return -1;
    }

    theSpectrum = new JonswapSpectrum(tag, data[0], data[1], data[2], data[3], data[4]);    
  }
  else if (strcmp(type,"narrowband") == 0) {
    if (OPS_GetNumRemainingInputArgs() < 3) {
      opserr << "ERROR: insufficient arguments for narrowband spectrum" << endln;
      return -1;
    }
    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata,data) < 0) {
      opserr << "ERROR: invalid double data for spectrum with tag " << tag << endln;
      return -1;
    }

    theSpectrum = new NarrowBandSpectrum(tag, data[0], data[1], data[2]);
  }
  else if (strcmp(type,"points") == 0) {
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 4) {
      opserr << "ERROR: insufficient arguments for points spectrum, need at least two points" << endln;
      return -1;
    }    

    if (numArgs % 2 == 1) {
      opserr << "Warning - points spectrum, odd number of values entered, ignoring final input" << endln;
    }
    int numPoints = numArgs / 2;
    Vector frequencies(numPoints);
    Vector amplitudes(numPoints);
    numdata = 1;
    for (int i = 0; i < numPoints; i++) {
      double data;
      if (OPS_GetDoubleInput(&numdata,&data) < 0) {
	opserr << "ERROR: invalid double data for spectrum with tag " << tag << endln;
	return -1;
      }
      frequencies(i) = data;
      if (OPS_GetDoubleInput(&numdata,&data) < 0) {
	opserr << "ERROR: invalid double data for spectrum with tag " << tag << endln;
	return -1;
      }
      amplitudes(i) = data;      
    }

    theSpectrum = new PointsSpectrum(tag, frequencies, amplitudes);    
  } 
  else {
    opserr << "Unknown spectrum type: " << type << endln;
    return -1;
  }

  if (theSpectrum == 0) {
    opserr << "ERROR: ran out of memory creating spectrum \n";
    opserr << type << ' ' << tag << endln;
    return -1;
  }
  
  // ADD THE OBJECT TO THE DOMAIN
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  if (theReliabilityDomain->addSpectrum(theSpectrum) == false) {
    opserr << "ERROR: failed to add spectrum to the domain\n";
    opserr << type << ' ' << tag << endln;
    delete theSpectrum; // otherwise memory leak
    return -1;
  }
  
  return 0;
}

int OPS_modulatingFunction() {
  ModulatingFunction *theModulatingFunction = 0;
  int tag;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "ERROR: Invalid number of arguments to modulating function "
      "command: modulatingFunction tag type arg1 arg2 ..." << endln;
    return -1;
  }

  int numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "ERROR: invalid tag for modulatingFunction command" << endln;
    return -1;
  }

  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  
  const char *type = OPS_GetString();
  if (strcmp(type,"gamma") == 0) {

    if (OPS_GetNumRemainingInputArgs() < 4) {
      opserr << "ERROR: insufficient input for gamma modulatingFunction" << endln;
      return -1;
    }
    
    int filterTag;
    double abc[3];
    if (OPS_GetIntInput(&numdata, &filterTag) < 0) {
      opserr << "ERROR: invalid filter tag for modulating function" << endln;
      return -1;
    }

    Filter *theFilter = theReliabilityDomain->getFilterPtr(filterTag);
    if (theFilter == 0) {
      opserr << "ERROR: could not find filter with tag " << filterTag
	     << " for modulating function" << endln;
      return -1;
    }
    
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, abc) < 0) {
      opserr << "ERROR: invalid double data for modulating function" << endln;
      return -1;
    }    

    theModulatingFunction = new GammaModulatingFunction(tag, theFilter,
							abc[0], abc[1], abc[2]);
  }
  else if (strcmp(type,"constant") == 0) {

    if (OPS_GetNumRemainingInputArgs() < 2) {
      opserr << "ERROR: insufficient input for constant modulatingFunction" << endln;
      return -1;
    }
    
    int filterTag;
    double amplitude;
    if (OPS_GetIntInput(&numdata, &filterTag) < 0) {
      opserr << "ERROR: invalid filter tag for modulating function" << endln;
      return -1;
    }

    Filter *theFilter = theReliabilityDomain->getFilterPtr(filterTag);
    if (theFilter == 0) {
      opserr << "ERROR: could not find filter with tag " << filterTag
	     << " for modulating function" << endln;
      return -1;
    }
    
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &amplitude) < 0) {
      opserr << "ERROR: invalid amplitude for modulating function" << endln;
      return -1;
    }    

    theModulatingFunction = new ConstantModulatingFunction(tag, theFilter, amplitude);
  }
  else if (strcmp(type,"trapezoidal") == 0) {

    if (OPS_GetNumRemainingInputArgs() < 6) {
      opserr << "ERROR: insufficient input for trapezoidal modulatingFunction" << endln;
      return -1;
    }
    
    int filterTag;
    double ddata[5];
    if (OPS_GetIntInput(&numdata, &filterTag) < 0) {
      opserr << "ERROR: invalid filter tag for modulating function" << endln;
      return -1;
    }

    Filter *theFilter = theReliabilityDomain->getFilterPtr(filterTag);
    if (theFilter == 0) {
      opserr << "ERROR: could not find filter with tag " << filterTag
	     << " for modulating function" << endln;
      return -1;
    }
    
    numdata = 5;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
      opserr << "ERROR: invalid double data for modulating function" << endln;
      return -1;
    }    

    theModulatingFunction = new TrapezoidalModulatingFunction(tag, theFilter,
							      ddata[0], ddata[1], ddata[2],
							      ddata[3], ddata[4]);
  }  
  else if (strcmp(type,"Koo") == 0) {

    if (OPS_GetNumRemainingInputArgs() < 3) {
      opserr << "ERROR: insufficient input for Koo modulatingFunction" << endln;
      return -1;
    }
    
    int filterTag;
    double ddata[2];
    if (OPS_GetIntInput(&numdata, &filterTag) < 0) {
      opserr << "ERROR: invalid filter tag for modulating function" << endln;
      return -1;
    }

    Filter *theFilter = theReliabilityDomain->getFilterPtr(filterTag);
    if (theFilter == 0) {
      opserr << "ERROR: could not find filter with tag " << filterTag
	     << " for modulating function" << endln;
      return -1;
    }
    
    numdata = 2;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
      opserr << "ERROR: invalid double data for modulating function" << endln;
      return -1;
    }    

    theModulatingFunction = new KooModulatingFunction(tag, theFilter, ddata[0], ddata[1]);
  }
  else {
    opserr << "Unknown modulatingFunction type: " << type << endln;
    return -1;
  }

  if (theModulatingFunction == 0) {
    opserr << "ERROR: ran out of memory creating modulating function \n";
    opserr << type << ' ' << tag << endln;
    return -1;
  }
  
  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addModulatingFunction(theModulatingFunction) == false) {
    opserr << "ERROR: failed to add modulating function to the domain\n";
    opserr << type << ' ' << tag << endln;
    delete theModulatingFunction; // otherwise memory leak
    return -1;
  }
  
  return 0;
}

int OPS_randomVariable() {
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
        opserr
            << "ERROR: invalid number of arguments to randomVariable "
               "command : randomVariable tag dist -mean mean -stdv stdv "
               "-startPoint startPoint -parameters pram1 pram2 ...\n";
        return -1;
    }

    // GET TAG NUMBER
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "ERROR: invalid input: tag \n";
        return -1;
    }

    // GET DISTRIBUTION
    const char *dist = OPS_GetString();

    char *filename = 0;
    char *functionname = 0;

    // read options
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();

        if (strcmp(arg, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr
                    << "WARNING not enough args, need -file filename??\n";
                opserr << "for random variable: " << tag << "\n";
                return -1;
            }
            filename = (char *)OPS_GetString();
        }
        if (strcmp(arg, "-function") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING not enough args, need -function "
                          "functionname??\n";
                opserr << "for random variable: " << tag << "\n";
                return -1;
            }
            functionname = (char *)OPS_GetString();
        }

        // user specified mean directly
        if (strcmp(arg, "-mean") == 0) {
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
        else if (strcmp(arg, "-stdv") == 0) {
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
        else if (strcmp(arg, "-startPoint") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING not enough args, need -startPoint "
                          "startPt??\n";
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
        else if (strcmp(arg, "-parameters") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING not enough args, need -parameters "
                          "param1 ...??\n";
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
    if (strcmp(dist, "normal") == 0) {
        if (param_indx > 0)
            theRandomVariable = new NormalRV(tag, parameters);
        else
            theRandomVariable = new NormalRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "lognormal") == 0) {
        if (param_indx > 0)
            theRandomVariable = new LognormalRV(tag, parameters);
        else
            theRandomVariable = new LognormalRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "gamma") == 0) {
        if (param_indx > 0)
            theRandomVariable = new GammaRV(tag, parameters);
        else
            theRandomVariable = new GammaRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "shiftedExponential") == 0) {
        if (param_indx > 0)
            theRandomVariable = new ShiftedExponentialRV(tag, parameters);
        else
            theRandomVariable = new ShiftedExponentialRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "shiftedRayleigh") == 0) {
        if (param_indx > 0)
            theRandomVariable = new ShiftedRayleighRV(tag, parameters);
        else
            theRandomVariable = new ShiftedRayleighRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "exponential") == 0) {
        if (param_indx > 0)
            theRandomVariable = new ExponentialRV(tag, parameters);
        else
            theRandomVariable = new ExponentialRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "rayleigh") == 0) {
        if (param_indx > 0)
            theRandomVariable = new RayleighRV(tag, parameters);
        else {
            opserr << "Rayleigh random variable with tag " << tag
                   << " cannot be created with only mean/stdv."
                   << "\n";
            return TCL_ERROR;
        }
    }

    else if (strcmp(dist, "uniform") == 0) {
        if (param_indx > 0)
            theRandomVariable = new UniformRV(tag, parameters);
        else
            theRandomVariable = new UniformRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "beta") == 0) {
        if (param_indx > 0)
            theRandomVariable = new BetaRV(tag, parameters);
        else {
            opserr << "Beta random variable with tag " << tag
                   << " cannot be created with only mean/stdv."
                   << "\n";
            return TCL_ERROR;
        }
    }

    else if (strcmp(dist, "type1LargestValue") == 0) {
        if (param_indx > 0)
            theRandomVariable = new Type1LargestValueRV(tag, parameters);
        else
            theRandomVariable = new Type1LargestValueRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "type1SmallestValue") == 0) {
        if (param_indx > 0)
            theRandomVariable = new Type1SmallestValueRV(tag, parameters);
        else
            theRandomVariable = new Type1SmallestValueRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "type2LargestValue") == 0) {
        if (param_indx > 0)
            theRandomVariable = new Type2LargestValueRV(tag, parameters);
        else
            theRandomVariable = new Type2LargestValueRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "type3SmallestValue") == 0) {
        if (param_indx > 0)
            theRandomVariable = new Type3SmallestValueRV(tag, parameters);
        else {
            opserr << "T3S random variable with tag " << tag
                   << " cannot be created with only mean/stdv."
                   << "\n";
            return TCL_ERROR;
        }
    }

    else if (strcmp(dist, "chiSquare") == 0) {
        if (param_indx > 0)
            theRandomVariable = new ChiSquareRV(tag, parameters);
        else
            theRandomVariable = new ChiSquareRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "gumbel") == 0) {
        if (param_indx > 0)
            theRandomVariable = new GumbelRV(tag, parameters);
        else
            theRandomVariable = new GumbelRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "weibull") == 0) {
        if (param_indx > 0)
            theRandomVariable = new WeibullRV(tag, parameters);
        else
            theRandomVariable = new WeibullRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "laplace") == 0) {
        if (param_indx > 0)
            theRandomVariable = new LaplaceRV(tag, parameters);
        else
            theRandomVariable = new LaplaceRV(tag, mean, stdv);
    }

    else if (strcmp(dist, "pareto") == 0) {
        if (param_indx > 0)
            theRandomVariable = new ParetoRV(tag, parameters);
        else {
            opserr << "Pareto random variable with tag " << tag
                   << " cannot be created with only mean/stdv."
                   << "\n";
            return TCL_ERROR;
        }
    }

    else if (strcmp(dist, "userdefined") == 0) {
        // note userdefined is a special case and will not have any input
        // read from the command line yet unless user defined mean and
        // standard deviation for some reason, which will break this input
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
        // 	    opserr << "ERROR: Invalid x point to user-defined
        // random variable." << "\n"; 	    return TCL_ERROR;
        // 	}
        // 	if (Tcl_GetDouble(interp, argv[5+2*i], &pdf) != TCL_OK) {
        // 	    opserr << "ERROR: Invalid PDF value point to
        // user-defined random variable." << "\n"; 	    return
        // TCL_ERROR;
        // 	}
        // 	if (i>0 && x<=x_old) {
        // 	    opserr << "ERROR: x-points to user-defined random
        // variable must be consecutive!" << "\n"; 	    return
        // TCL_ERROR;
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
        // 	opserr << "File " << argv[4] << " could not be opened. " <<
        // "\n"; 	return TCL_ERROR;
        //     }

        //     // Loop through file to see how many entries there are
        //     double dummy;
        //     numPoints = 0;
        //     while (inputFile >> dummy) {
        // 	inputFile >> dummy;
        // 	numPoints++;
        //     }
        //     if (numPoints == 0) {
        // 	opserr << "ERROR: No entries in the direction file read by
        // " << "\n"
        // 	       << "user-defined random variable, number " << tag <<
        // "\n"; 	return TCL_ERROR;
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
        //     opserr << "ERROR: Invalid argument to user-defined random
        //     variable, number " << tag << "\n"; return TCL_ERROR;
        // }

        // theRandomVariable = new UserDefinedRV(tag, xPoints, PDFpoints);

    }

    else if (strcmp(dist, "python") == 0) {
#ifdef PY_VERSION
        if (filename == 0 || functionname == 0) {
            opserr
                << "ERROR: PythonRV filename or functionname not specified"
                << endln;
            return -1;
        }
        if (param_indx > 0)
            theRandomVariable =
                new PythonRV(tag, parameters, filename, functionname);
        else
            theRandomVariable =
                new PythonRV(tag, mean, stdv, filename, functionname);
#endif
    }

    else {
        opserr << "ERROR: unknown random variable type: " << dist
               << " provided. Must be one of "
               << "\n";
        return -1;
    }

    if (theRandomVariable == 0) {
        opserr << "ERROR: could not create random variable number " << tag
               << "\n";
        return -1;
    }

    // set start point on object if user provided
    if (use_start_pt == 1) {
        theRandomVariable->setStartValue(startPt);
        theRandomVariable->setCurrentValue(startPt);
    } else {
        theRandomVariable->setStartValue(theRandomVariable->getMean());
        theRandomVariable->setCurrentValue(theRandomVariable->getMean());
    }

    // Add the random variable to the domain
    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain->addRandomVariable(theRandomVariable) ==
        false) {
        opserr << "ERROR: failed to add random variable to the domain "
                  "(wrong number of arguments?)\n";
        opserr << "random variable: " << tag << "\n";
        delete theRandomVariable;  // otherwise memory leak
        return -1;
    }

    // RVParameter *theRVParam = new RVParameter(tag, theRandomVariable);
    // theStructuralDomain->addParameter(theRVParam);

    return 0;
}

int OPS_getRVTags() {
    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) return -1;

    std::vector<int> rvTags;
    RandomVariable *theRV;
    RandomVariableIter &rvIter =
        theReliabilityDomain->getRandomVariables();
    while ((theRV = rvIter()) != 0) rvTags.push_back(theRV->getTag());

    int size = 0;
    int *data = 0;
    if (!rvTags.empty()) {
        size = (int)rvTags.size();
        data = &rvTags[0];
    }

    if (OPS_SetIntOutput(&size, data, false) < 0) {
        opserr << "ERROR: failed to set outputs in getRVTags" << endln;
        return -1;
    }

    return 0;
}

int OPS_getRVParamTag() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: need getRVParamTag rvTag\n";
        return -1;
    }

    // random variable tag input
    int rvTag;
    int num = 1;
    if (OPS_GetIntInput(&num, &rvTag) < 0) {
        opserr << "ERROR: failed to get rvTag\n";
        return -1;
    }

    // reliability domain
    auto *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "ERROR: reliability domain is null\n";
        return -1;
    }

    // get random variable pointer
    auto *rv = theReliabilityDomain->getRandomVariablePtr(rvTag);
    if (rv == 0) {
        opserr << "ERROR: random variable with tag " << rvTag
               << " not found\n";
        return -1;
    }

    // get random variable index
    int rvIndex = theReliabilityDomain->getRandomVariableIndex(rvTag);

    // get parameter index
    int paramIndex =
        theReliabilityDomain->getParameterIndexFromRandomVariableIndex(
            rvIndex);
    if (paramIndex < 0) {
        opserr
            << "ERROR: failed to get parameter index for random variable"
            << rvTag << "\n";
        return -1;
    }

    // get parameter pointer
    auto *domain = cmds->getStructuralDomain();
    if (domain == 0) {
        opserr << "ERROR: domain is null\n";
        return -1;
    }

    auto *param = domain->getParameterFromIndex(paramIndex);
    if (param == 0) {
        opserr << "ERROR: failed to get parameter for random variable"
               << rvTag << "\n";
        return -1;
    }

    int paramTag = param->getTag();
    if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
        opserr << "ERROR: failed to set paramTag output\n";
        return -1;
    }

    return 0;
}

int OPS_getRVValue() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: need getRVValue rvTag\n";
        return -1;
    }

    // random variable tag input
    int rvTag;
    int num = 1;
    if (OPS_GetIntInput(&num, &rvTag) < 0) {
        opserr << "ERROR: failed to get rvTag\n";
        return -1;
    }

    // reliability domain
    auto *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "ERROR: reliability domain is null\n";
        return -1;
    }

    // get random variable pointer
    auto *rv = theReliabilityDomain->getRandomVariablePtr(rvTag);
    if (rv == 0) {
        opserr << "ERROR: random variable with tag " << rvTag
               << " not found\n";
        return -1;
    }

    // get value
    double value = rv->getCurrentValue();

    if (OPS_SetDoubleOutput(&num, &value, true) < 0) {
        opserr << "ERROR: failed to set value output\n";
        return -1;
    }

    return 0;
}

int OPS_getRVMean() {
    // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: invalid number of arguments to getMean command "
                  ": getMean rvTag\n";
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
        opserr << "ERROR: getMean - random variable with tag " << rvTag
               << " not found" << endln;
        return -1;
    }

    double mean = rv->getMean();
    if (OPS_SetDoubleOutput(&numData, &mean, true) < 0) {
        opserr << "ERROR: getMean - failed to set double output\n";
        return -1;
    }

    return 0;
}

int OPS_getRVStdv() {
    // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: invalid number of arguments to getStdv command "
                  ": getStdv rvTag\n";
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
        opserr << "ERROR: getStdv - random variable with tag " << rvTag
               << " not found" << endln;
        return -1;
    }

    double stdv = rv->getStdv();
    if (OPS_SetDoubleOutput(&numData, &stdv, true) < 0) {
        opserr << "ERROR: getStdv - failed to set double output\n";
        return -1;
    }

    return 0;
}

int OPS_getRVPDF() {
    // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "ERROR: invalid number of arguments to getPDF command : "
                  "getPDF rvTag X\n";
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
        opserr << "ERROR: getPDF - random variable with tag " << rvTag
               << " not found" << endln;
        return -1;
    }

    double pdf = rv->getPDFvalue(x);
    if (OPS_SetDoubleOutput(&numData, &pdf, true) < 0) {
        opserr << "ERROR: getPDF - failed to set double output\n";
        return -1;
    }

    return 0;
}

int OPS_getRVCDF() {
    // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "ERROR: invalid number of arguments to getCDF command : "
                  "getCDF rvTag X\n";
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
        opserr << "ERROR: getCDF - random variable with tag " << rvTag
               << " not found" << endln;
        return -1;
    }

    double cdf = rv->getCDFvalue(x);
    if (OPS_SetDoubleOutput(&numData, &cdf, true) < 0) {
        opserr << "ERROR: getCDF - failed to set double output\n";
        return -1;
    }

    return 0;
}

int OPS_getRVInverseCDF() {
    // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "ERROR: invalid number of arguments to getInverseCDF "
                  "command : getInverseCDF rvTag p\n";
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
        opserr << "ERROR: getInverseCDF - random variable with tag "
               << rvTag << " not found" << endln;
        return -1;
    }

    double invcdf = rv->getInverseCDFvalue(p);
    if (OPS_SetDoubleOutput(&numData, &invcdf, true) < 0) {
        opserr << "ERROR: getInverseCDF - failed to set double output\n";
        return -1;
    }

    return 0;
}

int OPS_getLSFTags() {
    if (cmds == 0) {
        opserr << "WARNING: reliability cmds not defined\n";
        return -1;
    }
    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    LimitStateFunction *theLSF;
    LimitStateFunctionIter &lsfIter =
        theReliabilityDomain->getLimitStateFunctions();

    std::vector<int> tags;
    while ((theLSF = lsfIter()) != 0) {
        tags.push_back(theLSF->getTag());
    }

    int size = 0;
    int *data = 0;
    if (!tags.empty()) {
        size = (int)tags.size();
        data = &tags[0];
    }

    if (OPS_SetIntOutput(&size, data, false) < 0) {
        opserr << "WARNING failed to set outputs\n";
        return -1;
    }

    return 0;
}

int OPS_addCorrelate() {
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "ERROR: Wrong number of arguments to correlate command"
               << endln;
        return -1;
    }

    int rvTag[2];
    double correlationValue;

    int numData = 2;
    if (OPS_GetIntInput(&numData, rvTag) < 0) {
        opserr << "ERROR: invalid input to correlate: tag" << endln;
        ;
        return -1;
    }

    numData = 1;
    if (OPS_GetDoubleInput(&numData, &correlationValue) < 0) {
        opserr << "ERROR: invalid input to correlate: value" << endln;
        ;
        return -1;
    }

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    int tag = theReliabilityDomain->getNumberOfCorrelationCoefficients();
    CorrelationCoefficient *theCorrelationCoefficient =
        new CorrelationCoefficient(tag + 1, rvTag[0], rvTag[1],
                                   correlationValue);
    if (theCorrelationCoefficient == 0) {
        opserr << "ERROR: failed to add correlation coefficient to domain"
               << endln;
        return -1;
    }

    if (theReliabilityDomain->addCorrelationCoefficient(
            theCorrelationCoefficient) == false) {
        opserr
            << "ERROR: failed to add correlation coefficient to domain\n";
        opserr << "tag, rv1, rv2: " << tag << ' ' << rvTag[0] << ' '
               << rvTag[1] << endln;
        return -1;
    }

    return 0;
}

void OpenSeesReliabilityCommands::setProbabilityTransformation(
    ProbabilityTransformation *transform) {
    if (theProbabilityTransformation != 0) {
        delete theProbabilityTransformation;
        theProbabilityTransformation = 0;
    }

    theProbabilityTransformation = transform;
    if (transform == 0) return;
}

void OpenSeesReliabilityCommands::setRandomNumberGenerator(
    RandomNumberGenerator *generator) {
    if (theRandomNumberGenerator != 0) {
        delete theRandomNumberGenerator;
        theRandomNumberGenerator = 0;
    }

    theRandomNumberGenerator = generator;
    if (generator == 0) return;
}

void OpenSeesReliabilityCommands::setReliabilityConvergenceCheck(
    ReliabilityConvergenceCheck *check) {
    if (theReliabilityConvergenceCheck != 0) {
        delete theReliabilityConvergenceCheck;
        theReliabilityConvergenceCheck = 0;
    }

    theReliabilityConvergenceCheck = check;
    if (check == 0) return;
}

void OpenSeesReliabilityCommands::setSearchDirection(
    SearchDirection *search) {
    if (theSearchDirection != 0) {
        delete theSearchDirection;
        theSearchDirection = 0;
    }

    theSearchDirection = search;
    if (search == 0) return;
}

void OpenSeesReliabilityCommands::setMeritFunctionCheck(
    MeritFunctionCheck *merit) {
    if (theMeritFunctionCheck != 0) {
        delete theMeritFunctionCheck;
        theMeritFunctionCheck = 0;
    }

    theMeritFunctionCheck = merit;
    if (merit == 0) return;
}

void OpenSeesReliabilityCommands::setStepSizeRule(StepSizeRule *rule) {
    if (theStepSizeRule != 0) {
        delete theStepSizeRule;
        theStepSizeRule = 0;
    }

    theStepSizeRule = rule;
    if (rule == 0) return;
}

void OpenSeesReliabilityCommands::setRootFinding(RootFinding *root) {
    if (theRootFinding != 0) {
        delete theRootFinding;
        theRootFinding = 0;
    }

    theRootFinding = root;
    if (root == 0) return;
}

void OpenSeesReliabilityCommands::setFindDesignPointAlgorithm(
    FindDesignPointAlgorithm *algo) {
    if (theFindDesignPointAlgorithm != 0) {
        delete theFindDesignPointAlgorithm;
        theFindDesignPointAlgorithm = 0;
    }
    theFindDesignPointAlgorithm = algo;
}

void OpenSeesReliabilityCommands::setFindCurvatures(
    FindCurvatures *algo) {
    if (theFindCurvatures != 0) {
        delete theFindCurvatures;
        theFindCurvatures = 0;
    }
    theFindCurvatures = algo;
}

void OpenSeesReliabilityCommands::setGradientEvaluator(
    GradientEvaluator *eval) {
    if (theGradientEvaluator != 0) {
        delete theGradientEvaluator;
        theGradientEvaluator = 0;
    }

    theGradientEvaluator = eval;
    if (eval == 0) return;
}

void OpenSeesReliabilityCommands::setFunctionEvaluator(
    FunctionEvaluator *eval) {
    if (theFunctionEvaluator != 0) {
        delete theFunctionEvaluator;
        theFunctionEvaluator = 0;
    }

    theFunctionEvaluator = eval;
    if (eval == 0) {
        return;
    }
}

void OpenSeesReliabilityCommands::setFOSMAnalysis(FOSMAnalysis *analysis) {
    if (theFOSMAnalysis != 0) {
        delete theFOSMAnalysis;
        theFOSMAnalysis = 0;
    }
    theFOSMAnalysis = analysis;
}

void OpenSeesReliabilityCommands::setFORMAnalysis(FORMAnalysis *analysis) {
    if (theFORMAnalysis != 0) {
        delete theFORMAnalysis;
        theFORMAnalysis = 0;
    }
    theFORMAnalysis = analysis;
}

void OpenSeesReliabilityCommands::setSORMAnalysis(SORMAnalysis *analysis) {
    if (theSORMAnalysis != 0) {
        delete theSORMAnalysis;
        theSORMAnalysis = 0;
    }
    theSORMAnalysis = analysis;
}

void OpenSeesReliabilityCommands::setImportanceSamplingAnalysis(
    ImportanceSamplingAnalysis *analysis) {
    if (theImportanceSamplingAnalysis != 0) {
        delete theImportanceSamplingAnalysis;
        theImportanceSamplingAnalysis = 0;
    }
    theImportanceSamplingAnalysis = analysis;
}

void OpenSeesReliabilityCommands::setMonteCarloAnalysis(
    MonteCarloResponseAnalysis *analysis) {
    if (theMonteCarloAnalysis != 0) {
        delete theMonteCarloAnalysis;
        theMonteCarloAnalysis = 0;
    }
    theMonteCarloAnalysis = analysis;
}

int OPS_startPoint() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to startPoint"
               << endln;
        return -1;
    }

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    int nrv = theReliabilityDomain->getNumberOfRandomVariables();
    RandomVariable *aRandomVariable;

    // Get the type of start point
    const char *type = OPS_GetString();
    int meanOrZero = -1;
    if (strcmp(type, "Mean") == 0) meanOrZero = 1;
    if (strcmp(type, "Zero") == 0) meanOrZero = 0;
    if (strcmp(type, "Origin") == 0) meanOrZero = 0;

    RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
    while ((aRandomVariable = rvIter()) != 0) {
        if (meanOrZero == 0) {
            aRandomVariable->setStartValue(0.0);
        }
        if (meanOrZero == 1) {
            double mean = aRandomVariable->getMean();
            aRandomVariable->setStartValue(mean);
        }
    }

    if (meanOrZero >= 0) {
        return 0;
    }

    // file type
    if (strcmp(type, "-file") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WARNING: need file name which is space "
                      "delimited and contains a starting point\n";
            return -1;
        }

        const char *filename = OPS_GetString();
        std::ifstream inputFile(filename, std::ios::in);
        if (inputFile.fail()) {
            opserr << "File " << filename
                   << " could not be opened for startPoint.\n";
            return -1;
        }

        // Loop through file to see how many entries there are
        double dummy;
        int numEntries = 0;
        while (inputFile >> dummy) {
            numEntries++;
        }

        if (numEntries == 0) {
            opserr << "ERROR: No entries in the file read by "
                      "startPoint!\n";
            return -1;
        }
        if (numEntries != nrv) {
            opserr << "ERROR: Wrong number of entries in the file "
                      "read by startPoint.\n";
            return -1;
        }

        // rewind the file and pass values to the RVs
        inputFile.seekg(0, ios::beg);
        for (int i = 0; i < nrv; i++) {
            aRandomVariable =
                theReliabilityDomain->getRandomVariablePtrFromIndex(i);
            inputFile >> dummy;
            aRandomVariable->setStartValue(dummy);
        }

        inputFile.close();
        return 0;
    }

    opserr << "ERROR: Invalid type of start point is given.\n";
    return -1;
}

int OPS_randomNumberGenerator() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to "
                  "randomNumberGenerator"
               << endln;
        return -1;
    }

    // Get the type of generator
    const char *type = OPS_GetString();
    if (strcmp(type, "CStdLib") != 0) {
        opserr << "ERROR: unrecognized type of RandomNumberGenerator "
               << type << endln;
        return -1;
    }

    RandomNumberGenerator *theGenerator = new CStdLibRandGenerator();
    if (theGenerator == 0) {
        opserr << "ERROR: could not create randomNumberGenerator" << endln;
        return -1;
    }
    if (cmds != 0) {
        cmds->setRandomNumberGenerator(theGenerator);
    }

    return 0;
}

int OPS_reliabilityConvergenceCheck() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to "
                  "reliabilityConvergenceCheck"
               << endln;
        return -1;
    }

    // Get the type of convergence check
    const char *type = OPS_GetString();

    double e1 = 1.0e-3;
    double e2 = 1.0e-3;
    double scaleValue = 0.0;
    int print = 1;

    // Standard and Optimality convergence checks have the same
    // optional arguments
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        int numdata = 1;
        if (strcmp(arg, "-e1") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetDoubleInput(&numdata, &e1) < 0) {
                opserr << "ERROR: unable to read -e1 value for "
                          "reliability convergence check"
                       << endln;
                return -1;
            }
        }
        if (strcmp(arg, "-e2") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetDoubleInput(&numdata, &e2) < 0) {
                opserr << "ERROR: unable to read -e2 value for "
                          "reliability convergence check"
                       << endln;
                return -1;
            }
        }
        if (strcmp(arg, "-scaleValue") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetDoubleInput(&numdata, &scaleValue) < 0) {
                opserr << "ERROR: unable to read -scaleValue value for "
                          "reliability convergence check"
                       << endln;
                return -1;
            }
        }
        if (strcmp(arg, "-print") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numdata, &print) < 0) {
                opserr << "ERROR: unable to read -print value for "
                          "reliability convergence check"
                       << endln;
                return -1;
            }
        }
    }

    ReliabilityConvergenceCheck *theCheck = 0;
    if (strcmp(type, "Standard") == 0) {
        theCheck = new StandardReliabilityConvergenceCheck(
            e1, e2, scaleValue, print);
    } else if (strcmp(type, "OptimalityCondition") == 0) {
        theCheck = new OptimalityConditionReliabilityConvergenceCheck(
            e1, e2, scaleValue, print);
    } else {
        opserr << "ERROR: unrecognized type of "
                  "reliabilityConvergenceCheck "
               << type << endln;
        return -1;
    }

    if (theCheck == 0) {
        opserr << "ERROR: could not create reliabilityConvergenceCheck"
               << endln;
        return -1;
    } else {
        if (cmds != 0) cmds->setReliabilityConvergenceCheck(theCheck);
    }

    return 0;
}

int OPS_searchDirection() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to searchDirection"
               << endln;
        return -1;
    }

    // Get the type of search direction
    const char *type = OPS_GetString();

    SearchDirection *theSearch = 0;
    if (strcmp(type, "iHLRF") == 0) {
        theSearch = new HLRFSearchDirection();

    } else if (strcmp(type, "PolakHe") == 0) {
        double gamma = 1.0;
        double delta = 1.0;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;
            if (strcmp(arg, "-gamma") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &gamma) < 0) {
                    opserr << "ERROR: unable to read -gamma value for "
                              "PolakHe search direction"
                           << endln;
                    return -1;
                }
            } else if (strcmp(arg, "-delta") == 0 &&
                       OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &delta) < 0) {
                    opserr << "ERROR: unable to read -delta value for "
                              "PolakHe search direction"
                           << endln;
                    return -1;
                }
            } else {
                opserr << "ERROR: Invalid input to Polak-He algorithm.\n";
                return -1;
            }
        }
        cmds->setPolakHeDualPurpose(
            new PolakHeSearchDirectionAndMeritFunction(gamma, delta));
        theSearch = cmds->getPolakHeDualPurpose();

    } else if (strcmp(type, "GradientProjection") == 0) {
        // Check that a step size rule has been created
        StepSizeRule *theStepSizeRule = cmds->getStepSizeRule();
        if (theStepSizeRule == 0) {
            opserr << "Need theStepSizeRule before a "
                      "GradientProjectionSearchDirection can be "
                      "created\n";
            return -1;
        }

        // Check that a transformation has been created
        ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
        ProbabilityTransformation *theProbabilityTransformation =
            cmds->getProbabilityTransformation();
        if (theProbabilityTransformation == 0) {
            opserr << "Assume all RV's are independent" << endln;
            theProbabilityTransformation =
                new AllIndependentTransformation(theReliabilityDomain, 0);
            cmds->setProbabilityTransformation(
                theProbabilityTransformation);
        }

        // Check that a gfun evaluator has been created
        FunctionEvaluator *theFunctionEvaluator =
            cmds->getFunctionEvaluator();
        if (theFunctionEvaluator == 0) {
            opserr << "Need theGFunEvaluator before a "
                      "GradientProjectionSearchDirection can be "
                      "created\n";
            return -1;
        }

        // Check that a root-finding algorithm has been created
        RootFinding *theRootFindingAlgorithm = cmds->getRootFinding();
        if (theRootFindingAlgorithm == 0) {
            opserr << "Need theRootFindingAlgorithm before a "
                      "GradientProjectionSearchDirection can be "
                      "created\n";
            return -1;
        }

        theSearch = new GradientProjectionSearchDirection(
            theStepSizeRule, theProbabilityTransformation,
            theFunctionEvaluator, theRootFindingAlgorithm);

    } else if (strcmp(type, "SQP") == 0) {
        double c_bar = 200.0;
        double e_bar = 0.5;

        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;

            if (strcmp(arg, "-c_bar") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &c_bar) < 0) {
                    opserr << "ERROR: invalid input: c_bar for SQP "
                              "algorithm\n";
                    return -1;
                }

            } else if (strcmp(arg, "-e_bar") == 0 &&
                       OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &e_bar) < 0) {
                    opserr << "ERROR: invalid input: e_bar for SQP "
                              "algorithm\n";
                    return -1;
                }

            } else {
                opserr << "ERROR: Invalid input to SQP algorithm.\n";
                return -1;
            }
        }

        cmds->setSQPtriplePurpose(
            new SQPsearchDirectionMeritFunctionAndHessian(c_bar, e_bar));
        theSearch = cmds->getSQPtriplePurpose();

    } else {
        opserr << "ERROR: unrecognized type of searchDirection " << type
               << endln;
        return -1;
    }

    if (theSearch == 0) {
        opserr << "ERROR: could not create searchDirection" << endln;
        return -1;
    } else {
        if (cmds != 0) cmds->setSearchDirection(theSearch);
    }

    return 0;
}

int OPS_meritFunctionCheck() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to "
                  "meritFunctionCheck"
               << endln;
        return -1;
    }

    // Get the type of merit function check
    const char *type = OPS_GetString();

    MeritFunctionCheck *theFunction = 0;
    if (strcmp(type, "AdkZhang") == 0) {
        double multi = 2.0;
        double add = 10.0;
        double factor = 0.0;

        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;
            if (strcmp(arg, "-multi") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &multi) < 0) {
                    opserr << "ERROR: unable to read -multi value for "
                              "AdkZhang merit function check"
                           << endln;
                    return -1;
                }
            } else if (strcmp(arg, "-add") == 0 &&
                       OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &add) < 0) {
                    opserr << "ERROR: unable to read -add value for "
                              "AdkZhang merit function check"
                           << endln;
                    return -1;
                }
            } else if (strcmp(arg, "-factor") == 0 &&
                       OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &factor) < 0) {
                    opserr << "ERROR: unable to read -factor value for "
                              "AdkZhang merit function check"
                           << endln;
                    return -1;
                }
            } else {
                opserr << "ERROR: Invalid input to AdkZhang merit "
                          "function check. \n";
                return -1;
            }
        }

        // Do a quick input check
        if (multi < 1.0 || add < 0.0) {
            opserr << "ERROR: Invalid values of multi/add parameters "
                      "to AdkZhang merit function check"
                   << endln;
            return -1;
        }

        theFunction = new AdkZhangMeritFunctionCheck(multi, add, factor);

    } else if (strcmp(type, "PolakHe") == 0) {
        PolakHeSearchDirectionAndMeritFunction *thePolakHeDualPurpose =
            cmds->getPolakHeDualPurpose();
        if (thePolakHeDualPurpose == 0) {
            opserr << "Need thePolakHeSearchDirection before a "
                      "PolakHe merit function can be created\n";
            return -1;
        }

        double factor = 0.5;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;
            if (strcmp(arg, "-factor") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &factor) < 0) {
                    opserr << "ERROR: invalid input: factor \n";
                    return -1;
                }
            } else {
                opserr << "ERROR: Invalid input to Polak He merit "
                          "function check.\n";
                return -1;
            }
        }

        thePolakHeDualPurpose->setAlpha(factor);
        theFunction = thePolakHeDualPurpose;

    } else if (strcmp(type, "SQP") == 0) {
        // Check that the SQP search direction is already created
        SQPsearchDirectionMeritFunctionAndHessian *theSQPtriplePurpose =
            cmds->getSQPtriplePurpose();
        if (theSQPtriplePurpose == 0) {
            opserr << "Need theSQPSearchDirection before a SQP merit "
                      "function can be created\n";
            return -1;
        }

        double factor = 0.5;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;
            if (strcmp(arg, "-factor") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &factor) < 0) {
                    opserr << "ERROR: invalid input: factor \n";
                    return -1;
                }
            } else {
                opserr << "ERROR: Invalid input to SQP merit "
                          "function check.\n";
                return -1;
            }
        }

        theSQPtriplePurpose->setAlpha(factor);
        theFunction = theSQPtriplePurpose;

    } else {
        opserr << "ERROR: unrecognized type of meritFunctionCheck " << type
               << endln;
        return -1;
    }

    if (theFunction == 0) {
        opserr << "ERROR: could not create meritFunctionCheck" << endln;
        return -1;
    } else {
        if (cmds != 0) cmds->setMeritFunctionCheck(theFunction);
    }

    return 0;
}

int OPS_stepSizeRule() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to stepSizeRule"
               << endln;
        return -1;
    }

    // Get the type of step size rule
    const char *type = OPS_GetString();

    StepSizeRule *theRule = 0;
    if (strcmp(type, "Armijo") == 0) {
        double base = 0.5;
        int maxNumReductions = 10;
        double b0 = 1.0;
        int numberOfShortSteps = 2;
        double radius = 50.0;
        double surfaceDistance = 0.1;
        double evolution = 0.5;
        int printFlag = 0;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;
            if (strcmp(arg, "-print") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &printFlag) < 0) {
                    opserr << "ERROR: unable to read -print value for "
                              "Armijo step size rule"
                           << endln;
                    return -1;
                }
            }
            if (strcmp(arg, "-maxNum") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &maxNumReductions) < 0) {
                    opserr << "ERROR: unable to read -maxNum value for "
                              "Armijo step size rule"
                           << endln;
                    return -1;
                }
            }
            if (strcmp(arg, "-base") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &base) < 0) {
                    opserr << "ERROR: unable to read -base value for "
                              "Armijo step size rule"
                           << endln;
                    return -1;
                }
            }
            if (strcmp(arg, "-initial") == 0 &&
                OPS_GetNumRemainingInputArgs() > 1) {
                if (OPS_GetDoubleInput(&numdata, &b0) < 0) {
                    opserr << "ERROR: unable to read -initial b0 value "
                              "for Armijo step size rule"
                           << endln;
                    return -1;
                }
                if (OPS_GetIntInput(&numdata, &numberOfShortSteps) < 0) {
                    opserr << "ERROR: unable to read -initial "
                              "numberOfShortSteps value for Armijo step "
                              "size rule"
                           << endln;
                    return -1;
                }
            }
            if (strcmp(arg, "-sphere") == 0 &&
                OPS_GetNumRemainingInputArgs() > 2) {
                if (OPS_GetDoubleInput(&numdata, &radius) < 0) {
                    opserr << "ERROR: unable to read -sphere radius value "
                              "for Armijo step size rule"
                           << endln;
                    return -1;
                }
                if (OPS_GetDoubleInput(&numdata, &surfaceDistance) < 0) {
                    opserr
                        << "ERROR: unable to read -sphere surfaceDistance "
                           "value for Armijo step size rule"
                        << endln;
                    return -1;
                }
                if (OPS_GetDoubleInput(&numdata, &evolution) < 0) {
                    opserr << "ERROR: unable to read -sphere evolution "
                              "value for Armijo step size rule"
                           << endln;
                    return -1;
                }
            }
        }

        ReliabilityDomain *theRelDomain = cmds->getDomain();
        RootFinding *theAlgorithm = cmds->getRootFinding();
        ProbabilityTransformation *theTransformation =
            cmds->getProbabilityTransformation();
        if (theTransformation == 0) {
            opserr << "Assume all RV's are independent" << endln;
            theTransformation =
                new AllIndependentTransformation(theRelDomain, 0);
            cmds->setProbabilityTransformation(theTransformation);
        }
        FunctionEvaluator *theEvaluator = cmds->getFunctionEvaluator();
        if (theEvaluator == 0) {
            opserr << "Function evaluator must be defined before "
                      "ArmijoStepSize rule"
                   << endln;
            return -1;
        }
        MeritFunctionCheck *theMeritFunction =
            cmds->getMeritFunctionCheck();
        if (theMeritFunction == 0) {
            opserr << "Merit function check must be defined before "
                      "ArmijoStepSize rule"
                   << endln;
            return -1;
        }

        theRule = new ArmijoStepSizeRule(
            theRelDomain, theEvaluator, theTransformation,
            theMeritFunction, theAlgorithm, base, maxNumReductions, b0,
            numberOfShortSteps, radius, surfaceDistance, evolution,
            printFlag);
    } else if (strcmp(type, "Fixed") == 0) {
        double stepSize = 1.0;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;
            if (strcmp(arg, "-stepSize") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &stepSize) < 0) {
                    opserr << "ERROR: unable to read -stepSize value for "
                              "Fixed step size rule"
                           << endln;
                    return -1;
                }
            }
        }
        theRule = new FixedStepSizeRule(stepSize);
    } else {
        opserr << "ERROR: unrecognized type of stepSizeRule " << type
               << endln;
        return -1;
    }

    if (theRule == 0) {
        opserr << "ERROR: could not create stepSizeRule" << endln;
        return -1;
    } else {
        if (cmds != 0) cmds->setStepSizeRule(theRule);
    }

    return 0;
}

int OPS_rootFinding() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to rootFinding"
               << endln;
        return -1;
    }

    // Get the type of rootFinding
    const char *type = OPS_GetString();

    int maxIter = 50;
    double tol = 1.0e-3;
    double maxStepLength = 1.0;
    RootFinding *theFinding = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        int numdata = 1;
        if (strcmp(arg, "-maxIter") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numdata, &maxIter) < 0) {
                opserr << "ERROR: unable to read -maxIter value for "
                       << type << " root finding" << endln;
                return -1;
            }
        }
        if (strcmp(arg, "-tol") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
                opserr << "ERROR: unable to read -tol value for " << type
                       << " root finding" << endln;
                return -1;
            }
        }
        if (strcmp(arg, "-maxStepLength") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetDoubleInput(&numdata, &maxStepLength) < 0) {
                opserr << "ERROR: unable to read -maxStepLength value for "
                       << type << " root finding" << endln;
                return -1;
            }
        }
    }

    if (strcmp(type, "Secant") == 0) {
        ReliabilityDomain *theRelDomain = cmds->getDomain();
        ProbabilityTransformation *theTransformation =
            cmds->getProbabilityTransformation();
        if (theTransformation == 0) {
            opserr << "Assume all RV's are independent" << endln;
            theTransformation =
                new AllIndependentTransformation(theRelDomain, 0);
            cmds->setProbabilityTransformation(theTransformation);
        }
        FunctionEvaluator *theEvaluator = cmds->getFunctionEvaluator();
        if (theEvaluator == 0) {
            opserr << "Function evaluator must be defined before "
                      "ArmijoStepSize rule"
                   << endln;
            return -1;
        }
        theFinding = new SecantRootFinding(theRelDomain, theTransformation,
                                           theEvaluator, maxIter, tol,
                                           maxStepLength);
    } else {
        opserr << "ERROR: unrecognized type of rootFinding: " << type
               << endln;
        return -1;
    }

    if (theFinding == 0) {
        opserr << "ERROR: could not create rootFinding" << endln;
        return -1;
    } else {
        if (cmds != 0) cmds->setRootFinding(theFinding);
    }

    return 0;
}

int OPS_findCurvatures() {
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "ERROR: wrong number of arguments to findCurvatures" << endln;
    return -1;
  }
  if (cmds == 0) {
    opserr << "WARNING: reliability cmds not defined\n";
    return -1;
  }

  // type
  const char *type = OPS_GetString();

  // Check that the necessary ingredients are present
  ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
  if (theReliabilityDomain == 0) {
    opserr << "Need theReliabilityDomain before find curvatures" << endln;
    return -1;
  }

  FunctionEvaluator *theFunctionEvaluator = cmds->getFunctionEvaluator();
  if (theFunctionEvaluator == 0) {
    opserr << "Need theFunctionEvaluator before find curvatures" << endln;    
    return -1;
  }  

  FORMAnalysis *theFORMAnalysis = cmds->getFORMAnalysis();
  if (theFORMAnalysis == 0) {
    opserr << "FORMAnalysis must be performed prior to find curvatures" << endln;    
    return -1;
  }
  
  // 
  FindCurvatures *theFindCurvatures = 0;
  if (strcmp(type, "firstPrincipal") == 0) {
    theFindCurvatures = new FirstPrincipalCurvature(theReliabilityDomain,
						    theFunctionEvaluator,
						    theFORMAnalysis);
  }
  else if (strcmp(type, "bySearchAlgorithm") == 0) {
    int numCurvatures = -1;
    int numData = 1;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      if (OPS_GetIntInput(&numData, &numCurvatures) < 0) {
	opserr << "ERROR: unable to read numCurvatures for " << type
	       << " curvature finding" << endln;
	return -1;
      }
    }

    theFindCurvatures = new CurvaturesBySearchAlgorithm(theReliabilityDomain,
							theFunctionEvaluator,
							theFORMAnalysis, numCurvatures);    
  }
  else if (strcmp(type, "curvatureFitting") == 0) {
    opserr << "curvatureFitting not in Python interpreter yet, also need to add Hessian" << endln;
    return -1;
  }
  else {
    opserr << "ERROR: unrecognized type of FindCurvatures strategy" << endln;
    return -1;
  }

  if (theFindCurvatures == 0) {
    opserr << "ERROR: could not create theFindCurvatures" << endln;
    return -1;
  }
  
  cmds->setFindCurvatures(theFindCurvatures);
  
  return 0;
}

int OPS_findDesignPoint() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to "
                  "findDesignPoint\n";
        return -1;
    }
    if (cmds == 0) {
        opserr << "WARNING: reliability cmds not defined\n";
        return -1;
    }

    // args
    int printFlag = 0;
    char fileNamePrint[256];
    strcpy(fileNamePrint, "initialized");
    int maxNumIter = 100;

    // type
    const char *type = OPS_GetString();

    // Check that the necessary ingredients are present
    FunctionEvaluator *theFunctionEvaluator = cmds->getFunctionEvaluator();
    if (theFunctionEvaluator == 0) {
        opserr << "Need theGFunEvaluator before a "
                  "FindDesignPointAlgorithm can be created\n";
        return -1;
    }
    GradientEvaluator *theGradientEvaluator = cmds->getGradientEvaluator();
    if (theGradientEvaluator == 0) {
        opserr << "Need theGradientEvaluator before a "
                  "FindDesignPointAlgorithm can be created\n";
        return -1;
    }
    StepSizeRule *theStepSizeRule = cmds->getStepSizeRule();
    if (theStepSizeRule == 0) {
        opserr << "Need theStepSizeRule before a "
                  "FindDesignPointAlgorithm can be created\n";
        return -1;
    }
    SearchDirection *theSearchDirection = cmds->getSearchDirection();
    if (theSearchDirection == 0) {
        opserr << "Need theSearchDirection before a "
                  "FindDesignPointAlgorithm can be created\n";
        return -1;
    }

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "Need theReliabilityDomain before a "
                  "FindDesignPointAlgorithm can be created\n";
        return -1;
    }

    Domain *theStructuralDomain = cmds->getStructuralDomain();
    if (theReliabilityDomain == 0) {
        opserr << "Need theStructuralDomain before a "
                  "FindDesignPointAlgorithm can be created\n";
        return -1;
    }

    ProbabilityTransformation *theProbabilityTransformation =
        cmds->getProbabilityTransformation();
    if (theProbabilityTransformation == 0) {
        opserr << "Assume all RV's are independent" << endln;
        theProbabilityTransformation =
            new AllIndependentTransformation(theReliabilityDomain, 0);
	cmds->setProbabilityTransformation(theProbabilityTransformation);
    }

    ReliabilityConvergenceCheck *theReliabilityConvergenceCheck =
        cmds->getReliabilityConvergenceCheck();
    if (theReliabilityConvergenceCheck == 0) {
        opserr << "Need theReliabilityConvergenceCheck before a "
                  "FindDesignPointAlgorithm can be created\n";
        return -1;
    }

    // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
    int numdata = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *subtype = OPS_GetString();

        if (strcmp(subtype, "-maxNumIter") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numdata, &maxNumIter) < 0) {
                opserr << "ERROR: invalid input: maxNumIter \n";
                return -1;
            }
        } else if (strcmp(subtype, "-printAllPointsX") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            printFlag = 1;
            const char *name = OPS_GetString();
            strcpy(fileNamePrint, name);
        } else if (strcmp(subtype, "-printAllPointsY") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            printFlag = 2;
            const char *name = OPS_GetString();
            strcpy(fileNamePrint, name);
        } else if (strcmp(subtype, "-printDesignPointX") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            printFlag = 3;
            const char *name = OPS_GetString();
            strcpy(fileNamePrint, name);
        } else if (strcmp(subtype, "-printDesignPointY") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            printFlag = 4;
            const char *name = OPS_GetString();
            strcpy(fileNamePrint, name);
        } else if (strcmp(subtype, "-printCurrentPointX") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            printFlag = 5;
            const char *name = OPS_GetString();
            strcpy(fileNamePrint, name);
        } else if (strcmp(subtype, "-printCurrentPointY") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            printFlag = 6;
            const char *name = OPS_GetString();
            strcpy(fileNamePrint, name);
        } else {
            opserr << "ERROR: Invalid input to "
                      "SearchWithStepSizeAndStepDirection. \n";
            return -1;
        }
    }

    FindDesignPointAlgorithm *theFindDesignPointAlgorithm = 0;
    if (strcmp(type, "StepSearch") == 0) {
        theFindDesignPointAlgorithm =
            new SearchWithStepSizeAndStepDirection(
                maxNumIter, theReliabilityDomain, theStructuralDomain,
                theFunctionEvaluator, theGradientEvaluator,
                theStepSizeRule, theSearchDirection,
                theProbabilityTransformation,
                theReliabilityConvergenceCheck, printFlag, fileNamePrint);

    } else {
        opserr << "ERROR: unrecognized type of "
                  "FindDesignPointAlgorithm Algorithm \n";
        return -1;
    }

    if (theFindDesignPointAlgorithm == 0) {
        opserr << "ERROR: could not create "
                  "theFindDesignPointAlgorithm \n";
        return -1;
    }

    cmds->setFindDesignPointAlgorithm(theFindDesignPointAlgorithm);

    return 0;
}

int OPS_functionEvaluator() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to functionEvaluator"
               << endln;
        return -1;
    }

    if (cmds == 0) {
        opserr << "WARNING: Reliability is not initialized\n";
        return -1;
    }
    if (cmds->getStructuralDomain() == 0) {
        opserr << "WARNING: Reliability has no structural domain\n";
        return -1;
    }
    if (cmds->getDomain() == 0) {
        opserr << "WARNING: Reliability has no domain\n";
        return -1;
    }

    FunctionEvaluator *theEval = 0;

    // Get the type of functionEvaluator
    const char *type = OPS_GetString();
    const char *filename = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        filename = OPS_GetString();
        if (strcmp(filename, "-file") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            filename = OPS_GetString();
        }
    }
    if (strcmp(type, "Matlab") == 0) {
        opserr << "ERROR: Matlab function evaluator not implemented"
               << endln;
        return -1;
    } else if (strcmp(type, "Tcl") == 0) {
        opserr << "ERROR: Tcl function evaluator not implemented" << endln;
        return -1;
    } else if (strcmp(type, "Python") == 0) {
#ifdef _PYTHON3
        if (filename == 0) {
            theEval = new PythonEvaluator(cmds->getDomain(),
                                          cmds->getStructuralDomain());
        } else {
            theEval = new PythonEvaluator(
                cmds->getDomain(), cmds->getStructuralDomain(), filename);
        }
#else
        opserr << "ERROR: Python function evaluator not implemented"
               << endln;
        return -1;
#endif
    } else {
        opserr << "ERROR: unrecognized type of function evaluator: "
               << type << endln;
        return -1;
    }

    if (theEval == 0) {
        opserr << "ERROR: could not create function evaluator" << endln;
        return -1;
    } else {
        cmds->setFunctionEvaluator(theEval);
    }

    return 0;
}

int OPS_gradientEvaluator() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to "
                  "gradientEvaluator"
               << endln;
        return -1;
    }

    GradientEvaluator *theEval = 0;

    // Get the type of gradientEvaluator
    const char *type = OPS_GetString();
    if (strcmp(type, "FiniteDifference") == 0) {
        double perturbationFactor = 1000.0;
        // bool doGradientCheck = false;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            int numdata = 1;
            if (strcmp(arg, "-pert") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &perturbationFactor) <
                    0) {
                    opserr << "ERROR: unable to read -pert value for "
                           << type << " gradient evaluator" << endln;
                    return -1;
                }
            }
            if (strcmp(arg, "-check") == 0) {
                // doGradientCheck = true;
            }
        }

        ReliabilityDomain *theRelDomain = cmds->getDomain();
        Domain *theStrDomain = cmds->getStructuralDomain();
        FunctionEvaluator *theEvaluator = cmds->getFunctionEvaluator();
        if (theEvaluator == 0) {
            opserr << "Function evaluator must be defined before "
                      "gradient evaluator"
                   << endln;
            return -1;
        }

        theEval = new FiniteDifferenceGradient(theEvaluator, theRelDomain,
                                               theStrDomain);
    } else if (strcmp(type, "OpenSees") == 0 ||
               strcmp(type, "Implicit") == 0) {
        // bool doGradientCheck = false;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *arg = OPS_GetString();
            if (strcmp(arg, "-check") == 0) {
                // doGradientCheck = true;
            }
        }

        ReliabilityDomain *theRelDomain = cmds->getDomain();
        Domain *theStrDomain = cmds->getStructuralDomain();
        FunctionEvaluator *theEvaluator = cmds->getFunctionEvaluator();
        if (theEvaluator == 0) {
            opserr << "Function evaluator must be defined before "
                      "gradient evaluator"
                   << endln;
            return -1;
        }
        Integrator *sensAlgo = cmds->getSensitivityAlgorithm();
        if (sensAlgo == 0 && strcmp(type,"OpenSees") == 0) {
            opserr << "WARNING: OpenSees integrator must be defined before "
                      "gradient evaluator\n";
            return -1;
        }
        theEval = new ImplicitGradient(theEvaluator, theRelDomain,
                                       theStrDomain, sensAlgo);

    } else {
        opserr << "ERROR: unrecognized type of gradient evaluator: "
               << type << endln;
        return -1;
    }

    if (theEval == 0) {
        opserr << "ERROR: could not create function evaluator" << endln;
        return -1;
    } else {
        if (cmds != 0) cmds->setGradientEvaluator(theEval);
    }

    return 0;
}

int OPS_probabilityTransformation() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "ERROR: wrong number of arguments to "
                  "probabilityTransformation"
               << endln;
        return -1;
    }

    // Get transformation type
    const char *type = OPS_GetString();

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    ProbabilityTransformation *theTransf = 0;

    // Nataf and AllIndependent have the same optional arguments
    int print = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        int numdata = 1;
        if (strcmp(arg, "-print") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numdata, &print) < 0) {
                opserr << "ERROR: unable to read -print value for "
                          "probability transformation"
                       << endln;
                return -1;
            }
        }
    }

    if (strcmp(type, "Nataf") == 0) {
        theTransf = new NatafProbabilityTransformation(
            theReliabilityDomain, print);
    }

    else if (strcmp(type, "AllIndependent") == 0) {
        theTransf =
            new AllIndependentTransformation(theReliabilityDomain, print);
    }

    else {
        opserr << "ERROR: unrecognized type of "
                  "probabilityTransformation "
               << type << endln;
        return -1;
    }

    if (theTransf == 0) {
        opserr << "ERROR: could not create probabilityTransformation"
               << endln;
        return -1;
    } else {
        if (cmds != 0) cmds->setProbabilityTransformation(theTransf);
    }

    return 0;
}

int OPS_transformUtoX() {
    ProbabilityTransformation *theTransf =
        cmds->getProbabilityTransformation();
    if (theTransf == 0) {
        opserr << "ERROR: probability transformation has not been set"
               << endln;
        return -1;
    }

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    int nrv = theReliabilityDomain->getNumberOfRandomVariables();

    if (OPS_GetNumRemainingInputArgs() < nrv) {
        opserr << "ERROR: transformUtoX insufficient # args" << endln;
        return -1;
    }
    if (OPS_GetNumRemainingInputArgs() > nrv &&
        OPS_GetNumRemainingInputArgs() < 2 * nrv) {
        opserr << "ERROR: transformUtoX insufficient # rv tags" << endln;
        return -1;
    }

    int numData = 1;
    double val;
    Vector u(nrv);
    int loc = 0;
    // Read in u-values
    while (loc < nrv && OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &val) < 0) {
            OPS_ResetCurrentInputArg(-1);
            break;
        }
        u(loc) = val;
        loc++;
    }

    ID rvIndex(nrv);
    // Initialize the index to be sequential (default)
    for (int i = 0; i < nrv; i++) rvIndex(i) = i;

    int rvTag;
    loc = 0;
    // Now read in rv tags and get their indices
    while (loc < nrv && OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetIntInput(&numData, &rvTag) < 0) {
            OPS_ResetCurrentInputArg(-1);
            break;
        }
        rvIndex(loc) = theReliabilityDomain->getRandomVariableIndex(rvTag);
        loc++;
    }

    // Map in
    Vector uSorted(nrv);
    for (int i = 0; i < nrv; i++) uSorted(rvIndex(i)) = u(i);

    Vector x(nrv);
    theTransf->transform_u_to_x(uSorted, x);

    // Map out
    Vector xSorted(nrv);
    for (int i = 0; i < nrv; i++) xSorted(i) = x(rvIndex(i));

    if (OPS_SetDoubleOutput(&nrv, &xSorted[0], false) < 0) {
        opserr << "ERROR: failed to set output in transformUtoX" << endln;
        return -1;
    }

    return 0;
}

int OPS_runFOSMAnalysis() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: Wrong number of input parameter to FOSM "
                  "analysis\n";
        return -1;
    }

    // get file name
    const char *filename = OPS_GetString();

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "FOSMAnalysis -- ReliabilityDomain is not defined\n";
        return -1;
    }

    // Check for essential ingredients
    FunctionEvaluator *theFunctionEvaluator = cmds->getFunctionEvaluator();
    if (theFunctionEvaluator == 0) {
        opserr << "Need theGFunEvaluator before a FOSMAnalysis can "
                  "be created\n";
        return -1;
    }

    GradientEvaluator *theGradientEvaluator = cmds->getGradientEvaluator();
    if (theGradientEvaluator == 0) {
        opserr << "Need theGradientEvaluator before a FOSMAnalysis "
                  "can be created\n";
        return -1;
    }

    Domain *theStructuralDomain = cmds->getStructuralDomain();
    if (theStructuralDomain == 0) {
        opserr << "Structural Domain is not defined\n";
        return -1;
    }

    FOSMAnalysis *theFOSMAnalysis = new FOSMAnalysis(
        theReliabilityDomain, theStructuralDomain, theFunctionEvaluator,
        theGradientEvaluator, 0, filename);

    if (theFOSMAnalysis == 0) {
      opserr << "Unable to create FOSM analysis" << endln;
      return -1;
    }

    cmds->setFOSMAnalysis(theFOSMAnalysis);
    
    // Now run the analysis
    if (theFOSMAnalysis->analyze() < 0) {
        opserr << "WARNING: the FOSM analysis failed\n";
        return -1;
    }

    return 0;
}

int OPS_runFORMAnalysis() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: Wrong number of input parameter to FORM "
                  "analysis\n";
        return -1;
    }

    // get file name
    const char *filename = OPS_GetString();

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "FORMAnalysis -- ReliabilityDomain is not defined\n";
        return -1;
    }
    
    // Check for essential ingredients
    FunctionEvaluator *theFunctionEvaluator = cmds->getFunctionEvaluator();
    if (theFunctionEvaluator == 0) {
        opserr << "Need theGFunEvaluator before a FORMAnalysis can "
                  "be created\n";
        return -1;
    }

    FindDesignPointAlgorithm *theFindDesignPointAlgorithm =
        cmds->getFindDesignPointAlgorithm();
    if (theFindDesignPointAlgorithm == 0) {
        opserr << "Need theFindDesignPointAlgorithm before a "
                  "FORMAnalysis "
                  "can be created\n";
        return -1;
    }

    ProbabilityTransformation *theProbabilityTransformation =
        cmds->getProbabilityTransformation();
    if (theProbabilityTransformation == 0) {
      opserr << "FORMAnalysis - probability transformation not defined - ";
      opserr << "assuming all independent transformation" << endln;
	theProbabilityTransformation =
            new AllIndependentTransformation(theReliabilityDomain, 0);
	cmds->setProbabilityTransformation(theProbabilityTransformation);      
    }

    Domain *theStructuralDomain = cmds->getStructuralDomain();
    if (theStructuralDomain == 0) {
        opserr << "Structural Domain is not defined\n";
        return -1;
    }

    // Read input parameter(s)
    int relSensTag = 0;
    if (OPS_GetNumRemainingInputArgs() > 1) {
        const char *type = OPS_GetString();
        if (strcmp(type, "-relSens") == 0) {
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &relSensTag) < 0) {
                opserr << "ERROR: invalid input: relSensTag \n";
                return -1;
            }
        } else {
            opserr << "ERROR: Invalid input to FORMAnalysis.\n";
            return -1;
        }
    }

    // Create the analysis object
    FORMAnalysis *theFORMAnalysis = new FORMAnalysis(
        theReliabilityDomain, theFindDesignPointAlgorithm,
        theFunctionEvaluator, theProbabilityTransformation, filename,
        relSensTag);

    if (theFORMAnalysis == 0) {
      opserr << "Unable to create FORM analysis" << endln;
      return -1;
    }
    
    cmds->setFORMAnalysis(theFORMAnalysis);
    
    // Now run the analysis
    if (theFORMAnalysis->analyze() < 0) {
        opserr << "WARNING: the FORM analysis failed\n";
        return -1;
    }

    return 0;
}

int OPS_runSORMAnalysis() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: Wrong number of input parameter to SORM "
	       << "analysis" << endln;;
        return -1;
    }

    // get file name
    const char *filename = OPS_GetString();

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
      opserr << "FORMAnalysis -- ReliabilityDomain is not defined" << endln;
      return -1;
    }
    
    // Check for essential ingredients
    FindCurvatures *theFindCurvatures = cmds->getFindCurvatures();
    if (theFindCurvatures == 0) {
        opserr << "Need FindCurvature approach before a SORMAnalysis can "
	       << "be created" << endln;
        return -1;
    }

    FORMAnalysis *theFORMAnalysis =
        cmds->getFORMAnalysis();
    if (theFORMAnalysis == 0) {
        opserr << "Need to run a FORM analysis before "
	       << "SORMAnalysis can be created" << endln;
        return -1;
    }

    // Check for essential ingredients
    FunctionEvaluator *theFunctionEvaluator = cmds->getFunctionEvaluator();
    if (theFunctionEvaluator == 0) {
        opserr << "Need theGFunEvaluator before a SORMAnalysis can "
	       << "be created" << endln;
        return -1;
    }
    
    // Create the analysis object
    SORMAnalysis *theSORMAnalysis = new SORMAnalysis(
        theReliabilityDomain, theFunctionEvaluator,
	theFORMAnalysis, theFindCurvatures, filename);

    if (theSORMAnalysis == 0) {
      opserr << "Unable to create SORM analysis" << endln;
      return -1;
    }
    
    cmds->setSORMAnalysis(theSORMAnalysis);
    
    // Now run the analysis
    if (theSORMAnalysis->analyze() < 0) {
        opserr << "WARNING: the SORM analysis failed\n";
        return -1;
    }

    return 0;
}

int OPS_runImportanceSamplingAnalysis() {
    if (cmds == 0) {
        return -1;
    }

    // filename
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: need filename\n";
        return -1;
    }
    const char *filename = OPS_GetString();

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "Need theReliabilityDomain before a "
                  "ImportanceSamplingAnalysis can be created\n";
        return -1;
    }

    // Check for essential tools
    ProbabilityTransformation *theProbabilityTransformation =
        cmds->getProbabilityTransformation();
    if (theProbabilityTransformation == 0) {
      opserr << "ImportanceSampingAnalysis - probability transformation not defined - ";
      opserr << "assuming all independent transformation" << endln;
	theProbabilityTransformation =
            new AllIndependentTransformation(theReliabilityDomain, 0);
	cmds->setProbabilityTransformation(theProbabilityTransformation);
    }
    FunctionEvaluator *theFunctionEvaluator = cmds->getFunctionEvaluator();
    if (theFunctionEvaluator == 0) {
        opserr << "Need theGFunEvaluator before a "
                  "ImportanceSamplingAnalysis "
                  "can be created\n";
        return -1;
    }
    RandomNumberGenerator *theRandomNumberGenerator =
        cmds->getRandomNumberGenerator();
    if (theRandomNumberGenerator == 0) {
      // Don't really need to say this - just do it - bc CStdLib is the only option -- MHS
      //opserr << "ImportanceSampingAnalysis - random number generator not defined - ";
      //opserr << "assuming CStdLib generator" << endln;      
	theRandomNumberGenerator = new CStdLibRandGenerator();
	cmds->setRandomNumberGenerator(theRandomNumberGenerator);
    }

    Domain *theStructuralDomain = cmds->getStructuralDomain();
    if (theStructuralDomain == 0) {
        opserr << "Need theStructuralDomain before a "
                  "ImportanceSamplingAnalysis can be created\n";
        return -1;
    }

    // The following switches are available (default values are
    // provided) (The sampling is performed around theStartPoint,
    // except for response statistics sampling; then the mean is
    // used together with unit sampling variance.)
    //
    //     -type  failureProbability (1)......... this is the
    //     default -type  responseStatistics (2) -type  saveGvalues
    //     (3)
    //
    //     -variance 1.0  ....................... this is the
    //     default
    //
    //     -maxNum 1000  ........................ this is the
    //     default
    //
    //     -targetCOV 0.05  ..................... this is the
    //     default
    //
    //     -print 0   (print nothing) ........... this is the
    //     default -print 1   (print to screen) -print 2   (print
    //     to restart file)
    //

    // Declaration of input parameters
    long int numberOfSimulations = 1000;
    double targetCOV = 0.05;
    double samplingVariance = 1.0;
    int printFlag = 0;
    int analysisTypeTag = 1;

    while (OPS_GetNumRemainingInputArgs() > 1) {
        const char *type = OPS_GetString();

        if (strcmp(type, "-type") == 0) {
            const char *subtype = OPS_GetString();
            if (strcmp(subtype, "failureProbability") == 0) {
                analysisTypeTag = 1;

            } else if (strcmp(subtype, "outCrossingFailureProbability") ==
                       0) {
                analysisTypeTag = 4;

            } else if ((strcmp(subtype, "responseStatistics") == 0) ||
                       (strcmp(subtype, "saveGvalues") == 0)) {
                if (strcmp(subtype, "responseStatistics") == 0) {
                    analysisTypeTag = 2;
                } else {
                    analysisTypeTag = 3;
                }
                if (samplingVariance != 1.0) {
                    opserr
                        << "ERROR:: sampling variance must be 1.0 for \n"
                        << " response statistics sampling.\n";
                    return -1;
                }

            } else {
                opserr << "ERROR: invalid input: type \n";
                return -1;
            }
        } else if (strcmp(type, "-variance") == 0) {
            // GET INPUT PARAMETER (double)
            int numdata = 1;
            if (OPS_GetDoubleInput(&numdata, &samplingVariance) < 0) {
                opserr << "ERROR: invalid input: samplingVariance \n";
                return -1;
            }
            if (analysisTypeTag == 2 && samplingVariance != 1.0) {
                opserr << "ERROR:: sampling variance must be 1.0 for \n"
                       << " response statistics sampling.\n";
                return -1;
            }
        } else if (strcmp(type, "-maxNum") == 0) {
            int numdata = 1;
            double data = 0.0;
            if (OPS_GetDoubleInput(&numdata, &data) < 0) {
                opserr << "ERROR: invalid input: maxNum \n";
                return -1;
            }
            numberOfSimulations = long(data);

        } else if (strcmp(type, "-targetCOV") == 0) {
            int numdata = 1;
            if (OPS_GetDoubleInput(&numdata, &targetCOV) < 0) {
                opserr << "ERROR: invalid input: targetCOV \n";
                return -1;
            }

        } else if (strcmp(type, "-print") == 0) {
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &printFlag) < 0) {
                opserr << "ERROR: invalid input: printFlag \n";
                return -1;
            }

        } else {
            opserr << "ERROR: invalid input to sampling analysis. \n";
            return -1;
        }
    }

    // Warn about illegal combinations
    if (analysisTypeTag == 2 && printFlag == 2) {
        opserr << "ERROR:: The restart option of the sampling "
                  "analysis cannot be \n"
               << " used together with the response statistics "
                  "option. \n";
        return -1;
    }

    ImportanceSamplingAnalysis *theImportanceSamplingAnalysis
      = new ImportanceSamplingAnalysis(
            theReliabilityDomain, theStructuralDomain,
            theProbabilityTransformation, theFunctionEvaluator,
            theRandomNumberGenerator, 0, numberOfSimulations, targetCOV,
            samplingVariance, printFlag, filename, analysisTypeTag);

    if (theImportanceSamplingAnalysis == 0) {
      opserr << "Unable to create ImportanceSampling analysis" << endln;
      return -1;
    }
    
    cmds->setImportanceSamplingAnalysis(theImportanceSamplingAnalysis);

    if (theImportanceSamplingAnalysis->analyze() < 0) {
        opserr << "WARNING: failed to run ImportanceSamplingAnalysis\n";
        return -1;
    }

    return 0;
}

int OPS_runMonteCarloAnalysis() {
    if (cmds == 0) {
        return -1;
    }

    // filename
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: need filename\n";
        return -1;
    }
    const char *filename = OPS_GetString();

    ReliabilityDomain *theReliabilityDomain = cmds->getDomain();
    if (theReliabilityDomain == 0) {
        opserr << "Need theReliabilityDomain before a "
                  "ImportanceSamplingAnalysis can be created\n";
        return -1;
    }

    // Check for essential tools
    ProbabilityTransformation *theProbabilityTransformation =
        cmds->getProbabilityTransformation();
    if (theProbabilityTransformation == 0) {
      opserr << "ImportanceSampingAnalysis - probability transformation not defined - ";
      opserr << "assuming all independent transformation" << endln;
	theProbabilityTransformation =
            new AllIndependentTransformation(theReliabilityDomain, 0);
	cmds->setProbabilityTransformation(theProbabilityTransformation);
    }
    FunctionEvaluator *theFunctionEvaluator = cmds->getFunctionEvaluator();
    if (theFunctionEvaluator == 0) {
        opserr << "Need theGFunEvaluator before a "
                  "ImportanceSamplingAnalysis "
                  "can be created\n";
        return -1;
    }
    RandomNumberGenerator *theRandomNumberGenerator =
        cmds->getRandomNumberGenerator();
    if (theRandomNumberGenerator == 0) {
      // Don't really need to say this - just do it - bc CStdLib is the only option -- MHS
      //opserr << "ImportanceSampingAnalysis - random number generator not defined - ";
      //opserr << "assuming CStdLib generator" << endln;      
	theRandomNumberGenerator = new CStdLibRandGenerator();
	cmds->setRandomNumberGenerator(theRandomNumberGenerator);
    }

    // Declaration of input parameters
    long int numberOfSimulations = 1000;
    double targetCOV = 0.05;
    int seed = 0;
    int printFlag = 0;
    char outputFile[80] = "";
    char *tclFileName = 0;

    while (OPS_GetNumRemainingInputArgs() > 1) {
        const char *type = OPS_GetString();

        if (strcmp(type, "-maxNum") == 0 || strcmp(type,"maxNum") == 0) {
            int numdata = 1;
	    double data = 0.0;
            if (OPS_GetDoubleInput(&numdata, &data) < 0) {
                opserr << "ERROR: invalid input: numberOfSimulations \n";
                return -1;
            }
            numberOfSimulations = long(data);

        } else if (strcmp(type, "-seed") == 0) {
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &seed) < 0) {
                opserr << "ERROR: invalid input: seed \n";
                return -1;
            }

        } else if (strcmp(type, "-print") == 0 || strcmp(type, "-printFlag") == 0) {
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &printFlag) < 0) {
                opserr << "ERROR: invalid input: printFlag \n";
                return -1;
            }

        } else {
            opserr << "ERROR: invalid input to Monte Carlo analysis. \n";
            return -1;
        }
    }

    
    MonteCarloResponseAnalysis *theMonteCarloAnalysis = 0;
    /*
      theMonteCarloAnalysis
      = new MonteCarloResponseAnalysis(theReliabilityDomain, 0,
				       theProbabilityTransformation, //theFunctionEvaluator,
				       theRandomNumberGenerator, numberOfSimulations,
				       printFlag, outputFile, tclFileName, seed);

    if (theMonteCarloAnalysis == 0) {
      opserr << "Unable to create MonteCarlo analysis" << endln;
      return -1;
    }
    
    cmds->setMonteCarloAnalysis(theMonteCarloAnalysis);

    if (theMonteCarloAnalysis->analyze() < 0) {
        opserr << "WARNING: failed to run MonteCarloAnalysis\n";
        return -1;
    }
*/

    return 0;
}
