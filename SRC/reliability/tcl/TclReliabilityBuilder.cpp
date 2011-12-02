/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.12 $
// $Date: 2003-04-28 20:51:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/tcl/TclReliabilityBuilder.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>

#include <Domain.h>

#include <ReliabilityDomain.h>
#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <ParameterPositioner.h>
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
#include <LaplaceRV.h>
#include <ParetoRV.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <ProbabilityTransformation.h>
#include <NatafProbabilityTransformation.h>
#include <FindDesignPointAlgorithm.h>
#include <ReliabilityAnalysis.h>
#include <HLRFSearchDirection.h>
#include <ArmijoStepSizeRule.h>
#include <FixedStepSizeRule.h>
#include <OpenSeesGFunEvaluator.h>
#include <OpenSeesGradGEvaluator.h>
#include <BasicGFunEvaluator.h>
#include <TclGFunEvaluator.h>
#include <FiniteDifferenceGradGEvaluator.h>
#include <SearchWithStepSizeAndStepDirection.h>
#include <FORMAnalysis.h>
#include <MVFOSMAnalysis.h>
#include <FragilityAnalysis.h>
#include <GFunVisualizationAnalysis.h>
#include <OutCrossingAnalysis.h>
#include <SamplingAnalysis.h>
#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>
#include <FindCurvatures.h>
#include <FirstPrincipalCurvature.h>
#include <CurvaturesBySearchAlgorithm.h>
#include <SORMAnalysis.h>
#include <SystemAnalysis.h>
#include <Filter.h>
#include <StandardLinearOscillatorDisplacementFilter.h>
#include <StandardLinearOscillatorVelocityFilter.h>
#include <StandardLinearOscillatorAccelerationFilter.h>
#include <ModulatingFunction.h>
#include <GammaModulatingFunction.h>
#include <ConstantModulatingFunction.h>
#include <TrapezoidalModulatingFunction.h>
#include <Spectrum.h>
#include <JonswapSpectrum.h>
#include <NarrowBandSpectrum.h>
#include <PointsSpectrum.h>
#include <SensitivityAlgorithm.h>
#include <ReliabilityConvergenceCheck.h>
#include <StandardReliabilityConvergenceCheck.h>
#include <OptimalityConditionReliabilityConvergenceCheck.h>
#include <MeritFunctionCheck.h>
#include <AdkZhangMeritFunctionCheck.h>
#include <CriteriaReductionMeritFunctionCheck.h>
#include <PolakHeSearchDirectionAndMeritFunction.h>
#include <SQPsearchDirectionMeritFunctionAndHessian.h>
#include <HessianApproximation.h>
#include <GradientProjectionSearchDirection.h>
#include <RootFinding.h>
#include <SecantRootFinding.h>


#include <TclReliabilityBuilder.h>
extern SensitivityAlgorithm *theSensitivityAlgorithm;

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//
ReliabilityDomain *theReliabilityDomain = 0;
static Domain *theStructuralDomain = 0;

static GFunEvaluator *theGFunEvaluator = 0;
static GradGEvaluator *theGradGEvaluator = 0;
static StepSizeRule *theStepSizeRule = 0;
static SearchDirection *theSearchDirection = 0;
static HessianApproximation *theHessianApproximation = 0;
static MeritFunctionCheck *theMeritFunctionCheck = 0;
static PolakHeSearchDirectionAndMeritFunction *thePolakHeDualPurpose = 0;
static SQPsearchDirectionMeritFunctionAndHessian *theSQPtriplePurpose = 0;
static ProbabilityTransformation *theProbabilityTransformation = 0;
static ReliabilityConvergenceCheck *theReliabilityConvergenceCheck = 0;
static Vector *theStartPoint = 0;
static RootFinding *theRootFindingAlgorithm = 0;
RandomNumberGenerator *theRandomNumberGenerator = 0;
static FindDesignPointAlgorithm *theFindDesignPointAlgorithm = 0;
static FindCurvatures *theFindCurvatures = 0;
static GFunVisualizationAnalysis *theGFunVisualizationAnalysis = 0;
static FORMAnalysis *theFORMAnalysis = 0;
static MVFOSMAnalysis *theMVFOSMAnalysis = 0;
static FragilityAnalysis *theFragilityAnalysis = 0;
static OutCrossingAnalysis *theOutCrossingAnalysis = 0;
static SORMAnalysis *theSORMAnalysis = 0;
static SamplingAnalysis *theSamplingAnalysis = 0;
static SystemAnalysis *theSystemAnalysis = 0;

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//
int TclReliabilityModelBuilder_addRandomVariable(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addCorrelate(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_correlateGroup(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_correlationStructure(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addLimitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addRandomVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addParameterPositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addModulatingFunction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFilter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addSpectrum(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addProbabilityTransformation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addStartPoint(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addRootFinding(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addRandomNumberGenerator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addSearchDirection(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addHessianApproximation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addMeritFunctionCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addReliabilityConvergenceCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addStepSizeRule(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addgFunEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addGradGEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFindDesignPointAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFindCurvatures(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runFORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runMVFOSMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runFragilityAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runGFunVisualizationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runOutCrossingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSystemAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSamplingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_tempCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_inputCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getMean(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getStdv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);



//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclReliabilityBuilder::TclReliabilityBuilder(Domain &passedDomain, Tcl_Interp *interp)
{
  // Set the interpreter (the destructor needs it to delete commands)
  // Well... not any more. 
  theInterp = interp;

  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "randomVariable",	TclReliabilityModelBuilder_addRandomVariable,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "correlate", TclReliabilityModelBuilder_addCorrelate,(ClientData)NULL, NULL); 
  Tcl_CreateCommand(interp, "correlateGroup", TclReliabilityModelBuilder_correlateGroup,(ClientData)NULL, NULL); 
  Tcl_CreateCommand(interp, "correlationStructure", TclReliabilityModelBuilder_correlationStructure,(ClientData)NULL, NULL); 
  Tcl_CreateCommand(interp, "performanceFunction", TclReliabilityModelBuilder_addLimitState,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "randomVariablePositioner",TclReliabilityModelBuilder_addRandomVariablePositioner,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "parameterPositioner",TclReliabilityModelBuilder_addParameterPositioner,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "modulatingFunction",TclReliabilityModelBuilder_addModulatingFunction,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "filter",TclReliabilityModelBuilder_addFilter,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "spectrum",TclReliabilityModelBuilder_addSpectrum,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "findDesignPoint",	TclReliabilityModelBuilder_addFindDesignPointAlgorithm,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "startPoint",	TclReliabilityModelBuilder_addStartPoint,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "rootFinding",	TclReliabilityModelBuilder_addRootFinding,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "gFunEvaluator",	TclReliabilityModelBuilder_addgFunEvaluator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "gradGEvaluator",TclReliabilityModelBuilder_addGradGEvaluator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "stepSizeRule",TclReliabilityModelBuilder_addStepSizeRule,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "searchDirection",	TclReliabilityModelBuilder_addSearchDirection,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HessianApproximation",	TclReliabilityModelBuilder_addHessianApproximation,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "meritFunctionCheck",	TclReliabilityModelBuilder_addMeritFunctionCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "reliabilityConvergenceCheck",	TclReliabilityModelBuilder_addReliabilityConvergenceCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "probabilityTransformation",	TclReliabilityModelBuilder_addProbabilityTransformation,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "findCurvatures",	TclReliabilityModelBuilder_addFindCurvatures,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "randomNumberGenerator",TclReliabilityModelBuilder_addRandomNumberGenerator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runFORMAnalysis",TclReliabilityModelBuilder_runFORMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runMVFOSMAnalysis",TclReliabilityModelBuilder_runMVFOSMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runFragilityAnalysis",TclReliabilityModelBuilder_runFragilityAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runGFunVizAnalysis",TclReliabilityModelBuilder_runGFunVisualizationAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runOutCrossingAnalysis",TclReliabilityModelBuilder_runOutCrossingAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runSORMAnalysis",TclReliabilityModelBuilder_runSORMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runSystemAnalysis",TclReliabilityModelBuilder_runSystemAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runSamplingAnalysis",TclReliabilityModelBuilder_runSamplingAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "tempCommand",TclReliabilityModelBuilder_tempCommand,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "inputCheck",TclReliabilityModelBuilder_inputCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getMean",TclReliabilityModelBuilder_getMean,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getStdv",TclReliabilityModelBuilder_getStdv,(ClientData)NULL, NULL);


  // set the static pointers in this file
  theStructuralDomain	= &passedDomain;
  theReliabilityDomain	= new ReliabilityDomain();


}

TclReliabilityBuilder::~TclReliabilityBuilder()
{
	// Delete objects
	if (theReliabilityDomain != 0)
		delete theReliabilityDomain;
	if (theGFunEvaluator != 0)
		delete theGFunEvaluator;
	if (theGradGEvaluator != 0)
		delete theGradGEvaluator;
	if (theStepSizeRule != 0)
		delete theStepSizeRule;
	if (theSearchDirection != 0)
		delete theSearchDirection;
	if (theHessianApproximation != 0)
		delete theHessianApproximation;
	if (thePolakHeDualPurpose != 0)
		delete thePolakHeDualPurpose;
	if (theSQPtriplePurpose != 0)
		delete theSQPtriplePurpose;
	if (theMeritFunctionCheck != 0)
		delete theMeritFunctionCheck;
	if (theReliabilityConvergenceCheck != 0)
		delete theReliabilityConvergenceCheck;
	if (theProbabilityTransformation != 0)
		delete theProbabilityTransformation;
	if (theStartPoint != 0)
		delete theStartPoint;
	if (theRootFindingAlgorithm != 0)
		delete theRootFindingAlgorithm;
	if (theRandomNumberGenerator != 0)
		delete theRandomNumberGenerator;
	if (theFindDesignPointAlgorithm != 0)
		delete theFindDesignPointAlgorithm;
	if (theFindCurvatures != 0)
		delete theFindCurvatures;
	if (theFORMAnalysis != 0)
		delete theFORMAnalysis;
	if (theMVFOSMAnalysis != 0)
		delete theMVFOSMAnalysis;
	if (theFragilityAnalysis != 0)
		delete theFragilityAnalysis;
	if (theSORMAnalysis != 0)
		delete theSORMAnalysis;
	if (theSamplingAnalysis != 0)
		delete theSamplingAnalysis;
	if (theSystemAnalysis != 0)
		delete theSystemAnalysis;

	// Delete commands
	Tcl_DeleteCommand(theInterp, "randomVariable");
	Tcl_DeleteCommand(theInterp, "correlate");
	Tcl_DeleteCommand(theInterp, "correlateGroup");
	Tcl_DeleteCommand(theInterp, "correlationStructure");
	Tcl_DeleteCommand(theInterp, "limitState");
	Tcl_DeleteCommand(theInterp, "randomVariablePositioner");
	Tcl_DeleteCommand(theInterp, "positionerPositioner");
	Tcl_DeleteCommand(theInterp, "modulatingFunction");
	Tcl_DeleteCommand(theInterp, "filter");
	Tcl_DeleteCommand(theInterp, "spectrum");
	Tcl_DeleteCommand(theInterp, "findDesignPoint");
	Tcl_DeleteCommand(theInterp, "gFunEvaluator");
	Tcl_DeleteCommand(theInterp, "GradGEvaluator");
	Tcl_DeleteCommand(theInterp, "stepSizeRule");
	Tcl_DeleteCommand(theInterp, "searchDirection");
	Tcl_DeleteCommand(theInterp, "HessianApproximation");
	Tcl_DeleteCommand(theInterp, "meritFunctionCheck");
	Tcl_DeleteCommand(theInterp, "reliabilityConvergenceCheck");
	Tcl_DeleteCommand(theInterp, "ProbabilityTransformation");
	Tcl_DeleteCommand(theInterp, "startPoint");
	Tcl_DeleteCommand(theInterp, "rootFinding");
	Tcl_DeleteCommand(theInterp, "findCurvatures");
	Tcl_DeleteCommand(theInterp, "randomNumberGenerator");
	Tcl_DeleteCommand(theInterp, "runFORMAnalysis");
	Tcl_DeleteCommand(theInterp, "runMVFOSMAnalysis");
	Tcl_DeleteCommand(theInterp, "runFragilityAnalysis");
	Tcl_DeleteCommand(theInterp, "runGFunVizAnalysis");
	Tcl_DeleteCommand(theInterp, "runSORMAnalysis");
	Tcl_DeleteCommand(theInterp, "runSystemAnalysis");
	Tcl_DeleteCommand(theInterp, "runSamplingAnalysis");
	Tcl_DeleteCommand(theInterp, "tempCommand");
	Tcl_DeleteCommand(theInterp, "inputCheck");
	Tcl_DeleteCommand(theInterp, "getMean");
	Tcl_DeleteCommand(theInterp, "getStdv");
}


//
// CLASS METHODS
//


ReliabilityDomain *
TclReliabilityBuilder::getReliabilityDomain()
{
	return theReliabilityDomain;
}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addRandomVariable(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{
  RandomVariable *theRandomVariable = 0;
  int tag;
  double mean;
  double stdv;
  double startPt;
  double parameter1;
  double parameter2;
  double parameter3;
  double parameter4;
  int numberOfArguments = argc;


  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 5) {
		opserr << "ERROR: invalid number of arguments to randomVariable command \n";
		return TCL_ERROR;
  }


  // CHECK THAT THE USER HAS PROVIDED A TYPE
  if ((strcmp(argv[2],"beta")				== 0) ||
	  (strcmp(argv[2],"chiSquare")			== 0) ||
	  (strcmp(argv[2],"exponential")		== 0) ||
	  (strcmp(argv[2],"gamma")				== 0) ||
	  (strcmp(argv[2],"gumbel")				== 0) ||
	  (strcmp(argv[2],"laplace")			== 0) ||
	  (strcmp(argv[2],"lognormal")			== 0) ||
	  (strcmp(argv[2],"normal")				== 0) ||
      (strcmp(argv[2],"pareto")				== 0) || 
	  (strcmp(argv[2],"rayleigh")			== 0) ||	  
	  (strcmp(argv[2],"shiftedExponential") == 0) ||
	  (strcmp(argv[2],"shiftedRayleigh")	== 0) ||
	  (strcmp(argv[2],"type1LargestValue")	== 0) ||
	  (strcmp(argv[2],"type1SmallestValue") == 0) ||
	  (strcmp(argv[2],"type2LargestValue")	== 0) ||
	  (strcmp(argv[2],"type3SmallestValue") == 0) ||
	  (strcmp(argv[2],"uniform")			== 0) ||
	  (strcmp(argv[2],"weibull")			== 0) 
	 ) 
  {
  }
  else {
	  opserr << "ERROR: A correct type has not been provided for a random variable." << endln
		  << " (Available types: normal, lognormal, uniform, etc.)" << endln
		  << " Syntax: randomVariable tag? type mean? stdv? <startPt?>" << endln
		  << "     or: randomVariable tag? type par1? par2? par3? par4? <startPt?>" << endln;
	  return TCL_ERROR;
  }


  // GET TAG NUMBER
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }


  // NOW START CREATING THE RANDOM VARIBLE OBJECT
  if (numberOfArguments==5)  {   // (Use mean/stdv WITHOUT startPt)

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[3], &mean) != TCL_OK) {
		opserr << "ERROR: invalid input: mean \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[4], &stdv) != TCL_OK) {
		opserr << "ERROR: invalid input: stdv \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	  if (strcmp(argv[2],"normal") == 0) {
		  if (stdv <= 0.0) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new NormalRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"lognormal") == 0) {
		  if (mean == 0.0 || stdv <= 0.0) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			  theRandomVariable = new LognormalRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"gamma") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			  theRandomVariable = new GammaRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedExponential") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ShiftedExponentialRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedRayleigh") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ShiftedRayleighRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"exponential") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ExponentialRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"rayleigh") == 0) {
		  opserr << "Random variable with tag " << tag << "cannot be created with only mean/stdv." << endln;
		  return TCL_ERROR;
	  }
	  else if (strcmp(argv[2],"uniform") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new UniformRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"beta") == 0) {
		  opserr << "ERROR:: 'Beta' type random variable: use parameters to create!\n";
		  return TCL_ERROR;
	  }
	  else if (strcmp(argv[2],"type1LargestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type1LargestValueRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"type1SmallestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type1SmallestValueRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"type2LargestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type2LargestValueRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"type3SmallestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type3SmallestValueRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"chiSquare") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ChiSquareRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"gumbel") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new GumbelRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"weibull") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new WeibullRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"laplace") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new LaplaceRV(tag, mean, stdv);
		  }
	  }
	  else if (strcmp(argv[2],"pareto") == 0) {
		  opserr << "Random variable with tag " << tag << "cannot be created with only mean/stdv." << endln;
		  return TCL_ERROR;
	  }
	  else {
		opserr << "ERROR: unrecognized type of random variable number " << tag << endln;
		return TCL_ERROR;
	  }

	  if (theRandomVariable == 0) {
		opserr << "ERROR: could not create random variable number " << tag << endln;
		return TCL_ERROR;
	  }
  }
  else if (numberOfArguments==6)  {   // (Use mean/stdv AND startPt)

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[3], &mean) != TCL_OK) {
		opserr << "ERROR: invalid input: mean \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[4], &stdv) != TCL_OK) {
		opserr << "ERROR: invalid input: stdv \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[5], &startPt) != TCL_OK) {
		opserr << "ERROR: invalid input: startPt \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	  if (strcmp(argv[2],"normal") == 0) {
		  if (stdv <= 0.0) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new NormalRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"lognormal") == 0) {
		  if (mean == 0.0 || stdv <= 0.0) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			  theRandomVariable = new LognormalRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"gamma") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			  theRandomVariable = new GammaRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedExponential") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ShiftedExponentialRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedRayleigh") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ShiftedRayleighRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"exponential") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ExponentialRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"rayleigh") == 0) {
		  opserr << "Random variable with tag " << tag << "cannot be created with only mean/stdv." << endln;
	  }
	  else if (strcmp(argv[2],"uniform") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new UniformRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"beta") == 0) {
		  opserr << "ERROR:: 'Beta' type random variable: use parameters to create!\n";
		  return TCL_ERROR;
	  }
	  else if (strcmp(argv[2],"type1LargestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type1LargestValueRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"type1SmallestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type1SmallestValueRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"type2LargestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type2LargestValueRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"type3SmallestValue") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new Type3SmallestValueRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"chiSquare") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new ChiSquareRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"gumbel") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new GumbelRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"weibull") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new WeibullRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"laplace") == 0) {
		  if ( stdv <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else {
			theRandomVariable = new LaplaceRV(tag, mean, stdv, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"pareto") == 0) {
		  opserr << "Random variable with tag " << tag << "cannot be created with only mean/stdv." << endln;
	  }
	  else {
		opserr << "ERROR: unrecognized type of random variable number " << tag << endln;
		return TCL_ERROR;
	  }

	  if (theRandomVariable == 0) {
		opserr << "ERROR: could not create random variable number " << tag << endln;
		return TCL_ERROR;
	  }
  }
  else if (numberOfArguments==7)  {  // (Use parameters WITHOUT startPt)

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[3], &parameter1) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter1 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[4], &parameter2) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter2 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[5], &parameter3) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter3 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[6], &parameter4) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter4 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	  if (strcmp(argv[2],"normal") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new NormalRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"lognormal") == 0) {
//		  if ( parameter2 <= 0.0 ) { Now assume that this indicates negative lognormal 
//			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
//			  return TCL_ERROR;
//		  }
//		  else  {
			theRandomVariable = new LognormalRV(tag, parameter1, parameter2, parameter3, parameter4);
//		  }
	  }
	  else if (strcmp(argv[2],"gamma") == 0) {
		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new GammaRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedExponential") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ShiftedExponentialRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedRayleigh") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ShiftedRayleighRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"exponential") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ExponentialRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"rayleigh") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new RayleighRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"uniform") == 0) {
		  if ( parameter1 >= parameter2 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new UniformRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"beta") == 0) {
		  if ( parameter1 >= parameter2  ||  parameter3 <= 0.0  || parameter4 <= 0.0  ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new BetaRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"type1LargestValue") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type1LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"type1SmallestValue") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type1SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"type2LargestValue") == 0) {
		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type2LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"type3SmallestValue") == 0) {
		  if ( parameter2 <= 0.0  ||  parameter3 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type3SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"chiSquare") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ChiSquareRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"gumbel") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new GumbelRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"weibull") == 0) {
		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new WeibullRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"laplace") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new LaplaceRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else if (strcmp(argv[2],"pareto") == 0) {
		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ParetoRV(tag, parameter1, parameter2, parameter3, parameter4);
		  }
	  }
	  else {
		opserr << "ERROR: unrecognized type of random variable \n";
		opserr << "random variable: " << tag << endln;
		return TCL_ERROR;
	  }

	  if (theRandomVariable == 0) {
		opserr << "ERROR: unrecognized type of random variable number " << tag << endln;
		return TCL_ERROR;
	  }
  }
  else if (numberOfArguments==8)  {  // (Use parameters AND startPt)

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[3], &parameter1) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter1 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[4], &parameter2) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter2 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[5], &parameter3) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter3 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[6], &parameter4) != TCL_OK) {
		opserr << "ERROR: invalid input: parameter4 \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (double)
	  if (Tcl_GetDouble(interp, argv[7], &startPt) != TCL_OK) {
		opserr << "ERROR: invalid input: startPt \n";
		return TCL_ERROR;
	  }

	  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	  if (strcmp(argv[2],"normal") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new NormalRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"lognormal") == 0) {
//		  if ( parameter2 <= 0.0 ) { Now assume that this indicates negative lognormal
//			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
//			  return TCL_ERROR;
//		  }
//		  else  {
			theRandomVariable = new LognormalRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
//		  }
	  }
	  else if (strcmp(argv[2],"gamma") == 0) {
		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new GammaRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedExponential") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ShiftedExponentialRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"shiftedRayleigh") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ShiftedRayleighRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"exponential") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ExponentialRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"rayleigh") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new RayleighRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"uniform") == 0) {
		  if ( parameter1 >= parameter2 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new UniformRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"beta") == 0) {
		  if ( parameter1 >= parameter2  ||  parameter3 <= 0.0  || parameter4 <= 0.0  ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new BetaRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"type1LargestValue") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type1LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"type1SmallestValue") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type1SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"type2LargestValue") == 0) {
		  if ( parameter1 <= 0.0  || parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type2LargestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"type3SmallestValue") == 0) {
		  if ( parameter2 <= 0.0  ||  parameter3 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new Type3SmallestValueRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"chiSquare") == 0) {
		  if ( parameter1 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ChiSquareRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"gumbel") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new GumbelRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"weibull") == 0) {
		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new WeibullRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"laplace") == 0) {
		  if ( parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new LaplaceRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else if (strcmp(argv[2],"pareto") == 0) {
		  if ( parameter1 <= 0.0  ||  parameter2 <= 0.0 ) {
			  opserr << "ERROR: Invalid parameter input to random variable number " << tag << endln;
			  return TCL_ERROR;
		  }
		  else  {
			theRandomVariable = new ParetoRV(tag, parameter1, parameter2, parameter3, parameter4, startPt);
		  }
	  }
	  else {
		opserr << "ERROR: unrecognized type of random variable number " << tag << endln;
		return TCL_ERROR;
	  }

	  if (theRandomVariable == 0) {
		opserr << "ERROR: could not create random variable number " << tag << endln;
		return TCL_ERROR;
	  }
  }


  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addRandomVariable(theRandomVariable) == false) {
	opserr << "ERROR: failed to add random variable to the domain (wrong number of arguments?)\n";
	opserr << "random variable: " << tag << endln;
	delete theRandomVariable; // otherwise memory leak
	return TCL_ERROR;
  }

  return TCL_OK;
}
					   



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_getMean(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int tag;
	RandomVariable *rv;
	if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
		opserr << "ERROR: Invalid random variable number tag to getMean command." << endln;
		return TCL_ERROR;
	}
	rv = theReliabilityDomain->getRandomVariablePtr(tag);
	if (rv == 0) {
		opserr << "ERROR: Invalid tag number to getMean command. " << endln;
		return TCL_ERROR;
	}
	opserr << "Mean of random variable number " << tag << ": " << rv->getMean() << endln;

	return TCL_OK;
}




//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_getStdv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int tag;
	RandomVariable *rv;
	if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
		opserr << "ERROR: Invalid random variable number tag to getStdv command." << endln;
		return TCL_ERROR;
	}
	rv = theReliabilityDomain->getRandomVariablePtr(tag);
	if (rv == 0) {
		opserr << "ERROR: Invalid tag number to getStdv command. " << endln;
		return TCL_ERROR;
	}
	opserr << "Standard deviation of random variable number " << tag << ": " << rv->getStdv() << endln;
	
	return TCL_OK;
}




//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addCorrelate(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (argc != 4) {
		opserr << "ERROR: Wrong number of arguments to correlate command." << endln;
		return TCL_ERROR;
	}
	
  CorrelationCoefficient *theCorrelationCoefficient = 0;
  int tag;
  int rv1;
  int rv2;
  double correlationValue;


  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[1], &rv1) != TCL_OK) {
	opserr << "ERROR: invalid input: rv1 \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[2], &rv2) != TCL_OK) {
	opserr << "ERROR: invalid input: rv2 \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (double)
  if (Tcl_GetDouble(interp, argv[3], &correlationValue) != TCL_OK) {
	opserr << "ERROR: invalid input: correlationValue \n";
	return TCL_ERROR;
  }

  // CREATE THE OBJECT
  tag = theReliabilityDomain->getNumberOfCorrelationCoefficients();
  theCorrelationCoefficient = new CorrelationCoefficient(tag+1, rv1, rv2, correlationValue);

  if (theCorrelationCoefficient == 0) {
	opserr << "ERROR: ran out of memory creating correlation coefficient \n";
	opserr << "correlation coefficient: " << tag << endln;
	return TCL_ERROR;
  }

  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addCorrelationCoefficient(theCorrelationCoefficient) == false) {
	opserr << "ERROR: failed to add correlation coefficient to the domain\n";
	opserr << "correlation coefficient: " << tag << endln;
	delete theCorrelationCoefficient; // otherwise memory leak
	return TCL_ERROR;
  }

  return TCL_OK;
}




//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_correlateGroup(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int firstRV, lastRV;
	double correlationValue;

	// GET INPUT PARAMETER (integer)
	if (Tcl_GetInt(interp, argv[1], &firstRV) != TCL_OK) {
		opserr << "ERROR: invalid input: firstRV \n";
		return TCL_ERROR;
	}

	// GET INPUT PARAMETER (integer)
	if (Tcl_GetInt(interp, argv[2], &lastRV) != TCL_OK) {
		opserr << "ERROR: invalid input: lastRV \n";
		return TCL_ERROR;
	}

	// GET INPUT PARAMETER (double)
	if (Tcl_GetDouble(interp, argv[3], &correlationValue) != TCL_OK) {
		opserr << "ERROR: invalid input: correlationValue \n";
		return TCL_ERROR;
	}

	// Assume that previos corr. coeffs. have been added in order
	char theCorrelateCommand[50];
	for (int i=firstRV; i<=lastRV; i++) {
		for (int j=i+1; j<=lastRV; j++) {
			sprintf(theCorrelateCommand,"correlate %d %d %10.5f",i,j,correlationValue);
			Tcl_Eval(interp, theCorrelateCommand );
		}
	}

	return TCL_OK;
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_correlationStructure(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int firstRV, lastRV, i;
	double theta, correlationValue;
	char theCorrelateCommand[50];

	// GET INPUT PARAMETER (integer)
	if (Tcl_GetInt(interp, argv[2], &firstRV) != TCL_OK) {
		opserr << "ERROR: invalid input: firstRV \n";
		return TCL_ERROR;
	}

	// GET INPUT PARAMETER (integer)
	if (Tcl_GetInt(interp, argv[3], &lastRV) != TCL_OK) {
		opserr << "ERROR: invalid input: lastRV \n";
		return TCL_ERROR;
	}

	// GET INPUT PARAMETER (double)
	if (Tcl_GetDouble(interp, argv[4], &theta) != TCL_OK) {
		opserr << "ERROR: invalid input: theta \n";
		return TCL_ERROR;
	}


	// Create appropriate correlation coefficients
	if (strcmp(argv[1],"homogeneous1") == 0) {
		for (int i=firstRV; i<=lastRV; i++) {
			for (int j=i+1; j<=lastRV; j++) {
				correlationValue = exp(-fabs(i-j)/theta);
				sprintf(theCorrelateCommand,"correlate %d %d %10.5f",i,j,correlationValue);
				Tcl_Eval(interp, theCorrelateCommand );
			}
		}
	}
	else if (strcmp(argv[1],"homogeneous2") == 0) {
		for (int i=firstRV; i<=lastRV; i++) {
			for (int j=i+1; j<=lastRV; j++) {
				correlationValue = exp(-pow((i-j)/theta,2.0));
				sprintf(theCorrelateCommand,"correlate %d %d %10.5f",i,j,correlationValue);
				Tcl_Eval(interp, theCorrelateCommand );
			}
		}
	}
	else if (strcmp(argv[1],"homogeneous3") == 0) {
		for (int i=firstRV; i<=lastRV; i++) {
			for (int j=i+1; j<=lastRV; j++) {
				correlationValue = 1.0/(1.0*theta*pow(i-j,2.0));
				sprintf(theCorrelateCommand,"correlate %d %d %10.5f",i,j,correlationValue);
				Tcl_Eval(interp, theCorrelateCommand );
			}
		}
	}
	else if (strcmp(argv[1],"homogeneous4") == 0) {
		for (int i=firstRV; i<=lastRV; i++) {
			for (int j=i+1; j<=lastRV; j++) {
				if (fabs(i-j)<theta) {
					correlationValue = 1.0-(fabs(i-j)/theta);
					sprintf(theCorrelateCommand,"correlate %d %d %10.5f",i,j,correlationValue);
					Tcl_Eval(interp, theCorrelateCommand );
				}
			}
		}
	}
	else if (strcmp(argv[1],"vectorProduct") == 0) {

		// Open file where the vector is given
		ifstream inputFile( "correlationVector.txt", ios::in );
		if (inputFile.fail()) {
			opserr << "File correlationVector.txt could not be opened. " << endln;
			return TCL_ERROR;
		}

		// Loop through file to see how many entries there are
		double dummy;
		int numEntries = 0;
		while (inputFile >> dummy) {
			numEntries++;
		}
		if (numEntries == 0) {
			opserr << "ERROR: No entries in the correlationVector.txt file!" << endln;
			return TCL_ERROR;
		}

		// Give a warning if the number of elements of the vector is
		// different from the number of random variables being correlated
		if (numEntries != lastRV-firstRV+1) {
			opserr << "WARNING: The number of entries in the correlationVector.txt file " << endln
				<< " is not equal to the number of random variables that are being correlated." << endln;
		}

		// Close the file
		inputFile.close();

		// Open it again, now being ready to store the results in a vector
		ifstream inputFile2( "correlationVector.txt", ios::in );
		if (inputFile.fail()) {
			opserr << "File correlationVector.txt could not be opened. " << endln;
			return TCL_ERROR;
		}

		// Store the vector
		Vector theVector(numEntries);
		for (i=0; i<numEntries; i++) {
			inputFile2 >> theVector(i);
		}
		inputFile2.close();

		for (int i=firstRV; i<=lastRV; i++) {
			for (int j=i+1; j<=lastRV; j++) {
				correlationValue = theta * theVector(i-firstRV+1) * theVector(j-firstRV+1);
				sprintf(theCorrelateCommand,"correlate %d %d %10.5f",i,j,correlationValue);
				Tcl_Eval(interp, theCorrelateCommand );
			}
		}
	}
	else {
		opserr << "ERROR: Invalid type of correlation structure. " << endln;
		return TCL_ERROR;
	}

	return TCL_OK;
}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addLimitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  LimitStateFunction *theLimitStateFunction = 0;
  int tag;

  	if (theGFunEvaluator != 0 ) {
		opserr << "ERROR: A limit-state function should not be created after the GFunEvaluator has been instantiated." << endln;
		return TCL_ERROR;
	}

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }
  
  // CREATE THE OBJECT (passing on argv[2])
  theLimitStateFunction = new LimitStateFunction(tag, argv[2]);
  if (theLimitStateFunction == 0) {
	opserr << "ERROR: ran out of memory creating limit-state function \n";
	opserr << "limit-state function: " << tag << endln;
	return TCL_ERROR;
  }

  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addLimitStateFunction(theLimitStateFunction) == false) {
	opserr << "ERROR: failed to add limit-state function to the domain\n";
	opserr << "limit-state function: " << tag << endln;
	delete theLimitStateFunction; // otherwise memory leak
	return TCL_ERROR;
  }

  return TCL_OK;




/*

  LimitStateFunction *theLimitStateFunction = 0;
  int tag;
  int node;
  int dof;
  double displacementLimit;

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[3], &node) != TCL_OK) {
	opserr << "ERROR: invalid input: node \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[4], &dof) != TCL_OK) {
	opserr << "ERROR: invalid input: dof \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (double)
  if (Tcl_GetDouble(interp, argv[5], &displacementLimit) != TCL_OK) {
	opserr << "ERROR: invalid input: displacementLimit \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
  if (strcmp(argv[1],"disp") == 0) {
	  theLimitStateFunction = new LimitStateFunction(tag, node, dof, displacementLimit);
  }
  else {
	opserr << "ERROR: unrecognized type of limit-state function \n";
	opserr << "limit-state function: " << tag << endln;
  }

  if (theLimitStateFunction == 0) {
	opserr << "ERROR: ran out of memory creating limit-state function \n";
	opserr << "limit-state function: " << tag << endln;
	return TCL_ERROR;
  }

  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addLimitStateFunction(theLimitStateFunction) == false) {
	opserr << "ERROR: failed to add limit-state function to the domain\n";
	opserr << "limit-state function: " << tag << endln;
	delete theLimitStateFunction; // otherwise memory leak
	return TCL_ERROR;
  }

  return TCL_OK;
*/
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addRandomVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	RandomVariablePositioner *theRandomVariablePositioner = 0;
	int tag;
	int rvNumber;
	int tagOfObject;
	DomainComponent *theObject;
	int argvCounter = 1;


	// READ THE TAG NUMBER
	if (Tcl_GetInt(interp, argv[argvCounter++], &tag) != TCL_OK) {
		opserr << "ERROR: Invalid input tag to random variable positioner." << endln;
		return TCL_ERROR;
	}


	// CHECK IF THE USER WANTS TO CREATE THE RANDOM VARIABLE HERE
	if (strcmp(argv[argvCounter],"-createRV3") == 0) {
		argvCounter++;

		char theTclCommand[100];
		TCL_Char *rvType;
		double mean,stdv;

		rvNumber = tag;
		rvType = argv[argvCounter++];

		// READ MEAN
		if (Tcl_GetDouble(interp, argv[argvCounter], &mean) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv mean \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ STDV
		if (Tcl_GetDouble(interp, argv[argvCounter], &stdv) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv stdv \n";
			return TCL_ERROR;
		}
		argvCounter++;

		// LET TCL CREATE THE RANDOM VARIABLE
		sprintf(theTclCommand,"randomVariable %d %s %15.10e %15.10e",rvNumber,rvType,mean,stdv);
		Tcl_Eval( interp, theTclCommand );

	}
	else if (strcmp(argv[argvCounter],"-createRV4") == 0) {
		argvCounter++;

		char theTclCommand[100];
		TCL_Char *rvType;
		double mean,stdv,startPt;

		rvNumber = tag;
		rvType = argv[argvCounter++];

		// READ MEAN
		if (Tcl_GetDouble(interp, argv[argvCounter], &mean) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv mean \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ STDV
		if (Tcl_GetDouble(interp, argv[argvCounter], &stdv) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv stdv \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ STARTVALUE
		if (Tcl_GetDouble(interp, argv[argvCounter], &startPt) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv startPt \n";
			return TCL_ERROR;
		}
		argvCounter++;

		// LET TCL CREATE THE RANDOM VARIABLE
		sprintf(theTclCommand,"randomVariable %d %s %15.10e %15.10e %15.10e",rvNumber,rvType,mean,stdv,startPt);
		Tcl_Eval( interp, theTclCommand );
	}
	else if (strcmp(argv[argvCounter],"-createRV5") == 0) {
		argvCounter++;

		char theTclCommand[100];
		TCL_Char *rvType;
		double par1, par2, par3, par4;

		rvNumber = tag;
		rvType = argv[argvCounter++];

		// READ PARAMETER 1
		if (Tcl_GetDouble(interp, argv[argvCounter], &par1) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 1 \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ PARAMETER 2
		if (Tcl_GetDouble(interp, argv[argvCounter], &par2) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 2 \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ PARAMETER 3
		if (Tcl_GetDouble(interp, argv[argvCounter], &par3) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 3 \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ PARAMETER 4
		if (Tcl_GetDouble(interp, argv[argvCounter], &par4) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 4 \n";
			return TCL_ERROR;
		}
		argvCounter++;

		// LET TCL CREATE THE RANDOM VARIABLE
		sprintf(theTclCommand,"randomVariable %d %s %15.10e %15.10e %15.10e %15.10e",rvNumber,rvType,par1,par2,par3,par4);
		Tcl_Eval( interp, theTclCommand );
	}
	else if (strcmp(argv[argvCounter],"-createRV6") == 0) {
		argvCounter++;

		char theTclCommand[100];
		TCL_Char *rvType;
		double par1, par2, par3, par4, startPt;

		rvNumber = tag;
		rvType = argv[argvCounter++];

		// READ PARAMETER 1
		if (Tcl_GetDouble(interp, argv[argvCounter], &par1) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 1 \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ PARAMETER 2
		if (Tcl_GetDouble(interp, argv[argvCounter], &par2) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 2 \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ PARAMETER 3
		if (Tcl_GetDouble(interp, argv[argvCounter], &par3) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 3 \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ PARAMETER 4
		if (Tcl_GetDouble(interp, argv[argvCounter], &par4) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv parameter 4 \n";
			return TCL_ERROR;
		}
		argvCounter++;
		// READ START VALUE
		if (Tcl_GetDouble(interp, argv[argvCounter], &startPt) != TCL_OK) {
			opserr << "ERROR: invalid input in positioner: rv startPt \n";
			return TCL_ERROR;
		}
		argvCounter++;

		// LET TCL CREATE THE RANDOM VARIABLE
		sprintf(theTclCommand,"randomVariable %d %s %15.10e %15.10e %15.10e %15.10e %15.10e",rvNumber,rvType,par1,par2,par3,par4,startPt);
		Tcl_Eval( interp, theTclCommand );
	}
	else if (strcmp(argv[argvCounter],"-rvNum") == 0) {
		argvCounter++;
		
		// READ THE RANDOM VARIABLE NUMBER
		if (Tcl_GetInt(interp, argv[argvCounter++], &rvNumber) != TCL_OK) {
			opserr << "ERROR: invalid input: rvNumber \n";
			return TCL_ERROR;
		}

		// CHECK THAT THE RANDOM VARIABLE ACTUALLY EXISTS
		RandomVariable *theRandomVariable = 0;
		theRandomVariable = theReliabilityDomain->getRandomVariablePtr(rvNumber);
		if (theRandomVariable == 0){
			opserr << "ERROR:: A non-existing random variable number " << rvNumber << " is being positioned in the model " << endln;
			return TCL_ERROR;
		}
	}
	else {
		opserr << "ERROR: Illegal random variable specification in random " << endln
			<< " variable positioner command. " << endln;
		return TCL_ERROR;
	}
	

	const char **data = new const char *[argc-argvCounter-2];
	int ii,jj;
	for (ii=argvCounter+2, jj=0; ii<argc; ii++, jj++)
	  data[jj] = argv[ii];

	// IF UNCERTAIN *ELEMENT* PROPERTY
	if (strcmp(argv[argvCounter],"-element") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			argvCounter++;
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}

		theObject = (DomainComponent *)theStructuralDomain->getElement(tagOfObject);

		theRandomVariablePositioner = new RandomVariablePositioner(tag,
									   rvNumber,
									   theObject,
									   data,
									   argc-argvCounter);

		int rvnumber = theRandomVariablePositioner->getRvNumber();
	}

	// IF UNCERTAIN *LOAD*
	else if (strcmp(argv[argvCounter],"-loadPattern") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getLoadPattern(tagOfObject);


//		if (argc > 8) {
//
//			// GET INPUT PARAMETER (double)
//			double factor = 1.0;
//			if (Tcl_GetDouble(interp, argv[8], &factor) != TCL_OK) {
//				opserr << "ERROR: invalid input: factor \n";
//				return TCL_ERROR;
//			}
//			theRandomVariablePositioner = new RandomVariablePositioner(tag,rvNumber,theObject,&argv[5],argc-5,factor);
//		}
//		else {
//			theRandomVariablePositioner = new RandomVariablePositioner(tag,rvNumber,theObject,&argv[5],argc-5);
//		}
		theRandomVariablePositioner = new RandomVariablePositioner(tag,
									   rvNumber,
									   theObject,
									   data,
									   argc-argvCounter);
	}

	// IF UNCERTAIN *NODE* PROPERTY
	else if (strcmp(argv[argvCounter],"-node") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getNode(tagOfObject);

		theRandomVariablePositioner = new RandomVariablePositioner(tag,
									   rvNumber,
									   theObject,
									   data,
									   argc-argvCounter);
	}
	else {
		opserr << "ERROR: Unknown parameter in randomVariablePositioner" << endln;
		return TCL_ERROR;
	}

	delete [] data;

	// ADD THE RANDOMVARIABLEPOSITIONER TO THE DOMAIN
	if (theReliabilityDomain->addRandomVariablePositioner(theRandomVariablePositioner) == false) {
		opserr << "ERROR: failed to add random variable positioner number " << tag << " to the domain." << endln;
		delete theRandomVariablePositioner; // otherwise memory leak
		return TCL_ERROR;
	}

	return TCL_OK;







/*
// THE OLD VERSION OF THIS FUNCTION:

  RandomVariablePositioner *theRandomVariablePositioner = 0;
  int tag;
  int rvNumber;
  int typeOfObject;
  int tagOfObject;
  int typeOfParameterInObject;

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[2], &rvNumber) != TCL_OK) {
	opserr << "ERROR: invalid input: rvNumber \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[3], &typeOfObject) != TCL_OK) {
	opserr << "ERROR: invalid input: typeOfObject \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[4], &tagOfObject) != TCL_OK) {
	opserr << "ERROR: invalid input: tagOfObject \n";
	return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[5], &typeOfParameterInObject) != TCL_OK) {
	opserr << "ERROR: invalid input: typeOfParameterInObject \n";
	return TCL_ERROR;
  }


  // CREATE THE OBJECT
  theRandomVariablePositioner = new RandomVariablePositioner(tag, rvNumber, typeOfObject, tagOfObject, typeOfParameterInObject);

  if (theRandomVariablePositioner == 0) {
	opserr << "ERROR: ran out of memory creating random variable identificator \n";
	opserr << "randomVariableID: " << tag << endln;
	return TCL_ERROR;
  }

  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addRandomVariablePositioner(theRandomVariablePositioner) == false) {
	opserr << "ERROR: failed to add random variable identificator to the domain\n";
	opserr << "randomvariableID: " << tag << endln;
	delete theRandomVariablePositioner; // otherwise memory leak
	return TCL_ERROR;
  }

  return TCL_OK;
*/
}





//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addParameterPositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	ParameterPositioner *theParameterPositioner = 0;
	int tag;
	int tagOfObject;
	DomainComponent *theObject;
	int argvCounter = 1;


	// READ THE TAG NUMBER
	if (Tcl_GetInt(interp, argv[argvCounter++], &tag) != TCL_OK) {
		opserr << "ERROR: Invalid tag given to parameterPositioner. " << endln;
		return TCL_ERROR;
	}


	// IF UNCERTAIN *LOAD*
	if (strcmp(argv[argvCounter],"-loadPattern") == 0) {
		argvCounter++;
		
		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
		  opserr << "ERROR: invalid input: tagOfObject \n";
		  return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getLoadPattern(tagOfObject);

		const char **data = new const char *[argc-argvCounter];
		int ii,jj;
		for (ii=argvCounter, jj=0; ii<argc; ii++, jj++)
		  data[jj] = argv[ii];		
		

		theParameterPositioner = new ParameterPositioner(tag,
								 theObject,
								 data,
								 argc-argvCounter);
		
		delete [] data;
	}
	else {
		opserr << "ERROR: Unknown parameter in parameterPositioner" << endln;
		return TCL_ERROR;
	}


	// ADD THE PARAMETERPOSITIONER TO THE DOMAIN
	if (theReliabilityDomain->addParameterPositioner(theParameterPositioner) == false) {
		opserr << "ERROR: failed to add parameter positioner number " << tag << " to the domain." << endln;
		delete theParameterPositioner; // otherwise memory leak
		return TCL_ERROR;
	}

	return TCL_OK;

}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addModulatingFunction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	ModulatingFunction *theModulatingFunction = 0;

	if (strcmp(argv[2],"gamma") == 0) {

		int thisTag, filterTag;
		double b,c;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[1], &thisTag) != TCL_OK) {
			opserr << "ERROR: invalid input: tag \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[3], &filterTag) != TCL_OK) {
			opserr << "ERROR: invalid input: filterTag \n";
			return TCL_ERROR;
		}

		// CHECK THAT THE FILTER EXISTS
		Filter *theFilter = 0;
		theFilter = theReliabilityDomain->getFilter(filterTag);
		if (theFilter == 0) {
			opserr << "ERROR: Could not find the filter with tag " << filterTag << endln;
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[4], &b) != TCL_OK) {
			opserr << "ERROR: invalid input: b \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[5], &c) != TCL_OK) {
			opserr << "ERROR: invalid input: c \n";
			return TCL_ERROR;
		}

		// CREATE THE OBJECT
		theModulatingFunction = new GammaModulatingFunction(thisTag,theFilter,b,c);

		if (theModulatingFunction == 0) {
			opserr << "ERROR: ran out of memory creating modulating function \n";
			opserr << "modulating function: " << thisTag << endln;
			return TCL_ERROR;
		}

		// ADD THE OBJECT TO THE DOMAIN
		if (theReliabilityDomain->addModulatingFunction(theModulatingFunction) == false) {
			opserr << "ERROR: failed to add modulating function to the domain\n";
			opserr << "modulating function: " << thisTag << endln;
			delete theModulatingFunction; // otherwise memory leak
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[2],"constant") == 0) {
	
		int thisTag, filterTag;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[1], &thisTag) != TCL_OK) {
			opserr << "ERROR: invalid input: tag \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[3], &filterTag) != TCL_OK) {
			opserr << "ERROR: invalid input: filterTag \n";
			return TCL_ERROR;
		}

		// CHECK THAT THE FILTER EXISTS
		Filter *theFilter = 0;
		theFilter = theReliabilityDomain->getFilter(filterTag);
		if (theFilter == 0) {
			opserr << "ERROR: Could not find the filter with tag " << filterTag << endln;
			return TCL_ERROR;
		}

		// CREATE THE OBJECT
		theModulatingFunction = new ConstantModulatingFunction(thisTag,theFilter);

		if (theModulatingFunction == 0) {
			opserr << "ERROR: ran out of memory creating modulating function \n";
			opserr << "modulating function: " << thisTag << endln;
			return TCL_ERROR;
		}

		// ADD THE OBJECT TO THE DOMAIN
		if (theReliabilityDomain->addModulatingFunction(theModulatingFunction) == false) {
			opserr << "ERROR: failed to add modulating function to the domain\n";
			opserr << "modulating function: " << thisTag << endln;
			delete theModulatingFunction; // otherwise memory leak
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[2],"trapezoidal") == 0) {
	
		int thisTag, filterTag;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[1], &thisTag) != TCL_OK) {
			opserr << "ERROR: invalid input: tag \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[3], &filterTag) != TCL_OK) {
			opserr << "ERROR: invalid input: filterTag \n";
			return TCL_ERROR;
		}

		// CHECK THAT THE FILTER EXISTS
		Filter *theFilter = 0;
		theFilter = theReliabilityDomain->getFilter(filterTag);
		if (theFilter == 0) {
			opserr << "ERROR: Could not find the filter with tag " << filterTag << endln;
			return TCL_ERROR;
		}

		double t1, t2, t3, t4;

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[4], &t1) != TCL_OK) {
			opserr << "ERROR: invalid input to modulating function: t1 \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[5], &t2) != TCL_OK) {
			opserr << "ERROR: invalid input to modulating function: t2 \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[6], &t3) != TCL_OK) {
			opserr << "ERROR: invalid input to modulating function: t3 \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[7], &t4) != TCL_OK) {
			opserr << "ERROR: invalid input to modulating function: t4 \n";
			return TCL_ERROR;
		}

		// CREATE THE OBJECT
		theModulatingFunction = new TrapezoidalModulatingFunction(thisTag,theFilter,t1,t2,t3,t4);

		if (theModulatingFunction == 0) {
			opserr << "ERROR: ran out of memory creating modulating function \n";
			opserr << "modulating function: " << thisTag << endln;
			return TCL_ERROR;
		}

		// ADD THE OBJECT TO THE DOMAIN
		if (theReliabilityDomain->addModulatingFunction(theModulatingFunction) == false) {
			opserr << "ERROR: failed to add modulating function to the domain\n";
			opserr << "modulating function: " << thisTag << endln;
			delete theModulatingFunction; // otherwise memory leak
			return TCL_ERROR;
		}
	}
	else {
		opserr << "ERROR:: Unknown type of modulating function. " << endln;
		return TCL_ERROR;
	}

	return TCL_OK;
}
//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addFilter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	Filter *theFilter = 0;

	if (argc != 5) {
		opserr << "ERROR: Wrong number of arguments to filter command." << endln;
		return TCL_ERROR;
	}


	int tag;
	double period_Tn, damping;

	// GET INPUT PARAMETER (integer)
	if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
		opserr << "ERROR: invalid input: tag \n";
		return TCL_ERROR;
	}

	// GET INPUT PARAMETER (double)
	if (Tcl_GetDouble(interp, argv[3], &period_Tn) != TCL_OK) {
		opserr << "ERROR: invalid input: freq_wn \n";
		return TCL_ERROR;
	}

	// GET INPUT PARAMETER (double)
	if (Tcl_GetDouble(interp, argv[4], &damping) != TCL_OK) {
		opserr << "ERROR: invalid input: damping \n";
		return TCL_ERROR;
	}



	if (strcmp(argv[2],"standardDisplacement") == 0) {

		theFilter = new StandardLinearOscillatorDisplacementFilter(tag,period_Tn,damping);
	}
	else if (strcmp(argv[2],"standardVelocity") == 0) {

		theFilter = new StandardLinearOscillatorVelocityFilter(tag,period_Tn,damping);
	}
	else if (strcmp(argv[2],"standardAcceleration") == 0) {

		theFilter = new StandardLinearOscillatorAccelerationFilter(tag,period_Tn,damping);
	}
	else {
		opserr << "ERROR:: Unknown type of filter. " << endln;
		return TCL_ERROR;
	}


	if (theFilter == 0) {
		opserr << "ERROR: ran out of memory creating filter \n";
		opserr << "filter: " << tag << endln;
		return TCL_ERROR;
	}

	// ADD THE OBJECT TO THE DOMAIN
	if (theReliabilityDomain->addFilter(theFilter) == false) {
		opserr << "ERROR: failed to add filter to the domain\n";
		opserr << "filter: " << tag << endln;
		delete theFilter; // otherwise memory leak
		return TCL_ERROR;
	}


	return TCL_OK;
}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addSpectrum(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	Spectrum *theSpectrum = 0;

	if (strcmp(argv[2],"jonswap") == 0) {

		int tag;
		double minFreq, maxFreq, alpha, wp, gamma;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
			opserr << "ERROR: invalid input: tag \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[3], &minFreq) != TCL_OK) {
			opserr << "ERROR: invalid input: minFreq \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[4], &maxFreq) != TCL_OK) {
			opserr << "ERROR: invalid input: maxFreq \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[5], &alpha) != TCL_OK) {
			opserr << "ERROR: invalid input: alpha \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[6], &wp) != TCL_OK) {
			opserr << "ERROR: invalid input: wp \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[7], &gamma) != TCL_OK) {
			opserr << "ERROR: invalid input: gamma \n";
			return TCL_ERROR;
		}

		// CREATE THE OBJECT 
		theSpectrum = new JonswapSpectrum(tag, minFreq, maxFreq, alpha, wp, gamma);

		if (theSpectrum == 0) {
			opserr << "ERROR: ran out of memory creating spectrum \n";
			opserr << "spectrum: " << tag << endln;
			return TCL_ERROR;
		}

		// ADD THE OBJECT TO THE DOMAIN
		if (theReliabilityDomain->addSpectrum(theSpectrum) == false) {
			opserr << "ERROR: failed to add spectrum to the domain\n";
			opserr << "spectrum: " << tag << endln;
			delete theSpectrum; // otherwise memory leak
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[2],"narrowband") == 0) {

		int tag;
		double minFreq, maxFreq, amplitude;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
			opserr << "ERROR: invalid input: tag \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[3], &minFreq) != TCL_OK) {
			opserr << "ERROR: invalid input: minFreq \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[4], &maxFreq) != TCL_OK) {
			opserr << "ERROR: invalid input: maxFreq \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[5], &amplitude) != TCL_OK) {
			opserr << "ERROR: invalid input: amplitude \n";
			return TCL_ERROR;
		}


		// CREATE THE OBJECT 
		theSpectrum = new NarrowBandSpectrum(tag, minFreq, maxFreq, amplitude);

		if (theSpectrum == 0) {
			opserr << "ERROR: ran out of memory creating spectrum \n";
			opserr << "spectrum: " << tag << endln;
			return TCL_ERROR;
		}

		// ADD THE OBJECT TO THE DOMAIN
		if (theReliabilityDomain->addSpectrum(theSpectrum) == false) {
			opserr << "ERROR: failed to add spectrum to the domain\n";
			opserr << "spectrum: " << tag << endln;
			delete theSpectrum; // otherwise memory leak
			return TCL_ERROR;
		}
	}
	else if (strcmp(argv[2],"points") == 0) {

		int tag;
		double frequency, amplitude;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
			opserr << "ERROR: invalid input: tag \n";
			return TCL_ERROR;
		}

		if ( fmod((argc-3),2.0) ) {
			opserr << "ERROR: Inconsistent number of points given to spectrum " << tag << endln;
			return TCL_ERROR;
		}

		
		int numPoints = (int)((argc-3)/2.0);

		Vector frequencies(numPoints);
		Vector amplitudes(numPoints);
		for (int iii=1; iii<=numPoints; iii++) {

			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[(iii-1)*2+3], &frequency) != TCL_OK) {
				opserr << "ERROR: invalid input: frequency \n";
				return TCL_ERROR;
			}

			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[(iii-1)*2+4], &amplitude) != TCL_OK) {
				opserr << "ERROR: invalid input: amplitude \n";
				return TCL_ERROR;
			}

			frequencies(iii-1) = frequency;
			amplitudes(iii-1) = amplitude;
		}



		// CREATE THE OBJECT 
		theSpectrum = new PointsSpectrum(tag, frequencies, amplitudes);

		if (theSpectrum == 0) {
			opserr << "ERROR: ran out of memory creating spectrum \n";
			opserr << "spectrum: " << tag << endln;
			return TCL_ERROR;
		}

		// ADD THE OBJECT TO THE DOMAIN
		if (theReliabilityDomain->addSpectrum(theSpectrum) == false) {
			opserr << "ERROR: failed to add spectrum to the domain\n";
			opserr << "spectrum: " << tag << endln;
			delete theSpectrum; // otherwise memory leak
			return TCL_ERROR;
		}
	}
	else {
		opserr << "ERROR:: Unknown type of spectrum. " << endln;
		return TCL_ERROR;
	}

	return TCL_OK;
}


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addRandomNumberGenerator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theRandomNumberGenerator != 0) {
		delete theRandomNumberGenerator;
		theRandomNumberGenerator = 0;
	}


  // GET INPUT PARAMETER (string) AND CREATE THE OBJECT
  if (strcmp(argv[1],"CStdLib") == 0) {
	  theRandomNumberGenerator = new CStdLibRandGenerator();
  }
  else {
	opserr << "ERROR: unrecognized type of RandomNumberGenerator \n";
	return TCL_ERROR;
  }

  if (theRandomNumberGenerator == 0) {
	opserr << "ERROR: could not create theRandomNumberGenerator \n";
	return TCL_ERROR;
  }
  return TCL_OK;
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addProbabilityTransformation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theProbabilityTransformation != 0) {
		delete theProbabilityTransformation;
		theProbabilityTransformation = 0;
	}


	// Check number of arguments
	if (argc!= 2 && argc!= 4) {
		opserr << "ERROR: Wrong number of arguments to probability transformation." << endln;
		return TCL_ERROR;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"Nataf") == 0) {

		int printFlag = 0; 
		
		if (argc > 2) {
			if (strcmp(argv[2],"-print") == 0) {

				if (Tcl_GetInt(interp, argv[3], &printFlag) != TCL_OK) {
					opserr << "ERROR: invalid input: printFlag to Nataf transformation \n";
					return TCL_ERROR;
				}
			}
		}

		theProbabilityTransformation = new NatafProbabilityTransformation(theReliabilityDomain,printFlag);
  }
  else {
	opserr << "ERROR: unrecognized type of ProbabilityTransformation \n";
	return TCL_ERROR;
  }

  if (theProbabilityTransformation == 0) {
	opserr << "ERROR: could not create theProbabilityTransformation \n";
	return TCL_ERROR;
  }
  return TCL_OK;
}




//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addSearchDirection(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theSearchDirection != 0) {
		delete theSearchDirection;
		theSearchDirection = 0;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"iHLRF") == 0) {

		if (argc != 2) {
			opserr << "ERROR: Wrong number of arguments to iHLRF search direction. " << endln;
			return TCL_ERROR;
		}

		theSearchDirection = new HLRFSearchDirection();

		if (theSearchDirection == 0) {
			opserr << "ERROR: could not create theSearchDirection \n";
			return TCL_ERROR;
		}


	}
	else if (strcmp(argv[1],"PolakHe") == 0) {

		if (argc != 2) {
			opserr << "ERROR: Wrong number of arguments to PolakHe search direction. " << endln;
			return TCL_ERROR;
		}
		thePolakHeDualPurpose = new PolakHeSearchDirectionAndMeritFunction();
		theSearchDirection = thePolakHeDualPurpose;
	}
	else if (strcmp(argv[1],"GradientProjection") == 0) {

		if (argc != 2) {
			opserr << "ERROR: Wrong number of arguments to GradientProjection search direction. " << endln;
			return TCL_ERROR;
		}


		// Check that a step size rule has been created
		if (theStepSizeRule == 0 ) {
			opserr << "Need theStepSizeRule before a GradientProjectionSearchDirection can be created" << endln;
			return TCL_ERROR;
		}

		// Check that a transformation has been created
		if (theProbabilityTransformation == 0 ) {
			opserr << "Need theProbabilityTransformation before a GradientProjectionSearchDirection can be created" << endln;
			return TCL_ERROR;
		}

		// Check that a gfun evaluator has been created
		if (theGFunEvaluator == 0 ) {
			opserr << "Need theGFunEvaluator before a GradientProjectionSearchDirection can be created" << endln;
			return TCL_ERROR;
		}

		// Check that a root-finding algorithm has been created
		if (theRootFindingAlgorithm == 0 ) {
			opserr << "Need theRootFindingAlgorithm before a GradientProjectionSearchDirection can be created" << endln;
			return TCL_ERROR;
		}


		theSearchDirection = new GradientProjectionSearchDirection(theStepSizeRule,
																   theProbabilityTransformation,
																   theGFunEvaluator,
																   theRootFindingAlgorithm);
	}
	else if (strcmp(argv[1],"SQP") == 0) {

		double c_bar = 200.0;
		double e_bar = 0.5;

		int argvCounter = 2;
		while (argc > argvCounter) {
			if (strcmp(argv[argvCounter],"-c_bar") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &c_bar) != TCL_OK) {
					opserr << "ERROR: invalid input: c_bar for algorithm" << endln;
					return TCL_ERROR;
				}
				argvCounter++;

			}
			else if (strcmp(argv[argvCounter],"-e_bar") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &e_bar) != TCL_OK) {
					opserr << "ERROR: invalid input: e_bar for SQP algorithm" << endln;
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to SQP algorithm." << endln;
				return TCL_ERROR;
			}
		}


		theSQPtriplePurpose = new SQPsearchDirectionMeritFunctionAndHessian(c_bar,e_bar);
		theSearchDirection = theSQPtriplePurpose;

		// Set default Hessian approximation in case user forgets
		theHessianApproximation = theSQPtriplePurpose;

		// Set the Hessian approximation in the search direction
		theSQPtriplePurpose->setHessianApproximation(theHessianApproximation);

	
	}
	else {
		opserr << "ERROR: unrecognized type of SearchDirection \n";
		return TCL_ERROR;
	}

	if (theSearchDirection == 0) {
		opserr << "ERROR: could not create theSearchDirection \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}






//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addHessianApproximation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theHessianApproximation != 0) {
		delete theHessianApproximation;
		theHessianApproximation = 0;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"SQP_BFGS") == 0) {

		// Check that the SQP search direction is already created
		if (theSQPtriplePurpose == 0 ) {
			opserr << "Need theSQPSearchDirection before a SQP Hessian Approximation can be created" << endln;
			return TCL_ERROR;
		}

		theHessianApproximation = theSQPtriplePurpose;

		// Set the Hessian approximation in the search direction
		// (this needs to be changed for generatlity; new method of search direction)
		theSQPtriplePurpose->setHessianApproximation(theHessianApproximation);

	}
	else {
		opserr << "ERROR: unrecognized type of HessianApproximation \n";
		return TCL_ERROR;
	}

	if (theHessianApproximation == 0) {
		opserr << "ERROR: could not create theHessianApproximation \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}





//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addMeritFunctionCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theMeritFunctionCheck != 0) {
		delete theMeritFunctionCheck;
		theMeritFunctionCheck = 0;
	}

	int argvCounter = 1;

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[argvCounter],"AdkZhang") == 0) {
		argvCounter++;

		double multi = 2.0;
		double add = 10.0;
		double factor = 0.0;
		
		while (argvCounter < argc) {

			if (strcmp(argv[argvCounter],"-multi") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &multi) != TCL_OK) {
					opserr << "ERROR: invalid input: multi \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-add") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &add) != TCL_OK) {
					opserr << "ERROR: invalid input: add \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-factor") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &factor) != TCL_OK) {
					opserr << "ERROR: invalid input: factor \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to AdkZhang merit function check. " << endln;
				return TCL_ERROR;
			}
		}

		// Do a quick input check
		if (multi<1.0 || add<0.0) {
			opserr << "ERROR: Invalid values of multi/add parameters to AdkZhang merit function check." << endln;
			return TCL_ERROR;
		}

		theMeritFunctionCheck = new AdkZhangMeritFunctionCheck(multi,add,factor);
	}
	else if (strcmp(argv[argvCounter],"criteriaReduction") == 0) {
		argvCounter++;

		// Check for necessary ingredients
		if (theReliabilityConvergenceCheck == 0) {
			opserr << "The TerjeMeritFunctionCheck needs a ReliabilityConvergenceCheck!" << endln;
			return TCL_ERROR;
		}

		theMeritFunctionCheck = new CriteriaReductionMeritFunctionCheck(theReliabilityConvergenceCheck);

	}
	else if (strcmp(argv[argvCounter],"PolakHe") == 0) {
		argvCounter++;

		// Check that the PolakHe search direction is already created
		if (thePolakHeDualPurpose == 0 ) {
			opserr << "Need thePolakHeSearchDirection before a PolakHe merit function can be created" << endln;
			return TCL_ERROR;
		}
		double factor = 0.5;
		while (argvCounter < argc) {

			if (strcmp(argv[argvCounter],"-factor") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &factor) != TCL_OK) {
					opserr << "ERROR: invalid input: factor \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to Polak He merit function check. " << endln;
				return TCL_ERROR;
			}
		}

		thePolakHeDualPurpose->setAlpha(factor);
		theMeritFunctionCheck = thePolakHeDualPurpose;

	}
	else if (strcmp(argv[argvCounter],"SQP") == 0) {
		argvCounter++;

		// Check that the SQP search direction is already created
		if (theSQPtriplePurpose == 0 ) {
			opserr << "Need theSQPSearchDirection before a SQP merit function can be created" << endln;
			return TCL_ERROR;
		}
		
		double factor = 0.5;
		while (argvCounter < argc) {

			if (strcmp(argv[argvCounter],"-factor") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &factor) != TCL_OK) {
					opserr << "ERROR: invalid input: factor \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to SQP merit function check. " << endln;
				return TCL_ERROR;
			}
		}

		theSQPtriplePurpose->setAlpha(factor);
		theMeritFunctionCheck = theSQPtriplePurpose;

	}
	else {
		opserr << "ERROR: unrecognized type of MeritFunctionCheck \n";
		return TCL_ERROR;
	}

	if (theMeritFunctionCheck == 0) {
		opserr << "ERROR: could not create theMeritFunctionCheck \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}








//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addReliabilityConvergenceCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theReliabilityConvergenceCheck != 0) {
		delete theReliabilityConvergenceCheck;
		theReliabilityConvergenceCheck = 0;
	}


	if (argc < 2) {
		opserr << "ERROR: Wrong number of arguments to reliability convergence check." << endln;
		return TCL_ERROR;
	}

	// Initial declarations
	int argvCounter = 1;


	if (strcmp(argv[argvCounter],"Standard") == 0) {
		argvCounter++;

		double e1 = 1.0e-3;
		double e2 = 1.0e-3;
		int print=1;

		while (argvCounter < argc) {

			if (strcmp(argv[argvCounter],"-e1") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &e1) != TCL_OK) {
					opserr << "ERROR: invalid input: e1 \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-e2") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &e2) != TCL_OK) {
					opserr << "ERROR: invalid input: e2 \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-print") == 0) {
				argvCounter++;

				if (Tcl_GetInt(interp, argv[argvCounter], &print) != TCL_OK) {
					opserr << "ERROR: invalid input: print \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to standard reliability convergence check. " << endln;
				return TCL_ERROR;
			}
		}
			theReliabilityConvergenceCheck = new StandardReliabilityConvergenceCheck(e1,e2,print);
	}
	else if (strcmp(argv[argvCounter],"OptimalityCondition") == 0) {
		argvCounter++;

		double e1 = 1.0e-3;
		double e2 = 1.0e-3;
		int print = 1;

		while (argvCounter < argc) {

			if (strcmp(argv[argvCounter],"-e1") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &e1) != TCL_OK) {
					opserr << "ERROR: invalid input: e1 \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-e2") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &e2) != TCL_OK) {
					opserr << "ERROR: invalid input: e2 \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-print") == 0) {
				argvCounter++;

				if (Tcl_GetInt(interp, argv[argvCounter], &print) != TCL_OK) {
					opserr << "ERROR: invalid input: print \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to standard reliability convergence check. " << endln;
				return TCL_ERROR;
			}
		}
		theReliabilityConvergenceCheck = new OptimalityConditionReliabilityConvergenceCheck(e1,e2,print);
	}
	else {
		opserr << "ERROR: unrecognized type of ReliabilityConvergenceCheck \n";
		return TCL_ERROR;
	}

	if (theReliabilityConvergenceCheck == 0) {
		opserr << "ERROR: could not create theReliabilityConvergenceCheck \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addStepSizeRule(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theStepSizeRule != 0) {
		delete theStepSizeRule;
		theStepSizeRule = 0;
	}


	// Initial declarations
	int argvCounter = 1;


	if (strcmp(argv[argvCounter],"Armijo") == 0) {
		argvCounter++;

		// Check that the necessary ingredients are present
		if (theGFunEvaluator == 0 ) {
			opserr << "Need theGFunEvaluator before an ArmijoStepSizeRule can be created" << endln;
			return TCL_ERROR;
		}
		if (theProbabilityTransformation == 0 ) {
			opserr << "Need theProbabilityTransformation before a ArmijoStepSizeRule can be created" << endln;
			return TCL_ERROR;
		}
		if (theGradGEvaluator == 0 ) {
			opserr << "Need theGradGEvaluator before a ArmijoStepSizeRule can be created" << endln;
			return TCL_ERROR;
		}
		if (theMeritFunctionCheck == 0 ) {
			opserr << "Need theMeritFunctionCheck before a ArmijoStepSizeRule can be created" << endln;
			return TCL_ERROR;
		}

		// Get input parameters
		double base = 0.5;
		int maxNumReductions = 10;

		double b0 = 1.0;
		int numberOfShortSteps = 2;
		
		double radius = 10.0;
		double surfaceDistance = 0.1;
		double evolution = 0.5;

		int printFlag = 0;
		
		while (argvCounter < argc) {


			if (strcmp(argv[argvCounter],"-print") == 0) {
				argvCounter++;

				if (Tcl_GetInt(interp, argv[argvCounter], &printFlag) != TCL_OK) {
					opserr << "ERROR: invalid input: printFlag \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-maxNum") == 0) {
				argvCounter++;

				if (Tcl_GetInt(interp, argv[argvCounter], &maxNumReductions) != TCL_OK) {
					opserr << "ERROR: invalid input: maxNumReductions \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-base") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &base) != TCL_OK) {
					opserr << "ERROR: invalid input: base \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-initial") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &b0) != TCL_OK) {
					opserr << "ERROR: invalid input: b0 \n";
					return TCL_ERROR;
				}
				argvCounter++;

				if (Tcl_GetInt(interp, argv[argvCounter], &numberOfShortSteps) != TCL_OK) {
					opserr << "ERROR: invalid input: numberOfShortSteps \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-sphere") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &radius) != TCL_OK) {
					opserr << "ERROR: invalid input: radius \n";
					return TCL_ERROR;
				}
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &surfaceDistance) != TCL_OK) {
					opserr << "ERROR: invalid input: surfaceDistance \n";
					return TCL_ERROR;
				}
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &evolution) != TCL_OK) {
					opserr << "ERROR: invalid input: evolution \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to Armijo rule. " << endln;
				return TCL_ERROR;
			}
		}

		theStepSizeRule = new ArmijoStepSizeRule(theGFunEvaluator,
												 theGradGEvaluator,
												 theProbabilityTransformation,
												 theMeritFunctionCheck,
												 theRootFindingAlgorithm, 
												 base,
												 maxNumReductions,
												 b0,
												 numberOfShortSteps,
												 radius,
												 surfaceDistance,
												 evolution,
												 printFlag);


	}
	else if (strcmp(argv[argvCounter],"Fixed") == 0) {
		argvCounter++;

		double stepSize = 1.0;

		while (argvCounter < argc) {

			if (strcmp(argv[argvCounter],"-stepSize") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &stepSize) != TCL_OK) {
					opserr << "ERROR: Invalid step size input to Fixed step size rule." << endln;
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to Fixed step size rule. " << endln;
				return TCL_ERROR;
			}

		}


		theStepSizeRule = new FixedStepSizeRule(stepSize);
	}
	else {
		opserr << "ERROR: unrecognized type of StepSizeRule \n";
		return TCL_ERROR;
	}

	if (theStepSizeRule == 0) {
		opserr << "ERROR: could not create theStepSizeRule \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addgFunEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theGFunEvaluator != 0) {
		delete theGFunEvaluator;
		theGFunEvaluator = 0;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"Matlab") == 0) {
		opserr << "ERROR: The Matlab g-function evaluator is not implemented in this " << endln
			<< " version of your OpenSees executable file. Please contact the " << endln
			<< " developer for more information." << endln;
		return TCL_ERROR;
	}
	else if (strcmp(argv[1],"Tcl") == 0) {

		if (argc != 3) {
			opserr << "ERROR: Wrong number of arguments to Tcl g-function evaluator." << endln;
			return TCL_ERROR;
		}
		theGFunEvaluator = new TclGFunEvaluator(interp, theReliabilityDomain, argv[2]);

	}
	else if (strcmp(argv[1],"OpenSees") == 0) {

		// There are several alternatives for this command:
		// gFunEvaluator  OpenSees  -file <filename>
		// gFunEvaluator  OpenSees  -runToMaxTimeInGFun
		// gFunEvaluator  OpenSees  -analyze <numSteps> <dt(optional)>

		if (argc < 3) {
			opserr << "ERROR: Too few arguments to gFunEvaluator" << endln;
			return TCL_ERROR;
		}

		if (strcmp(argv[2],"-file") == 0) {

			// Try to open the file to make sure it exists
			ifstream inputFile( argv[3], ios::in );
			if (inputFile.fail()) {
				opserr << "File " << *argv[3] << " could not be opened. " << endln;
				return TCL_ERROR;
			}
			inputFile.close();

			theGFunEvaluator = new OpenSeesGFunEvaluator(
				interp, theReliabilityDomain, argv[3]);
		}
		else if (strcmp(argv[2],"-analyze") == 0) {

			// Get number of steps
			int nsteps;
			if (Tcl_GetInt(interp, argv[3], &nsteps) != TCL_OK) {
				opserr << "ERROR: invalid input: numSteps for OpenSees GFunEvaluator \n";
				return TCL_ERROR;
			}
			// Get optional time increment
			double dt = 0.0;
			if (argc == 5) {
				if (Tcl_GetDouble(interp, argv[4], &dt) != TCL_OK) {
					opserr << "ERROR: invalid input: dt for OpenSees GFunEvaluator \n";
					return TCL_ERROR;
				}
			}
			theGFunEvaluator = new OpenSeesGFunEvaluator(interp, theReliabilityDomain, nsteps, dt);

		}
		else {
			opserr << "ERROR: unrecognized parameter in OpenSees GFunEvaluator \n";
			return TCL_ERROR;
		}
	
	}
	else if (strcmp(argv[1],"Basic") == 0) {
		theGFunEvaluator = new BasicGFunEvaluator(interp, theReliabilityDomain);
	}
	else {
		opserr << "ERROR: unrecognized type of GFunEvaluator \n";
		return TCL_ERROR;
	}
	
	if (theGFunEvaluator == 0) {
		opserr << "ERROR: could not create the theGFunEvaluator \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}





//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addGradGEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theGradGEvaluator != 0) {
		delete theGradGEvaluator;
		theGradGEvaluator = 0;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"FiniteDifference") == 0) {

		double perturbationFactor = 1000.0;
		bool doGradientCheck = false;

		// Check that the necessary ingredients are present
		if (theGFunEvaluator == 0 ) {
			opserr << "Need theGFunEvaluator before a FiniteDifferenceGradGEvaluator can be created" << endln;
			return TCL_ERROR;
		}

		// Possibly read perturbation factor
		if (argc>2) {
			int numExtras;
			if (argc==3 || argc==4) {
				numExtras = 1;
			}
			else if (argc==5) {
				numExtras = 2;
			}
			else {
				opserr << "ERROR: Wrong number of arguments to gradGEvaluator. " << endln;
				return TCL_ERROR;
			}

			for (int i=1; i<=numExtras; i++) {

				if (strcmp(argv[2],"-pert") == 0) {

					if (Tcl_GetDouble(interp, argv[3], &perturbationFactor) != TCL_OK) {
						opserr << "ERROR: invalid input: perturbationFactor \n";
						return TCL_ERROR;
					}
				}
				else if (strcmp(argv[2],"-check") == 0) {
					doGradientCheck = true;
				}
				else {
					opserr << "ERROR: Error in input to gradGEvaluator. " << endln;
					return TCL_ERROR;
				}
			}
		}


		theGradGEvaluator = new FiniteDifferenceGradGEvaluator(theGFunEvaluator, theReliabilityDomain,interp, perturbationFactor,doGradientCheck);
	}
	else if (strcmp(argv[1],"OpenSees") == 0) {

		bool doGradientCheck = false;

		if (theSensitivityAlgorithm == 0) {
			opserr << "Need a DDM sensitivity algorithm before a OpenSees sensitivity evaluator can be created" << endln;
			return TCL_ERROR;
		}
		if (theGFunEvaluator == 0) {
			opserr << "Need theGFunEvaluator before a OpenSees sensitivity evaluator can be created" << endln;
			return TCL_ERROR;
		}

		if (argc==2) {
			// Do nothing
		}
		else if (argc==3) {
			if (strcmp(argv[2],"-check") == 0) {
				doGradientCheck = true;
			}
		}
		else {
			opserr << "ERROR: Wrong number of arguments to gradGEvaluator. " << endln;
			return TCL_ERROR;
		}

		theGradGEvaluator = new OpenSeesGradGEvaluator(theGFunEvaluator, interp, theReliabilityDomain,doGradientCheck);
	}
	else {
		opserr << "ERROR: unrecognized type of GradGEvaluator \n";
		return TCL_ERROR;
	}

	if (theGradGEvaluator == 0) {
		opserr << "ERROR: could not create theGradGEvaluator \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}


//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addFindCurvatures(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theFindCurvatures != 0) {
		delete theFindCurvatures;
		theFindCurvatures = 0;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"firstPrincipal") == 0) {

		theFindCurvatures = new FirstPrincipalCurvature();
	}
	else if (strcmp(argv[1],"bySearchAlgorithm") == 0) {

		// Check that the necessary ingredients are present
		if (theFindDesignPointAlgorithm == 0 ) {
			opserr << "Need theFindDesignPointAlgorithm before a CurvaturesBySearchAlgorithm can be created" << endln;
			return TCL_ERROR;
		}

		int numberOfCurvatures;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[2], &numberOfCurvatures) != TCL_OK) {
			opserr << "ERROR: invalid input: numberOfCurvatures \n";
			return TCL_ERROR;
		}

		theFindCurvatures = new CurvaturesBySearchAlgorithm(numberOfCurvatures,theFindDesignPointAlgorithm);
	}
	else {
		opserr << "ERROR: unrecognized type of FindCurvatures \n";
		return TCL_ERROR;
	}

	if (theFindCurvatures == 0) {
		opserr << "ERROR: could not create theFindCurvatures \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addFindDesignPointAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theFindDesignPointAlgorithm != 0) {
		delete theFindDesignPointAlgorithm;
		theFindDesignPointAlgorithm = 0;
	}

	if (argc < 2) {
		opserr << "ERROR: Wrong number of arguments to find design point algorithm." << endln;
		return TCL_ERROR;
	}

	int argvCounter = 1;

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[argvCounter],"StepSearch") == 0) {
		argvCounter++;

		// Check that the necessary ingredients are present
		if (theGFunEvaluator == 0 ) {
			opserr << "Need theGFunEvaluator before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theGradGEvaluator == 0 ) {
			opserr << "Need theGradGEvaluator before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theStepSizeRule == 0 ) {
			opserr << "Need theStepSizeRule before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theSearchDirection == 0 ) {
			opserr << "Need theSearchDirection before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theProbabilityTransformation == 0 ) {
			opserr << "Need theProbabilityTransformation before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theStartPoint == 0 ) {
			opserr << "Need theStartPoint before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theReliabilityConvergenceCheck == 0 ) {
			opserr << "Need theReliabilityConvergenceCheck before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}

		int printFlag=0;
		char *fileNamePrint;
		fileNamePrint = new char[256];
		strcpy(fileNamePrint,"initialized");


		int maxNumIter = 100;
		while (argvCounter < argc) {

			if (strcmp(argv[argvCounter],"-maxNumIter") == 0) {
				argvCounter++;

				if (Tcl_GetInt(interp, argv[argvCounter], &maxNumIter) != TCL_OK) {
					opserr << "ERROR: invalid input: maxNumIter \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printAllPointsX") == 0) {
				argvCounter++;
				printFlag = 1;
				strcpy(fileNamePrint,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printAllPointsY") == 0) {
				argvCounter++;
				printFlag = 2;
				strcpy(fileNamePrint,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printDesignPointX") == 0) {
				argvCounter++;
				printFlag = 3;
				strcpy(fileNamePrint,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printDesignPointY") == 0) {
				argvCounter++;
				printFlag = 4;
				strcpy(fileNamePrint,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printCurrentPointX") == 0) {
				argvCounter++;
				printFlag = 5;
				strcpy(fileNamePrint,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printCurrentPointY") == 0) {
				argvCounter++;
				printFlag = 6;
				strcpy(fileNamePrint,argv[argvCounter]);
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to SearchWithStepSizeAndStepDirection. " << endln;
				return TCL_ERROR;
			}
		}
		
		theFindDesignPointAlgorithm = new SearchWithStepSizeAndStepDirection(
					maxNumIter, 
					theGFunEvaluator,
					theGradGEvaluator,
					theStepSizeRule,
					theSearchDirection,
					theProbabilityTransformation,
					theHessianApproximation,
					theReliabilityConvergenceCheck,
					printFlag,
					fileNamePrint,
					theStartPoint);

		delete [] fileNamePrint;
		
	}
	else {
		opserr << "ERROR: unrecognized type of FindDesignPointAlgorithm Algorithm \n";
		return TCL_ERROR;
	}

	if (theFindDesignPointAlgorithm == 0) {
		opserr << "ERROR: could not create theFindDesignPointAlgorithm \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}




//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addStartPoint(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theStartPoint != 0) {
		delete theStartPoint;
		theStartPoint = 0;
	}


	// Check that there are enough arguments
	if (argc<2) {
		opserr << "ERROR: Not enough arguments to theStartPoint. " << endln;
		return TCL_ERROR;
	}

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	RandomVariable *aRandomVariable;


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"Mean") == 0) {

		theStartPoint = new Vector(nrv);

		for ( int i=1; i<=nrv; i++ )
		{
			aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i);
			if (aRandomVariable == 0) {
				opserr << "ERROR: when creating theStartPoint - could not find" << endln
					<< " random variable with tag #" << i << "." << endln;
				return TCL_ERROR;
			}
			(*theStartPoint)(i-1) = aRandomVariable->getMean();
		}
	}
	else if (strcmp(argv[1],"Origin") == 0) {
		// This is the default option (at least from now on...)
		// Do nothing; theStartPoint==0 is the indication of this case
	}
	else if (strcmp(argv[1],"Given") == 0) {

		theStartPoint = new Vector(nrv);

		for ( int i=1; i<=nrv; i++ )
		{
			aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i);
			if (aRandomVariable == 0) {
				opserr << "ERROR: when creating theStartPoint - could not find" << endln
					<< " random variable with tag #" << i << "." << endln;
				return TCL_ERROR;
			}
			(*theStartPoint)(i-1) = aRandomVariable->getStartValue();
		}
	}
	else if (strcmp(argv[1],"-file") == 0) {

		theStartPoint = new Vector(nrv);

		ifstream inputFile( argv[2], ios::in );
		if (inputFile.fail()) {
			opserr << "File " << *argv[2] << " could not be opened. " << endln;
			return TCL_ERROR;
		}

		// Loop through file to see how many entries there are
		int numRVs = theReliabilityDomain->getNumberOfRandomVariables();
		double dummy;
		int numEntries = 0;
		while (inputFile >> dummy) {
			numEntries++;
		}
		if (numEntries == 0) {
			opserr << "ERROR: No entries in the file read by startPoint!" << endln;
			return TCL_ERROR;
		}
		if (numEntries != numRVs) {
			opserr << "ERROR: Wrong number of entries in the file read by startPoint." << endln;
			return TCL_ERROR;
		}

		// Close the file
		inputFile.close();

		// Open it again, now being ready to store the results
		ifstream inputFile2( argv[2], ios::in );
		for (int i=0; i<nrv; i++){
			inputFile2 >> (*theStartPoint)(i);
		}
		inputFile2.close();

	}
	else {
		opserr << "ERROR: Invalid type of start point is given. " << endln;
		return TCL_ERROR;
	}

	// Check that the vector is of correct size
	if (theStartPoint==0) {
		opserr << "ERROR: Could not create the start point. " << endln;
		return TCL_ERROR;
	}
	else {
		if (theStartPoint->Size() != nrv) {
			opserr << "ERROR: The size of the start point vector is NOT equal " << endln
				<< " to the number of random variables in the model! " << endln;
			return TCL_ERROR;
		}
	}

	return TCL_OK;
}







//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addRootFinding(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theRootFindingAlgorithm != 0) {
		delete theRootFindingAlgorithm;
		theRootFindingAlgorithm = 0;
	}


	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a root-finding algorithm can be created" << endln;
		return TCL_ERROR;
	}

	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a root-finding algorithm can be created" << endln;
		return TCL_ERROR;
	}


	int maxIter = 50;
	double tol = 1.0e-3;
	double maxStepLength = 1.0;

	int argvCounter = 2;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-maxIter") == 0) {
			argvCounter++;

			if (Tcl_GetInt(interp, argv[argvCounter], &maxIter) != TCL_OK) {
				opserr << "ERROR: invalid input: maxIter for projection" << endln;
				return TCL_ERROR;
			}
			argvCounter++;

		}
		else if (strcmp(argv[argvCounter],"-tol") == 0) {
			argvCounter++;

			if (Tcl_GetDouble(interp, argv[argvCounter], &tol) != TCL_OK) {
				opserr << "ERROR: invalid input: tol factor \n";
				return TCL_ERROR;
			}		
			argvCounter++;

		}
		else if (strcmp(argv[argvCounter],"-maxStepLength") == 0) {
			argvCounter++;


			if (Tcl_GetDouble(interp, argv[argvCounter], &maxStepLength) != TCL_OK) {
				opserr << "ERROR: invalid input: maxStepLength factor \n";
				return TCL_ERROR;
			}		
			argvCounter++;
		}
		else {
			opserr << "ERROR: Invalid input to projection algorithm. " << endln;
			return TCL_ERROR;
		}

	}


	if (strcmp(argv[1],"Secant") == 0) {

		theRootFindingAlgorithm = new SecantRootFinding(
			theReliabilityDomain,
			theProbabilityTransformation,
			theGFunEvaluator,
			maxIter,
			tol,
			maxStepLength);
		
	}
	else {
		opserr << "ERROR: unrecognized type of root-finding algorithm \n";
		return TCL_ERROR;
	}

	if (theRootFindingAlgorithm == 0) {
		opserr << "ERROR: could not create root-finding algorithm \n";
		return TCL_ERROR;
	}


	return TCL_OK;
}




//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runFORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theFORMAnalysis != 0) {
		delete theFORMAnalysis;
		theFORMAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	// Check number of arguments
	if ( (argc!=2) && (argc!=4))  {
		opserr << "ERROR: Wrong number of input parameter to FORM analysis" << endln;
		return TCL_ERROR;
	}


	// Check for essential tools
	if (theFindDesignPointAlgorithm == 0 ) {
		opserr << "Need theFindDesignPointAlgorithm before a FORMAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a FORMAnalysis can be created" << endln;
		return TCL_ERROR;
	}


	// Read input parameter(s)
	int relSensTag = 0;
	if (argc == 4) {
		if (strcmp(argv[2],"-relSens") == 0) {
			if (Tcl_GetInt(interp, argv[3], &relSensTag) != TCL_OK) {
				opserr << "ERROR: invalid input: relSensTag \n";
				return TCL_ERROR;
			}
		}
		else {
			opserr << "ERROR: Invalid input to FORMAnalysis." << endln;
			return TCL_ERROR;
		}
	}


	// Create the analysis object
	theFORMAnalysis 
		= new FORMAnalysis( theReliabilityDomain, 
							theFindDesignPointAlgorithm, 
							theProbabilityTransformation, 
							argv[1] , 
							relSensTag);


	// Check that it really was created
	if (theFORMAnalysis == 0) {
		opserr << "ERROR: could not create theFORMAnalysis \n";
		return TCL_ERROR;
	}


	// Now run the analysis
	theFORMAnalysis->analyze();

	return TCL_OK;
}





//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runMVFOSMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theMVFOSMAnalysis != 0) {
		delete theMVFOSMAnalysis;
		theMVFOSMAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	// Check number of arguments
	if (argc != 2)  {
		opserr << "ERROR: Wrong number of input parameter to MVFOSM analysis" << endln;
		return TCL_ERROR;
	}


	// Check for essential ingredients
	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a MVFOSMAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGradGEvaluator == 0 ) {
		opserr << "Need theGradGEvaluator before a MVFOSMAnalysis can be created" << endln;
		return TCL_ERROR;
	}


	theMVFOSMAnalysis = new MVFOSMAnalysis( theReliabilityDomain,
											theGFunEvaluator,
											theGradGEvaluator,
											interp,
											argv[1]);

	if (theMVFOSMAnalysis == 0) {
		opserr << "ERROR: could not create theMVFOSMAnalysis \n";
		return TCL_ERROR;
	}

	// Now run the analysis
	theMVFOSMAnalysis->analyze();

	return TCL_OK;
}






//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runFragilityAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theFragilityAnalysis != 0) {
		delete theFragilityAnalysis;
		theFragilityAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	// Check number of arguments
	if (argc != 9)  {
		opserr << "ERROR: Wrong number of input parameter to Fragility analysis" << endln;
		return TCL_ERROR;
	}


	// Check for essential ingredients
	if (theFindDesignPointAlgorithm == 0 ) {
		opserr << "Need theFindDesignPointAlgorithm before a FragilityAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGradGEvaluator == 0 ) {
		opserr << "Need theGradGEvaluator before a FragilityAnalysis can be created" << endln;
		return TCL_ERROR;
	}


	// Read input
	bool parGiven = false;
	bool rangeGiven = false; 
	bool numIntGiven = false;
	int parameterNumber;
	double first;
	double last;
	int numIntervals;
	int counter = 3;
	for (int i=1; i<=3; i++) {

		if (strcmp(argv[counter-1],"-par") == 0) {
			// GET INPUT PARAMETER (int)
			if (Tcl_GetInt(interp, argv[counter], &parameterNumber) != TCL_OK) {
				opserr << "ERROR: invalid input: parameterNumber \n";
				return TCL_ERROR;
			}
			counter++; counter++;
			parGiven = true;
		}
		else if (strcmp(argv[counter-1],"-range") == 0) {
			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[counter], &first) != TCL_OK) {
				opserr << "ERROR: invalid input: first bound to range \n";
				return TCL_ERROR;
			}
			counter++; 
			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[counter], &last) != TCL_OK) {
				opserr << "ERROR: invalid input: last bound to range \n";
				return TCL_ERROR;
			}
			counter++; counter++;
			rangeGiven = true;
		}
		else if (strcmp(argv[counter-1],"-numInt") == 0) {
			// GET INPUT PARAMETER (int)
			if (Tcl_GetInt(interp, argv[counter], &numIntervals) != TCL_OK) {
				opserr << "ERROR: invalid input: number of intervals \n";
				return TCL_ERROR;
			}
			counter++; counter++;
			numIntGiven = true;
		}
		else {
			opserr << "ERROR: invalid input to Fragility analysis " << endln;
			return TCL_ERROR;
		}

	}

	if (parGiven && rangeGiven && numIntGiven) {
		theFragilityAnalysis = new FragilityAnalysis( theReliabilityDomain,
												  theFindDesignPointAlgorithm,
												  theGradGEvaluator,
												  parameterNumber,
												  first,
												  last,
												  numIntervals,
												  argv[1],
												  interp);
	}
	else {
		opserr << "ERROR:: some input to theFragilityAnalysis was not provided" << endln;
		return TCL_ERROR;
	}

	if (theFragilityAnalysis == 0) {
		opserr << "ERROR: could not create theFragilityAnalysis \n";
		return TCL_ERROR;
	}

	// Now run the analysis
	theFragilityAnalysis->analyze();

	return TCL_OK;
}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runSORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theSORMAnalysis != 0) {
		delete theSORMAnalysis;
		theSORMAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	if (theFindCurvatures == 0 ) {
		opserr << "Need theFindCurvatures before a SORMAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theFORMAnalysis == 0 ) {
		opserr << "ERROR: The current SORM implementation requires a FORM analysis" << endln
			<< " to have been executed previously in the same session." << endln;
		return TCL_ERROR;
	}

	if (argc != 2)  {
		opserr << "ERROR: Wrong number of arguments to SORM analysis" << endln;
		return TCL_ERROR;
	}

	theSORMAnalysis 
		= new SORMAnalysis(theReliabilityDomain, theFindCurvatures , argv[1]);

	if (theSORMAnalysis == 0) {
		opserr << "ERROR: could not create theSORMAnalysis \n";
		return TCL_ERROR;
	}

	// Now run the analysis
	theSORMAnalysis->analyze();

	return TCL_OK;
}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runSystemAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theSystemAnalysis != 0) {
		delete theSystemAnalysis;
		theSystemAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	if (argc != 3)  {
		opserr << "ERROR: Wrong number of arguments to System Reliability analysis" << endln;
		return TCL_ERROR;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[2],"allInSeries") == 0) {

		theSystemAnalysis = new SystemAnalysis(theReliabilityDomain, argv[1]);

	}
	else {
		opserr << "ERROR: Invalid input to system reliability analysis" << endln;
		return TCL_ERROR;
	}


	if (theSystemAnalysis == 0) {
		opserr << "ERROR: Could not create theSystemAnalysis. " << endln;
		return TCL_ERROR;
	}


	// Now run the analysis
	theSystemAnalysis->analyze();

	return TCL_OK;
}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runSamplingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theSamplingAnalysis != 0) {
		delete theSamplingAnalysis;
		theSamplingAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	// Check for essential tools
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a SimulationAnalyis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a SimulationAnalyis can be created" << endln;
		return TCL_ERROR;
	}
	if (theRandomNumberGenerator == 0 ) {
		opserr << "Need theRandomNumberGenerator before a SimulationAnalyis can be created" << endln;
		return TCL_ERROR;
	}

	
	// The following switches are available (default values are provided)
	// (The sampling is performed around theStartPoint, except 
	// for response statistics sampling; then the mean is used together
	// with unit sampling variance.)
	//
	//     -type  failureProbability (1)......... this is the default
	//     -type  responseStatistics (2)
	//
	//     -variance 1.0  ....................... this is the default
	//
	//     -maxNum 1000  ........................ this is the default
	//
	//     -targetCOV 0.05  ..................... this is the default
	//
	//     -print 0   (print nothing) ........... this is the default
	//     -print 1   (print to screen)
	//     -print 2   (print to restart file)
	//

	if (argc!=2 && argc!=4 && argc!=6 && argc!=8 && argc!=10 && argc!=12) {
		opserr << "ERROR: Wrong number of arguments to Sampling analysis" << endln;
		return TCL_ERROR;
	}


	// Declaration of input parameters
	int numberOfSimulations	= 1000;
	double targetCOV		= 0.05;
	double samplingVariance	= 1.0;
	int printFlag			= 0;
	int analysisTypeTag		= 1;


	for (int i=2; i<argc; i=i+2) {

		if (strcmp(argv[i],"-type") == 0) {

			if (strcmp(argv[i+1],"failureProbability") == 0) {

				analysisTypeTag = 1;

				if (theStartPoint == 0 ) {
					opserr << "Need theStartPoint before a SimulationAnalyis can be created" << endln;
					return TCL_ERROR;
				}
			}
			else if (strcmp(argv[i+1],"responseStatistics") == 0) {
				analysisTypeTag = 2;
				if (samplingVariance != 1.0) {
					opserr << "ERROR:: sampling variance must be 1.0 for " << endln
						<< " response statistics sampling." << endln;
					return TCL_ERROR;
				}
				// Make sure that the mean point is the sampling center
				int nrv = theReliabilityDomain->getNumberOfRandomVariables();
				RandomVariable *aRandomVariable;
				if (theStartPoint == 0) {
					theStartPoint = new Vector(nrv);
				}
				for ( int i=1; i<=nrv; i++ )
				{
					aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i);
					if (aRandomVariable == 0) {
						opserr << "ERROR: when creating theStartPoint - could not find" << endln
							<< " random variable with tag #" << i << "." << endln;
						return TCL_ERROR;
					}
					(*theStartPoint)(i-1) = aRandomVariable->getMean();
				}
				opserr << "NOTE: The startPoint is set to the Mean due to the choice " << endln
					<< "      of responseStatistics as sampling analysis type." << endln;
			}
			else {
				opserr << "ERROR: invalid input: type \n";
				return TCL_ERROR;
			}
		}
		else if (strcmp(argv[i],"-variance") == 0) {
			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[i+1], &samplingVariance) != TCL_OK) {
				opserr << "ERROR: invalid input: samplingVariance \n";
				return TCL_ERROR;
			}
			if (analysisTypeTag == 2 && samplingVariance != 1.0) {
				opserr << "ERROR:: sampling variance must be 1.0 for " << endln
					<< " response statistics sampling." << endln;
				return TCL_ERROR;
			}
		}
		else if (strcmp(argv[i],"-maxNum") == 0) {
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[i+1], &numberOfSimulations) != TCL_OK) {
				opserr << "ERROR: invalid input: numberOfSimulations \n";
				return TCL_ERROR;
			}
		}
		else if (strcmp(argv[i],"-targetCOV") == 0) {
			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[i+1], &targetCOV) != TCL_OK) {
				opserr << "ERROR: invalid input: targetCOV \n";
				return TCL_ERROR;
			}
		}
		else if (strcmp(argv[i],"-print") == 0) {
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[i+1], &printFlag) != TCL_OK) {
				opserr << "ERROR: invalid input: printFlag \n";
				return TCL_ERROR;
			}
		}
		else {
			opserr << "ERROR: invalid input to sampling analysis. " << endln;
			return TCL_ERROR;
		}
	}

	// Warn about illegal combinations
	if (analysisTypeTag==2 && printFlag==2) {
		opserr << "ERROR:: The restart option of the sampling analysis cannot be " << endln
			<< " used together with the response statistics option. " << endln;
		return TCL_ERROR;
	}
	
	
	theSamplingAnalysis 
			= new SamplingAnalysis(theReliabilityDomain, 
									theProbabilityTransformation, 
									theGFunEvaluator, 
									theRandomNumberGenerator, 
									numberOfSimulations, 
									targetCOV,
									samplingVariance,
									printFlag,
									argv[1],
									theStartPoint,
									analysisTypeTag);

	if (theSamplingAnalysis == 0) {
		opserr << "ERROR: could not create theSamplingAnalysis \n";
		return TCL_ERROR;
	}

	// Now run analysis
	theSamplingAnalysis->analyze();

	return TCL_OK;

}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runOutCrossingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theOutCrossingAnalysis != 0) {
		delete theOutCrossingAnalysis;
		theOutCrossingAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	if (theFindDesignPointAlgorithm == 0 ) {
		opserr << "Need theFindDesignPointAlgorithm before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGradGEvaluator == 0 ) {
		opserr << "Need theGradGEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}

	int sampleFreq = 1;
	double littleDt = 0.01;
	int analysisType = 1;

	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-resultFreq") == 0) {
			argvCounter++;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &sampleFreq) != TCL_OK) {
				opserr << "ERROR: invalid input sampleFreq to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-littleDt") == 0) {
			argvCounter++;

			  // GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &littleDt) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-Koo") == 0) {
			argvCounter++;
			analysisType = 2;
		}
		else if (strcmp(argv[argvCounter],"-twoSearches") == 0) {
			argvCounter++;
			analysisType = 1;
		}
		else {
			opserr << "ERROR: Invalid input to theOutCrossingAnalysis." << endln;
			argvCounter++;
		}
	}

	theOutCrossingAnalysis 
			= new OutCrossingAnalysis(
				theReliabilityDomain,
				theGFunEvaluator,
				theGradGEvaluator,
				theFindDesignPointAlgorithm,
				analysisType,
				sampleFreq,
				littleDt,
				argv[1]);

	if (theOutCrossingAnalysis == 0) {
		opserr << "ERROR: could not create theOutCrossingAnalysis \n";
		return TCL_ERROR;
	}

	// Now run analysis
	theOutCrossingAnalysis->analyze();

	return TCL_OK;

}


//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runGFunVisualizationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theGFunVisualizationAnalysis != 0) {
		delete theGFunVisualizationAnalysis;
		theGFunVisualizationAnalysis = 0;
	}


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a GFunVisualizationAnalysis can be created" << endln;
		return TCL_ERROR;
	}

	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a GFunVisualizationAnalysis can be created" << endln;
		return TCL_ERROR;
	}


	// Initial declarations
	int rv1 = 0;
	int rv2 = 0;
	int numPts1 = 0;
	int numPts2 = 0;
	double from1 = 0.0;
	double to1 = 0.0;
	double from2 = 0.0;
	double to2 = 0.0;

	int rvDir;
	Matrix theMatrix;
	int numLinePts;
	Vector theDirectionVector;
	Vector axesVector;
	int convFileArgv = 0;

	// Tags to keep track of which options the users chooses
	// (and to check which ones have been given)
	int convResults = 0;
	int space = 0;
	int funSurf = 0;
	int dir = 0;
	int axes = 0;


	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-convResults") == 0) {
			argvCounter++;

			convFileArgv = argvCounter;
			if ((argc-1)<argvCounter) {
				opserr << "ERROR: No file name found for visualization of convergence results. " << endln;
				return TCL_ERROR;
			}
			argvCounter++;

			convResults = 1;		
		}
		else if (strcmp(argv[argvCounter],"-space") == 0) {
			argvCounter++;

			if (strcmp(argv[argvCounter],"X") == 0 || strcmp(argv[argvCounter],"x") == 0) {
				space = 1;
			}
			else if (strcmp(argv[argvCounter],"Y") == 0 || strcmp(argv[argvCounter],"y") == 0) {
				space = 2;
			}
			else {
				opserr << "ERROR: Invalid input to visualization analysis. " << endln;
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-funSurf") == 0) {
			argvCounter++;

			if (strcmp(argv[argvCounter],"function") == 0 ) {
				funSurf = 1;
			}
			else if (strcmp(argv[argvCounter],"surface") == 0 ) {
				funSurf = 2;
			}
			else {
				opserr << "ERROR: Invalid input to visualization analysis. " << endln;
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-dir") == 0) {
			argvCounter++;

			if (strcmp(argv[argvCounter],"rv") == 0 ) {
				argvCounter++;

				dir = 1;

				// GET INPUT PARAMETER (integer)
				if (Tcl_GetInt(interp, argv[argvCounter], &rvDir) != TCL_OK) {
					opserr << "ERROR: invalid input: rvDir  in theGFunVisualizationAnalysis \n";
					return TCL_ERROR;
				}
				argvCounter++;

			}
			else if (strcmp(argv[argvCounter],"file") == 0 ) {
				argvCounter++;

				dir = 2;

				// Open file where the vectors are given
				ifstream inputFile( argv[argvCounter], ios::in );
				if (inputFile.fail()) {
					opserr << "File " << *argv[argvCounter] << " could not be opened. " << endln;
					return TCL_ERROR;
				}

				// Loop through file to see how many entries there are
				int numRVs = theReliabilityDomain->getNumberOfRandomVariables();
				double dummy;
				int numEntries = 0;
				while (inputFile >> dummy) {
					numEntries++;
				}
				if (numEntries == 0) {
					opserr << "ERROR: No entries in the direction file read by visualization analysis!" << endln;
					return TCL_ERROR;
				}

				// Check that the number of points are ok
				if (numEntries != numRVs) {
					opserr << "ERROR: Wrong number of entries in the the file " << argv[argvCounter] << endln;
					return TCL_ERROR;
				}

				// Close the file
				inputFile.close();

				// Open it again, now being ready to store the results in a matrix
				ifstream inputFile2( argv[argvCounter], ios::in );
				if (inputFile2.fail()) {
					opserr << "File " << *argv[argvCounter] << " could not be opened. " << endln;
					return TCL_ERROR;
				}
				argvCounter++;

				// Store the vector
				Vector dummyDirectionVector(numRVs);
				for (int i=0; i<numRVs; i++) {
						inputFile2 >> dummyDirectionVector(i);
				}
				inputFile2.close();

				theDirectionVector = dummyDirectionVector;

				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to visualization analysis. " << endln;
				return TCL_ERROR;
			}
		}
		else if (strcmp(argv[argvCounter],"-coords1") == 0) {

			axes = 1;
			argvCounter++;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &rv1) != TCL_OK) {
				opserr << "ERROR: invalid input: rv1  in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			  // GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &from1) != TCL_OK) {
				opserr << "ERROR: invalid input: from1 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &to1) != TCL_OK) {
				opserr << "ERROR: invalid input: to1 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (double)
			if (Tcl_GetInt(interp, argv[argvCounter], &numPts1) != TCL_OK) {
				opserr << "ERROR: invalid input: numPts1 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			Vector dummy(4);
			dummy(0) = (double)rv1;
			dummy(1) = from1;
			dummy(2) = to1;
			dummy(3) = (double)numPts1;
			axesVector = dummy;

		}
		else if (strcmp(argv[argvCounter],"-coords2") == 0) {

			axes = 2;
			argvCounter++;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &rv1) != TCL_OK) {
				opserr << "ERROR: invalid input: rv1  in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			  // GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &from1) != TCL_OK) {
				opserr << "ERROR: invalid input: from1 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &to1) != TCL_OK) {
				opserr << "ERROR: invalid input: to1 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (double)
			if (Tcl_GetInt(interp, argv[argvCounter], &numPts1) != TCL_OK) {
				opserr << "ERROR: invalid input: numPts1 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &rv2) != TCL_OK) {
				opserr << "ERROR: invalid input: rv2  in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			  // GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &from2) != TCL_OK) {
				opserr << "ERROR: invalid input: from2 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &to2) != TCL_OK) {
				opserr << "ERROR: invalid input: to2 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (double)
			if (Tcl_GetInt(interp, argv[argvCounter], &numPts2) != TCL_OK) {
				opserr << "ERROR: invalid input: numPts2 in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			Vector dummy(8);
			dummy(0) = (double)rv1;
			dummy(1) = from1;
			dummy(2) = to1;
			dummy(3) = (double)numPts1;
			dummy(4) = (double)rv2;
			dummy(5) = from2;
			dummy(6) = to2;
			dummy(7) = (double)numPts2;
			axesVector = dummy;

		}
		else if (strcmp(argv[argvCounter],"-file") == 0) {

			axes = 3;
			argvCounter++;

			// Open file where the vectors are given
			ifstream inputFile( argv[argvCounter], ios::in );
			if (inputFile.fail()) {
				opserr << "File " << *argv[argvCounter] << " could not be opened. " << endln;
				return TCL_ERROR;
			}

			// Loop through file to see how many entries there are
			int numRVs = theReliabilityDomain->getNumberOfRandomVariables();
			int numVectors;
			double dummy;
			int numEntries = 0;
			while (inputFile >> dummy) {
				numEntries++;
			}
			if (numEntries == 0) {
				opserr << "ERROR: No entries in the file read by visualization analysis!" << endln;
				return TCL_ERROR;
			}

			// Check that the number of points are ok
			if (fmod(numEntries,numRVs)!=0.0) {
				opserr << "ERROR: Wrong number of entries in the the file " << argv[argvCounter] << endln;
				return TCL_ERROR;
			}
			numVectors = (int)(numEntries/numRVs);

			// Close the file
			inputFile.close();

			// Open it again, now being ready to store the results in a matrix
			ifstream inputFile2( argv[argvCounter], ios::in );
			if (inputFile2.fail()) {
				opserr << "File " << *argv[argvCounter] << " could not be opened. " << endln;
				return TCL_ERROR;
			}
			argvCounter++;

			// Store the vectors in a matrix
			Matrix dummyMatrix(numRVs,numVectors);
			for (int i=0; i<numVectors; i++) {
				for (int j=0; j<numRVs; j++ ) {
					inputFile2 >> dummyMatrix(j,i);
				}
			}
			inputFile2.close();

			theMatrix = dummyMatrix;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &numLinePts) != TCL_OK) {
				opserr << "ERROR: invalid input: numPts  in theGFunVisualizationAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else {
			opserr << "ERROR: invalid input to theGFunVisualizationAnalysis." << endln;
			return TCL_ERROR;
		}
	}


	// Check that the input was more or less reasonable
	// convResults [ 0:no,                   1:yes                 ]
	// space       [ 0:error,                1:X,        2:Y       ]
	// funSurf     [ 0:error,                1:function, 2:surface ]
	// dir         [ 0:(needed for surface), 1:rv,       2:file    ] (pass rvDir or theDirectionVector)
	// axes        [ 0:error,   1:coords1,   2:coords2,  3:file    ] (pass axesVector... or theMatrix+numLinePts)

	if (space==0 || funSurf==0 || axes==0) {
		opserr << "ERROR: Some input is missing to the visualization analysis." << endln;
		return TCL_ERROR;
	}
	if (dir==0 && funSurf==2) {
		opserr << "A direction is needed for visualization of the limit-state surface." << endln;
		return TCL_ERROR;
	}

	theGFunVisualizationAnalysis = new GFunVisualizationAnalysis(
											theReliabilityDomain, 
											theGFunEvaluator, 
											theProbabilityTransformation, 
											argv[1],
											argv[convFileArgv],
											convResults,
											space,
											funSurf,
											axes,
											dir);


	// Pass stuff to the analysis object
	if (dir == 1) {
		theGFunVisualizationAnalysis->setDirection(rvDir);
	}
	else if (dir == 2) {
		theGFunVisualizationAnalysis->setDirection(theDirectionVector);
	}

	if (axes == 1 || axes == 2) {
		theGFunVisualizationAnalysis->setAxes(axesVector);
	}
	else if (axes == 3) {
		theGFunVisualizationAnalysis->setAxes(theMatrix);
		theGFunVisualizationAnalysis->setNumLinePts(numLinePts)	;	
	}
	
	if (axes == 1 || axes == 2) {

		if (theStartPoint == 0 ) {
			opserr << "Need theStartPoint before this GFunVisualizationAnalysis can be created" << endln;
			return TCL_ERROR;
		}
		
		theGFunVisualizationAnalysis->setStartPoint(theStartPoint);
	}

	if (convResults == 1) {

		if (theGradGEvaluator == 0 ) {
			opserr << "Need theGradGEvaluator before this GFunVisualizationAnalysis can be created" << endln;
			return TCL_ERROR;
		}

		if (theMeritFunctionCheck == 0 ) {
			opserr << "Need theMeritFunctionCheck before this GFunVisualizationAnalysis can be created" << endln;
			return TCL_ERROR;
		}

		if (theReliabilityConvergenceCheck == 0 ) {
			opserr << "Need theReliabilityConvergenceCheck before this GFunVisualizationAnalysis can be created" << endln;
			return TCL_ERROR;
		}

		theGFunVisualizationAnalysis->setGradGEvaluator(theGradGEvaluator);
		theGFunVisualizationAnalysis->setMeritFunctionCheck(theMeritFunctionCheck);
		theGFunVisualizationAnalysis->setReliabilityConvergenceCheck(theReliabilityConvergenceCheck);
	}

	if (funSurf == 2) {

		if (theRootFindingAlgorithm == 0 ) {
			opserr << "Need theRootFindingAlgorithm before this GFunVisualizationAnalysis can be created" << endln;
			return TCL_ERROR;
		}

		theGFunVisualizationAnalysis->setRootFindingAlgorithm(theRootFindingAlgorithm);
	}


	// It is chosen to have only one constructor. 
	// Hence, pass some stuff via methods, depending on analysis options. 


	if (theGFunVisualizationAnalysis == 0) {
		opserr << "ERROR: could not create theGFunVisualizationAnalysis \n";
		return TCL_ERROR;
	}

	// Now run the analysis
	theGFunVisualizationAnalysis->analyze();

	return TCL_OK;
}






//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_inputCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// Check that tagged objects are consequtive
	int i, num;

	num = theReliabilityDomain->getNumberOfRandomVariables();
	ReliabilityDomainComponent *component;
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getRandomVariablePtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive random variable list." << endln;
			return TCL_ERROR;
		}
	}
	
	num = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive random variable positioner list." << endln;
			return TCL_ERROR;
		}
	}
	
	num = theReliabilityDomain->getNumberOfCorrelationCoefficients();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getCorrelationCoefficientPtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive correlation coefficient list." << endln;
			return TCL_ERROR;
		}
	}
	
	num = theReliabilityDomain->getNumberOfFilters();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getFilter(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive filter list." << endln;
			return TCL_ERROR;
		}
	}
	
	num = theReliabilityDomain->getNumberOfLimitStateFunctions();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getLimitStateFunctionPtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive limit-state (performance) function list." << endln;
			return TCL_ERROR;
		}
	}
	
	num = theReliabilityDomain->getNumberOfModulatingFunctions();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getModulatingFunction(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive modulating function list." << endln;
			return TCL_ERROR;
		}
	}
	
	num = theReliabilityDomain->getNumberOfSpectra();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getSpectrum(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive spectrum list." << endln;
			return TCL_ERROR;
		}
	}



	// Check that the correlation matrix is positive definite
	// theCorrelationMatrix


	return TCL_OK;
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_tempCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	opserr << "The temp command does nothing now!" << endln;
	return TCL_OK;
}


