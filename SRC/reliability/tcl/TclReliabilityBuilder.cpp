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
                                                                        
// $Revision: 1.57 $
// $Date: 2010-09-13 21:40:25 $
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
#include <RandomVariableIter.h>
#include <CorrelationCoefficient.h>
#include <Cutset.h>
#include <CutsetIter.h>
#include <PerformanceFunction.h>
#include <LimitStateFunction.h>
#include <LimitStateFunctionIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <RVParameter.h>

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

#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <HessianEvaluator.h>
#include <TclEvaluator.h>
#include <ImplicitGradient.h>
#include <FiniteDifferenceGradient.h>
#include <FiniteDifferenceHessian.h>

#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <ProbabilityTransformation.h>
#include <NatafProbabilityTransformation.h>
#include <FindDesignPointAlgorithm.h>
#include <HLRFSearchDirection.h>
#include <ArmijoStepSizeRule.h>
#include <FixedStepSizeRule.h>
#include <SearchWithStepSizeAndStepDirection.h>
#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>
#include <FindCurvatures.h>
#include <FirstPrincipalCurvature.h>
#include <CurvaturesBySearchAlgorithm.h>
#include <CurvatureFitting.h>
#include <ReliabilityConvergenceCheck.h>
#include <StandardReliabilityConvergenceCheck.h>
#include <OptimalityConditionReliabilityConvergenceCheck.h>
#include <MeritFunctionCheck.h>
#include <AdkZhangMeritFunctionCheck.h>
#include <PolakHeSearchDirectionAndMeritFunction.h>
#include <SQPsearchDirectionMeritFunctionAndHessian.h>
#include <GradientProjectionSearchDirection.h>

#include <ReliabilityAnalysis.h>
//#include <OpenSeesGFunEvaluator.h>
//#include <OpenSeesGradGEvaluator.h>
#include <FORMAnalysis.h>
#include <FOSMAnalysis.h>
//#include <ParametricReliabilityAnalysis.h>
#include <GFunVisualizationAnalysis.h>
#include <OutCrossingAnalysis.h>
#include <ImportanceSamplingAnalysis.h>
#include <SORMAnalysis.h>
#include <SystemAnalysis.h>
#include <PCM.h>
#include <IPCM.h>
#include <SCIS.h>
#include <MVNcdf.h>

#include <Filter.h>
#include <KooFilter.h>
#include <StandardLinearOscillatorDisplacementFilter.h>
#include <StandardLinearOscillatorVelocityFilter.h>
#include <StandardLinearOscillatorAccelerationFilter.h>
#include <DeltaFilter.h>
#include <ModulatingFunction.h>
#include <GammaModulatingFunction.h>
#include <ConstantModulatingFunction.h>
#include <TrapezoidalModulatingFunction.h>
#include <KooModulatingFunction.h>
#include <Spectrum.h>
#include <JonswapSpectrum.h>
#include <NarrowBandSpectrum.h>
#include <PointsSpectrum.h>

//#include <SensitivityAlgorithm.h>
#include <RootFinding.h>
#include <SecantRootFinding.h>

//Quan---
#include <MonteCarloResponseAnalysis.h>
#include <OrthogonalPlaneSamplingAnalysis.h>
#include <MultiDimVisPrincPlane.h>
#include <DP_RSM_Sim.h>
#include <DP_RSM_Sim_TimeVariant.h>
//---Quan

#include <TclReliabilityBuilder.h>
/////////////////////////////////////////////////////////
///S added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////
#include <Analyzer.h>
#include <StaticAnalyzer.h>
#include <DynamicAnalyzer.h>
#include <InitialStaticAnalysis.h>
#include <SelectLoadInitialStaticAnalysis.h>
#include <AnalyzerGFunEvaluator.h>
#include <AnalyzerGradGEvaluator.h>
#include <NewWhitenoiseFilter.h>
#include <NewStandardLinearOscillatorAccelerationFilter.h>
#include <NewSearchWithStepSizeAndStepDirection.h>
#include <InitialPointBuilder.h>
#include <ThresholdIncInitialPointBuilder.h>
#include <CrossingRateAnalyzer.h>
#include <FOSeriesSimulation.h>
#include <FirstPassageAnalyzer.h>
#include <StatFirstPassageAnalyzer.h>
#include <NonStatFirstPassageAnalyzer.h>
#include <RandomVibrationSimulation.h>
#include <StatRandomVibrationSimulation.h>
#include <NonStatRandomVibrationSimulation.h>
#include <RandomVibrationAnalysis.h>
#include <AllIndependentTransformation.h>
#include<Integrator.h>//Abbas
/////////////////////////////////////////////////////////
/////E Modified by K Fujimura /////////////////////////////
/////////////////////////////////////////////////////////

//extern SensitivityAlgorithm *theSensitivityAlgorithm;
extern Integrator *theSensitivityAlgorithm;

/////////////////////////////////////////////////////////
/////S Modified by K Fujimura /////////////////////////////
/////////////////////////////////////////////////////////
extern ReliabilityStaticAnalysis* theReliabilityStaticAnalysis;
extern ReliabilityDirectIntegrationAnalysis* theReliabilityTransientAnalysis;
extern Integrator* theSensitivityIntegrator;
//extern SensitivityIntegrator* theSensitivityIntegrator; //Abbas



/////////////////////////////////////////////////////////
/////E Modified by K Fujimura /////////////////////////////
/////////////////////////////////////////////////////////

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

// Quan --
// ---------------- define global pointer for SNOPT ----------
MonteCarloResponseAnalysis * theMonteCarloResponseAnalysis=0;
SamplingAnalysis * theSamplingAnalysis =0;
// --- Quan

ReliabilityDomain *theReliabilityDomain = 0;
static Domain *theStructuralDomain = 0;

// base class static pointers
static FunctionEvaluator *theFunctionEvaluator = 0;
static GradientEvaluator *theGradientEvaluator = 0;
static StepSizeRule *theStepSizeRule = 0;
static SearchDirection *theSearchDirection = 0;
static HessianEvaluator *theHessianEvaluator = 0;
static MeritFunctionCheck *theMeritFunctionCheck = 0;
static ProbabilityTransformation *theProbabilityTransformation = 0;
static ReliabilityConvergenceCheck *theReliabilityConvergenceCheck = 0;
static RootFinding *theRootFindingAlgorithm = 0;
static FindCurvatures *theFindCurvatures = 0;
static FindDesignPointAlgorithm *theFindDesignPointAlgorithm = 0;
RandomNumberGenerator *theRandomNumberGenerator = 0;

// mixed pointers
static PolakHeSearchDirectionAndMeritFunction *thePolakHeDualPurpose = 0;
static SQPsearchDirectionMeritFunctionAndHessian *theSQPtriplePurpose = 0;

// analysis base class pointers
static GFunVisualizationAnalysis *theGFunVisualizationAnalysis = 0;
static FORMAnalysis *theFORMAnalysis = 0;
static FOSMAnalysis *theFOSMAnalysis = 0;
//static ParametricReliabilityAnalysis *theParametricReliabilityAnalysis = 0;
static OutCrossingAnalysis *theOutCrossingAnalysis = 0;
static SORMAnalysis *theSORMAnalysis = 0;
static ImportanceSamplingAnalysis *theImportanceSamplingAnalysis = 0;
static SystemAnalysis *theSystemAnalysis = 0;

/////////////////////////////////////////////////////////
///S added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////
static Analyzer *theAnalyzer=0;
static InitialStaticAnalysis *theInitialStaticAnalysis=0;
static InitialPointBuilder *theInitialPointBuilder=0;
static CrossingRateAnalyzer *theCrossingRateAnalyzer=0;
static FOSeriesSimulation *theFOSeriesSimulation= 0;
static FirstPassageAnalyzer *theFirstPassageAnalyzer= 0;
static RandomVibrationSimulation *theRandomVibrationSimulation= 0;
static RandomVibrationAnalysis *theRandomVibrationAnalysis = 0;
/////////////////////////////////////////////////////////
///E added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//
int TclReliabilityModelBuilder_addRandomVariable(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addCorrelate(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addCutset(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_correlateGroup(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_correlationStructure(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addLimitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addGradLimitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addModulatingFunction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFilter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addSpectrum(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addProbabilityTransformation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addStartPoint(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addRootFinding(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addRandomNumberGenerator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addSearchDirection(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addHessianEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addMeritFunctionCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addReliabilityConvergenceCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addStepSizeRule(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFunctionEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addGradientEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFindDesignPointAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFindCurvatures(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runFORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runFOSMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
//int TclReliabilityModelBuilder_runParametricReliabilityAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runGFunVisualizationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runOutCrossingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSystemAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runImportanceSamplingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_printReliability(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getMean(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getStdv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
//int TclReliabilityModelBuilder_rvReduction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getBetaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getGammaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getAlphaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getPDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getStdNormalPDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getStdNormalCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getInverseCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getStdNormalInverseCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getRVTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getLSFTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getCutsetTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getCutsetComponents(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
/////////////////////////////////////////////////////////
///S added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////
int TclReliabilityModelBuilder_addAnalyzer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addInitialStaticAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addInitialPointBuilder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addCrossingRateAnalyzer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFOSeriesSimulation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addFirstPassageAnalyzer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addRandomVibrationSimulation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runRandomVibrationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
/////////////////////////////////////////////////////////
///E added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////

//--Quan --
int TclReliabilityModelBuilder_runMonteCarloResponseAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_updateParameterValue(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runOrthogonalPlaneSamplingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_computeHessian(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_MultiDimVisPrincPlane(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_transformXtoU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_transformUtoX(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runDP_RSM_SimTimeInvariantAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runDP_RSM_SimTimeVariantAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


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
  Tcl_CreateCommand(interp, "cutset", TclReliabilityModelBuilder_addCutset,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "performanceFunction", TclReliabilityModelBuilder_addLimitState,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "gradPerformanceFunction", TclReliabilityModelBuilder_addGradLimitState,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "modulatingFunction",TclReliabilityModelBuilder_addModulatingFunction,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "filter",TclReliabilityModelBuilder_addFilter,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "spectrum",TclReliabilityModelBuilder_addSpectrum,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "findDesignPoint",	TclReliabilityModelBuilder_addFindDesignPointAlgorithm,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "startPoint",	TclReliabilityModelBuilder_addStartPoint,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "rootFinding",	TclReliabilityModelBuilder_addRootFinding,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "functionEvaluator",	TclReliabilityModelBuilder_addFunctionEvaluator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "gradientEvaluator",TclReliabilityModelBuilder_addGradientEvaluator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "stepSizeRule",TclReliabilityModelBuilder_addStepSizeRule,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "searchDirection",	TclReliabilityModelBuilder_addSearchDirection,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "hessianEvaluator",	TclReliabilityModelBuilder_addHessianEvaluator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "meritFunctionCheck",	TclReliabilityModelBuilder_addMeritFunctionCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "reliabilityConvergenceCheck",	TclReliabilityModelBuilder_addReliabilityConvergenceCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "probabilityTransformation",	TclReliabilityModelBuilder_addProbabilityTransformation,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "findCurvatures",	TclReliabilityModelBuilder_addFindCurvatures,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "randomNumberGenerator",TclReliabilityModelBuilder_addRandomNumberGenerator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runFORMAnalysis",TclReliabilityModelBuilder_runFORMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runFOSMAnalysis",TclReliabilityModelBuilder_runFOSMAnalysis,(ClientData)NULL, NULL);
  //  Tcl_CreateCommand(interp, "runParametricReliabilityAnalysis",TclReliabilityModelBuilder_runParametricReliabilityAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runGFunVizAnalysis",TclReliabilityModelBuilder_runGFunVisualizationAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runOutCrossingAnalysis",TclReliabilityModelBuilder_runOutCrossingAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runSORMAnalysis",TclReliabilityModelBuilder_runSORMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runSystemAnalysis",TclReliabilityModelBuilder_runSystemAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runImportanceSamplingAnalysis",TclReliabilityModelBuilder_runImportanceSamplingAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "printReliability",TclReliabilityModelBuilder_printReliability,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getMean",TclReliabilityModelBuilder_getMean,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getStdv",TclReliabilityModelBuilder_getStdv,(ClientData)NULL, NULL);
  //  Tcl_CreateCommand(interp, "rvReduction",TclReliabilityModelBuilder_rvReduction,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "betaFORM",TclReliabilityModelBuilder_getBetaFORM,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "gammaFORM",TclReliabilityModelBuilder_getGammaFORM,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "alphaFORM",TclReliabilityModelBuilder_getAlphaFORM,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getPDF",TclReliabilityModelBuilder_getPDF,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getStdNormalPDF",TclReliabilityModelBuilder_getStdNormalPDF,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getCDF",TclReliabilityModelBuilder_getCDF,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getStdNormalCDF",TclReliabilityModelBuilder_getStdNormalCDF,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getInverseCDF",TclReliabilityModelBuilder_getInverseCDF,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getStdNormalInverseCDF",TclReliabilityModelBuilder_getStdNormalInverseCDF,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getRVTags",TclReliabilityModelBuilder_getRVTags,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getLSFTags",TclReliabilityModelBuilder_getLSFTags,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getCutsetTags",TclReliabilityModelBuilder_getCutsetTags,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getCutsetComponents",TclReliabilityModelBuilder_getCutsetComponents,(ClientData)NULL, NULL);
/////////////////////////////////////////////////////////
///S added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////
  Tcl_CreateCommand(interp, "analyzer",TclReliabilityModelBuilder_addAnalyzer,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "initialstaticanalysis",TclReliabilityModelBuilder_addInitialStaticAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "initialpoint",TclReliabilityModelBuilder_addInitialPointBuilder,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "crossingrateanalyzer",TclReliabilityModelBuilder_addCrossingRateAnalyzer,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "foseries",TclReliabilityModelBuilder_addFOSeriesSimulation,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "firstpassage",TclReliabilityModelBuilder_addFirstPassageAnalyzer,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "randomvibrationsimulation",TclReliabilityModelBuilder_addRandomVibrationSimulation,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "randomvibrationanalysis",TclReliabilityModelBuilder_runRandomVibrationAnalysis,(ClientData)NULL, NULL);
/////////////////////////////////////////////////////////
///E added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////



  // Quan --

   Tcl_CreateCommand(interp, "runMonteCarloResponseAnalysis", TclReliabilityModelBuilder_runMonteCarloResponseAnalysis,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "updateParameterValue", TclReliabilityModelBuilder_updateParameterValue,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "runOrthogonalPlaneSamplingAnalysis",TclReliabilityModelBuilder_runOrthogonalPlaneSamplingAnalysis,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "computeHessian",TclReliabilityModelBuilder_computeHessian,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "runMultiDimVisualPrinPlane",TclReliabilityModelBuilder_MultiDimVisPrincPlane,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "transformXtoU",TclReliabilityModelBuilder_transformXtoU,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "transformUtoX",TclReliabilityModelBuilder_transformUtoX,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "runDP_RSM_SimTimeInvariantAnalysis",TclReliabilityModelBuilder_runDP_RSM_SimTimeInvariantAnalysis,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "runDP_RSM_SimTimeVariantAnalysis",TclReliabilityModelBuilder_runDP_RSM_SimTimeVariantAnalysis,(ClientData)NULL, NULL);


   //--Quan
	
	// set the static pointers in this file
  theStructuralDomain	= &passedDomain;
  theReliabilityDomain	= new ReliabilityDomain(theStructuralDomain);


}

TclReliabilityBuilder::~TclReliabilityBuilder()
{

  // Delete objects
  if (theReliabilityDomain != 0)
    delete theReliabilityDomain;

  // base class pointers
  if (theFunctionEvaluator != 0)
    delete theFunctionEvaluator;
  if (theGradientEvaluator != 0)
    delete theGradientEvaluator;
  if (theStepSizeRule != 0)
    delete theStepSizeRule;
  if (theSearchDirection != 0)
    delete theSearchDirection;
  if (theHessianEvaluator != 0)
    delete theHessianEvaluator;
  if (theMeritFunctionCheck != 0)
    delete theMeritFunctionCheck;
  if (theReliabilityConvergenceCheck != 0)
    delete theReliabilityConvergenceCheck;
  if (theProbabilityTransformation != 0)
    delete theProbabilityTransformation;
  if (theRootFindingAlgorithm != 0)
    delete theRootFindingAlgorithm;
  if (theRandomNumberGenerator != 0)
    delete theRandomNumberGenerator;
  if (theFindDesignPointAlgorithm != 0)
    delete theFindDesignPointAlgorithm;
  if (theFindCurvatures != 0)
    delete theFindCurvatures;
	
  // mixed pointers
  if (thePolakHeDualPurpose != 0)
    delete thePolakHeDualPurpose;
  if (theSQPtriplePurpose != 0)
    delete theSQPtriplePurpose;
	
  // analysis pointers
  if (theFORMAnalysis != 0)
    delete theFORMAnalysis;
  if (theFOSMAnalysis != 0)
    delete theFOSMAnalysis;
  //  if (theParametricReliabilityAnalysis != 0)
  //delete theParametricReliabilityAnalysis;
  if (theSORMAnalysis != 0)
    delete theSORMAnalysis;
  if (theImportanceSamplingAnalysis != 0)
    delete theImportanceSamplingAnalysis;
  if (theSystemAnalysis != 0)
    delete theSystemAnalysis;
  if (theGFunVisualizationAnalysis != 0)
    delete theGFunVisualizationAnalysis;
  if (theOutCrossingAnalysis != 0)
    delete theOutCrossingAnalysis;
  
  /////S added by K Fujimura /////
  if (theAnalyzer != 0)
    delete theAnalyzer;
  if (theInitialStaticAnalysis != 0)
    delete theInitialStaticAnalysis;
  if (theInitialPointBuilder != 0)
	delete theInitialPointBuilder;
  if (theCrossingRateAnalyzer != 0)
	delete theCrossingRateAnalyzer;
  if (theFOSeriesSimulation !=0)
	delete theFOSeriesSimulation ;
  if (theFirstPassageAnalyzer !=0)
	delete theFirstPassageAnalyzer ;
  if (theRandomVibrationSimulation !=0)
	delete theRandomVibrationSimulation ;
  if (theRandomVibrationAnalysis !=0)
	delete theRandomVibrationAnalysis ;
  /////E added by K Fujimura /////

  // Quan ---
  if (theMonteCarloResponseAnalysis != 0)
    delete theMonteCarloResponseAnalysis;
  if (theSamplingAnalysis != 0)
    delete theSamplingAnalysis;
  // ---Quan 

  theReliabilityDomain = 0;
  theFunctionEvaluator = 0;
  theGradientEvaluator = 0;
  theStepSizeRule = 0;
  theSearchDirection = 0;
  theHessianEvaluator = 0;
  theMeritFunctionCheck = 0;
  theReliabilityConvergenceCheck = 0;
  theProbabilityTransformation = 0;
  theRootFindingAlgorithm = 0;
  theRandomNumberGenerator = 0;
  theFindDesignPointAlgorithm = 0;
  theFindCurvatures = 0;
  
  thePolakHeDualPurpose =0;
  theSQPtriplePurpose =0;
  
  theFORMAnalysis = 0;
  theFOSMAnalysis = 0;
  //  theParametricReliabilityAnalysis = 0;
  theSORMAnalysis = 0;
  theImportanceSamplingAnalysis = 0;
  theSystemAnalysis = 0;
  theGFunVisualizationAnalysis = 0;
  theOutCrossingAnalysis = 0;
  
 /////S added by K Fujimura /////
  theAnalyzer=0;
  theInitialStaticAnalysis=0;
  theInitialPointBuilder = 0;
  theCrossingRateAnalyzer=0;
  theFOSeriesSimulation= 0;
  theFirstPassageAnalyzer= 0;
  theRandomVibrationSimulation= 0;
  theRandomVibrationAnalysis = 0;
 /////E added by K Fujimura /////
 
  // Quan ---
  theMonteCarloResponseAnalysis = 0;
  theSamplingAnalysis = 0;
  // ---Quan

  // Delete commands
  Tcl_DeleteCommand(theInterp, "randomVariable");
  Tcl_DeleteCommand(theInterp, "correlate");
  Tcl_DeleteCommand(theInterp, "correlateGroup");
  Tcl_DeleteCommand(theInterp, "correlationStructure");
  Tcl_DeleteCommand(theInterp, "cutset");
  Tcl_DeleteCommand(theInterp, "performanceFunction");
  Tcl_DeleteCommand(theInterp, "gradPerformanceFunction");
  Tcl_DeleteCommand(theInterp, "modulatingFunction");
  Tcl_DeleteCommand(theInterp, "filter");
  Tcl_DeleteCommand(theInterp, "spectrum");
  Tcl_DeleteCommand(theInterp, "findDesignPoint");
  Tcl_DeleteCommand(theInterp, "startPoint");
  Tcl_DeleteCommand(theInterp, "rootFinding");
  Tcl_DeleteCommand(theInterp, "gFunEvaluator");
  Tcl_DeleteCommand(theInterp, "gradGEvaluator");
  Tcl_DeleteCommand(theInterp, "stepSizeRule");
  Tcl_DeleteCommand(theInterp, "searchDirection");
  Tcl_DeleteCommand(theInterp, "hessianEvaluator");
  Tcl_DeleteCommand(theInterp, "meritFunctionCheck");
  Tcl_DeleteCommand(theInterp, "reliabilityConvergenceCheck");
  Tcl_DeleteCommand(theInterp, "probabilityTransformation");
  Tcl_DeleteCommand(theInterp, "findCurvatures");
  Tcl_DeleteCommand(theInterp, "randomNumberGenerator");
  Tcl_DeleteCommand(theInterp, "runFORMAnalysis");
  Tcl_DeleteCommand(theInterp, "runFOSMAnalysis");
  //  Tcl_DeleteCommand(theInterp, "runParametricReliabilityAnalysis");
  Tcl_DeleteCommand(theInterp, "runGFunVizAnalysis");
  Tcl_DeleteCommand(theInterp, "runOutCrossingAnalysis");
  Tcl_DeleteCommand(theInterp, "runSORMAnalysis");
  Tcl_DeleteCommand(theInterp, "runSystemAnalysis");
  Tcl_DeleteCommand(theInterp, "runImportanceSamplingAnalysis");
  Tcl_DeleteCommand(theInterp, "printReliability");
  Tcl_DeleteCommand(theInterp, "getMean");
  Tcl_DeleteCommand(theInterp, "getStdv");
  //  Tcl_DeleteCommand(theInterp, "rvReduction");
  Tcl_DeleteCommand(theInterp, "betaFORM");
  Tcl_DeleteCommand(theInterp, "gammaFORM");
  Tcl_DeleteCommand(theInterp, "getPDF");
  Tcl_DeleteCommand(theInterp, "getStdNormalPDF");
  Tcl_DeleteCommand(theInterp, "getCDF");
  Tcl_DeleteCommand(theInterp, "getStdNormalCDF");
  Tcl_DeleteCommand(theInterp, "getInverseCDF");
  Tcl_DeleteCommand(theInterp, "getStdNormalInverseCDF");
  Tcl_DeleteCommand(theInterp, "getRVTags");
  Tcl_DeleteCommand(theInterp, "getLSFTags");
  Tcl_DeleteCommand(theInterp, "getCutsetTags");
  Tcl_DeleteCommand(theInterp, "getCutsetComponents");

  /////S added by K Fujimura /////
  Tcl_DeleteCommand(theInterp, "analyzer");
  Tcl_DeleteCommand(theInterp, "initialstaticanalysis");
  Tcl_DeleteCommand(theInterp, "initialpoint");
  Tcl_DeleteCommand(theInterp, "crossingrateanalyzer");
  Tcl_DeleteCommand(theInterp, "foseries");
  Tcl_DeleteCommand(theInterp, "firstpassage");
  Tcl_DeleteCommand(theInterp, "randomvibrationsimulation");
  Tcl_DeleteCommand(theInterp, "randomvibrationanalysis");
  /////E added by K Fujimura /////

  Tcl_DeleteCommand(theInterp, "runMonteCarloResponseAnalysis");
  Tcl_DeleteCommand(theInterp, "updateParameterValue");
  Tcl_DeleteCommand(theInterp, "runOrthogonalPlaneSamplingAnalysis");
  Tcl_DeleteCommand(theInterp, "computeHessian");
  Tcl_DeleteCommand(theInterp, "runMultiDimVisualPrinPlane");
  Tcl_DeleteCommand(theInterp, "transformXtoU");
  Tcl_DeleteCommand(theInterp, "transformUtoX");
  Tcl_DeleteCommand(theInterp, "runDP_RSM_SimTimeInvariantAnalysis");
  Tcl_DeleteCommand(theInterp, "runDP_RSM_SimTimeVariantAnalysis");
}


//
// CLASS METHODS
//


ReliabilityDomain *
TclReliabilityBuilder::getReliabilityDomain()
{
	return theReliabilityDomain;
}


int 
inputCheck()
{
	// Check that tagged objects are consecutive
	int i, num;
	ReliabilityDomainComponent *component;
    
	// Clear out old parameter positioners so we don't produce a memory leak
	/*
     theReliabilityDomain->removeAllParameterPositioners();
     
     ParameterIter &paramIter = theStructuralDomain->getParameters();
     Parameter *theParam;
     i = 1;
     while ((theParam = paramIter()) != 0) {
     ParameterPositioner *theParamPos = 
     new ParameterPositioner(i, *theParam);
     theParamPos->setGradNumber(i);
     if (theReliabilityDomain->addParameterPositioner(theParamPos) == false) {
     opserr << "ERROR: failed to add parameter positioner " << i << endln;
     delete theParamPos; // otherwise memory leak
     return TCL_ERROR;
     }
     i++;
     }
     */
	/*
     num = theReliabilityDomain->getNumberOfRandomVariablePositioners();
     for (i=1; i<=num; i++) {
     component = theReliabilityDomain->getRandomVariablePositionerPtr(i);
     if (component == 0) {
     opserr << "ERROR: Non-consecutive random variable positioner list." << endln;
     return TCL_ERROR;
     }
     }
     */

	num = theReliabilityDomain->getNumberOfFilters();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getFilter(i);
		if (component == 0) {
			opserr << "ERROR: Non-consecutive filter list." << endln;
			return -1;
		}
	}
	
	num = theReliabilityDomain->getNumberOfModulatingFunctions();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getModulatingFunction(i);
		if (component == 0) {
			opserr << "ERROR: Non-consecutive modulating function list." << endln;
			return -1;
		}
	}
	
	num = theReliabilityDomain->getNumberOfSpectra();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getSpectrum(i);
		if (component == 0) {
			opserr << "ERROR: Non-consecutive spectrum list." << endln;
			return -1;
		}
	}
    
	// Check that the correlation matrix is positive definite
	// theCorrelationMatrix
    
    // set defaults
    if (theProbabilityTransformation == 0) {
        opserr << "No probabilityTransformation specified, assuming AllIndependent" << endln;
        theProbabilityTransformation = new AllIndependentTransformation(theReliabilityDomain,0);
    }
    
    //reliabilityConvergenceCheck  Standard         -e1 1.0e-3    -e2 1.0e-3  -print 1
    //functionEvaluator            Tcl
    //gradientEvaluator            FiniteDifference -pert 1000
    
    if (theSearchDirection == 0) {
        opserr << "No searchDirectin specified, assuming Standard" << endln;
        theSearchDirection = new HLRFSearchDirection();
    }

   
    //meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0   -factor 0.5
    //stepSizeRule                 Armijo           -maxNum 5    -base 0.5   -initial 1.0 2  -print 0
    //startPoint                   Mean
    //findDesignPoint              StepSearch       -maxNumIter 30   -printDesignPointX CalRel_manual_1_output/1_designX.out
    //randomNumberGenerator        CStdLib
    
	return 0;
}


//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addRandomVariable(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{
	RandomVariable *theRandomVariable = 0;
	int tag;
	double mean = 0;
	double stdv = 0;
	double startPt = 0;
    int use_start_pt = 0;
	
	int param_indx = 0;
	double param = 0;
	Vector param_temp(4);

	// CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
	if (argc < 5) {
		opserr << "ERROR: invalid number of arguments to randomVariable command \n";
		return TCL_ERROR;
	}

	// GET TAG NUMBER
	if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
		opserr << "ERROR: invalid input: tag \n";
		return TCL_ERROR;
	}
	
	int argi = 3;
	while (argi < argc) {
		// user specified mean directly
		if (strcmp(argv[argi],"-mean") == 0) {
			if (argc < argi+2) {
				opserr << "WARNING not enough args, need -mean mean??\n";
				opserr << argv[1] << " for random variable: " << tag << endln;
				return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[argi+1], &mean) != TCL_OK) {
				opserr << "WARNING invalid mean\n";
				opserr << argv[1] << " for random variable: " << tag << endln;
				return TCL_ERROR;
			}
			argi += 2;
		}
		
		// user specified standard deviation directly
		else if (strcmp(argv[argi],"-stdv") == 0) {
			if (argc < argi+2) {
				opserr << "WARNING not enough args, need -stdv stdv??\n";
				opserr << argv[1] << " for random variable: " << tag << endln;
				return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[argi+1], &stdv) != TCL_OK) {
				opserr << "WARNING invalid standard deviation\n";
				opserr << argv[1] << " for random variable: " << tag << endln;
				return TCL_ERROR;
			}
			argi += 2;
		}
		
		// user specified starting point directly
		else if (strcmp(argv[argi],"-startPoint") == 0) {
			if (argc < argi+2) {
				opserr << "WARNING not enough args, need -startPoint startPt??\n";
				opserr << argv[1] << " for random variable: " << tag << endln;
				return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[argi+1], &startPt) != TCL_OK) {
				opserr << "WARNING invalid starting point\n";
				opserr << argv[1] << " for random variable: " << tag << endln;
				return TCL_ERROR;
			}
			argi += 2;
            use_start_pt = 1;
		}
		
		// user input distribution specific parameters directly
		else if (strcmp(argv[argi],"-parameters") == 0) {
			if (argc < argi+2) {
				opserr << "WARNING not enough args, need -parameters param1 ...??\n";
				opserr << argv[1] << " for random variable: " << tag << endln;
				return TCL_ERROR;
			}
			
			argi++;
			int err_break = 0;
			while (argi < argc && err_break == 0) {
				if (Tcl_GetDouble(interp, argv[argi], &param) != TCL_OK)
					err_break = 1;
				
				param_temp(param_indx) = param;
				param_indx++;
				argi++;
			}
		}
		
		// otherwise skip argument
		else
			argi++;
    }
    
    // resize parameter vector based on actual inputs
    Vector parameters;
    if (param_indx > 0) {
        parameters.resize(param_indx);
        for (int kl = 0; kl < param_indx; kl++)
            parameters(kl) = param_temp(kl);
    }
    

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[2],"normal") == 0) {
		if (param_indx > 0)
			theRandomVariable = new NormalRV(tag, parameters);
		else
			theRandomVariable = new NormalRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"lognormal") == 0) {
		if (param_indx > 0)
			theRandomVariable = new LognormalRV(tag, parameters);
		else
			theRandomVariable = new LognormalRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"gamma") == 0) {
		if (param_indx > 0)
			theRandomVariable = new GammaRV(tag, parameters);
		else
			theRandomVariable = new GammaRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"shiftedExponential") == 0) {
		if (param_indx > 0)
			theRandomVariable = new ShiftedExponentialRV(tag, parameters);
		else
			theRandomVariable = new ShiftedExponentialRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"shiftedRayleigh") == 0) {
		if (param_indx > 0)
			theRandomVariable = new ShiftedRayleighRV(tag, parameters);
		else
			theRandomVariable = new ShiftedRayleighRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"exponential") == 0) {
		if (param_indx > 0)
			theRandomVariable = new ExponentialRV(tag, parameters);
		else
			theRandomVariable = new ExponentialRV(tag, mean, stdv);
	}

	else if (strcmp(argv[2],"rayleigh") == 0) {
		if (param_indx > 0)
			theRandomVariable = new RayleighRV(tag, parameters);
		else {
			opserr << "Rayleigh random variable with tag " << tag << " cannot be created with only mean/stdv." << endln;
			return TCL_ERROR;
		}
	}
	
	else if (strcmp(argv[2],"uniform") == 0) {
		if (param_indx > 0)
			theRandomVariable = new UniformRV(tag, parameters);
		else
			theRandomVariable = new UniformRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"beta") == 0) {
		if (param_indx > 0)
			theRandomVariable = new BetaRV(tag, parameters);
		else {
			opserr << "Beta random variable with tag " << tag << " cannot be created with only mean/stdv." << endln;
			return TCL_ERROR;
		}
	}
	
	else if (strcmp(argv[2],"type1LargestValue") == 0) {
		if (param_indx > 0)
			theRandomVariable = new Type1LargestValueRV(tag, parameters);
		else
			theRandomVariable = new Type1LargestValueRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"type1SmallestValue") == 0) {
		if (param_indx > 0)
			theRandomVariable = new Type1SmallestValueRV(tag, parameters);
		else
			theRandomVariable = new Type1SmallestValueRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"type2LargestValue") == 0) {
		if (param_indx > 0)
			theRandomVariable = new Type2LargestValueRV(tag, parameters);
		else
			theRandomVariable = new Type2LargestValueRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"type3SmallestValue") == 0) {
		if (param_indx > 0)
			theRandomVariable = new Type3SmallestValueRV(tag, parameters);
		else {
			opserr << "T3S random variable with tag " << tag << " cannot be created with only mean/stdv." << endln;
			return TCL_ERROR;
		}
	}
		
	else if (strcmp(argv[2],"chiSquare") == 0) {
		if (param_indx > 0)
			theRandomVariable = new ChiSquareRV(tag, parameters);
		else
			theRandomVariable = new ChiSquareRV(tag, mean, stdv);
	}
		
	else if (strcmp(argv[2],"gumbel") == 0) {
		if (param_indx > 0)
			theRandomVariable = new GumbelRV(tag, parameters);
		else
			theRandomVariable = new GumbelRV(tag, mean, stdv);
	}

	else if (strcmp(argv[2],"weibull") == 0) {
		if (param_indx > 0)
			theRandomVariable = new WeibullRV(tag, parameters);
		else
			theRandomVariable = new WeibullRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"laplace") == 0) {
		if (param_indx > 0)
			theRandomVariable = new LaplaceRV(tag, parameters);
		else
			theRandomVariable = new LaplaceRV(tag, mean, stdv);
	}
	
	else if (strcmp(argv[2],"pareto") == 0) {
		if (param_indx > 0)
			theRandomVariable = new ParetoRV(tag, parameters);
		else {
			opserr << "Pareto random variable with tag " << tag << " cannot be created with only mean/stdv." << endln;
			return TCL_ERROR;
		}
	}

	else if (strcmp(argv[2],"userdefined") == 0) {
		// note userdefined is a special case and will not have any input read from the command line yet
		// unless user defined mean and standard deviation for some reason, which will break this input
		// because we assume argi starts at 3 here
        
        // KRM 4/22/2012 userdefined currently not implemented.
		Vector xPoints;
		Vector PDFpoints;
		int numPoints = 0;
		
		if (strcmp(argv[3],"-list") == 0) {
			
			numPoints = (argc-4) % 2;
			Vector temp_xPoints(numPoints);
			Vector temp_PDFpoints(numPoints);
			
			double x = 0.0;
			double pdf = 0.0;
			double x_old = 0.0;
			
			// Read the points
			for (int i=0; i < numPoints; i++) {
				if (Tcl_GetDouble(interp, argv[4+2*i], &x) != TCL_OK) {
					opserr << "ERROR: Invalid x point to user-defined random variable." << endln;
					return TCL_ERROR;
				}
				if (Tcl_GetDouble(interp, argv[5+2*i], &pdf) != TCL_OK) {
					opserr << "ERROR: Invalid PDF value point to user-defined random variable." << endln;
					return TCL_ERROR;
				}
				if (i>0 && x<=x_old) {
					opserr << "ERROR: x-points to user-defined random variable must be consecutive!" << endln;
					return TCL_ERROR;
				}
				temp_xPoints(i) = x;
				temp_PDFpoints(i) = pdf;
				x_old = x;
			}
			
			xPoints = temp_xPoints;
			PDFpoints = temp_PDFpoints;
			
		}
		else if (strcmp(argv[3],"-file") == 0) {
			
			// Open file where the vectors are given
			ifstream inputFile( argv[4], ios::in );
			if (inputFile.fail()) {
				opserr << "File " << argv[4] << " could not be opened. " << endln;
				return TCL_ERROR;
			}
			
			// Loop through file to see how many entries there are
			double dummy;
			numPoints = 0;
			while (inputFile >> dummy) {
				inputFile >> dummy;
				numPoints++;
			}
			if (numPoints == 0) {
				opserr << "ERROR: No entries in the direction file read by " << endln
					<< "user-defined random variable, number " << tag << endln;
				return TCL_ERROR;
			}
			
			// rewind
			inputFile.clear(); 
			inputFile.seekg(0); 
			
			// Allocate vectors of correct size
			Vector temp_xPoints(numPoints);
			Vector temp_PDFpoints(numPoints);
						
			// Store the vector
			for (int i=0; i<numPoints; i++) {
				inputFile >> temp_xPoints(i);
				inputFile >> temp_PDFpoints(i);
			}
			inputFile.close();
			
			xPoints = temp_xPoints;
			PDFpoints = temp_PDFpoints;
		}
		else {
			opserr << "ERROR: Invalid argument to user-defined random variable, number " << tag << endln;
			return TCL_ERROR;
		}
		
		//theRandomVariable = new UserDefinedRV(tag, xPoints, PDFpoints);
		
	}
	
	else {
		opserr << "ERROR: unknown random variable type: " << argv[2] << " provided. Must be one of " << endln;
		return TCL_ERROR;
	}

	if (theRandomVariable == 0) {
		opserr << "ERROR: could not create random variable number " << tag << endln;
		return TCL_ERROR;
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
	if (theReliabilityDomain->addRandomVariable(theRandomVariable) == false) {
		opserr << "ERROR: failed to add random variable to the domain (wrong number of arguments?)\n";
		opserr << "random variable: " << tag << endln;
		delete theRandomVariable; // otherwise memory leak
		return TCL_ERROR;
	}

	//RVParameter *theRVParam = new RVParameter(tag, theRandomVariable);
	//theStructuralDomain->addParameter(theRVParam);

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

	char buffer[40];
	sprintf(buffer,"%35.20f",rv->getMean());
	
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);

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
	
	char buffer[40];
	sprintf(buffer,"%35.20f",rv->getStdv());
	
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);

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
TclReliabilityModelBuilder_addCutset(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Cutset *theCutset = 0;
  int tag, set;

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }
  
  Vector components(argc-2);
  int argi = 2;
  while (argi < argc) {
     if (Tcl_GetInt(interp, argv[argi], &set) != TCL_OK) {
	    opserr << "ERROR: invalid input: cutset components \n";
	    return TCL_ERROR;
     }
	 
	 // check LSF exists in domain
	 if (theReliabilityDomain->getLimitStateFunctionIndex( abs(set) ) < 0) {
	    opserr << "ERROR: LSF does not exist \n";
	    return TCL_ERROR;
     }
	 
	 components(argi-2) = set;
	 argi++;
  }
  
  // CREATE THE OBJECT
  theCutset = new Cutset(tag, components);
  if (theCutset == 0) {
	opserr << "ERROR: ran out of memory creating cutset \n";
	opserr << "cutset: " << tag << endln;
	return TCL_ERROR;
  }

  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addCutset(theCutset) == false) {
	opserr << "ERROR: failed to add cutset to the domain\n";
	opserr << "cutset: " << tag << endln;
	delete theCutset; // otherwise memory leak
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

	// Assume that previous corr. coeffs. have been added in order
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
	// DO NOT TYPE STUFF INTO THE INTERPRETER!
	if (strcmp(argv[1],"homogeneous1") == 0) {
		for (int i=firstRV; i<=lastRV; i++) {
			for (int j=i+1; j<=lastRV; j++) {
				correlationValue = exp(-abs(i-j)/theta);
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
				correlationValue = 1.0/(theta*(i-j)*(i-j));
				sprintf(theCorrelateCommand,"correlate %d %d %10.5f",i,j,correlationValue);
				Tcl_Eval(interp, theCorrelateCommand );
			}
		}
	}
	else if (strcmp(argv[1],"homogeneous4") == 0) {
		for (int i=firstRV; i<=lastRV; i++) {
			for (int j=i+1; j<=lastRV; j++) {
				if (abs(i-j)<theta) {
					correlationValue = 1.0-(abs(i-j)/theta);
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

  if (theFunctionEvaluator != 0 && argc > 2) {
    opserr << "ERROR: A limit-state function should not be created after the GFunEvaluator has been instantiated." << endln;
    return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }

  if (argc == 2) {
    theReliabilityDomain->setTagOfActiveLimitStateFunction(tag);
    //theFunctionEvaluator->evaluateG(x);
    double g = theFunctionEvaluator->evaluateExpression();
    //opserr << "TclReliabilityBuilder: " << g << endln;
    char buffer[40];  
    sprintf(buffer,"%35.20f",g);
    
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
   
    return TCL_OK;
  }

  // CREATE THE OBJECT
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


}

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addGradLimitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  LimitStateFunction *theLimitStateFunction = 0;
  int lsfTag, rvTag;

  if (theFunctionEvaluator != 0 ) {
    //opserr << "ERROR: A limit-state function should not be created after the FunctionEvaluator has been instantiated." << endln;
    //return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[1], &lsfTag) != TCL_OK) {
    opserr << "ERROR: invalid input: lsfTag \n";
    return TCL_ERROR;
  }

  // GET INPUT PARAMETER (integer)
  if (Tcl_GetInt(interp, argv[2], &rvTag) != TCL_OK) {
    opserr << "ERROR: invalid input: rvTag \n";
    return TCL_ERROR;
  }
    
  // GET LSF pointer
  theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);
  if (theLimitStateFunction == 0) {
    opserr << "ERROR: limit state function with tag " << lsfTag
	   << " does not exist" << endln;
    return TCL_ERROR;
  }

  // ADD THE OBJECT TO THE LSF
  int ok = theLimitStateFunction->addGradientExpression(argv[3], rvTag);
  if (ok < 0) {
    opserr << "ERROR: could not add gradient of LSF " << lsfTag
	   << " for random variable " << rvTag << endln;
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

		if (argc!=7) {
			opserr << "ERROR: Incorrect number of arguments to gamma modulating function" << endln;
			return TCL_ERROR;
		}

		int thisTag, filterTag;
		double a,b,c;

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
		if (Tcl_GetDouble(interp, argv[4], &a) != TCL_OK) {
			opserr << "ERROR: invalid input: a \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
			opserr << "ERROR: invalid input: b \n";
			return TCL_ERROR;
		}

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[6], &c) != TCL_OK) {
			opserr << "ERROR: invalid input: c \n";
			return TCL_ERROR;
		}

		// CREATE THE OBJECT
		theModulatingFunction = new GammaModulatingFunction(thisTag,theFilter,a,b,c);

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
		double amplitude=0.0;

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
		if (Tcl_GetDouble(interp, argv[4], &amplitude) != TCL_OK) {
			opserr << "ERROR: invalid input: amplitude \n";
			return TCL_ERROR;
		}

		// CREATE THE OBJECT
		theModulatingFunction = new ConstantModulatingFunction(thisTag,theFilter,amplitude);

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

		double t1, t2, t3, t4, amplitude;

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

		// GET INPUT PARAMETER (double)
		if (Tcl_GetDouble(interp, argv[8], &amplitude) != TCL_OK) {
			opserr << "ERROR: invalid input to modulating function: amplitude \n";
			return TCL_ERROR;
		}

		// CREATE THE OBJECT
		theModulatingFunction = new TrapezoidalModulatingFunction(thisTag,theFilter,t1,t2,t3,t4,amplitude);

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
	}	else if (strcmp(argv[2],"Koo") == 0) {
	
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

		double t1, t2;

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


		// CREATE THE OBJECT
		theModulatingFunction = new KooModulatingFunction(thisTag,theFilter,t1,t2);

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
	int tag;
	double period_Tn, damping, dtpulse;

	if ( (strcmp(argv[2],"standard") == 0) || (strcmp(argv[2],"standardDisplacement") == 0) 
		|| ( strcmp(argv[2],"Koo") == 0 ) || (strcmp(argv[2],"standardVelocity")==0)
		|| (strcmp(argv[2],"standardAcceleration")==0)) {	
		if (argc != 5) {
			opserr << "ERROR: Wrong number of arguments to filter command." << endln;
			return TCL_ERROR;
		}
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
		if ( (strcmp(argv[2],"standard") == 0) || (strcmp(argv[2],"standardDisplacement") == 0) ) {
			theFilter = new StandardLinearOscillatorDisplacementFilter(tag,period_Tn,damping);
		}else if ( strcmp(argv[2],"Koo") == 0 ) {
			theFilter = new KooFilter(tag,period_Tn,damping);
		}else if (strcmp(argv[2],"standardVelocity") == 0) {
			theFilter = new StandardLinearOscillatorVelocityFilter(tag,period_Tn,damping);
		}else if (strcmp(argv[2],"standardAcceleration") == 0) {
			theFilter = new StandardLinearOscillatorAccelerationFilter(tag,period_Tn,damping);
		}else {
			opserr << "ERROR:: Unknown type of filter. " << endln;
			return TCL_ERROR;
		}
	}else if ((strcmp(argv[2],"whitenoise") == 0) 
		|| (strcmp(argv[2],"NewStandardLinearOscillatorAcceleration") == 0) ){
		if ( strcmp(argv[2],"whitenoise") == 0){
			if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
				opserr << "ERROR: invalid input: tag \n";
				return TCL_ERROR;
			}
			// GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[3], &period_Tn) != TCL_OK) {
				opserr << "ERROR: invalid input: freq_wn \n";
				return TCL_ERROR;
			}
			theFilter = new NewWhitenoiseFilter(tag,period_Tn);
		}else{
			if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
				opserr << "ERROR: invalid input: tag \n";
				return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[3], &period_Tn) != TCL_OK) {
				opserr << "ERROR: invalid input: freq_wn \n";
				return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[4], &damping) != TCL_OK) {
				opserr << "ERROR: invalid input: damping \n";
				return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[4], &dtpulse) != TCL_OK) {
				opserr << "ERROR: invalid input: damping \n";
				return TCL_ERROR;
			}
			theFilter = new NewStandardLinearOscillatorAccelerationFilter
				(tag, period_Tn, damping, dtpulse);
		}
	}else{
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

/*
int 
TclReliabilityModelBuilder_addFilter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	Filter *theFilter = 0;

	int tag;
	double period_Tn, damping;

	// GET INPUT PARAMETER (integer)
	if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
		opserr << "ERROR: invalid input: tag \n";
		return TCL_ERROR;
	}
// Quan and Michele
	if ( (strcmp(argv[2],"delta") == 0) || (strcmp(argv[2],"Delta") == 0) ) {

			theFilter = new DeltaFilter(tag);
			
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

	}    // if "constant"


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



	if ( (strcmp(argv[2],"standard") == 0) || (strcmp(argv[2],"standardDisplacement") == 0) ) {

		theFilter = new StandardLinearOscillatorDisplacementFilter(tag,period_Tn,damping);
	}
	else if ( strcmp(argv[2],"Koo") == 0 ) {

		theFilter = new KooFilter(tag,period_Tn,damping);
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
*/

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
    
	else if (strcmp(argv[1],"AllIndependent") == 0) {

		int printFlag = 0; 
		if (argc > 2) {
			if (strcmp(argv[2],"-print") == 0) {

				if (Tcl_GetInt(interp, argv[3], &printFlag) != TCL_OK) {
					opserr << "ERROR: invalid input: printFlag to AllIndependent transformation \n";
					return TCL_ERROR;
				}
			}
		}

		theProbabilityTransformation = new AllIndependentTransformation(theReliabilityDomain,printFlag);
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

		double gamma = 1.0;
		double delta = 1.0;

		int argvCounter = 2;
		while (argc > argvCounter) {
			if (strcmp(argv[argvCounter],"-gamma") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &gamma) != TCL_OK) {
					opserr << "ERROR: invalid input: gamma for Polak-He algorithm" << endln;
					return TCL_ERROR;
				}
				argvCounter++;

			}
			else if (strcmp(argv[argvCounter],"-delta") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &delta) != TCL_OK) {
					opserr << "ERROR: invalid input: delta for Polak-He algorithm" << endln;
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to Polak-He algorithm." << endln;
				return TCL_ERROR;
			}
		}


		thePolakHeDualPurpose = new PolakHeSearchDirectionAndMeritFunction(gamma,delta);
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
		  opserr << "Assume all RV's are independent" << endln;
		  theProbabilityTransformation = 
		    new AllIndependentTransformation(theReliabilityDomain,0);
		}

		// Check that a gfun evaluator has been created
		if (theFunctionEvaluator == 0 ) {
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
									   theFunctionEvaluator,
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
        // KRM 5-19-2012 all these mixed classes are going to need to change, hessian needs to go in 
        // hessian, merit function needs to go in merit function, etc.
		//theHessianEvaluator = theSQPtriplePurpose;

		// Set the Hessian approximation in the search direction
		//theSQPtriplePurpose->setHessianApproximation(theHessianEvaluator);

	
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
TclReliabilityModelBuilder_addHessianEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theHessianEvaluator != 0) {
		delete theHessianEvaluator;
		theHessianEvaluator = 0;
	}
    
    
	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"FiniteDifference") == 0) {
        
		double perturbationFactor = 1000.0;
		bool doGradientCheck = false;
        
		// Check that the necessary ingredients are present
		if (theFunctionEvaluator == 0 ) {
			opserr << "Need FunctionEvaluator before a FiniteDifferenceHessian can be created" << endln;
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
				opserr << "ERROR: Wrong number of arguments to FiniteDifferenceGradient. " << endln;
				return TCL_ERROR;
			}
            
			int counter = 2;
            
			for (int i=1; i<=numExtras; i++) {
                
				if (strcmp(argv[counter],"-pert") == 0) {
					counter ++;
                    
					if (Tcl_GetDouble(interp, argv[counter], &perturbationFactor) != TCL_OK) {
						opserr << "ERROR: invalid input: perturbationFactor \n";
						return TCL_ERROR;
					}
					counter++;
				}
				else if (strcmp(argv[counter],"-check") == 0) {
					counter++;
					doGradientCheck = true;
				}
				else {
					opserr << "ERROR: Error in input to FiniteDifferenceGradient. " << endln;
					return TCL_ERROR;
				}
			}
		}
        
		theHessianEvaluator = new FiniteDifferenceHessian(theFunctionEvaluator, theReliabilityDomain, 
                                                            theStructuralDomain);
	}

	else if (strcmp(argv[1],"SQP_BFGS") == 0) {

		// Check that the SQP search direction is already created
		if (theSQPtriplePurpose == 0 ) {
			opserr << "Need theSQPSearchDirection before a SQP Hessian Approximation can be created" << endln;
			return TCL_ERROR;
		}

        // KRM 5-19-2012 all these mixed classes are going to need to change, hessian needs to go in 
        // hessian, merit function needs to go in merit function, etc.
		//theHessianEvaluator = theSQPtriplePurpose;

		// Set the Hessian approximation in the search direction
		// (this needs to be changed for generatlity; new method of search direction)
		//theSQPtriplePurpose->setHessianApproximation(theHessianEvaluator);

	}
    
	else {
		opserr << "ERROR: unrecognized type of HessianEvaluator \n";
		return TCL_ERROR;
	}

	if (theHessianEvaluator == 0) {
		opserr << "ERROR: could not create theHessianEvaluator \n";
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
		double scaleValue = 0.0;
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
			else if (strcmp(argv[argvCounter],"-scaleValue") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &scaleValue) != TCL_OK) {
					opserr << "ERROR: invalid input: scaleValue \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to standard reliability convergence check. " << endln;
				return TCL_ERROR;
			}
		}
			theReliabilityConvergenceCheck = new StandardReliabilityConvergenceCheck(e1,e2,scaleValue,print);
	}
	
    else if (strcmp(argv[argvCounter],"OptimalityCondition") == 0) {
		argvCounter++;

		double e1 = 1.0e-3;
		double e2 = 1.0e-3;
		double scaleValue = 0.0;
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
			else if (strcmp(argv[argvCounter],"-scaleValue") == 0) {
				argvCounter++;

				if (Tcl_GetDouble(interp, argv[argvCounter], &scaleValue) != TCL_OK) {
					opserr << "ERROR: invalid input: scaleValue \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to standard reliability convergence check. " << endln;
				return TCL_ERROR;
			}
		}
		theReliabilityConvergenceCheck = new OptimalityConditionReliabilityConvergenceCheck(e1,e2,scaleValue,print);
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
		if (theFunctionEvaluator == 0 ) {
			opserr << "Need theGFunEvaluator before an ArmijoStepSizeRule can be created" << endln;
			return TCL_ERROR;
		}
		if (theProbabilityTransformation == 0 ) {
//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
			opserr << "Assume all RV's are independent" << endln;
			theProbabilityTransformation = 
			new AllIndependentTransformation(theReliabilityDomain,0);
//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
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
		
		double radius = 50.0;
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

		theStepSizeRule = new ArmijoStepSizeRule(theReliabilityDomain, theFunctionEvaluator,
							 theProbabilityTransformation, theMeritFunctionCheck,
							 theRootFindingAlgorithm, 
							 base, maxNumReductions, b0,
							 numberOfShortSteps, radius, surfaceDistance,
							 evolution, printFlag);
		

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
TclReliabilityModelBuilder_addFunctionEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theFunctionEvaluator != 0) {
		delete theFunctionEvaluator;
		theFunctionEvaluator = 0;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"Matlab") == 0) {
		opserr << "ERROR: The Matlab function evaluator is not implemented in this " << endln
			<< " version of your OpenSees executable file. Please contact the " << endln
			<< " developer for more information." << endln;
		return TCL_ERROR;
	}
	else if (strcmp(argv[1],"Tcl") == 0) {
        if (argc == 2)
            theFunctionEvaluator = new TclEvaluator(interp, theReliabilityDomain, theStructuralDomain);

		else if (argc == 4) {
            if (strcmp(argv[2],"-file") != 0 && strcmp(argv[2],"-command") != 0) {
                opserr << "ERROR: Tcl function evaluator only takes -file and -command arguments." << endln;
                return TCL_ERROR;    
            } else {
                theFunctionEvaluator = new TclEvaluator(interp, theReliabilityDomain, 
							theStructuralDomain, argv[3]);
            }
        }
        
        else {
            opserr << "ERROR: Wrong input to Tcl function evaluator." << endln;
			return TCL_ERROR;
		}
	}

/////////////////////////////////////////
////////S modified by K Fujimura 10/10/2004
/////////////////////////////////////////
	else if (strcmp(argv[1],"Analyzer") == 0) {

		// There are several alternatives for this command:
		// gFunEvaluator  OpenSees  -file <filename>
		// gFunEvaluator  OpenSees  -runToMaxTimeInGFun
		// gFunEvaluator  OpenSees  -analyze <numSteps> <dt(optional)>

		if (theAnalyzer == 0 ){
			opserr << "Fatalerror \n";
			opserr << "Analyzer must be defined before \n";
			opserr << "AnalyzerGfunEvaluator \n";
			return TCL_ERROR;
		}
		if (theReliabilityDomain == 0 ){
			opserr << "Fatalerror \n";
			opserr << "theReliabilityDomain must be defined before \n";
			opserr << "AnalyzerGfunEvaluator \n";
			return TCL_ERROR;
		}
		if (theReliabilityDomain == 0 ){
			opserr << "Fatalerror \n";
			opserr << "theStructuralDomain must be defined before \n";
			opserr << "AnalyzerGfunEvaluator \n";
			return TCL_ERROR;
		}
		//theFunctionEvaluator = new AnalyzerGFunEvaluator(interp, theReliabilityDomain,
		//theStructuralDomain,theAnalyzer);
	}
/////////////////////////////////////////
////////E modified by K Fujimura 10/10/2004
/////////////////////////////////////////
	
	else {
		opserr << "ERROR: unrecognized type of FunctionEvaluator \n";
		return TCL_ERROR;
	}
	
	if (theFunctionEvaluator == 0) {
		opserr << "ERROR: could not create the theFunctionEvaluator \n";
		return TCL_ERROR;
	}
	return TCL_OK;
}



//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addGradientEvaluator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theGradientEvaluator != 0) {
		delete theGradientEvaluator;
		theGradientEvaluator = 0;
	}


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"FiniteDifference") == 0) {

		double perturbationFactor = 1000.0;
		bool doGradientCheck = false;

		// Check that the necessary ingredients are present
		if (theFunctionEvaluator == 0 ) {
			opserr << "Need FunctionEvaluator before a FiniteDifferenceGradient can be created" << endln;
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
				opserr << "ERROR: Wrong number of arguments to FiniteDifferenceGradient. " << endln;
				return TCL_ERROR;
			}

			int counter = 2;

			for (int i=1; i<=numExtras; i++) {

				if (strcmp(argv[counter],"-pert") == 0) {
					counter ++;

					if (Tcl_GetDouble(interp, argv[counter], &perturbationFactor) != TCL_OK) {
						opserr << "ERROR: invalid input: perturbationFactor \n";
						return TCL_ERROR;
					}
					counter++;
				}
				else if (strcmp(argv[counter],"-check") == 0) {
					counter++;
					doGradientCheck = true;
				}
				else {
					opserr << "ERROR: Error in input to FiniteDifferenceGradient. " << endln;
					return TCL_ERROR;
				}
			}
		}

		theGradientEvaluator = new FiniteDifferenceGradient(theFunctionEvaluator, theReliabilityDomain, 
								    theStructuralDomain);
	}

	else if (strcmp(argv[1],"OpenSees") == 0 || strcmp(argv[1],"Implicit") == 0) {

		bool doGradientCheck = false;

	//Quan Apr. 2006	
	//	if (theSensitivityAlgorithm == 0) {
	//		opserr << "Warning:Need a DDM sensitivity algorithm before a OpenSees sensitivity evaluator can be created" << endln;
	//		return TCL_ERROR;
	//	}

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

		theGradientEvaluator = new ImplicitGradient(theFunctionEvaluator, 
					theReliabilityDomain, theStructuralDomain, theSensitivityAlgorithm);
	
	}
	////////////////////////////////////////
	//////S modified by K Fujimura 10/10/2004
	////////////////////////////////////////
	else if (strcmp(argv[1],"Analyzer") == 0) {

		bool doGradientCheck = false;

		if (theAnalyzer == 0) {
			opserr << "Need Analyzer before a Analyzer sensitivity evaluator can be created" << endln;
			return TCL_ERROR;
		}
		if (theSensitivityAlgorithm == 0) {
			opserr << "Need a DDM sensitivity algorithm before a Analyzer sensitivity evaluator can be created" << endln;
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
		/*
		theGradientEvaluator = new AnalyzerGradGEvaluator(interp, theReliabilityDomain,
							 theStructuralDomain, theFunctionEvaluator,
							 doGradientCheck);
		*/
	}
	////////////////////////////////////////
	//////E modified by K Fujimura 10/10/2004
	////////////////////////////////////////


	else {
		opserr << "ERROR: unrecognized type of gradientEvaluator \n";
		return TCL_ERROR;
	}

	if (theGradientEvaluator == 0) {
		opserr << "ERROR: could not create theGradientEvaluator \n";
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
    
    if (theFunctionEvaluator == 0 ) {
        opserr << "Need theFunctionEvaluator before a findCurvatures can be created" << endln;
        return TCL_ERROR;
    }
    if (theFORMAnalysis == 0 ) {
		opserr << "Need theFORMAnalysis before a findCurvatures can be created" << endln;
        return TCL_ERROR;
	}

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"firstPrincipal") == 0) {
        // KRM 5-22-2012
        // combining former firstPrincipal and bySearchAlgorithm because they are doing same thing
		theFindCurvatures = new FirstPrincipalCurvature(theReliabilityDomain, theFunctionEvaluator,
                                                        theFORMAnalysis);
		
	}
	else if (strcmp(argv[1],"bySearchAlgorithm") == 0) {

		int numberOfCurvatures;

		// GET INPUT PARAMETER (integer)
		if (Tcl_GetInt(interp, argv[2], &numberOfCurvatures) != TCL_OK) {
			opserr << "ERROR: invalid input: numberOfCurvatures \n";
			return TCL_ERROR;
		}

		theFindCurvatures = new CurvaturesBySearchAlgorithm(theReliabilityDomain, theFunctionEvaluator,
                                                            theFORMAnalysis, numberOfCurvatures);
	}
    else if (strcmp(argv[1],"curvatureFitting") == 0) {
        // needs Hessian
        if (theHessianEvaluator == 0 ) {
            opserr << "Need theHessianEvaluator before curvatureFitting can be created" << endln;
            return TCL_ERROR;
        }
        
        // needs probability transformation
        if (theProbabilityTransformation == 0 ) {
            opserr << "Need theProbabilityTransformation before curvatureFitting can be created" << endln;
            return TCL_ERROR;
        }
        
		theFindCurvatures = new CurvatureFitting(theReliabilityDomain, theStructuralDomain, theFunctionEvaluator,
                                                 theFORMAnalysis, theHessianEvaluator, theProbabilityTransformation);
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

	int printFlag=0;
	char fileNamePrint[256];
	strcpy(fileNamePrint,"initialized");
	int maxNumIter = 100;
	int argvCounter = 1;

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[argvCounter],"StepSearch") == 0) {
		argvCounter++;

		// Check that the necessary ingredients are present
		if (theFunctionEvaluator == 0 ) {
			opserr << "Need theGFunEvaluator before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theGradientEvaluator == 0 ) {
			opserr << "Need theGradientEvaluator before a FindDesignPointAlgorithm can be created" << endln;
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
//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
			opserr << "Assume all RV's are independent" << endln;
			theProbabilityTransformation = 
			new AllIndependentTransformation(theReliabilityDomain,0);
//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
			opserr << "Need theProbabilityTransformation before a FindDesignPointAlgorithm can be created" << endln;
//			return TCL_ERROR;
		}
		if (theReliabilityConvergenceCheck == 0 ) {
			opserr << "Need theReliabilityConvergenceCheck before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}

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
					maxNumIter, theReliabilityDomain, theStructuralDomain,
					theFunctionEvaluator,
					theGradientEvaluator,
					theStepSizeRule,
					theSearchDirection,
					theProbabilityTransformation,
					theReliabilityConvergenceCheck,
					printFlag, fileNamePrint);
		
	}   //if StepSearch

	else if (strcmp(argv[argvCounter],"NewStepSearch") == 0) {
		
		argvCounter++;

		// Check that the necessary ingredients are present
		if (theFunctionEvaluator == 0 ) {
			opserr << "Need theGFunEvaluator before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
		if (theGradientEvaluator == 0 ) {
			opserr << "Need theGradientEvaluator before a FindDesignPointAlgorithm can be created" << endln;
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
//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
			opserr << "Assume all RV's are independent" << endln;
			theProbabilityTransformation = 
			new AllIndependentTransformation(theReliabilityDomain,0);
//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//			return TCL_ERROR;
		}
		if (theReliabilityConvergenceCheck == 0 ) {
			opserr << "Need theReliabilityConvergenceCheck before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}

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
		/*
		theFindDesignPointAlgorithm = new NewSearchWithStepSizeAndStepDirection(
					maxNumIter, theReliabilityDomain, 
					theFunctionEvaluator,
					theGradientEvaluator,
					theStepSizeRule,
					theSearchDirection,
					theProbabilityTransformation,
					theHessianApproximation,
					theReliabilityConvergenceCheck,
					startAtOrigin,
					printFlag, fileNamePrint);
		*/
		
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

	// Check that there are enough arguments
	if (argc<2) {
		opserr << "ERROR: Not enough arguments to theStartPoint. " << endln;
		return TCL_ERROR;
	}

	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	RandomVariable *aRandomVariable;


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[1],"Mean") == 0) {

		RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
		while ((aRandomVariable = rvIter()) != 0) {
            //int tag = aRandomVariable->getTag();
            double mean = aRandomVariable->getMean();
            aRandomVariable->setStartValue(mean);
		}
	}
    
	else if (strcmp(argv[1],"Origin") == 0) {
		
        RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
        while ((aRandomVariable = rvIter()) != 0) {
            //int tag = aRandomVariable->getTag();
            aRandomVariable->setStartValue(0.0);
        }
	}
    
	else if (strcmp(argv[1],"-file") == 0) {
        // space delimited file containing a starting point
		ifstream inputFile( argv[2], ios::in );
		if (inputFile.fail()) {
			opserr << "File " << argv[2] << " could not be opened for startPoint. " << endln;
			return TCL_ERROR;
		}

		// Loop through file to see how many entries there are
		double dummy;
		int numEntries = 0;
		while (inputFile >> dummy)
			numEntries++;
		
		if (numEntries == 0) {
			opserr << "ERROR: No entries in the file read by startPoint!" << endln;
			return TCL_ERROR;
		}
		if (numEntries != nrv) {
			opserr << "ERROR: Wrong number of entries in the file read by startPoint." << endln;
			return TCL_ERROR;
		}

		// rewind the file and pass values to the RVs
        inputFile.seekg(0, ios::beg);
        for (int i = 0; i < nrv; i++) {
            aRandomVariable = theReliabilityDomain->getRandomVariablePtrFromIndex(i);
			inputFile >> dummy;
            aRandomVariable->setStartValue(dummy);
		}
        
		inputFile.close();
	}
    
	else {
        opserr << "ERROR: Invalid type of start point is given. " << endln;
        return TCL_ERROR;
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


	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a root-finding algorithm can be created" << endln;
		return TCL_ERROR;
	}

	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a root-finding algorithm can be created" << endln;
		//////////////////////////////////////////////////////////////////////////////////
/////////////S Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
		opserr << "Assume all RV's are independent" << endln;
		theProbabilityTransformation = 
		new AllIndependentTransformation(theReliabilityDomain,0);
//		return TCL_ERROR;
//////////////////////////////////////////////////////////////////////////////////
/////////////E Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
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
			theFunctionEvaluator,
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


	// Check number of arguments
	if ( (argc!=2) && (argc!=4))  {
		opserr << "ERROR: Wrong number of input parameter to FORM analysis" << endln;
		return TCL_ERROR;
	}

    // Do input check
	inputCheck();

	// Check for essential tools
	if (theFindDesignPointAlgorithm == 0 ) {
		opserr << "Need theFindDesignPointAlgorithm before a FORMAnalysis can be created" << endln;
		return TCL_ERROR;
	}
    if (theFunctionEvaluator == 0 ) {
		opserr << "Need theFunctionEvaluator before a FORMAnalysis can be created" << endln;
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
                            theFunctionEvaluator,
							theProbabilityTransformation, 
							argv[1], relSensTag);


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
TclReliabilityModelBuilder_runFOSMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theFOSMAnalysis != 0) {
		delete theFOSMAnalysis;
		theFOSMAnalysis = 0;
	}

    
	// Check number of arguments
	if (argc != 2)  {
		opserr << "ERROR: Wrong number of input parameter to FOSM analysis" << endln;
		return TCL_ERROR;
	}

    // Do input check
	inputCheck();
    
	// Check for essential ingredients
	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a FOSMAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGradientEvaluator == 0 ) {
		opserr << "Need theGradientEvaluator before a FOSMAnalysis can be created" << endln;
		return TCL_ERROR;
	}


	theFOSMAnalysis = new FOSMAnalysis( theReliabilityDomain, theStructuralDomain,
											theFunctionEvaluator, theGradientEvaluator,
											interp, argv[1]);

	if (theFOSMAnalysis == 0) {
		opserr << "ERROR: could not create theFOSMAnalysis \n";
		return TCL_ERROR;
	}

	// Now run the analysis
	theFOSMAnalysis->analyze();

	return TCL_OK;
}





/*
//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runParametricReliabilityAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theParametricReliabilityAnalysis != 0) {
		delete theParametricReliabilityAnalysis;
		theParametricReliabilityAnalysis = 0;
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
		opserr << "Need theFindDesignPointAlgorithm before a ParametricReliabilityAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGradientEvaluator == 0 ) {
		opserr << "Need theGradientEvaluator before a ParametricReliabilityAnalysis can be created" << endln;
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
		theParametricReliabilityAnalysis = new ParametricReliabilityAnalysis( theReliabilityDomain,
												  theFindDesignPointAlgorithm,
												  theGradientEvaluator,
												  parameterNumber,
												  first,
												  last,
												  numIntervals,
												  argv[1],
												  interp);
	}
	else {
		opserr << "ERROR:: some input to theParametricReliabilityAnalysis was not provided" << endln;
		return TCL_ERROR;
	}

	if (theParametricReliabilityAnalysis == 0) {
		opserr << "ERROR: could not create theParametricReliabilityAnalysis \n";
		return TCL_ERROR;
	}

	// Now run the analysis
	theParametricReliabilityAnalysis->analyze();

	return TCL_OK;
}
*/

//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_runSORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theSORMAnalysis != 0) {
		delete theSORMAnalysis;
		theSORMAnalysis = 0;
	}

    
    // check minimum arguments
    if (argc != 2)  {
		opserr << "ERROR: Wrong number of arguments to SORM analysis" << endln;
		return TCL_ERROR;
	}
    
	// Do input check
	inputCheck();

    // check for essential ingredients
	if (theFindCurvatures == 0 ) {
		opserr << "Need theFindCurvatures before a SORMAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theFORMAnalysis == 0 ) {
		opserr << "ERROR: The current SORM implementation requires a FORM analysis" << endln
			<< " to have been executed previously in the same session." << endln;
		return TCL_ERROR;
	}
	

	theSORMAnalysis 
		= new SORMAnalysis(theReliabilityDomain, theFunctionEvaluator,
                           theFORMAnalysis, theFindCurvatures , argv[1]);

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
    
    // Check for essential ingredients
	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theFunctionEvaluator before a SystemAnalysis can be created" << endln;
		return TCL_ERROR;
	}

    // Do input check
	inputCheck();
    
	
	int aType = 1;
	char betaFile[MAX_FILENAMELENGTH] = "";
	char rhoFile[MAX_FILENAMELENGTH] = "";
	long int nMax = 0;
	double tol = 0;

	if (argc < 4 || argc % 2 > 0 )  {
		opserr << "ERROR: Wrong number of arguments to System Reliability analysis" << endln;
		opserr << "Want: runSystemAnalysis fileName? analysisMethod? (allInParallel | allInSeries) <-Nmax val?> <-tol val?> <B_fileName R_fileName>" << endln;
		opserr << "analysisMethod options are: PCM, IPCM, MVN, and SCIS" << endln;
		return TCL_ERROR;
	} else {		
		int argi = 4;
		while (argi < argc) {
			if (strcmp(argv[argi],"-Nmax") == 0) {
				nMax = atol(argv[argi+1]);
				if ( nMax <= 0 ) {
					opserr << "WARNING invalid Nmax = " << argv[argi+1] << endln;
					return TCL_ERROR;
				}
			}
			else if (strcmp(argv[argi],"-tol") == 0) {
				if (Tcl_GetDouble(interp, argv[argi+1], &tol) != TCL_OK) {
					opserr << "WARNING invalid tol = " << argv[argi+1] << endln;
					return TCL_ERROR;
				}
			}
			else {
				// option of specifying files with beta and rho instead of using information in the reliability domain
				strcpy(betaFile,argv[argi]);
				strcpy(rhoFile,argv[argi+1]);
			}
			
			argi += 2;
		}
	}
	
	// GET INPUT PARAMETER (string)
	if (strcmp(argv[3],"allInParallel") == 0)
		aType = 0;
	else if (strcmp(argv[3],"allInSeries") == 0)
		aType = 1;
	else if (strcmp(argv[3],"cutsets") == 0)
		aType = 2;
	else {
		opserr << "ERROR: Invalid system reliability analysis type input:" << argv[3] << endln;
		return TCL_ERROR;
	}

	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT
	if (strcmp(argv[2],"PCM") == 0)
		theSystemAnalysis = new PCM(theReliabilityDomain, theFunctionEvaluator, 
                                    argv[1], aType, betaFile, rhoFile);
	else if (strcmp(argv[2],"IPCM") == 0)
		theSystemAnalysis = new IPCM(theReliabilityDomain, theFunctionEvaluator,
                                     argv[1], aType, betaFile, rhoFile);
	else if (strcmp(argv[2],"MVN") == 0)
		theSystemAnalysis = new MVNcdf(theReliabilityDomain, theFunctionEvaluator,
                                       argv[1], aType, betaFile, rhoFile, nMax, tol);
	else if (strcmp(argv[2],"SCIS") == 0)
		theSystemAnalysis = new SCIS(theReliabilityDomain, theFunctionEvaluator,
                                     argv[1], aType, betaFile, rhoFile, nMax, tol);
	else {
		opserr << "ERROR: Invalid system reliability analysis type input:" << argv[2] << endln;
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
TclReliabilityModelBuilder_runImportanceSamplingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theImportanceSamplingAnalysis != 0) {
		delete theImportanceSamplingAnalysis;
		theImportanceSamplingAnalysis = 0;
	}


    // Do input check
	inputCheck();

	// Check for essential tools
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a SimulationAnalyis can be created" << endln;
		return TCL_ERROR;
	}
	if (theFunctionEvaluator == 0 ) {
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
	//     -type  saveGvalues (3)
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
	long int numberOfSimulations	= 1000;
	double targetCOV		= 0.05;
	double samplingVariance	= 1.0;
	int printFlag			= 0;
	int analysisTypeTag		= 1;


	for (int i=2; i<argc; i=i+2) {

		if (strcmp(argv[i],"-type") == 0) {

			if (strcmp(argv[i+1],"failureProbability") == 0) {
				analysisTypeTag = 1;
			}

// Michele and Quan -------------------------
			else if (strcmp(argv[i+1],"outCrossingFailureProbability") == 0) {
				analysisTypeTag = 4;
			}

			else if ( (strcmp(argv[i+1],"responseStatistics") == 0) || (strcmp(argv[i+1],"saveGvalues") == 0) ) {

				if (strcmp(argv[i+1],"responseStatistics") == 0) {
					analysisTypeTag = 2;
				}
				else {
					analysisTypeTag = 3;
				}
				if (samplingVariance != 1.0) {
					opserr << "ERROR:: sampling variance must be 1.0 for " << endln
						<< " response statistics sampling." << endln;
					return TCL_ERROR;
				}
				
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
            numberOfSimulations = atol(argv[i+1]);
			//if (Tcl_GetInt(interp, argv[i+1], &numberOfSimulations) != TCL_OK) {
			//	opserr << "ERROR: invalid input: numberOfSimulations \n";
			//	return TCL_ERROR;
			//}
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
	
	
	theImportanceSamplingAnalysis 
			= new ImportanceSamplingAnalysis(theReliabilityDomain, theStructuralDomain, 
							 theProbabilityTransformation, 
							 theFunctionEvaluator, 
							 theRandomNumberGenerator, 
							 interp, 
							 numberOfSimulations, targetCOV, samplingVariance,
							 printFlag,
							 argv[1],
							 analysisTypeTag);

	if (theImportanceSamplingAnalysis == 0) {
		opserr << "ERROR: could not create theImportanceSamplingAnalysis \n";
		return TCL_ERROR;
	}

	// Now run analysis
	theImportanceSamplingAnalysis->analyze();

	return TCL_OK;

}

//////////////////////////////////////////////////////////////////


// Quan and Michele Feb 2006

// command "runOutCrossingAnalysis  filename?  -results stepsToStart?  stepsToEnd?  samplefreq? impulseFreq?   -littleDt dt? -analysisType 
// option for analysisType 1:   -twoSearches   <-integralTolerance  tol? -useFirstDesignPoint>
//            analysisType 2:    -Koo
//                                       
int 
TclReliabilityModelBuilder_runOutCrossingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theOutCrossingAnalysis != 0) {
		delete theOutCrossingAnalysis;
		theOutCrossingAnalysis = 0;
	}


    // Do input check
	inputCheck();

    // check for essential ingredients
	if (theFindDesignPointAlgorithm == 0 ) {
		opserr << "Need theFindDesignPointAlgorithm before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGradientEvaluator == 0 ) {
		opserr << "Need theGradientEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}

	int stepsToStart = 0;
	int stepsToEnd = 0;
	int sampleFreq = 1;
	double littleDt = 0.01;
	int analysisType = 1;
	
	int impulseFreq;

	double integralTolerance=1.e-10;
	bool useFirstDesignPt = false;



	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-results") == 0) {
			argvCounter++;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &stepsToStart) != TCL_OK) {
				opserr << "ERROR: invalid input stepsToStart to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &stepsToEnd) != TCL_OK) {
				opserr << "ERROR: invalid input stepsToEnd to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			// GET INPUT PARAMETER (integer)
			if (Tcl_GetInt(interp, argv[argvCounter], &sampleFreq) != TCL_OK) {
				opserr << "ERROR: invalid input sampleFreq to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;

			if (Tcl_GetInt(interp, argv[argvCounter], &impulseFreq) != TCL_OK) {
				opserr << "ERROR: invalid input impulseFreq to theOutCrossingAnalysis \n";
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
			if (strcmp(argv[argvCounter],"-integralTolerance") == 0) {
				argvCounter++;
				if (Tcl_GetDouble(interp, argv[argvCounter], &integralTolerance) != TCL_OK) {
					opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
					return TCL_ERROR;
				}
				argvCounter++;
			 }
			else if (strcmp(argv[argvCounter],"-useFirstDesignPoint") == 0) {
				argvCounter++;
				useFirstDesignPt =true;
				
			 }

		}
		else {
			opserr << "ERROR: Invalid input to theOutCrossingAnalysis." << endln;
			return TCL_ERROR;
			argvCounter++;
		}
	}

	theOutCrossingAnalysis 
			= new OutCrossingAnalysis(
				theReliabilityDomain,
				theFunctionEvaluator,
				theGradientEvaluator,
				theFindDesignPointAlgorithm,
				analysisType,
				stepsToStart,
				stepsToEnd,
				sampleFreq,
				impulseFreq,
				littleDt,
				argv[1],
				integralTolerance,
				useFirstDesignPt);

	if (theOutCrossingAnalysis == 0) {
		opserr << "ERROR: could not create theOutCrossingAnalysis \n";
		return TCL_ERROR;
	}

	// Now run analysis
	theOutCrossingAnalysis->analyze();

	return TCL_OK;

}


//////////////////////////////////////////////////////////////////
// Quan and Michele April 2006

// command "runOrthogonalPlaneSamplingAnalysis  -fileName filename?  -maxNum number?   -type  analysisType? -targetCOV cov? -print printFlag? 
// -funcTol tol1? -varTol tol2? -maxIter iter? -littleDt littleDt?....
// option for analysisType  "failureProbability"   ---1:   failure probability.
//            analysisType "outCrossing"           ---2:   upcrossing problem.
//                                       
// Not finish yet ????????????????????????????????????????????????????????????????????????????
int 
TclReliabilityModelBuilder_runOrthogonalPlaneSamplingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theSamplingAnalysis != 0) {
		delete theSamplingAnalysis;
		theSamplingAnalysis = 0;
	}


    // Do input check
	inputCheck();

	// Check for essential tools
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a SimulationAnalyis can be created" << endln;
		return TCL_ERROR;
	}
	if (theFunctionEvaluator == 0 ) {
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
	//     -type  outcrossing (2)
	//     -type  saveGvalues (3) .... not yet
    //
	//     -maxNum 1000  ........................ this is the default
	//
	//     -targetCOV 0.05  ..................... this is the default
	//
	//     -print 0   (print nothing) 
	//     -print 1   failure prob. and cov  .... this is the default
	//     -print 2   ............................ restart file  // 2007 Feb. Quan
	//     -print 5    recorder the surface ......this is for visualization
	//     -funcTol  1.e-5 ....................... this is the default
	//     -varTol   1.e-3........................ this is the default
	//     -maxIter   20  .........................this is the default

	if (argc<3) {
		opserr << "command: runOrthogonalPlaneSamplingAnalysis  -fileName filename?  -maxNum number?   -type  analysisType? -targetCOV cov? -print printFlag? ";
		opserr<<" -funcTol tol1? -varTol tol2? -maxIter iter?" << endln;
		return TCL_ERROR;
	}

	// Declaration of input parameters
	int numberOfSimulations	= 1000;
	double targetCOV		= 0.05;
	int printFlag			= 1;
	double funcTol = 1.e-5;
	double varTol = 1.e-3;
	int maxIter = 20;
	int analysisTypeTag		= 1;
	char name[50];
	Vector * theDesignPoint;
	double littleDt = 1.0e-3;

	/*
	if (theStartPoint == 0 ) {
		opserr << "orthogonalPlaneSamplingAnalysis can not run. Need StartPoint !" << endln;
		return TCL_ERROR;
	}
	else theDesignPoint = theStartPoint;
	*/
	//theReliabilityDomain->getStartPoint(*theDesignPoint);
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	for (int i = 0; i < nrv; i++) {
	  RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(i);
	  (*theDesignPoint)(i) = theRV->getStartValue();
	}

	int argvCounter = 1;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-fileName") == 0)||(strcmp(argv[argvCounter],"-filename") == 0)) {
			argvCounter++;
			strcpy(name,argv[argvCounter]);
			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-maxNum")==0) ||(strcmp(argv[argvCounter],"-maxnum") == 0)) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &numberOfSimulations) != TCL_OK) {
			opserr << "ERROR: invalid input: numberOfSimulations \n";
			return TCL_ERROR;
			}
			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-type") == 0) {
			argvCounter++;
			if (strcmp(argv[argvCounter],"failureProbability") == 0) {
				analysisTypeTag = 1;
			}
			else if (strcmp(argv[argvCounter],"outCrossing") == 0) {
				analysisTypeTag = 2;
			}
			argvCounter++;
		}

		else if (strcmp(argv[argvCounter],"-targetCOV") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter++], &targetCOV) != TCL_OK) {
				opserr << "ERROR: invalid input: targetCOV \n";
				return TCL_ERROR;
			}
		}
		else if (strcmp(argv[argvCounter],"-print") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter++], &printFlag) != TCL_OK) {
				opserr << "ERROR: invalid input: printFlag \n";
				return TCL_ERROR;
			}
		}
 
		else if (strcmp(argv[argvCounter],"-funcTol") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter++], &funcTol) != TCL_OK) {
				opserr << "ERROR: invalid input: funcTol \n";
				return TCL_ERROR;
			}
		}
				
		else if (strcmp(argv[argvCounter],"-varTol") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter++], &varTol) != TCL_OK) {
				opserr << "ERROR: invalid input: varTol \n";
				return TCL_ERROR;
			}
		}
		
		else if (strcmp(argv[argvCounter],"-maxIter") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter++], &maxIter) != TCL_OK) {
				opserr << "ERROR: invalid input: maxIter \n";
				return TCL_ERROR;
			}
		}
		else if (strcmp(argv[argvCounter],"-littleDt") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter++], &littleDt) != TCL_OK) {
				opserr << "ERROR: invalid input: littleDt \n";
				return TCL_ERROR;
			}
		}
		
		else {
			opserr << "ERROR: invalid input to sampling analysis. " << endln;
			return TCL_ERROR;
		}
	}


	/*
		OrthogonalPlaneSamplingAnalysis(   Tcl_Interp *interp,
		                ReliabilityDomain *passedReliabilityDomain,
						ProbabilityTransformation *passedProbabilityTransformation,
						GFunEvaluator *passedGFunEvaluator,
						RandomNumberGenerator *passedRandomNumberGenerator,
						int passedNumberOfSimulations,
						int passedMaxNumOfIterations,
						double passedTargetCOV,
						double samplingStdv,
						int printFlag,
						TCL_Char *fileName,
						Vector * pDesignPoint,
						int analysisTypeTag,
						int zeroFindingType);
	
	*/

	
	// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx
	theSamplingAnalysis 
			= new OrthogonalPlaneSamplingAnalysis(interp, 
									theReliabilityDomain, 
									theProbabilityTransformation, 
									theFunctionEvaluator, 
									theRandomNumberGenerator, 
									numberOfSimulations,
									maxIter,
									targetCOV,
     								printFlag,
									name,
									theDesignPoint,
									analysisTypeTag,
									1,
									funcTol,
									varTol,
									maxIter,
									littleDt);

	if (theSamplingAnalysis == 0) {
		opserr << "ERROR: could not create theOrthogonalSamplingAnalysis \n";
		return TCL_ERROR;
	}

	// Now run analysis
	theSamplingAnalysis->analyze();

 
	return TCL_OK;

}


//////////////////////////////////////////////////////////////////


// Quan & Michele: add command for visualization of another zerofinding algorithm
//  command: runGFunVizAnalysis outputfile -space y    -funSurf surface     -dir file designpoint.out -file filename numPts? -zeroFindingAlgorithm safeguardedZeroFinding



int 
TclReliabilityModelBuilder_runGFunVisualizationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theGFunVisualizationAnalysis != 0) {
		delete theGFunVisualizationAnalysis;
		theGFunVisualizationAnalysis = 0;
	}


    // Do input check
	inputCheck();

    // check for essential ingredients
	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a GFunVisualizationAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a GFunVisualizationAnalysis can be created" << endln;
		return TCL_ERROR;		
	}


	// Initial declarations
	int zeroFindingAlg =0;

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
					opserr << "File " << argv[argvCounter] << " could not be opened. " << endln;
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
					opserr << "File " << argv[argvCounter] << " could not be opened. " << endln;
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

			//	argvCounter++;   -- wrong  Quan
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
				opserr << "File " << argv[argvCounter] << " could not be opened. " << endln;
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
			if ((numEntries % numRVs) !=0.0) {
				opserr << "ERROR: Wrong number of entries in the the file " << argv[argvCounter] << endln;
				return TCL_ERROR;
			}
			numVectors = (int)(numEntries/numRVs);

			// Close the file
			inputFile.close();

			// Open it again, now being ready to store the results in a matrix
			ifstream inputFile2( argv[argvCounter], ios::in );
			if (inputFile2.fail()) {
				opserr << "File " << argv[argvCounter] << " could not be opened. " << endln;
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
		else if (strcmp(argv[argvCounter],"-zeroFindingAlgorithm") == 0) {
			argvCounter++;
			if (strcmp(argv[argvCounter],"safeguardedZeroFinding") == 0){
				zeroFindingAlg = 1;
				argvCounter++;
			}

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

	if (zeroFindingAlg ==0) 
	   theGFunVisualizationAnalysis = 
	     new GFunVisualizationAnalysis(theReliabilityDomain, 
					   theFunctionEvaluator, 
					   theProbabilityTransformation, 
					   argv[1],
					   argv[convFileArgv],
					   convResults,
					   space,
					   funSurf,
					   axes,
					   dir);
/*	else if (zeroFindingAlg ==1) 
	   theGFunVisualizationAnalysis = new GFunVisualizationSamplingAnalysis(
											theReliabilityDomain, 
											theFunctionEvaluator, 
											theProbabilityTransformation, 
											argv[1],
											argv[convFileArgv],
											convResults,
											space,
											funSurf,
											axes,
											dir,
											1);

*/
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

//		if (theStartPoint == 0 ) {
//			opserr << "Need theStartPoint before this GFunVisualizationAnalysis can be created" << endln;
//			return TCL_ERROR;
//		}
		
		//theGFunVisualizationAnalysis->setStartPoint(theStartPoint);
	}

	if (convResults == 1) {

		if (theGradientEvaluator == 0 ) {
			opserr << "Need theGradientEvaluator before this GFunVisualizationAnalysis can be created" << endln;
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

		theGFunVisualizationAnalysis->setGradGEvaluator(theGradientEvaluator);
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
TclReliabilityModelBuilder_printReliability(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc > 1)
    theReliabilityDomain->Print(opserr, 1);
  else
    theReliabilityDomain->Print(opserr);

  return TCL_OK;
}


// ---------- Quan Gu ------------------------


///Command:  runMonteCarloResponseAnalysis  -outPutFile  m.out -maxNum 1000 -print 1 -tclFileToRun test.tcl <-seed 1>
int 
TclReliabilityModelBuilder_runMonteCarloResponseAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theMonteCarloResponseAnalysis != 0) {
		delete theMonteCarloResponseAnalysis;
		theMonteCarloResponseAnalysis = 0;
	}
	
	int seed=1;

	// Do input check
	inputCheck();

	// Check for essential tools
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a SimulationAnalyis can be created" << endln;
		return TCL_ERROR;
	}

	if (theRandomNumberGenerator == 0 ) {
		opserr << "Need theRandomNumberGenerator before a SimulationAnalyis can be created" << endln;
		return TCL_ERROR;
	}

  int numberOfArguments = argc;
  if (numberOfArguments < 4) {
		opserr << "ERROR: invalid number of arguments to designVariable command "<<endln;
	    opserr <<"command: runMonteCarloResponseAnalysis  -outPutFile  m.out -maxNum 1000 -print 1 -tclFileToRun test.tcl"<<endln;
		return TCL_ERROR;
  }	


	// Declaration of input parameters
	int numberOfSimulations	= 1000;
	int printFlag			= 0;
	char outPutFile[25]="";
	char * tclFileName = 0;

	int argvCounter = 1;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-maxNum") == 0)||(strcmp(argv[argvCounter],"-maxnum") == 0)) {

			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &numberOfSimulations) != TCL_OK) {
			opserr << "ERROR: invalid input: numberOfSimulations \n";
			return TCL_ERROR;
			}

			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-print") == 0)||(strcmp(argv[argvCounter],"-printFlag") == 0)) {

			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &printFlag) != TCL_OK) {
			opserr << "ERROR: invalid input: printFlag \n";
			return TCL_ERROR;
			}

			argvCounter++;
		}// if
		else if ((strcmp(argv[argvCounter],"-outPutFile") == 0)||(strcmp(argv[argvCounter],"-outputfile") == 0)) {
			
			argvCounter++;
			
			strcpy(outPutFile,argv[argvCounter]);
			
			argvCounter++;
		}// else if
		
		else if ((strcmp(argv[argvCounter],"-tclFileToRun") == 0)||(strcmp(argv[argvCounter],"-TclFileToRun") == 0)) {
			argvCounter++;
			tclFileName = new char[25];
			strcpy(tclFileName,argv[argvCounter]);
			argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-seed") == 0) {
			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &seed) != TCL_OK) {
			opserr << "ERROR: invalid input: seed \n";
			return TCL_ERROR;
			}
			argvCounter++;
		}// else if

		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		} //else

	};  // while


	
	theMonteCarloResponseAnalysis
			= new MonteCarloResponseAnalysis(theReliabilityDomain,
						interp,
						theProbabilityTransformation,
						theRandomNumberGenerator,
						numberOfSimulations,
						printFlag,
						outPutFile,
						tclFileName,
						seed);
			
			
	if (theMonteCarloResponseAnalysis == 0) {
		opserr << "ERROR: could not create theMonteCarloResponseAnalysis \n";
		return TCL_ERROR;
	}

	if (tclFileName !=0) delete [] tclFileName;

	// Now run analysis
	theMonteCarloResponseAnalysis->analyze();
	
	return TCL_OK;

}


///////
///  Command:  updateParameterValue  -rv 1 -value 20.0  or updateParameter   -startPoint 
///            updateParameterValue  -dv 3 -value 20.0
int 
TclReliabilityModelBuilder_updateParameterValue(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  /*
	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	// Declaration of input parameters
	int numberOfRVDVPositioners;

	int dvrv;
	double value;
	bool usingStartPt = false;

	int dvrvMark = 0; // by default, is randomvariable 0
	int argvCounter = 1;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-dv") == 0)||(strcmp(argv[argvCounter],"-dvNum") == 0)) {

			argvCounter++;
	
			if (Tcl_GetInt(interp, argv[argvCounter], &dvrv) != TCL_OK) {
				opserr << "ERROR: invalid input: dv number \n";
				return TCL_ERROR;
			}
			dvrvMark=1; // dv
			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-rv") == 0)||(strcmp(argv[argvCounter],"-rvNum") == 0)) {

			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &dvrv) != TCL_OK) {
				opserr << "ERROR: invalid input: rv Number \n";
				return TCL_ERROR;
			}
			dvrvMark=0; // rv
			argvCounter++;
		}// if
		else if ((strcmp(argv[argvCounter],"-value") == 0)||(strcmp(argv[argvCounter],"-VALUE") == 0)) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &value) != TCL_OK) {
				opserr << "ERROR: invalid input: value \n";
				return TCL_ERROR;
			}
//			opserr<<"value is:"<<value<<endln;
			
			argvCounter++;
		}

		else if ((strcmp(argv[argvCounter],"-startPoint") == 0)||(strcmp(argv[argvCounter],"-startpoint") == 0)) {
			argvCounter++;
			usingStartPt = true;
		}

		else {
			opserr<<"warning: unknown command: updateparameter" <<argv[argvCounter]<<endln;
			argvCounter++;
		
		} //else

	};  // while

	if (! usingStartPt) {
		if (dvrvMark==0) { // rv

			RandomVariablePositioner * theRandomVariablePositioner =0;
			numberOfRVDVPositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
			if (numberOfRVDVPositioners==0) {opserr<<"warnning: updateParameter no randomVariablePositioner"<<endln;  }
			else {
				int rvNumber; 
				for (int i=1 ; i<=numberOfRVDVPositioners ; i++ )  {
					theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
					rvNumber = theRandomVariablePositioner->getRvNumber();
					if (rvNumber == dvrv )  theRandomVariablePositioner->update(value);
				}	
			}
		} 
		else if (dvrvMark==1) { // dv
			DesignVariable * theDesignVariable = theReliabilityDomain->getDesignVariablePtr(dvrv);
			theDesignVariable->update(value);

		}
	}
	
	else if (usingStartPt){ // only rv is possible
	
			RandomVariablePositioner * theRandomVariablePositioner =0;
			int numberOfRVDVPositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
			for (int i=1 ; i<=numberOfRVDVPositioners ; i++ )  {
				theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
				int rvNumber = theRandomVariablePositioner->getRvNumber();
				value = (*theStartPoint)(rvNumber-1);
				theRandomVariablePositioner->update(value);
			}
	
				
	}	

*/ 

	return TCL_OK;

}


///////
///  Command:  computeHessian -FDM -file $filename1 -designPoint $fileName -perturbation $pTol
//////
int 
TclReliabilityModelBuilder_computeHessian(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	char fileName[20]="HessianByFDM.out";
	char designPointFile[20] = "theDesignPoint.out";
	int argvCounter = 1;
	bool FDM = false;
	double pTol = 1.0e-5; // default
	
	while (argc > argvCounter) {

		if ((strcmp(argv[argvCounter],"-FDM") == 0)||(strcmp(argv[argvCounter],"-FFD") == 0)) {

			argvCounter++;
			FDM =true;
	
		}// if

		else if ((strcmp(argv[argvCounter],"-file") == 0)||(strcmp(argv[argvCounter],"-File") == 0)) {
			argvCounter++;
			strcpy(fileName,argv[argvCounter]);
			argvCounter++;
		}// else if

		else if ((strcmp(argv[argvCounter],"-designPoint") == 0)||(strcmp(argv[argvCounter],"-designpoint") == 0)) {
			argvCounter++;
			strcpy(designPointFile,argv[argvCounter]);
			argvCounter++;
		}// else if
		else if ((strcmp(argv[argvCounter],"-perturbation") == 0)||(strcmp(argv[argvCounter],"-perturbationTolerance") == 0)) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &pTol) != TCL_OK) {
				opserr << "ERROR: invalid input: perturbationTolerance \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}// else if
		else {
		
			opserr<< "unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}
	}

	int size = theReliabilityDomain->getNumberOfRandomVariables();
	Vector * designPoint = new Vector(size);

	ofstream resultsOutputFile( fileName, ios::out);
    ifstream inputFile( designPointFile, ios::in);

	if (inputFile.good()){
		int ii=0;
		double tmp;
		while(!inputFile.eof() && ii<size){ 
			inputFile >> tmp;
			(*designPoint)(ii)=tmp;
			ii++;
		}
	} 
	else {
		opserr<<"ERROR: designpoint can not read from file "<<designPointFile<<endln;
		exit(-1);
	}
	
    // KRM 5-19-2012
    // this needs to be implemented properly using the reliability base classes (FunctionEvaluator, 
    // GradientEvaluator, HessianEvaluator, etc.)
	//Hessian * theHessian = new Hessian(size,theReliabilityDomain,theProbabilityTransformation,theFunctionEvaluator,theGradientEvaluator,pTol);
 
	if (FDM) { 
		//theHessian->formReducedHessian(designPoint);
		
		
		//Matrix hessian=	theHessian->getHessianApproximation();

//		resultsOutputFile <<"Hessian in U space: \n";
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				//resultsOutputFile <<hessian(i,j)<<"   ";
			}
			resultsOutputFile <<"\n";
		}
	
	//	opserr<<"------ theHessian: ------- \n"<<hessian<<endln;

	//	opserr<<"\n---------DesignPoint:   ----------\n "<<*designPoint<<endln;


	//	Matrix reducedHessian = theHessian->getReducedHessian();
	//	opserr<<"\n -------- theReducedHessian: ---------\n"<<reducedHessian<<endln;


		delete designPoint;	
	} 

	return 0;

}






// MultiDimVisualPrinPlane -funSurf function -designPt dp.out  -ndir $n -output vis.out <-gridInfo {0  minY  maxY nPts0 1 minX1  maxX1 nPts1 2 minX2  maxX2 nPts2 ...}  -timeVariant -littleDt 0.001> <-saveHessian $filename>

int 
TclReliabilityModelBuilder_MultiDimVisPrincPlane(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int argvCounter = 1;
	int type =1;
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	int numPPlane = 1;
	Vector * designPoint = new Vector(nrv);
	Vector * gridInfo = 0;
    Matrix * passedHessian =0;
	char outputFile[25];
	char * hessianFileName=0;
	int analysisType = 0;  // 
	double littleDt  =1.0e-3;
 
	while (argc > argvCounter) {

		if ((strcmp(argv[argvCounter],"-funSurf") == 0)||(strcmp(argv[argvCounter],"-FunSurf") == 0)) {
			argvCounter++;
			if (strcmp(argv[argvCounter],"function") == 0)
				type = 1;
			else if (strcmp(argv[argvCounter],"surface") == 0)
				type =0;
			argvCounter++;
			
	
		}// if

		else if ((strcmp(argv[argvCounter],"-designPt") == 0)||(strcmp(argv[argvCounter],"-designPoint") == 0)) {
			char fileName[20];
			argvCounter++;
			strcpy(fileName,argv[argvCounter]);
			argvCounter++;
			ifstream inputFile( fileName, ios::in);

			if (inputFile.good()){
				int ii=0;
				double tmp;
				while(!inputFile.eof() && ii<nrv){ 
					inputFile >> tmp;
					(*designPoint)(ii)=tmp;
					ii++;
				}
				inputFile.close();
			} 
			else {
				opserr<<"ERROR: designpoint can not read from file "<<fileName<<endln;
				exit(-1);
			}
						
		}// else if

		else if (strcmp(argv[argvCounter],"-ndir") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &numPPlane) != TCL_OK) {
				opserr << "ERROR: invalid input: numPPlane \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if
		else if (strcmp(argv[argvCounter],"-output") == 0) {
			argvCounter++;
			strcpy(outputFile,argv[argvCounter]);
			argvCounter++;
		}// else if
		else if (strcmp(argv[argvCounter],"-gridInfo") == 0) {
			argvCounter++;

			int pathSize;
		    TCL_Char **pathStrings;
		  
		    if (Tcl_SplitList(interp, argv[argvCounter], 
					&pathSize, &pathStrings) != TCL_OK) {
				  
			  opserr << "WARNING problem splitting path list in gridInfo\n";
			  return 0;
			}
		  
		    gridInfo = new Vector(pathSize);

		    for (int i=0; i<pathSize; i++) {
			  double value;
			  if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
			    opserr << "WARNING problem reading path data value " << pathStrings[i] << "\n";

			    return 0;
			  }
			  (*gridInfo)(i) = value;
		  }  //for
		  // free up the array of pathsStrings .. see tcl man pages as to why
//		  cleanup(pathStrings);

			argvCounter++;
		}// else if
		else if (strcmp(argv[argvCounter],"-saveHessian") == 0){
			argvCounter++;	
			hessianFileName = new char[30];
			strcpy(hessianFileName,argv[argvCounter]);
			argvCounter++;
			
			ifstream inputFile( hessianFileName, ios::in);

			if (inputFile.good()){
				passedHessian = new Matrix(nrv,nrv);
				int ii=0;
				int jj=0;
				double tmp;
				while(ii<nrv){ 
					for (jj=0; jj<nrv;jj++){
						if (!inputFile.eof()){
							inputFile >> tmp;
							(*passedHessian)(ii,jj)=tmp;
						}
						else{
							opserr<<"-saveHessian, size of Hessian in file: "<<hessianFileName<< "is wrong"<<endln;
							exit(-1);
						}
					}	
					
					ii++;
				}// while 
				inputFile.close();

			} 
			
			else {
				opserr<<"warning: no data in hessianfrom file: "<<hessianFileName<<endln;
				
			}
		
		} //else if hessian
		else if (strcmp(argv[argvCounter],"-timeVariant") == 0) {
			argvCounter++;
			analysisType =1;
		}// else if
		else if (strcmp(argv[argvCounter],"-littleDt") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &littleDt) != TCL_OK) {
				opserr << "ERROR: invalid input: littleDt \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if
		else {
		
			opserr<< "unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}
	} //while


    MultiDimVisPrincPlane * theMultiDimVisPrincPlane = new MultiDimVisPrincPlane(
		            theReliabilityDomain,
					theFunctionEvaluator,
					theProbabilityTransformation,
					outputFile,
					theGradientEvaluator,
					designPoint, numPPlane, type, gridInfo,interp, passedHessian, hessianFileName, analysisType, littleDt);

	if (passedHessian !=0) delete passedHessian;
	if (hessianFileName !=0) delete hessianFileName;
	theMultiDimVisPrincPlane->analyze();


    delete designPoint;
	return 0;

}




///////
///  Command:  transformXtoU  -fileX pointx.out   -fileU pointu.out   
///            
// KRM - I see no reason for this code to be in TclReliabilityBuilder, will delete soon 
// not being used by any other classes

int 
TclReliabilityModelBuilder_transformXtoU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theProbabilityTransformation ==0){
    opserr<<"Fatal: theProbabilityTransformation does not exist!"<<endln;
    exit(-1);
  }

  int nrv = theReliabilityDomain->getNumberOfRandomVariables();
  Vector x(nrv);

  int numInput = 0;
  TCL_Char **rvInput;
  // Version using a Tcl list
  if (argc == 2) {
    if (Tcl_SplitList(interp, argv[1], &numInput, &rvInput) != TCL_OK) {
      opserr << "transformXtoU -- error splitting input list of X realizations" << endln;
      return 0;
    }
    if (numInput != nrv) {
      opserr << "transformXtoU -- num X realizations not equal to nrv" << endln;
      return 0;
    }
    for (int i = 0; i < nrv; i++) {
      double value;
      if ( Tcl_GetDouble(interp, rvInput[i], &value) != TCL_OK) {
	opserr << "WARNING problem reading X realization " << rvInput[i] << "\n";
	
	return 0;
      }
      x(i) = value;
    }
  }

	char filenameX[30]="pointx.out";
	char filenameU[30]="pointu.out";

	int argvCounter = 1;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-fileX") == 0) {

			argvCounter++;
			strcpy(filenameX,argv[argvCounter]);
			argvCounter++;
		}// if

		else if (strcmp(argv[argvCounter],"-fileU") == 0) {

			argvCounter++;
			strcpy(filenameU,argv[argvCounter]);
			argvCounter++;
		}// if


		else {
			opserr<<"warning: unknown command: updateparameter" <<argv[argvCounter]<<endln;
			argvCounter++;
		
		} //else

	};  // while
	

	ifstream inputFile( filenameX, ios::in);

	Vector pointX(nrv);

	if (inputFile.good()){
		int ii=0;
		double tmp;
		while(ii<nrv){ 
			if (!inputFile.eof()){
				inputFile >> tmp;
				pointX(ii)=tmp;
			}
			else{
				opserr<<"file" <<filenameX<<" has different size than numRV"<<endln;
				exit(-1);
			}
			ii++;
		}// while 
		inputFile.close();

	} 
	else {
		opserr<<"warning: no data in file: "<<filenameX<<endln;
		exit(-1);			
	}

	Vector pointU;
	theProbabilityTransformation->transform_x_to_u(pointU);

	ofstream outputFile( filenameU, ios::out);
	outputFile.precision(16);

	for (int ii =0; ii<nrv; ii++){

		outputFile << pointU(ii)<<endln;
	}// for
 
	outputFile.close();

	return TCL_OK;

}

///////
///  Command:  transformUtoX  -fileX pointx.out   -fileU pointu.out   
///            
// KRM - I see no reason for this code to be in TclReliabilityBuilder, will delete soon 
// not being used by any other classes

int 
TclReliabilityModelBuilder_transformUtoX(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theProbabilityTransformation ==0){
    opserr<<"Fatal: theProbabilityTransformation does not exist!"<<endln;
    exit(-1);
  }

  int nrv = theReliabilityDomain->getNumberOfRandomVariables();

  // Version using a Tcl list
  int numInput = 0;
  TCL_Char **rvInput;
  if (argc == 2) {
    if (Tcl_SplitList(interp, argv[1], &numInput, &rvInput) != TCL_OK) {
      opserr << "transformUtoX -- error splitting input list of U realizations" << endln;
      return 0;
    }
    if (numInput != nrv) {
      opserr << "transformUtoX -- num U realizations not equal to nrv" << endln;
      return 0;
    }
    Vector u(nrv);
    for (int i = 0; i < nrv; i++) {
      double value;
      if ( Tcl_GetDouble(interp, rvInput[i], &value) != TCL_OK) {
	opserr << "WARNING problem reading U realization " << rvInput[i] << "\n";
	
	return 0;
      }
      u(i) = value;
    }

    Vector x(nrv);
    theProbabilityTransformation->transform_u_to_x(u, x);

    char buffer[20];
    for (int i = 0; i < nrv; i++) {
      sprintf(buffer, "%f ", x(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
    
    return TCL_OK;
  }




	
//	double value;
	char filenameX[30]="pointx.out";
	char filenameU[30]="pointu.out";

	int argvCounter = 1;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-fileX") == 0) {

			argvCounter++;
			strcpy(filenameX,argv[argvCounter]);
			argvCounter++;
		}// if

		else if (strcmp(argv[argvCounter],"-fileU") == 0) {

			argvCounter++;
			strcpy(filenameU,argv[argvCounter]);
			argvCounter++;
		}// if


		else {
			opserr<<"warning: unknown command: updateparameter" <<argv[argvCounter]<<endln;
			argvCounter++;
		
		} //else

	};  // while
	

	ifstream inputFile( filenameU, ios::in);

	Vector pointU(nrv);

	if (inputFile.good()){
		int ii=0;
		double tmp;
		while(ii<nrv){ 
			if (!inputFile.eof()){
				inputFile >> tmp;
				pointU(ii)=tmp;
			}
			else{
				opserr<<"file" <<filenameX<<" has different size than numRV"<<endln;
				exit(-1);
			}
			ii++;
		}// while 
		inputFile.close();

	} 
	else {
		opserr<<"warning: no data in file: "<<filenameX<<endln;
		exit(-1);			
	}
	
    Vector pointX;
	theProbabilityTransformation->transform_u_to_x(pointU, pointX);


	ofstream outputFile( filenameX, ios::out);
	outputFile.precision(16);

	for (int ii =0; ii<nrv; ii++){

		outputFile << pointX(ii)<<endln;
	}// for
 
	outputFile.close();



	return TCL_OK;

}





// command: runDP_RSM_SimTimeInvariantAnalysis -designPt dp.out  -output results.out  -ndir $n <-experimentalPointRule Uniform -gridInfo {-1  minY  maxY nPts 0  minY  maxY nPts0 1 minX1  maxX1 nPts1 2 minX2  maxX2 nPts2 ...}> 
//  -saveHessian hession.out <-surfaceDesign UnivariateDecomposition -simulation ImportanceSampling -tarCOV 0.1 -numSimulation 100000>

int 
TclReliabilityModelBuilder_runDP_RSM_SimTimeInvariantAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
 

	if ((theRandomNumberGenerator ==0)||(theFunctionEvaluator==0)||(theProbabilityTransformation==0) ||(theGradientEvaluator==0)) {
	
		opserr<<"theRandomNumberGenerator ==0)||(theGFunEvaluator)||(theProbabilityTransformation==0) ||(theGradientEvaluator==0)"<<endln;
		exit(-1);
	
	}
	int argvCounter = 1;
	
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	
	int numPPlane = 1;
	Vector * designPoint = new Vector(nrv);
	Vector * gridInfo = 0;
    Matrix * passedHessian =0;
	char outputFile[25];
	char * hessianFileName=0;
	char surfDesign[50] = "UnivariateDecomposition";
	char typeRespSurfaceSimulation[50]= "ImportanceSampling";
	double tarCov = 0.1;
	int numberOfSimulations = 1000000;
	char experimentalPointRule[40]="Uniform";




	while (argc > argvCounter) {

		if ((strcmp(argv[argvCounter],"-designPt") == 0)||(strcmp(argv[argvCounter],"-designPoint") == 0)) {
			char fileName[20];
			argvCounter++;
			strcpy(fileName,argv[argvCounter]);
			argvCounter++;
			ifstream inputFile( fileName, ios::in);

			if (inputFile.good()){
				int ii=0;
				double tmp;
				while(!inputFile.eof() && ii<nrv){ 
					inputFile >> tmp;
					(*designPoint)(ii)=tmp;
					ii++;
				}
				inputFile.close();
			} 
			else {
				opserr<<"ERROR: designpoint can not read from file "<<fileName<<endln;
				exit(-1);
			}
						
		}// else if

		else if (strcmp(argv[argvCounter],"-ndir") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &numPPlane) != TCL_OK) {
				opserr << "ERROR: invalid input: numPPlane \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if

		//====  <-surfaceDesign UnivariateDecomposition -simulation ImportanceSampling -tarCOV 0.1 -numSimulation 100000>
		else if (strcmp(argv[argvCounter],"-surfaceDesign") == 0) {
			argvCounter++;
			strcpy(surfDesign,argv[argvCounter]);
			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-experimentalPointRule") == 0) {
			argvCounter++;
			strcpy(experimentalPointRule,argv[argvCounter]);
			argvCounter++;
		}// else if


		else if (strcmp(argv[argvCounter],"-simulation") == 0) {
			argvCounter++;
			strcpy(typeRespSurfaceSimulation,argv[argvCounter]);
			argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-tarCOV") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &tarCov) != TCL_OK) {
				opserr << "ERROR: invalid input: tarCov \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-numSimulation") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &numberOfSimulations) != TCL_OK) {
				opserr << "ERROR: invalid input: numOfSimulations \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if

 

// ---------
		else if (strcmp(argv[argvCounter],"-output") == 0) {
			argvCounter++;
			strcpy(outputFile,argv[argvCounter]);
			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-gridInfo") == 0) {
			argvCounter++;

			int pathSize;
		    TCL_Char **pathStrings;
		  
		    if (Tcl_SplitList(interp, argv[argvCounter], 
					&pathSize, &pathStrings) != TCL_OK) {
				  
			  opserr << "WARNING problem splitting path list in gridInfo\n";
			  return 0;
			}
		  
		    gridInfo = new Vector(pathSize);

		    for (int i=0; i<pathSize; i++) {
			  double value;
			  if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
			    opserr << "WARNING problem reading path data value " << pathStrings[i] << "\n";

			    return 0;
			  }
			  (*gridInfo)(i) = value;
		  }  //for
		  // free up the array of pathsStrings .. see tcl man pages as to why
//		  cleanup(pathStrings);

			argvCounter++;
		}// else if
		else if (strcmp(argv[argvCounter],"-saveHessian") == 0){
			argvCounter++;	
			hessianFileName = new char[30];
			strcpy(hessianFileName,argv[argvCounter]);
			argvCounter++;
			
			ifstream inputFile( hessianFileName, ios::in);

			if (inputFile.good()){
				passedHessian = new Matrix(nrv,nrv);
				int ii=0;
				int jj=0;
				double tmp;
				while(ii<nrv){ 
					for (jj=0; jj<nrv;jj++){
						if (!inputFile.eof()){
							inputFile >> tmp;
							(*passedHessian)(jj,ii)=tmp;
						}
						else{
							opserr<<"-saveHessian, size of Hessian in file: "<<hessianFileName<< "is wrong"<<endln;
							exit(-1);
						}
					}	
					
					ii++;
				}// while 
				inputFile.close();

			} 
			else {
				opserr<<"warning: no data in hessianfrom file: "<<hessianFileName<<endln;
				
			}
		
		}
		else {
		
			opserr<< "unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}
	} //while

	
 


    DP_RSM_Sim * theDP_RSM_Sim = new DP_RSM_Sim(theReliabilityDomain,
					theFunctionEvaluator,
					theProbabilityTransformation,
					outputFile,
					theGradientEvaluator, 
					designPoint, 
					numPPlane, 
					experimentalPointRule,
					interp, 
					passedHessian, 
					hessianFileName, 
					surfDesign, 
					typeRespSurfaceSimulation, 
					gridInfo,
					theRandomNumberGenerator,
					tarCov,
					numberOfSimulations);

	if (passedHessian !=0) delete passedHessian;
	if (hessianFileName !=0) delete hessianFileName;
	theDP_RSM_Sim->analyze();

	if (gridInfo !=0) delete gridInfo;

    delete designPoint;
	return 0;

}






// command: runDP_RSM_SimTimeVariantAnalysis -designPt dp.out  -output results.out  -ndir $n <-experimentalPointRule Uniform -gridInfo {-1  minY  maxY nPts 0  minY  maxY nPts0 1 minX1  maxX1 nPts1 ..}> 
//  -saveHessian hession.out <-surfaceDesign UnivariateDecomposition -simulation ImportanceSampling -tarCOV 0.1 -numSimulation 100000 -littleDt dt -ImpulseInterval Dt>

int 
TclReliabilityModelBuilder_runDP_RSM_SimTimeVariantAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
 

	if ((theRandomNumberGenerator ==0)||(theFunctionEvaluator==0)||(theProbabilityTransformation==0) ||(theGradientEvaluator==0)) {
	
		opserr<<"theRandomNumberGenerator ==0)||(theGFunEvaluator)||(theProbabilityTransformation==0) ||(theGradientEvaluator==0)"<<endln;
		exit(-1);
	
	}
	int argvCounter = 1;
	
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	
	int numPPlane = 1;
	Vector * designPoint = new Vector(nrv);
	Vector * gridInfo = 0;
    Matrix * passedHessian =0;
	char outputFile[25];
	char * hessianFileName=0;
	char surfDesign[50] = "UnivariateDecomposition";
	char typeRespSurfaceSimulation[50]= "ImportanceSampling";
	double tarCov = 0.1;
	int numberOfSimulations = 1000000;
	char experimentalPointRule[40]="Uniform";

	double littleDt = 0.0;
	double ImpulseInterval =0.0;


	while (argc > argvCounter) {

		if ((strcmp(argv[argvCounter],"-designPt") == 0)||(strcmp(argv[argvCounter],"-designPoint") == 0)) {
			char fileName[20];
			argvCounter++;
			strcpy(fileName,argv[argvCounter]);
			argvCounter++;
			ifstream inputFile( fileName, ios::in);

			if (inputFile.good()){
				int ii=0;
				double tmp;
				while(!inputFile.eof() && ii<nrv){ 
					inputFile >> tmp;
					(*designPoint)(ii)=tmp;
					ii++;
				}
				inputFile.close();
			} 
			else {
				opserr<<"ERROR: designpoint can not read from file "<<fileName<<endln;
				exit(-1);
			}
						
		}// else if

		else if (strcmp(argv[argvCounter],"-ndir") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &numPPlane) != TCL_OK) {
				opserr << "ERROR: invalid input: numPPlane \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if

		//====  <-surfaceDesign UnivariateDecomposition -simulation ImportanceSampling -tarCOV 0.1 -numSimulation 100000>
		else if (strcmp(argv[argvCounter],"-surfaceDesign") == 0) {
			argvCounter++;
			strcpy(surfDesign,argv[argvCounter]);
			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-experimentalPointRule") == 0) {
			argvCounter++;
			strcpy(experimentalPointRule,argv[argvCounter]);
			argvCounter++;
		}// else if


		else if (strcmp(argv[argvCounter],"-simulation") == 0) {
			argvCounter++;
			strcpy(typeRespSurfaceSimulation,argv[argvCounter]);
			argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-tarCOV") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &tarCov) != TCL_OK) {
				opserr << "ERROR: invalid input: tarCov \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-numSimulation") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &numberOfSimulations) != TCL_OK) {
				opserr << "ERROR: invalid input: numOfSimulations \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if

 

// ---------
		else if (strcmp(argv[argvCounter],"-output") == 0) {
			argvCounter++;
			strcpy(outputFile,argv[argvCounter]);
			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-gridInfo") == 0) {
			argvCounter++;

			int pathSize;
		    TCL_Char **pathStrings;
		  
		    if (Tcl_SplitList(interp, argv[argvCounter], 
					&pathSize, &pathStrings) != TCL_OK) {
				  
			  opserr << "WARNING problem splitting path list in gridInfo\n";
			  return 0;
			}
		  
		    gridInfo = new Vector(pathSize);

		    for (int i=0; i<pathSize; i++) {
			  double value;
			  if ( Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
			    opserr << "WARNING problem reading path data value " << pathStrings[i] << "\n";

			    return 0;
			  }
			  (*gridInfo)(i) = value;
		  }  //for
		  // free up the array of pathsStrings .. see tcl man pages as to why
//		  cleanup(pathStrings);

			argvCounter++;
		}// else if
		else if (strcmp(argv[argvCounter],"-saveHessian") == 0){
			argvCounter++;	
			hessianFileName = new char[30];
			strcpy(hessianFileName,argv[argvCounter]);
			argvCounter++;
			
			ifstream inputFile( hessianFileName, ios::in);

			if (inputFile.good()){
				passedHessian = new Matrix(nrv,nrv);
				int ii=0;
				int jj=0;
				double tmp;
				while(ii<nrv){ 
					for (jj=0; jj<nrv;jj++){
						if (!inputFile.eof()){
							inputFile >> tmp;
							(*passedHessian)(jj,ii)=tmp;
						}
						else{
							opserr<<"-saveHessian, size of Hessian in file: "<<hessianFileName<< "is wrong"<<endln;
							exit(-1);
						}
					}	
					
					ii++;
				}// while 
				inputFile.close();

			} 
			else {
				opserr<<"warning: no data in hessianfrom file: "<<hessianFileName<<endln;
				
			}
		
		}
		
		else if (strcmp(argv[argvCounter],"-littleDt") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &littleDt) != TCL_OK) {
				opserr << "ERROR: invalid input: littleDt \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-ImpulseInterval") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &ImpulseInterval) != TCL_OK) {
				opserr << "ERROR: invalid input: ImpulseInterval \n";
				return TCL_ERROR;
			}

			argvCounter++;
		}// else if
		else {
		
			opserr<< "unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}
	} //while

	
 
	if ((littleDt ==0) ||(ImpulseInterval ==0) ){
		opserr<< "not enough parameters. command: "<<endln;
		opserr<<"runDP_RSM_SimTimeVariantAnalysis -designPt dp.out  -output results.out  -ndir $n <-experimentalPointRule Uniform -gridInfo {-1  minY  maxY nPts 0  minY  maxY nPts0 1 minX1  maxX1 nPts1 ..}>"<<endln; 
        opserr<<"-saveHessian hession.out <-surfaceDesign UnivariateDecomposition -simulation ImportanceSampling -tarCOV 0.1 -numSimulation 100000> -littleDt dt -ImpulseInterval Dt" <<endln;
		exit(-1);
	}


    DP_RSM_Sim_TimeVariant * theDP_RSM_Sim_TimeVariant = new DP_RSM_Sim_TimeVariant(theReliabilityDomain,
					theFunctionEvaluator,
					theProbabilityTransformation,
					outputFile,
					theGradientEvaluator, 
					designPoint, 
					numPPlane, 
					experimentalPointRule,
					interp, 
					passedHessian, 
					hessianFileName, 
					surfDesign, 
					typeRespSurfaceSimulation, 
					gridInfo,
					theRandomNumberGenerator,
					tarCov,
					numberOfSimulations,
					littleDt,
					ImpulseInterval);

	if (passedHessian !=0) delete passedHessian;
	if (hessianFileName !=0) delete hessianFileName;
	theDP_RSM_Sim_TimeVariant->analyze();


    delete designPoint;
	if (gridInfo !=0) delete gridInfo;
	return 0;

}


int 
TclReliabilityModelBuilder_getBetaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "ERROR: Invalid number of arguments to getBetaFORM command." << endln;
    return TCL_ERROR;
  }

  if (theFunctionEvaluator == 0) {
    opserr << "WARNING betaFORM -- no function evaluator defined\n";
    return TCL_ERROR;	        
  }

  int lsfTag;
  if (Tcl_GetInt(interp, argv[1], &lsfTag) != TCL_OK) {
    opserr << "WARNING betaFORM lsfTag? - could not read lsfTag\n";
    return TCL_ERROR;	        
  }   

  double beta = theFunctionEvaluator->getResponseVariable("betaFORM",lsfTag);

  char buffer[40];
  sprintf(buffer,"%35.20f",beta);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int 
TclReliabilityModelBuilder_getGammaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "ERROR: Invalid number of arguments to getGammaFORM command." << endln;
    return TCL_ERROR;
  }

  if (theFunctionEvaluator == 0) {
    opserr << "WARNING gammaFORM -- no function evaluator defined\n";
    return TCL_ERROR;	        
  }

  int lsfTag;
  if (Tcl_GetInt(interp, argv[1], &lsfTag) != TCL_OK) {
    opserr << "WARNING gammaFORM lsfTag? rvTag? - could not read lsfTag\n";
    return TCL_ERROR;	        
  }

  int rvTag;
  if (Tcl_GetInt(interp, argv[2], &rvTag) != TCL_OK) {
    opserr << "WARNING gammaFORM lsfTag? rvTag? - could not read rvTag\n";
    return TCL_ERROR;	        
  }   
    
  RandomVariable *theRV =
    theReliabilityDomain->getRandomVariablePtr(rvTag);
    
  if (theRV == 0) {
    opserr << "WARNING gammaFORM RV with tag " << rvTag << " not found\n";
    return TCL_ERROR;	        
  }

  int index = theReliabilityDomain->getRandomVariableIndex(rvTag);
  
  double gamma = theFunctionEvaluator->getResponseVariable("gammaFORM",lsfTag,index);

  char buffer[40];  
  sprintf(buffer,"%35.20f",gamma);
    
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int 
TclReliabilityModelBuilder_getAlphaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "ERROR: Invalid number of arguments to getAlphaFORM command." << endln;
    return TCL_ERROR;
  }

  if (theFunctionEvaluator == 0) {
    opserr << "WARNING alphaFORM -- no function evaluator defined\n";
    return TCL_ERROR;	        
  }

  int lsfTag;
  if (Tcl_GetInt(interp, argv[1], &lsfTag) != TCL_OK) {
    opserr << "WARNING alphaFORM lsfTag? rvTag? - could not read lsfTag\n";
    return TCL_ERROR;	        
  }

  int rvTag;
  if (Tcl_GetInt(interp, argv[2], &rvTag) != TCL_OK) {
    opserr << "WARNING alphaFORM lsfTag? rvTag? - could not read rvTag\n";
    return TCL_ERROR;	        
  }   
    
  RandomVariable *theRV =
    theReliabilityDomain->getRandomVariablePtr(rvTag);
    
  if (theRV == 0) {
    opserr << "WARNING alphaFORM RV with tag " << rvTag << " not found\n";
    return TCL_ERROR;	        
  }

  int index = theReliabilityDomain->getRandomVariableIndex(rvTag);
  
  double alpha = theFunctionEvaluator->getResponseVariable("alphaFORM",lsfTag,index);

  char buffer[40];  
  sprintf(buffer,"%35.20f",alpha);
    
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int
TclReliabilityModelBuilder_getPDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "WARNING getPDF tag? x? -- insufficient number of arguments\n";
    return TCL_ERROR;
  }

  int rvTag;
  double x;
  if (Tcl_GetInt(interp, argv[1], &rvTag) != TCL_OK) {
    opserr << "WARNING getPDF tag? x? -- could not read tag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &x) != TCL_OK) {
    opserr << "WARNING getPDF tag? x? -- could not read x\n";
    return TCL_ERROR;	        
  }

  RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (theRV == 0) {
    opserr << "WARNING getPDF tag? x? -- random variable with tag "
	   << rvTag << " does not exist in model\n";
    return TCL_ERROR;
  }

  char buffer[40];

  sprintf(buffer,"%35.20f", theRV->getPDFvalue(x));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
TclReliabilityModelBuilder_getStdNormalPDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING getStdNormalPDF x? -- insufficient number of arguments\n";
    return TCL_ERROR;
  }

  double x;
  if (Tcl_GetDouble(interp, argv[1], &x) != TCL_OK) {
    opserr << "WARNING getStdNormalPDF x? -- could not read x\n";
    return TCL_ERROR;	        
  }

  NormalRV theRV(0,0,1);

  char buffer[40];

  sprintf(buffer,"%35.20f", theRV.getPDFvalue(x));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
TclReliabilityModelBuilder_getCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "WARNING getCDF tag? x? -- insufficient number of arguments\n";
    return TCL_ERROR;
  }

  int rvTag;
  double x;
  if (Tcl_GetInt(interp, argv[1], &rvTag) != TCL_OK) {
    opserr << "WARNING getCDF tag? x? -- could not read tag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &x) != TCL_OK) {
    opserr << "WARNING getCDF tag? x? -- could not read x\n";
    return TCL_ERROR;	        
  }

  RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (theRV == 0) {
    opserr << "WARNING getCDF tag? x? -- random variable with tag "
	   << rvTag << " does not exist in model\n";
    return TCL_ERROR;
  }

  char buffer[40];

  sprintf(buffer,"%35.20f", theRV->getCDFvalue(x));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
TclReliabilityModelBuilder_getStdNormalCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING getStdNormalCDF x? -- insufficient number of arguments\n";
    return TCL_ERROR;
  }

  double x;
  if (Tcl_GetDouble(interp, argv[1], &x) != TCL_OK) {
    opserr << "WARNING getStdNormalCDF x? -- could not read x\n";
    return TCL_ERROR;	        
  }

  NormalRV theRV(0,0,1);

  char buffer[40];

  sprintf(buffer,"%35.20f", theRV.getCDFvalue(x));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int
TclReliabilityModelBuilder_getInverseCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "WARNING getInverseCDF tag? p? -- insufficient number of arguments\n";
    return TCL_ERROR;
  }

  int rvTag;
  double x;
  if (Tcl_GetInt(interp, argv[1], &rvTag) != TCL_OK) {
    opserr << "WARNING getInverseCDF tag? p? -- could not read tag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &x) != TCL_OK) {
    opserr << "WARNING getInverseCDF tag? p? -- could not read p\n";
    return TCL_ERROR;	        
  }

  RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
  if (theRV == 0) {
    opserr << "WARNING getInverseCDF tag? p? -- random variable with tag "
	   << rvTag << " does not exist in model\n";
    return TCL_ERROR;
  }

  char buffer[40];

  sprintf(buffer,"%35.20f", theRV->getInverseCDFvalue(x));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
TclReliabilityModelBuilder_getStdNormalInverseCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING getStdNormalInverseCDF p? -- insufficient number of arguments\n";
    return TCL_ERROR;
  }

  double x;
  if (Tcl_GetDouble(interp, argv[1], &x) != TCL_OK) {
    opserr << "WARNING getStdNormalInverseCDF p? -- could not read p\n";
    return TCL_ERROR;	        
  }

  NormalRV theRV(0,0,1);

  char buffer[40];

  sprintf(buffer,"%35.20f", theRV.getInverseCDFvalue(x));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int
TclReliabilityModelBuilder_getRVTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  RandomVariable *theEle;
  RandomVariableIter &eleIter = theReliabilityDomain->getRandomVariables();
  
  char buffer[20];
  
  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }
  
  return TCL_OK;
}


int
TclReliabilityModelBuilder_getLSFTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  LimitStateFunction *theEle;
  LimitStateFunctionIter &eleIter = theReliabilityDomain->getLimitStateFunctions();
  
  char buffer[20];
  
  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }
  
  return TCL_OK;
}

int
TclReliabilityModelBuilder_getCutsetTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Cutset *theEle;
  CutsetIter &eleIter = theReliabilityDomain->getCutsets();
  
  char buffer[20];
  
  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }
  
  return TCL_OK;
}

int
TclReliabilityModelBuilder_getCutsetComponents(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING getCutsetComponents tag? -- insufficient number of arguments\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING getCutsetComponents tag? -- could not read tag\n";
    return TCL_ERROR;
  }

  Cutset *theCutset = theReliabilityDomain->getCutsetPtr(tag);
  if (theCutset == 0) {
    opserr << "WARNING getCutsetComponents tag? -- cutset with tag "
	   << tag << " does not exist in model\n";
    return TCL_ERROR;
  }
  
  int numComponents = theCutset->getNumberOfComponents();
  const Vector &theComponents = theCutset->getComponents();

  char buffer[20];
  for (int i = 0; i < numComponents; i++) {
    sprintf(buffer, "%d ", int(theComponents(i)));
    Tcl_AppendResult(interp, buffer, NULL);
  }
  
  return TCL_OK;
}


/////////////////////////////////////////////////////////
/// (from here to the end) added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addAnalyzer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int argvCounter;
	int nstep;
	double delta;
	bool print=false;
	bool defined=false;
	bool initialstat=false;

	if(theAnalyzer!=0){
		delete theAnalyzer;
		theAnalyzer = 0;
	}
	argvCounter=1;
	while(argvCounter<argc){
		if (strcmp(argv[argvCounter],"-definedabove") == 0) {
			defined=true;
			argvCounter++;
		}
		else if(strcmp(argv[argvCounter],"-print") == 0) {
			print=true;
			argvCounter++;
		}
		else if(strcmp(argv[argvCounter],"-initialstatic") == 0) {
			initialstat=true;
			argvCounter++;
		}
		else if(strcmp(argv[argvCounter],"-delta") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &delta) != TCL_OK) {
				opserr << "Invalid Input for ratio \n";
				opserr << "for FreeVibration in tclModelbuidler \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if(strcmp(argv[argvCounter],"-step") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &nstep) != TCL_OK) {
				opserr << "Invalid Input for ratio \n";
				opserr << "for FreeVibration in tclModelbuidler \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else{
			opserr << "ERROR: Invalid argument to FreeVibration " << endln;
			return TCL_ERROR;
		}
	}
	if(initialstat){
		if(theInitialStaticAnalysis==0){
			opserr << "Need to define InitialShapeAnalysis \n";
			opserr << "before defining Analyzer with -initialstat\n";
			return TCL_ERROR;
		}
	}
	if (theReliabilityDomain == 0 ) {
		opserr << "Need ReliabilityDomain before an Analyzer can be created" << endln;
		return TCL_ERROR;
	}
	if (theStructuralDomain == 0 ) {
		opserr << "Need StructuralDomain before an Analyzer can be created" << endln;
		return TCL_ERROR;
	}
	if (theReliabilityTransientAnalysis == 0 && theReliabilityStaticAnalysis==0  ) {
		opserr << "Need Analysis before an Analyzer can be created" << endln;
		return TCL_ERROR;
	}
//	if (theSensitivityAlgorithm== 0 ) {
//		opserr << "Need SensitivityAlgorithm before an Analyzer can be created" << endln;
//		return TCL_ERROR;
//	}
//	if (theSensitivityIntegrator== 0 ) {
//		opserr << "Need SensitivityIntegrator before an Analyzer can be created" << endln;
//		return TCL_ERROR;
//	}
	if(defined){
		int numLoadPatterns=0;
		int *LoadPatterns=0;
		if(theReliabilityTransientAnalysis !=0){
			theAnalyzer = new DynamicAnalyzer
						  (theReliabilityDomain,
						   theStructuralDomain,
						   theInitialStaticAnalysis,
						   theReliabilityTransientAnalysis,
						   theSensitivityAlgorithm,
						   theSensitivityIntegrator,
						   nstep,
						   delta,
						   numLoadPatterns,
						   LoadPatterns,
						   print);
		}
		else
		{
			theAnalyzer = new StaticAnalyzer
						  (theReliabilityDomain,
						   theStructuralDomain,
						   theInitialStaticAnalysis,
						   theReliabilityStaticAnalysis,
						   theSensitivityAlgorithm,
						   theSensitivityIntegrator,
						   nstep,
						   delta,
						   numLoadPatterns,
						   LoadPatterns,
						   print);
		}
	}
	else {
		opserr << "ERROR: Invalid argument to Analyzer " << endln;
		opserr << "-definedabove is required " << endln;
		return TCL_ERROR;
	}
	return TCL_OK;
}


int 
TclReliabilityModelBuilder_addInitialStaticAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int* StaticLoadPatterns=0;
	int* temploads=0;
	int numLoadPatterns=0;
	int nstep=0;
	int loadtag=0;
    LoadPattern *thePattern=0;
	bool loadfound= false;
	int argvCounter = 1;
	bool print = false;
	temploads= new int[100];

	if(theInitialStaticAnalysis!=0){
		delete theInitialStaticAnalysis;
		theInitialStaticAnalysis=0;
	}

	if (strcmp(argv[argvCounter],"-selectLoad") == 0) {  
		argvCounter++;
		while(argvCounter<argc){
			if (strcmp(argv[argvCounter],"-print") == 0) {
				print=true;
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-nstep") == 0) {
				argvCounter++;
				if (Tcl_GetInt(interp, argv[argvCounter], &nstep) != TCL_OK) {
				opserr << "ERROR: Invalid input";
				opserr << " nstep for initial static analysis" << endln;
				return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-loads") == 0) {
				argvCounter++;
				numLoadPatterns=0;
				while(argvCounter<argc){
					if(argv[argvCounter][0]!= '-'){		
						if (Tcl_GetInt(interp, argv[argvCounter], &loadtag) != TCL_OK) {
						opserr << "Error invalid input for";
						opserr << " LoadPattern ID for the initial static analysis";
						opserr << endln;
						return TCL_ERROR;
						}
						argvCounter++;
						LoadPatternIter& thePatterns = theStructuralDomain->getLoadPatterns();
						loadfound = false;
						while((thePattern = thePatterns()) != 0){
							int tag=thePattern->getTag();
							if( tag == loadtag ) {
								loadfound = true;
								break;
							}
						}
						if(loadfound){
							numLoadPatterns++;
							temploads[numLoadPatterns-1]=loadtag;
						}
					}
					else break;
				}
				if( numLoadPatterns != 0 ){
					StaticLoadPatterns = new int[numLoadPatterns];
					for (int i=0; i<numLoadPatterns; i++) {
						StaticLoadPatterns[i]=temploads[i];
					}
				}
			}
		}
		if (theStructuralDomain== 0 ) {
			opserr << "Need StructuralDomain before a InitialStaticAnalysis can be created" << endln;
			return TCL_ERROR;
		}
		if (theReliabilityDomain== 0 ) {
			opserr << "Need ReliabilityDomain before a InitialStaticAnalysis can be created" << endln;
			return TCL_ERROR;
		}
		theInitialStaticAnalysis = new SelectLoadInitialStaticAnalysis 
									(theReliabilityDomain,
									 theStructuralDomain,
									 nstep,
									 numLoadPatterns,
									 StaticLoadPatterns,
									 print);

	} else if (strcmp(argv[argvCounter],"-file") == 0) {  

		opserr << " FATAL error \n";
		opserr << " -file option for InitialShapeBuilder ";
		opserr << " is not yete implemented \n";
		exit(-1);
//		argvCounter++;
//		ifstream inputFile( argv[argvCounter], ios::in );
//		if (inputFile.fail()) {
//			opserr << "File " << *argv[2] << " could not be opened. " << endln;
//			return TCL_ERROR;
//		}
//		argvCounter++;
//		inputFile.close();
//		if(argvCounter<argc){
//			if (strcmp(argv[argvCounter],"-print") == 0) print=true;
//			else{
//				opserr << "Error invalid input";
//			}
//		}
//
//		theInitialStaticAnalysis = new tclStaticAnalysis
//									  (interp,
//									   theReliabilityDomain,
//									   theStructuralDomain,
//									   argv[argvCounter],
//									   print);
	} else {
		opserr << "ERROR: Invalid argument to Initial Static Analysis. " << endln;
		return TCL_ERROR;
	}
	delete [] StaticLoadPatterns; 
	StaticLoadPatterns=0;
	delete [] temploads;
	temploads=0;
	return TCL_OK;
}
/////////////////////////////////////////////////////////
/// added by K Fujimura for Random Vibration Analysis ///
/////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_addInitialPointBuilder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	bool print=false;
	double eps=0.01;
	int MaxLineSearch=20;

	if(theInitialPointBuilder!=0){
		delete theInitialPointBuilder;
		theInitialPointBuilder=0;
	}
	if (theReliabilityDomain == 0 ) {
		opserr << "Need ReliabilityDomain before an RandomVibrationAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if(theFunctionEvaluator==0){
		opserr << "Need to define GFunEvaluator \n";
		opserr << "before defining MirrorImageBuilder\n";
		return TCL_ERROR;
	}

 	if (strcmp(argv[1],"-mirror") == 0) {
		// mirror image initial point builder //
		opserr << "Invalid Input for Initialpoint \n";
		opserr << "-mirror can not be selected \n";
		return TCL_ERROR;
//		argvCounter=2;
//		while(argvCounter<argc){
//			if (strcmp(argv[argvCounter],"-print") == 0) {
//				print=true;
//				argvCounter++;
//			}
//			else if (strcmp(argv[argvCounter],"-eps") == 0) {
//				argvCounter++;
//				if (Tcl_GetDouble(interp, argv[argvCounter], &eps) != TCL_OK) {
//					opserr << "Invalid Input for threshold \n";
//					opserr << "for InitialPointBuilder in tclModelbuidler \n";
//					return TCL_ERROR;
//				}
//				argvCounter++;
//			}
//			else if (strcmp(argv[argvCounter],"-maxlinesearch") == 0) {
//				argvCounter++;
//				if (Tcl_GetInt(interp, argv[argvCounter], &MaxLineSearch) != TCL_OK) {
//					opserr << "Invalid Input for MaxLineSearch \n";
//					opserr << "for InitialPointBuilder in tclModelbuidler \n";
//					return TCL_ERROR;
//				}
//				argvCounter++;
//			}
//			else{
//				opserr << "ERROR: Invalid argument to InitialPointBuilder " << endln;
//				return TCL_ERROR;
//			}
//		}
//		theInitialPointBuilder = new MirrorImageInitialPointBuilder
//								(theReliabilityDomain,
//								 theFunctionEvaluator,
//								 eps,
//								 MaxLineSearch,
//								 print);
//
	} else if (strcmp(argv[1],"-threshold") == 0) {

		int maxDivide=10;
		bool start_mirror=true;
		double mirroreps=0.01;
		int argvCounter=2;
		while(argvCounter<argc){
			if (strcmp(argv[argvCounter],"-print") == 0) {
				print=true;
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-eps") == 0) {
				argvCounter++;
				if (Tcl_GetDouble(interp, argv[argvCounter], &eps) != TCL_OK) {
					opserr << "Invalid Input for threshold \n";
					opserr << "for InitialPointBuilder in tclModelbuidler \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-maxdivide") == 0) {
				argvCounter++;
				if (Tcl_GetInt(interp, argv[argvCounter], &maxDivide) != TCL_OK) {
					opserr << "Invalid Input for threshold \n";
					opserr << "for InitialPointBuilder in tclModelbuidler \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-startpoint") == 0) {
				argvCounter++;
				if(strcmp(argv[argvCounter],"none")==0){
					start_mirror=false;
				}else if(strcmp(argv[argvCounter],"mirror")==0){
					opserr << "Invalid Input for -stattpoint for initialpoint \n";
					opserr << "mirror can not be selected\n";
					return TCL_ERROR;
					start_mirror=true;
					argvCounter++;
					if (Tcl_GetDouble(interp, argv[argvCounter], &mirroreps) != TCL_OK) {
						opserr << "Invalid Input for threshold \n";
						opserr << "for InitialPointBuilder in tclModelbuidler \n";
						return TCL_ERROR;
					}
					argvCounter++;
					if (Tcl_GetInt(interp, argv[argvCounter], &MaxLineSearch) != TCL_OK) {
						opserr << "Invalid Input for threshold \n";
						opserr << "for InitialPointBuilder in tclModelbuidler \n";
						return TCL_ERROR;
					}
				}else{
					opserr << "Invalid Input for threshold \n";
					opserr << "for InitialPointBuilder in tclModelbuidler \n";
					return TCL_ERROR;
				}
				argvCounter++;
			}
			else{
				opserr << "ERROR: Invalid argument to InitialPointBuilder " << endln;
				return TCL_ERROR;
			}
		}

		if (theFindDesignPointAlgorithm == 0 ) {
			opserr << "Need theNewSearchWithStepSizeAndStepDirection before a ThresholdIncInitialPointBuilder can be created" << endln;
			return TCL_ERROR;
		}

		theInitialPointBuilder = new ThresholdIncInitialPointBuilder
								(theReliabilityDomain,
								 theFunctionEvaluator,
								 theFindDesignPointAlgorithm,
								 maxDivide,
								 eps,
								 start_mirror,
								 MaxLineSearch,
								 mirroreps, 
								 print);
	}else{
		opserr << "ERROR: Invalid argument to InitialPointBuilder " << endln;
		opserr << argv[1]<< endln;
		return TCL_ERROR;
	}
	return TCL_OK;
}


int 
TclReliabilityModelBuilder_addCrossingRateAnalyzer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (theCrossingRateAnalyzer != 0) {
		delete theCrossingRateAnalyzer;
		theCrossingRateAnalyzer = 0;
	}

	int analysisType=2;
	double littleDt=0.1;
	bool print=false;

	int argvCounter = 1;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-littledt") == 0) {
			argvCounter++;
			  // GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &littleDt) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-koo") == 0) {
			argvCounter++;
			analysisType = 2;
		}
		else if (strcmp(argv[argvCounter],"-twosearches") == 0) {
			argvCounter++;
			analysisType = 1;
		}
		else if (strcmp(argv[argvCounter],"-print") == 0) {
			argvCounter++;
			print=true;
		}
		else {
			opserr << "ERROR: Invalid input to theOutCrossingAnalysis." << endln;
			return TCL_ERROR;
		}
	}
	if(analysisType==1){
		if (theReliabilityDomain== 0 ) {
		opserr << "Need ReliabilityDomain before a CrossingRateAnalyzer can be created" << endln;
		return TCL_ERROR;
		}
		if (theFindDesignPointAlgorithm== 0 ) {
		opserr << "Need FindDesignPointAlgorithm before a CrossingRateAnalyzer can be created" << endln;
		return TCL_ERROR;
		}
		if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an CrossingRateAnalyzer can be created" << endln;
		return TCL_ERROR;
		}
		if (theGradientEvaluator == 0 ) {
		opserr << "Need theGradientEvaluator before an CrossingRateAnalyzer can be created" << endln;
		return TCL_ERROR;
		}
	}
	theCrossingRateAnalyzer = new CrossingRateAnalyzer
								(theReliabilityDomain,
								 theFindDesignPointAlgorithm,
								 theFunctionEvaluator,
								 theGradientEvaluator,
								 analysisType,
						         littleDt,
  								 print);

	if (theCrossingRateAnalyzer == 0) {
		opserr << "ERROR: could not create theCrossingRateAnalyzer\n";
		return TCL_ERROR;
	}
	return TCL_OK;
}


int 
TclReliabilityModelBuilder_addFOSeriesSimulation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (theFOSeriesSimulation != 0) {
		delete theFOSeriesSimulation;
		theFOSeriesSimulation = 0;
	}

	int MaxSim=10000;
	int Interval=100;
	double Eps=0.05;
	bool twoside=false;
	int analysis=-1;
	bool print=false;

	int argvCounter = 1;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-maxsim") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &MaxSim) != TCL_OK) {
				opserr << "ERROR: invalid input maxsim to FOSeriesSimulation \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-interval") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &Interval) != TCL_OK) {
				opserr << "ERROR: invalid input interval to FOSeriesSimulation \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-analysis") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &analysis) != TCL_OK) {
				opserr << "ERROR: invalid input analysis to FOSeriesSimulation \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-eps") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &Eps) != TCL_OK) {
				opserr << "ERROR: invalid input eps to FOSeriesSimulation \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-twoside") == 0) {
			argvCounter++;
			twoside=true;
		}
		else if (strcmp(argv[argvCounter],"-print") == 0) {
			argvCounter++;
			print=true;
		}
		else {
			opserr << "ERROR: Invalid input to theOutCrossingAnalysis." << endln;
			return TCL_ERROR;
		}
	}

	if(analysis!=0&&analysis!=1&&analysis!=2){
		opserr << "ERROR: analysisType must be either of 0, 1, or 2" << endln;
		return TCL_ERROR;
	}

	theFOSeriesSimulation= new FOSeriesSimulation(MaxSim,
												  Interval,
												  Eps,
												  twoside,
												  analysis,
												  print);


	if (theFOSeriesSimulation == 0) {
		opserr << "ERROR: could not create theFOSeriesSimulation\n";
		return TCL_ERROR;
	}
	return TCL_OK;
}


int 
TclReliabilityModelBuilder_addFirstPassageAnalyzer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (theFirstPassageAnalyzer != 0) {
		delete theFirstPassageAnalyzer;
		theFirstPassageAnalyzer= 0;
	}

	int analysis=1;
	int interpolation=1;
	bool print=false;
	bool twoside=true;
	int stationary=0;

	int argvCounter = 1;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-analysis") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &analysis) != TCL_OK) {
				opserr << "ERROR: invalid input analysis to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-interpolation") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &interpolation) != TCL_OK) {
				opserr << "ERROR: invalid input interpolation to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-stationary") == 0) {
			stationary=1;
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-nonstationary") == 0) {
			stationary=2;
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-twoside") == 0) {
			argvCounter++;
			int ind;
			if (Tcl_GetInt(interp, argv[argvCounter], &ind) != TCL_OK) {
				opserr << "ERROR: invalid input twoside to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			if(ind==0) twoside=false;
			else twoside=true;
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-print") == 0) {
			argvCounter++;
			print=true;
		}
		else {
			opserr << "ERROR: Invalid input to theFirstPassageAnalyzer" << endln;
			return TCL_ERROR;
		}
	}
	if(stationary==0){
			opserr << "ERROR: Need to specify stationary/nonstationary" << endln;
			opserr << "for FirstPassageAnalyzer" << endln;
			return TCL_ERROR;
	}

	if (theReliabilityDomain == 0) {
		opserr << "Need theReliabilityDomain before an FirstPassageAnalyzer can be created" << endln;
		return TCL_ERROR;
	}
	if (theFindDesignPointAlgorithm == 0 ) {
		opserr << "Need theGFunEvaluator before an FirstPassageAnalyzer can be created" << endln;
		return TCL_ERROR;
	}
	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an FirstPassageAnalyzer can be created" << endln;
		return TCL_ERROR;
	}

	if(stationary==1){
		theFirstPassageAnalyzer=new StatFirstPassageAnalyzer
									(theReliabilityDomain,
									 theFindDesignPointAlgorithm,
									 theFunctionEvaluator,
									 theFOSeriesSimulation,
									 analysis,
									 twoside,
									 print);
	}else{
		theFirstPassageAnalyzer=new NonStatFirstPassageAnalyzer
								   (theReliabilityDomain,
									theFindDesignPointAlgorithm,
									theFunctionEvaluator,
									theFOSeriesSimulation,
									analysis,
									interpolation,
									twoside,
									print);
	}

	if (theFirstPassageAnalyzer == 0) {
		opserr << "ERROR: could not create theFirstPassageAnalyzer\n";
		return TCL_ERROR;
	}
	return TCL_OK;
}


int 
TclReliabilityModelBuilder_addRandomVibrationSimulation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theRandomVibrationSimulation != 0) {
		delete theRandomVibrationSimulation;
		theRandomVibrationSimulation = 0;
	}
	if (theReliabilityDomain == 0 ) {
		opserr << "Need ReliabilityDomain before an RandomVibrationAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theStructuralDomain == 0 ) {
		opserr << "Need StructuralDomain before an RandomVibrationAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a Outcrossing can be created" << endln;
		///////////// Modified by K Fujimura /////////////////////////////////////////////
		opserr << "Assume all RV's are independent" << endln;
		theProbabilityTransformation = 
		new AllIndependentTransformation(theReliabilityDomain,0);
		///////////// Modified by K Fujimura /////////////////////////////////////////////
	}
	double StartTime=1.0;
	double EndTime=20.0;
	double TimeInterval=1.0;
	double FragMin=1.0;
	double FragInt=0.0;
	int nFrag=1;

	int maxSim=100000;
	int intervalSim=200;
	double eps=0.05;

	int instantaneous=0;
	int firstpassage=0;

	int istationary=0;
	bool stationary;
	double stationaryTime=0.0;
	char *fileBinary=new char[100];
	double sampleAmp=0.0;
	double sampleTime=0.0;
	bool print=false;
	bool twoside=false;
	bool system=false;

	// Loop through arguments
	int argvCounter = 2;

	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-starttime") == 0) {
			argvCounter++;
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetDouble(interp, argv[argvCounter], &StartTime) != TCL_OK) {
				opserr << "ERROR: invalid input stepsToStart to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-endtime") == 0) {
			argvCounter++;
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetDouble(interp, argv[argvCounter], &EndTime) != TCL_OK) {
				opserr << "ERROR: invalid input stepsToEnd to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-interval") == 0) {
			argvCounter++;
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetDouble(interp, argv[argvCounter], &TimeInterval) != TCL_OK) {
				opserr << "ERROR: invalid input sampleFreq to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-fragility") == 0) {
			argvCounter++;
			  // GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &FragMin) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &FragInt) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &nFrag) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-maxsim") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &maxSim) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-intervalsim") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &intervalSim) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-instantaneous") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &instantaneous) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-firstpassage") == 0) {
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &firstpassage) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-eps") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &eps) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-stationarytime") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &stationaryTime) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-print") == 0) {
			argvCounter++;
			print=true;
		}
		else if (strcmp(argv[argvCounter],"-sample") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &sampleAmp) != TCL_OK) {
				opserr << "ERROR: invalid input samplAmp theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &sampleTime) != TCL_OK) {
				opserr << "ERROR: invalid input samplAmp theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-twoside") == 0) {
			argvCounter++;
			twoside=true;
		}
		else if (strcmp(argv[argvCounter],"-system") == 0) {
			argvCounter++;
			system=true;
		}
		else if (strcmp(argv[argvCounter],"-stationary") == 0) {
			istationary=1;
			stationary=true;
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-nonstationary") == 0) {
			istationary=2;
			stationary=false;
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-binary") == 0) {
			argvCounter++;
			strcpy(fileBinary,argv[argvCounter]);
			argvCounter++;
		}
		else {
			opserr << "ERROR: Invalid input to theOutCrossingAnalysis." << endln;
			return TCL_ERROR;
		}
	}

	if(istationary==0){
		opserr << "ERROR: Need to specify stationary" << endln;
		opserr << "for RandomVibrationSimulation" << endln;
		return TCL_ERROR;
	}
	if(instantaneous==0&&firstpassage==0){
		opserr << "ERROR: Need to specify either instantaneous or firstpassage" << endln;
		opserr << "for RandomVibrationSimulation" << endln;
		return TCL_ERROR;
	}
	if(stationary&&stationaryTime<0.0){
		opserr << "ERROR: Need to specify stationaryTime" << endln;
		opserr << "for RandomVibrationSimulator" << endln;
		return TCL_ERROR;
	}

	opserr << "=========================================\n";
	opserr << "\n";
	opserr << "     RandomVibrationSimulation  \n";
	opserr << "\n";
	opserr << "=========================================\n";
	opserr << "\n";
	opserr << " StartTime.........................." << StartTime << "\n";
	opserr << " EndTime............................" << EndTime << "\n";
	opserr << " Interval..........................." << TimeInterval << "\n";
	opserr << "\n";
	opserr << " FragMin............................" << FragMin << "\n";
	opserr << " FragInt............................" << FragInt << "\n";
	opserr << " nFrag.............................." << nFrag << "\n";
	opserr << "\n";
	opserr << " instantaneous......................" << instantaneous<< "\n";
	opserr << " firstpassage......................." << firstpassage << "\n";
	opserr << "\n";
	opserr << " stationary........................." << stationary << "\n";
	opserr << " stationaryTime....................." << stationaryTime<< "\n";
	opserr << "\n";
	opserr << " maxSim............................." << maxSim<< "\n";
	opserr << " siminterval........................" << intervalSim<< "\n";
	opserr << " eps................................" << eps<< "\n";
//  check for analysis

	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}

	if(stationary==1){
	 	theRandomVibrationSimulation=new StatRandomVibrationSimulation
								(theReliabilityDomain,
								 theStructuralDomain,
						         theFunctionEvaluator,
							     theProbabilityTransformation,
						         StartTime,EndTime,TimeInterval,
						         FragMin,FragInt,nFrag,
                                 stationaryTime,
								 twoside,
								 system,
						         maxSim,intervalSim,eps,
						         instantaneous,
						         firstpassage,
	  				             argv[1],
						         fileBinary,
	                             interp,
								 print,
								 sampleAmp,
								 sampleTime);
//							 theStartPoint);
	}else{
		theRandomVibrationSimulation=new NonStatRandomVibrationSimulation
								(theReliabilityDomain,
								 theStructuralDomain,
						         theFunctionEvaluator,
							     theProbabilityTransformation,
						         StartTime,EndTime,TimeInterval,
						         FragMin,FragInt,nFrag,
								 twoside,
								 system,
						         maxSim,intervalSim,eps,
						         instantaneous,
						         firstpassage,
	  				             argv[1],
						         fileBinary,
	                             interp,
								 print);
	}

	if (theRandomVibrationSimulation == 0) {
		opserr << "ERROR: could not create theRandomVibrationSimulation\n";
		return TCL_ERROR;
	}
	theRandomVibrationSimulation->analyze();
	return TCL_OK;
}


int 
TclReliabilityModelBuilder_runRandomVibrationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// In case this is a replacement
	if (theRandomVibrationAnalysis != 0) {
		delete theRandomVibrationAnalysis;
		theRandomVibrationAnalysis = 0;
	}
	if (theReliabilityDomain == 0 ) {
		opserr << "Need ReliabilityDomain before an RandomVibrationAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theStructuralDomain == 0 ) {
		opserr << "Need StructuralDomain before an RandomVibrationAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theFunctionEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a Outcrossing can be created" << endln;
		///////////// Modified by K Fujimura /////////////////////////////////////////////
		opserr << "Assume all RV's are independent" << endln;
		theProbabilityTransformation = 
		new AllIndependentTransformation(theReliabilityDomain,0);
		///////////// Modified by K Fujimura /////////////////////////////////////////////
	}
	double StartTime=1.0;
	double EndTime=20.0;
	double TimeInterval=1.0;
//	double StartAnalysis=-999.9;
	double StartAnalysis=12.0;
	double FragMin=1.0;
	double FragInt=0.0;
	int nFrag=1;
	int designPoint=1;
	char *fileBinary=new char[100];
	
//  indicator for design-point analysis
//	=0 : No design-point Analysis
//	=1 : design-point Analysis and outcrossing rate analysis 
//	=2 : design-point Analysis and outcrossing rate analysis
//							   and FO system first passage		
	bool stationary = true;  // true - stationary problem
	bool print=false;
	bool mirrorimage = false;
	bool initialpoint = true;
	bool firstpassage = true;
//  indicator for the first-passage analysis
//	=0 : no first passage analysis;
//  =1 : first passage analysis with first passage analyzer 

	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if (strcmp(argv[argvCounter],"-starttime") == 0) {
			argvCounter++;
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetDouble(interp, argv[argvCounter], &StartTime) != TCL_OK) {
				opserr << "ERROR: invalid input stepsToStart to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-endtime") == 0) {
			argvCounter++;
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetDouble(interp, argv[argvCounter], &EndTime) != TCL_OK) {
				opserr << "ERROR: invalid input stepsToEnd to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-interval") == 0) {
			argvCounter++;
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetDouble(interp, argv[argvCounter], &TimeInterval) != TCL_OK) {
				opserr << "ERROR: invalid input sampleFreq to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-startanalysis") == 0) {
			argvCounter++;
			// GET INPUT PARAMETER (integer)
			if (Tcl_GetDouble(interp, argv[argvCounter], &StartAnalysis) != TCL_OK) {
				opserr << "ERROR: invalid input stepsToEnd to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-stationary") == 0) {
			argvCounter++;
			stationary=true;
		}
		else if (strcmp(argv[argvCounter],"-print") == 0) {
			argvCounter++;
			print=true;
		}
		else if (strcmp(argv[argvCounter],"-fragility") == 0) {
			argvCounter++;
			  // GET INPUT PARAMETER (double)
			if (Tcl_GetDouble(interp, argv[argvCounter], &FragMin) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &FragInt) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
			if (Tcl_GetInt(interp, argv[argvCounter], &nFrag) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-designpoint") == 0) {
			argvCounter++;
			  // GET INPUT PARAMETER (double)
			if (Tcl_GetInt(interp, argv[argvCounter], &designPoint) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-mirrorimage") == 0) {
			opserr << "Invalid Input for randomvibrationanalysis\n";
			opserr << "-mirrorimage can not be selected \n";
			return TCL_ERROR;
//			argvCounter++;
//			  // GET INPUT PARAMETER (double)
//			int ind;
//			if (Tcl_GetInt(interp, argv[argvCounter], &ind) != TCL_OK) {
//				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
//				return TCL_ERROR;
//			}
//			if(ind==0) mirrorimage=false;
//			else mirrorimage=true;
//			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-binary") == 0) {
			argvCounter++;
			strcpy(fileBinary,argv[argvCounter]);
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-initialpoint") == 0) {
			argvCounter++;
			  // GET INPUT PARAMETER (double)
			int ind;
			if (Tcl_GetInt(interp, argv[argvCounter], &ind) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			if(ind==0) initialpoint=false;
			else initialpoint=true;
			argvCounter++;
		}
		else if (strcmp(argv[argvCounter],"-firstpassage") == 0) {
			argvCounter++;
			  // GET INPUT PARAMETER (double)
			int ind;
			if (Tcl_GetInt(interp, argv[argvCounter], &ind) != TCL_OK) {
				opserr << "ERROR: invalid input littleDt to theOutCrossingAnalysis \n";
				return TCL_ERROR;
			}
			if(ind==0) firstpassage=false;
			else firstpassage=true;
			argvCounter++;
		}
		else {
			opserr << "ERROR: Invalid input to theOutCrossingAnalysis." << endln;
			return TCL_ERROR;
		}
	}
	if(designPoint!=0){
		if (theFindDesignPointAlgorithm == 0 ) {
		opserr << "Need theFindDesignPointAlgorithm before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
		}
		if (theGradientEvaluator == 0 ) {
		opserr << "Need theGradientEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
		}
	}

	opserr << "=========================================\n";
	opserr << "\n";
	opserr << "     OutCrossingAnalysis  \n";
	opserr << "\n";
	opserr << "=========================================\n";
	opserr << "\n";
	opserr << " StartTime.........................." << StartTime << "\n";
	opserr << " EndTime............................" << EndTime << "\n";
	opserr << " Interval..........................." << TimeInterval << "\n";
	opserr << "\n";
	opserr << " StartAnalysis......................" << StartAnalysis << "\n";
	opserr << "\n";
	opserr << " FragMin............................" << FragMin << "\n";
	opserr << " FragInt............................" << FragInt << "\n";
	opserr << " nFrag.............................." << nFrag << "\n";
	opserr << "\n";
	opserr << " designpoint........................" << designPoint << "\n";
	opserr << "\n";
	opserr << " stationary........................." << stationary << "\n";
	opserr << " mirrorimage........................" << mirrorimage << "\n";
	opserr << " initialpoint......................." << initialpoint << "\n";
	opserr << "\n";
//  check for analysis

	if(mirrorimage){
//		if (theMirrorImageBuilder== 0 ) {
//		opserr << "Need MirrorImageBuilder before a RANDOMVIBRATION can be created" << endln;
//		return TCL_ERROR;
//		}
	}
	if(initialpoint){
		if (theInitialPointBuilder== 0 ) {
		opserr << "InitialPointBuilder is not specified before randomvibration \n"; 
		opserr << "default object is initiated\n"; 
//		theInitialPointBuilder = new MirrorImageInitialPointBuilder
//									 (theReliabilityDomain,
//									  theFunctionEvaluator);
		}
	}
	if(abs(designPoint)>1){
		if (theCrossingRateAnalyzer== 0 ) {
		opserr << "Need CrossingRateAnalyzer before a RANDOMVIBRATION can be created" << endln;
		return TCL_ERROR;
		}
	}
	if(firstpassage){
		if (theFirstPassageAnalyzer== 0 ) {
		opserr << "Need FirstPassageAnalyzer before a RANDOMVIBRATION can be created" << endln;
		return TCL_ERROR;
		}
	}

	theRandomVibrationAnalysis = new RandomVibrationAnalysis
								(theReliabilityDomain,
								 theFindDesignPointAlgorithm,
								 theStructuralDomain,
//								 theMirrorImageBuilder,
								 theInitialPointBuilder,
								 theCrossingRateAnalyzer,
								 theFirstPassageAnalyzer,
								 theFunctionEvaluator,
								 theGradientEvaluator,
								 theReliabilityConvergenceCheck,
								 StartTime,EndTime,TimeInterval,StartAnalysis,	
								 FragMin,FragInt,nFrag,
								 designPoint,
								 stationary,
								 mirrorimage,
								 initialpoint,
								 firstpassage,
								 argv[1],
								 fileBinary,
								 interp,
								 print);


	if (theRandomVibrationAnalysis == 0) {
		opserr << "ERROR: could not create theRandomVibrationAnalysis\n";
		return TCL_ERROR;
	}
	theRandomVibrationAnalysis->analyze();

	return TCL_OK;
}
