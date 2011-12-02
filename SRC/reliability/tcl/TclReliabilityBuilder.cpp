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
                                                                        
// $Revision: 1.54 $
// $Date: 2010-06-10 20:16:16 $
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
#include <LimitStateFunction.h>
#include <LimitStateFunctionIter.h>
#include <RandomVariablePositioner.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <ParameterPositioner.h>
#include <ParameterPositionerIter.h>

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
#include <FOSMAnalysis.h>
#include <ParametricReliabilityAnalysis.h>
#include <GFunVisualizationAnalysis.h>
#include <OutCrossingAnalysis.h>
#include <ImportanceSamplingAnalysis.h>
#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>
#include <FindCurvatures.h>
#include <FirstPrincipalCurvature.h>
#include <CurvaturesBySearchAlgorithm.h>
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
#include <SensitivityAlgorithm.h>
#include <ReliabilityConvergenceCheck.h>
#include <StandardReliabilityConvergenceCheck.h>
#include <OptimalityConditionReliabilityConvergenceCheck.h>
#include <MeritFunctionCheck.h>
#include <AdkZhangMeritFunctionCheck.h>
#include <PolakHeSearchDirectionAndMeritFunction.h>
#include <SQPsearchDirectionMeritFunctionAndHessian.h>
#include <HessianApproximation.h>
#include <GradientProjectionSearchDirection.h>
#include <RootFinding.h>
#include <SecantRootFinding.h>

//Quan---
#ifdef _SNOPT
#include <SnoptProblem.h>  
#include <SnoptAnalysis.h>   
#endif

#include <DesignVariable.h>
#include <DesignVariablePositioner.h>
#include <ConstraintFunction.h>
#include <ObjectiveFunction.h>
#include <MonteCarloResponseAnalysis.h>
#include <OrthogonalPlaneSamplingAnalysis.h>

#include <Hessian.h>
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
/////////////////////////////////////////////////////////
/////E Modified by K Fujimura /////////////////////////////
/////////////////////////////////////////////////////////

extern SensitivityAlgorithm *theSensitivityAlgorithm;
/////////////////////////////////////////////////////////
/////S Modified by K Fujimura /////////////////////////////
/////////////////////////////////////////////////////////
extern ReliabilityStaticAnalysis* theReliabilityStaticAnalysis;
extern ReliabilityDirectIntegrationAnalysis* theReliabilityTransientAnalysis;
extern SensitivityIntegrator* theSensitivityIntegrator;
/////////////////////////////////////////////////////////
/////E Modified by K Fujimura /////////////////////////////
/////////////////////////////////////////////////////////

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

// Quan --
// ---------------- define global pointer for SNOPT ----------
#ifdef _SNOPT
SnoptProblem * theSNOPT=0;
SNOPTAnalysis * theSNOPTAnalysis=0;
#endif
MonteCarloResponseAnalysis * theMonteCarloResponseAnalysis=0;
SamplingAnalysis * theSamplingAnalysis =0;
// --- Quan

ReliabilityDomain *theReliabilityDomain = 0;
static Domain *theStructuralDomain = 0;

// base class static pointers
static GFunEvaluator *theGFunEvaluator = 0;
static GradGEvaluator *theGradGEvaluator = 0;
static StepSizeRule *theStepSizeRule = 0;
static SearchDirection *theSearchDirection = 0;
static HessianApproximation *theHessianApproximation = 0;
static MeritFunctionCheck *theMeritFunctionCheck = 0;
static ProbabilityTransformation *theProbabilityTransformation = 0;
static ReliabilityConvergenceCheck *theReliabilityConvergenceCheck = 0;
static Vector *theStartPoint = 0;
static bool startAtOrigin = false;
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
static ParametricReliabilityAnalysis *theParametricReliabilityAnalysis = 0;
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
int TclReliabilityModelBuilder_runFOSMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runParametricReliabilityAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runGFunVisualizationAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runOutCrossingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSORMAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSystemAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runImportanceSamplingAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_printReliability(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_inputCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getMean(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getStdv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_rvReduction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_getBetaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);//not in K.F.
int TclReliabilityModelBuilder_getGammaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);//not in K.F.
int TclReliabilityModelBuilder_getAlphaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);//not in K.F.
int TclReliabilityModelBuilder_invNormalCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);//not in K.F.
int TclReliabilityModelBuilder_getRVTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);//not in K.F.
int TclReliabilityModelBuilder_getLSFTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);//not in K.F.
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
int TclReliabilityModelBuilder_addDesignVariable(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv);
int TclReliabilityModelBuilder_addDesignVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_addObjectiveFunction(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv);
int TclReliabilityModelBuilder_addConstraintFunction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclReliabilityModelBuilder_runSNOPTAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
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
  Tcl_CreateCommand(interp, "hessianApproximation",	TclReliabilityModelBuilder_addHessianApproximation,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "meritFunctionCheck",	TclReliabilityModelBuilder_addMeritFunctionCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "reliabilityConvergenceCheck",	TclReliabilityModelBuilder_addReliabilityConvergenceCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "probabilityTransformation",	TclReliabilityModelBuilder_addProbabilityTransformation,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "findCurvatures",	TclReliabilityModelBuilder_addFindCurvatures,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "randomNumberGenerator",TclReliabilityModelBuilder_addRandomNumberGenerator,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runFORMAnalysis",TclReliabilityModelBuilder_runFORMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runFOSMAnalysis",TclReliabilityModelBuilder_runFOSMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runParametricReliabilityAnalysis",TclReliabilityModelBuilder_runParametricReliabilityAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runGFunVizAnalysis",TclReliabilityModelBuilder_runGFunVisualizationAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runOutCrossingAnalysis",TclReliabilityModelBuilder_runOutCrossingAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runSORMAnalysis",TclReliabilityModelBuilder_runSORMAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runSystemAnalysis",TclReliabilityModelBuilder_runSystemAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "runImportanceSamplingAnalysis",TclReliabilityModelBuilder_runImportanceSamplingAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "printReliability",TclReliabilityModelBuilder_printReliability,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "inputCheck",TclReliabilityModelBuilder_inputCheck,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getMean",TclReliabilityModelBuilder_getMean,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getStdv",TclReliabilityModelBuilder_getStdv,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "rvReduction",TclReliabilityModelBuilder_rvReduction,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "betaFORM",TclReliabilityModelBuilder_getBetaFORM,(ClientData)NULL, NULL); //not in K.F.
  Tcl_CreateCommand(interp, "gammaFORM",TclReliabilityModelBuilder_getGammaFORM,(ClientData)NULL, NULL);//not in K.F.
  Tcl_CreateCommand(interp, "alphaFORM",TclReliabilityModelBuilder_getAlphaFORM,(ClientData)NULL, NULL);//not in K.F.
  Tcl_CreateCommand(interp, "invNormalCDF",TclReliabilityModelBuilder_invNormalCDF,(ClientData)NULL, NULL);//not in K.F.
  Tcl_CreateCommand(interp, "getRVTags",TclReliabilityModelBuilder_getRVTags,(ClientData)NULL, NULL);//not in K.F.
  Tcl_CreateCommand(interp, "getLSFTags",TclReliabilityModelBuilder_getLSFTags,(ClientData)NULL, NULL);//not in K.F.
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
   Tcl_CreateCommand(interp, "designVariablePositioner",TclReliabilityModelBuilder_addDesignVariablePositioner,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "constraintFunction", TclReliabilityModelBuilder_addConstraintFunction,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "objectiveFunction", TclReliabilityModelBuilder_addObjectiveFunction,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "runSNOPTAnalysis", TclReliabilityModelBuilder_runSNOPTAnalysis,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "updateParameterValue", TclReliabilityModelBuilder_updateParameterValue,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "designVariable", TclReliabilityModelBuilder_addDesignVariable,(ClientData)NULL, NULL);
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
  theReliabilityDomain	= new ReliabilityDomain();


}

TclReliabilityBuilder::~TclReliabilityBuilder()
{

  // Delete objects
  if (theReliabilityDomain != 0)
    delete theReliabilityDomain;

  // base class pointers
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
  if (theParametricReliabilityAnalysis != 0)
    delete theParametricReliabilityAnalysis;
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
 #ifdef _SNOPT
  if (theSNOPTAnalysis != 0)
    delete theSNOPTAnalysis;
  if (theSNOPT != 0)
    delete theSNOPT;
#endif
  if (theMonteCarloResponseAnalysis != 0)
    delete theMonteCarloResponseAnalysis;
  if (theSamplingAnalysis != 0)
    delete theSamplingAnalysis;
  // ---Quan 

  theReliabilityDomain = 0;
  theGFunEvaluator = 0;
  theGradGEvaluator = 0;
  theStepSizeRule = 0;
  theSearchDirection = 0;
  theHessianApproximation = 0;
  theMeritFunctionCheck = 0;
  theReliabilityConvergenceCheck = 0;
  theProbabilityTransformation = 0;
  theStartPoint = 0;
  theRootFindingAlgorithm = 0;
  theRandomNumberGenerator = 0;
  theFindDesignPointAlgorithm = 0;
  theFindCurvatures = 0;
  
  thePolakHeDualPurpose =0;
  theSQPtriplePurpose =0;
  
  theFORMAnalysis = 0;
  theFOSMAnalysis = 0;
  theParametricReliabilityAnalysis = 0;
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
 #ifdef _SNOPT
  theSNOPTAnalysis = 0;
  theSNOPT = 0;
#endif
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
  Tcl_DeleteCommand(theInterp, "randomVariablePositioner");
  Tcl_DeleteCommand(theInterp, "positionerPositioner");
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
  Tcl_DeleteCommand(theInterp, "hessianApproximation");
  Tcl_DeleteCommand(theInterp, "meritFunctionCheck");
  Tcl_DeleteCommand(theInterp, "reliabilityConvergenceCheck");
  Tcl_DeleteCommand(theInterp, "probabilityTransformation");
  Tcl_DeleteCommand(theInterp, "findCurvatures");
  Tcl_DeleteCommand(theInterp, "randomNumberGenerator");
  Tcl_DeleteCommand(theInterp, "runFORMAnalysis");
  Tcl_DeleteCommand(theInterp, "runFOSMAnalysis");
  Tcl_DeleteCommand(theInterp, "runParametricReliabilityAnalysis");
  Tcl_DeleteCommand(theInterp, "runGFunVizAnalysis");
  Tcl_DeleteCommand(theInterp, "runOutCrossingAnalysis");
  Tcl_DeleteCommand(theInterp, "runSORMAnalysis");
  Tcl_DeleteCommand(theInterp, "runSystemAnalysis");
  Tcl_DeleteCommand(theInterp, "runImportanceSamplingAnalysis");
  Tcl_DeleteCommand(theInterp, "printReliability");
  Tcl_DeleteCommand(theInterp, "inputCheck");
  Tcl_DeleteCommand(theInterp, "getMean");
  Tcl_DeleteCommand(theInterp, "getStdv");
  Tcl_DeleteCommand(theInterp, "rvReduction");
  Tcl_DeleteCommand(theInterp, "betaFORM");
  Tcl_DeleteCommand(theInterp, "gammaFORM");
  Tcl_DeleteCommand(theInterp, "invNormalCDF");
  Tcl_DeleteCommand(theInterp, "getRVTags");
  Tcl_DeleteCommand(theInterp, "getLSFTags");

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
  Tcl_DeleteCommand(theInterp, "designVariablePositioner");
  Tcl_DeleteCommand(theInterp, "constraintFunction");
  Tcl_DeleteCommand(theInterp, "objectiveFunction");
  Tcl_DeleteCommand(theInterp, "runSNOPTAnalysis");
  Tcl_DeleteCommand(theInterp, "updateParameterValue");
  Tcl_DeleteCommand(theInterp, "designVariable");
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
	  (strcmp(argv[2],"userdefined")		== 0) ||
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


  // TREAT THE USER-DEFINED RANDOM VARIABLE AS A SPECIAL CASE
  if (strcmp(argv[2],"userdefined") == 0) {

		Vector xPoints;
		Vector PDFpoints;
		int numPoints = 0;
		int i;

		if (strcmp(argv[3],"-list") == 0) {

			numPoints = (argc-4) % 2;
			Vector temp_xPoints(numPoints);
			Vector temp_PDFpoints(numPoints);

			double x = 0.0;
			double pdf = 0.0;
			double x_old = 0.0;

			// Read the points
			for (i=0; i<numPoints; i++) {
			  if (Tcl_GetDouble(interp, argv[4+2*i], &x) != TCL_OK) {
				  opserr << "ERROR: Invalid x point to user-defined random variable." << endln;
				  return TCL_ERROR;
			  }
			  if (Tcl_GetDouble(interp, argv[5+2*i], &pdf) != TCL_OK) {
				  opserr << "ERROR: Invalid PDF value point to user-defined random variable." << endln;
				  return TCL_ERROR;
			  }
			  if (i>0 && x<=x_old) {
				  opserr << "ERROR: x-points to user-defined random variable must be consequtive!" << endln;
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

			// Close the file
			inputFile.close();

			// Allocate vectors of correct size
			Vector temp_xPoints(numPoints);
			Vector temp_PDFpoints(numPoints);

			// Open it again, now being ready to store the results in a matrix
			ifstream inputFile2( argv[4], ios::in );
			if (inputFile2.fail()) {
				opserr << "File " << argv[4] << " could not be opened. " << endln;
				return TCL_ERROR;
			}

			// Store the vector
			for (int i=0; i<numPoints; i++) {
					inputFile2 >> temp_xPoints(i);
					inputFile2 >> temp_PDFpoints(i);
			}
			inputFile2.close();

			xPoints = temp_xPoints;
			PDFpoints = temp_PDFpoints;
		}
		else {
			opserr << "ERROR: Invalid argument to user-defined random variable, number " << tag << endln;
			return TCL_ERROR;
		}
	  
	
	  // Normalize the PDF
	  double sum = 0.0;
	  for ( i=1; i<numPoints; i++ ) {
		  sum += 0.5 * (PDFpoints(i)+PDFpoints(i-1)) * (xPoints(i)-xPoints(i-1));
	  }
	  double diff;
	  if ( fabs(sum-1.0) < 1.0e-6) {
		  // It's normalized enough...!
	  }
	  else if (sum < 1.0) {
		  diff = (1.0-sum)/(xPoints(numPoints-1)-xPoints(0));
		  opserr << "WARNING: The PDF of random variable " << tag << " is normalized by " << endln
			     << "         uniformly increasing the PDF values by " << diff << endln;
		  for (int i=0; i<numPoints; i++) {
			  PDFpoints(i) = PDFpoints(i) + diff;
		  }
	  }
	  else {
		  diff = (sum-1.0)/(xPoints(numPoints-1)-xPoints(0));
		  opserr << "WARNING: The PDF of random variable " << tag << " is normalized by " << endln
			     << "         uniformly decreasing the PDF values by " << diff << endln;
		  for (int i=0; i<numPoints; i++) {
			  PDFpoints(i) = PDFpoints(i) - diff;
		  }
	  }

	  // Check if the PDF became negative somewhere
	  for (i=0; i<numPoints; i++) {
		  if ( PDFpoints(i) < 0.0 ) {
			  opserr << "ERROR: Some PDF values of random variable " << tag << endln
				  << "became negative after normalization. " << endln;
			  return TCL_ERROR;
		  }
	  }

	  // Check that it has been normalized
	  sum = 0.0;
	  for ( i=1; i<numPoints; i++ ) {
		  sum += 0.5 * (PDFpoints(i)+PDFpoints(i-1)) * (xPoints(i)-xPoints(i-1));
	  }
      if ( fabs(1.0-sum) > 1.0e-6 ) {
		  opserr << "ERROR: Did not succeed in normalizing user-defined distribution." << endln;
		  return TCL_ERROR;
	  }


//	  // Check if the user has given a start point too, and create the random variables
//	  if ( floor((argc-4)/2.0) != (argc-4)/2.0 ) {
//	  
//		  double startPt = 0.0;
//		  if (Tcl_GetDouble(interp, argv[5], &startPt) != TCL_OK) {
//			  opserr << "ERROR: Invalid start point to user-defined random variable." << endln;
//			  return TCL_ERROR;
//		  }
//
//	  	  theRandomVariable = new UserDefinedRV(tag, xPoints, PDFpoints, startPt);
//	  }
//	  else {

		  theRandomVariable = new UserDefinedRV(tag, xPoints, PDFpoints);
//	  }

	  // Add the random variable to the domain
		  if (theReliabilityDomain->addRandomVariable(theRandomVariable, startPt) == false) {
		opserr << "ERROR: failed to add random variable to the domain (wrong number of arguments?)\n";
		opserr << "random variable: " << tag << endln;
		delete theRandomVariable; // otherwise memory leak
		return TCL_ERROR;
	  }
	  else {
		  return TCL_OK;
	  }

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
  if (theReliabilityDomain->addRandomVariable(theRandomVariable, startPt) == false) {
	opserr << "ERROR: failed to add random variable to the domain (wrong number of arguments?)\n";
	opserr << "random variable: " << tag << endln;
	delete theRandomVariable; // otherwise memory leak
	return TCL_ERROR;
  }

  // this has all moved to the gFunEvaluator classes, also note you cannot use theRandomVariable once it has been deleted
  //char tclAssignment[80];
  //sprintf(tclAssignment , "set xrv(%d)  %15.5f", tag, theRandomVariable->getStartValue());
  //if (Tcl_Eval(interp, tclAssignment) == TCL_ERROR) {						
  //  opserr << "ERROR GFunEvaluator -- Tcl_Eval returned error in limit state function" << endln;
  //  opserr << interp->result << endln;
  //  return TCL_ERROR;
  //}

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
  theLimitStateFunction = new LimitStateFunction(tag, argv[2], interp);
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
TclReliabilityModelBuilder_addGradLimitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  LimitStateFunction *theLimitStateFunction = 0;
  int lsfTag, rvTag;

  if (theGFunEvaluator != 0 ) {
    opserr << "ERROR: A limit-state function should not be created after the GFunEvaluator has been instantiated." << endln;
    return TCL_ERROR;
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

		if (strcmp(argv[argvCounter],"userdefined") == 0) {
			opserr << "ERROR: Can't create a user-defined random variable like this." << endln;
			return TCL_ERROR;
		}

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

		if (strcmp(argv[argvCounter],"userdefined") == 0) {
			opserr << "ERROR: Can't create a user-defined random variable like this." << endln;
			return TCL_ERROR;
		}

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

		if (strcmp(argv[argvCounter],"userdefined") == 0) {
			opserr << "ERROR: Can't create a user-defined random variable like this." << endln;
			return TCL_ERROR;
		}

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

		if (strcmp(argv[argvCounter],"userdefined") == 0) {
			opserr << "ERROR: Can't create a user-defined random variable like this." << endln;
			return TCL_ERROR;
		}

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
	
	RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtr(rvNumber);
	if (theRV == 0){
	  opserr << "ERROR:: RandomVariable " << rvNumber << " not found in the reliability domain " << endln;
	  return TCL_ERROR;
	}
	//int rvIndex = theRV->getIndex();
	int rvIndex = theReliabilityDomain->getRandomVariableIndex(rvNumber);

	if (strcmp(argv[argvCounter],"-parameter") == 0) {
	  argvCounter++;
	  int paramTag;
	  if (Tcl_GetInt(interp, argv[argvCounter++], &paramTag) != TCL_OK) {
	    opserr << "ERROR: invalid input in positioner: parameter tag \n";
	    return TCL_ERROR;
	  }

	  Parameter *theParameter = theStructuralDomain->getParameter(paramTag);

	  if (theParameter == 0) {
	    opserr << "ERROR: parameter with tag " << paramTag 
		   << " not found in structural domain\n";
	    return TCL_ERROR;
	  }
	  else {
	    theRandomVariablePositioner =
	      new RandomVariablePositioner(tag, rvIndex, *theParameter, theRV);
	  }
	}

	// IF UNCERTAIN *ELEMENT* PROPERTY
	else if (strcmp(argv[argvCounter],"-element") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			argvCounter++;
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}

		theObject = (DomainComponent *)theStructuralDomain->getElement(tagOfObject);

		theRandomVariablePositioner =
		  new RandomVariablePositioner(tag,
					       rvIndex,
					       theObject,
					       &argv[argvCounter],
					       argc-argvCounter,
						   theRV);

		//int rvnumber = theRandomVariablePositioner->getRvNumber();
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
//			theRandomVariablePositioner = new RandomVariablePositioner(tag,rvIndex,theObject,&argv[5],argc-5,factor);
//		}
//		else {
//			theRandomVariablePositioner = new RandomVariablePositioner(tag,rvIndex,theObject,&argv[5],argc-5);
//		}
		theRandomVariablePositioner =
		  new RandomVariablePositioner(tag,
					       rvIndex,
					       theObject,
					       &argv[argvCounter],
					       argc-argvCounter,
						   theRV);
	}

	// IF UNCERTAIN *NODE* PROPERTY
	else if (strcmp(argv[argvCounter],"-node") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getNode(tagOfObject);

		theRandomVariablePositioner =
		  new RandomVariablePositioner(tag,
					       rvIndex,
					       theObject,
					       &argv[argvCounter],
					       argc-argvCounter,
						   theRV);
	}
	else {
		opserr << "ERROR: Unknown parameter in randomVariablePositioner" << endln;
		return TCL_ERROR;
	}

	// ADD THE RANDOMVARIABLEPOSITIONER TO THE DOMAIN
	if (theReliabilityDomain->addRandomVariablePositioner(theRandomVariablePositioner) == false) {
		opserr << "ERROR: failed to add random variable positioner number " << tag << " to the domain." << endln;
		delete theRandomVariablePositioner; // otherwise memory leak
		return TCL_ERROR;
	}

	return TCL_OK;


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

	if (strcmp(argv[argvCounter],"-parameter") == 0) {
	  argvCounter++;
	  int paramTag;
	  if (Tcl_GetInt(interp, argv[argvCounter++], &paramTag) != TCL_OK) {
	    opserr << "ERROR: invalid input in positioner: parameter tag \n";
	    return TCL_ERROR;
	  }

	  Parameter *theParameter = theStructuralDomain->getParameter(paramTag);

	  if (theParameter == 0) {
	    opserr << "ERROR: parameter with tag " << paramTag 
		   << " not found in structural domain\n";
	    return TCL_ERROR;
	  }
	  else {
	    theParameterPositioner =
	      new ParameterPositioner(tag, *theParameter);
	  }
	}

	// IF UNCERTAIN *LOAD*
	else if (strcmp(argv[argvCounter],"-loadPattern") == 0) {
		argvCounter++;
		
		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
		  opserr << "ERROR: invalid input: tagOfObject \n";
		  return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getLoadPattern(tagOfObject);

		theParameterPositioner =
		  new ParameterPositioner(tag,
					  theObject,
					  &argv[argvCounter],
					  argc-argvCounter);
	}

	// IF UNCERTAIN element property
	else if (strcmp(argv[argvCounter],"-element") == 0) {
		argvCounter++;
		
		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
		  opserr << "ERROR: invalid input: tagOfObject \n";
		  return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getElement(tagOfObject);

		theParameterPositioner =
		  new ParameterPositioner(tag,
					  theObject,
					  &argv[argvCounter],
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

		theParameterPositioner =
		  new ParameterPositioner(tag,
					  theObject,
					  &argv[argvCounter],
					  argc-argvCounter);
	}	else {
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
		if (theGFunEvaluator == 0 ) {
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

		theStepSizeRule = new ArmijoStepSizeRule(theGFunEvaluator,
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

		if (argc != 4) {
			opserr << "ERROR: Wrong number of arguments to Tcl g-function evaluator." << endln;
			return TCL_ERROR;
		}

		if (strcmp(argv[2],"-file") != 0) {
			opserr << "ERROR: Wrong input to Tcl g-function evaluator." << endln;
			return TCL_ERROR;
		}
		theGFunEvaluator = new TclGFunEvaluator(interp, theReliabilityDomain, 
											theStructuralDomain, argv[3]);

	}
	else if (strcmp(argv[1],"OpenSees") == 0) {

		// There are several alternatives for this command:
		// gFunEvaluator  OpenSees  -file <filename> <dT (optional for outcrossing only)> 
		// gFunEvaluator  OpenSees  -runToMaxTimeInGFun
		// gFunEvaluator  OpenSees  -analyze <numSteps> <dt(optional)>

		if (argc < 3) {
			opserr << "ERROR: Too few arguments to gFunEvaluator" << endln;
			return TCL_ERROR;
		}

		if (strcmp(argv[2],"-file") == 0) {

			// Try to open the file to make sure it exists
			double dt = 0.0;
			ifstream inputFile( argv[3], ios::in );
			if (inputFile.fail()) {
				opserr << "File " << argv[3] << " could not be opened. " << endln;
				return TCL_ERROR;
			}
			inputFile.close();

			if (argc ==5) {
				if (Tcl_GetDouble(interp, argv[4], &dt) != TCL_OK) {
					opserr << "ERROR: invalid input: dt for OpenSees GFunEvaluator \n";
					return TCL_ERROR;
				}
			}
			//				interp, theReliabilityDomain, argv[3], dt);
			theGFunEvaluator = new OpenSeesGFunEvaluator(
								     interp, theReliabilityDomain,
								     theStructuralDomain, argv[3]);
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
			theGFunEvaluator = new OpenSeesGFunEvaluator(
				interp, theReliabilityDomain,
				theStructuralDomain, nsteps, dt);

		}
		else {
			opserr << "ERROR: unrecognized parameter in OpenSees GFunEvaluator \n";
			return TCL_ERROR;
		}
	
	}
	else if (strcmp(argv[1],"Basic") == 0) {
		theGFunEvaluator = new BasicGFunEvaluator(interp, theReliabilityDomain);
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
		theGFunEvaluator = new AnalyzerGFunEvaluator
					       (interp, theReliabilityDomain,
						   theStructuralDomain,theAnalyzer);
	}
/////////////////////////////////////////
////////E modified by K Fujimura 10/10/2004
/////////////////////////////////////////
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
					opserr << "ERROR: Error in input to gradGEvaluator. " << endln;
					return TCL_ERROR;
				}
			}
		}

		theGradGEvaluator = new FiniteDifferenceGradGEvaluator(theGFunEvaluator, theReliabilityDomain, 
					interp, perturbationFactor,doGradientCheck, false);
	}
	else if (strcmp(argv[1],"OpenSees") == 0) {

		bool doGradientCheck = false;

	//Quan Apr. 2006	
		if (theSensitivityAlgorithm == 0) {
			opserr << "Warning:Need a DDM sensitivity algorithm before a OpenSees sensitivity evaluator can be created" << endln;
	//		return TCL_ERROR;
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

		theGradGEvaluator = new OpenSeesGradGEvaluator(interp, theGFunEvaluator, 
					theReliabilityDomain, theStructuralDomain, theSensitivityAlgorithm, doGradientCheck);
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

		theGradGEvaluator = new AnalyzerGradGEvaluator(interp, theReliabilityDomain,
							 theStructuralDomain, theGFunEvaluator,
							 doGradientCheck);
	}
	////////////////////////////////////////
	//////E modified by K Fujimura 10/10/2004
	////////////////////////////////////////


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


	// GET INPUT PARAMETER (string) AND CREATE THE OBJECT   // modified by Quan Gu July 2006
	if (strcmp(argv[1],"firstPrincipal") == 0) {

		theFindCurvatures = new FirstPrincipalCurvature();
		if (argc>=3){
			if (strcmp(argv[2],"-exe") == 0)
				theFindCurvatures->computeCurvatures(theReliabilityDomain);
		}
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

	int printFlag=0;
	char fileNamePrint[256];
	strcpy(fileNamePrint,"initialized");
	int maxNumIter = 100;
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
//		if (theStartPoint == 0 ) {
//			opserr << "Need theStartPoint before a FindDesignPointAlgorithm can be created" << endln;
//			return TCL_ERROR;
//		}
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
					maxNumIter, theReliabilityDomain, 
					theGFunEvaluator,
					theGradGEvaluator,
					theStepSizeRule,
					theSearchDirection,
					theProbabilityTransformation,
					theHessianApproximation,
					theReliabilityConvergenceCheck,
					startAtOrigin,
					printFlag, fileNamePrint);
		
	}   //if StepSearch

	// Quan SNOPT interface ----
#ifdef _SNOPT
	else if (strcmp(argv[argvCounter],"SNOPT") == 0) {
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
/*		if (theStepSizeRule == 0 ) {
			opserr << "Need theStepSizeRule before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		} 
		if (theSearchDirection == 0 ) {
			opserr << "Need theSearchDirection before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}  */
		if (theProbabilityTransformation == 0 ) {
			opserr << "Need theProbabilityTransformation before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
//		if (theStartPoint == 0 ) {
//			opserr << "Need theStartPoint before a FindDesignPointAlgorithm can be created" << endln;
//			return TCL_ERROR;
//		}
/*		if (theReliabilityConvergenceCheck == 0 ) {
			opserr << "Need theReliabilityConvergenceCheck before a FindDesignPointAlgorithm can be created" << endln;
			return TCL_ERROR;
		}
*/
		maxNumIter = 250;
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
				strcpy(fileNamePrint1,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printAllPointsY") == 0) {
				argvCounter++;
				printFlag = 2;
				strcpy(fileNamePrint1,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printDesignPointX") == 0) {
				argvCounter++;
				printFlag = 3;
				strcpy(fileNamePrint1,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printDesignPointY") == 0) {
				argvCounter++;
				printFlag = 4;
				strcpy(fileNamePrint1,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printCurrentPointX") == 0) {
				argvCounter++;
				printFlag = 5;
				strcpy(fileNamePrint1,argv[argvCounter]);
				argvCounter++;
			}
			else if (strcmp(argv[argvCounter],"-printCurrentPointY") == 0) {
				argvCounter++;
				printFlag = 6;
				strcpy(fileNamePrint1,argv[argvCounter]);
				argvCounter++;
			}
			else {
				opserr << "ERROR: Invalid input to SearchWithStepSizeAndStepDirection. " << endln;
				return TCL_ERROR;
			}
		}

		char ProbType[] = "reliability";
		
		theFindDesignPointAlgorithm = new SnoptProblem(
					maxNumIter, 
					theGFunEvaluator,
					theGradGEvaluator,
					theProbabilityTransformation,
					startAtOrigin,
					printFlag, fileNamePrint,
					ProbType,theReliabilityDomain);
/*
   snoptProblem::snoptProblem(int passedMaxNumberOfIterations, 
					GFunEvaluator *passedGFunEvaluator,
					GradGEvaluator *passedGradGEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					int pPrintFlag,
					char *pFileNamePrint,
					Vector *pStartPoint, char * probType):
   
   */
		
	//if SNOPT  ---- Quan
	}
#endif // _SNOPT
	else if (strcmp(argv[argvCounter],"NewStepSearch") == 0) {
		
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
//		if (theStartPoint == 0 ) {
//			opserr << "Need theStartPoint before a FindDesignPointAlgorithm can be created" << endln;
//			return TCL_ERROR;
//		}
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
	
		theFindDesignPointAlgorithm = new NewSearchWithStepSizeAndStepDirection(
					maxNumIter, theReliabilityDomain, 
					theGFunEvaluator,
					theGradGEvaluator,
					theStepSizeRule,
					theSearchDirection,
					theProbabilityTransformation,
					theHessianApproximation,
					theReliabilityConvergenceCheck,
					startAtOrigin,
					printFlag, fileNamePrint);
		
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

		RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
		while ((aRandomVariable = rvIter()) != 0) {
		  int tag = aRandomVariable->getTag();
		  double mean = aRandomVariable->getMean();
		  theReliabilityDomain->setStartPoint(tag, mean);
		}
		startAtOrigin = false;
	}
	else if (strcmp(argv[1],"Origin") == 0) {
		
	  /*
	  RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
	  while ((aRandomVariable = rvIter()) != 0) {
	    int tag = aRandomVariable->getTag();
	    theReliabilityDomain->setStartPoint(tag, 0.0);
	  }
	  */
	  startAtOrigin = true;
	}
	else if (strcmp(argv[1],"Given") == 0) {
	  // This is now the default option
	  startAtOrigin = false;
	}
	else if (strcmp(argv[1],"-file") == 0) {
	  opserr << "startPoint -file option is currently disabled" << endln;
	  return TCL_ERROR;
	  /*
		theStartPoint = new Vector(nrv);

		ifstream inputFile( argv[2], ios::in );
		if (inputFile.fail()) {
			opserr << "File " << argv[2] << " could not be opened. " << endln;
			return TCL_ERROR;
		}

		// Loop through file to see how many entries there are
		double dummy;
		int numEntries = 0;
		while (inputFile >> dummy) {
			numEntries++;
		}
		if (numEntries == 0) {
			opserr << "ERROR: No entries in the file read by startPoint!" << endln;
			return TCL_ERROR;
		}
		if (numEntries != nrv) {
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
	  */

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


	if (theGFunEvaluator == 0 ) {
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
							interp, argv[1], relSensTag);


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


	// Do input check
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	// Check number of arguments
	if (argc != 2)  {
		opserr << "ERROR: Wrong number of input parameter to FOSM analysis" << endln;
		return TCL_ERROR;
	}


	// Check for essential ingredients
	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a FOSMAnalysis can be created" << endln;
		return TCL_ERROR;
	}
	if (theGradGEvaluator == 0 ) {
		opserr << "Need theGradGEvaluator before a FOSMAnalysis can be created" << endln;
		return TCL_ERROR;
	}


	theFOSMAnalysis = new FOSMAnalysis( theReliabilityDomain,
											theGFunEvaluator,
											theGradGEvaluator,
											interp,
											argv[1]);

	if (theFOSMAnalysis == 0) {
		opserr << "ERROR: could not create theFOSMAnalysis \n";
		return TCL_ERROR;
	}

	// Now run the analysis
	theFOSMAnalysis->analyze();

	return TCL_OK;
}






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
	if (theGradGEvaluator == 0 ) {
		opserr << "Need theGradGEvaluator before a ParametricReliabilityAnalysis can be created" << endln;
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
												  theGradGEvaluator,
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
		theSystemAnalysis = new PCM(theReliabilityDomain, argv[1], aType, betaFile, rhoFile);
	else if (strcmp(argv[2],"IPCM") == 0)
		theSystemAnalysis = new IPCM(theReliabilityDomain, argv[1], aType, betaFile, rhoFile);
	else if (strcmp(argv[2],"MVN") == 0)
		theSystemAnalysis = new MVNcdf(theReliabilityDomain, argv[1], aType, betaFile, rhoFile, nMax, tol);
	else if (strcmp(argv[2],"SCIS") == 0)
		theSystemAnalysis = new SCIS(theReliabilityDomain, argv[1], aType, betaFile, rhoFile, nMax, tol);
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
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	// Check for essential tools
	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a SimulationAnalyis can be created" << endln;
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
	int numberOfSimulations	= 1000;
	double targetCOV		= 0.05;
	double samplingVariance	= 1.0;
	int printFlag			= 0;
	int analysisTypeTag		= 1;


	for (int i=2; i<argc; i=i+2) {

		if (strcmp(argv[i],"-type") == 0) {

			if (strcmp(argv[i+1],"failureProbability") == 0) {

				analysisTypeTag = 1;

//				if (theStartPoint == 0 ) {
//					opserr << "Need theStartPoint before a SimulationAnalyis can be created" << endln;
//					return TCL_ERROR;
//				}
			}

// Michele and Quan -------------------------
			else if (strcmp(argv[i+1],"outCrossingFailureProbability") == 0) {

				analysisTypeTag = 4;

//				if (theStartPoint == 0 ) {
//					opserr << "Need theStartPoint before a SimulationAnalyis can be created" << endln;
//					return TCL_ERROR;
//				}
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
				// Make sure that the mean point is the sampling center
				int nrv = theReliabilityDomain->getNumberOfRandomVariables();
				RandomVariable *aRandomVariable;
				RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
				while ((aRandomVariable = rvIter()) != 0) {
				  int tag = aRandomVariable->getTag();
				  double mean = aRandomVariable->getMean();
				  theReliabilityDomain->setStartPoint(tag, mean);
				}
				opserr << "NOTE: The startPoint is set to the Mean due to the selected sampling analysis type." << endln;
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
	
	
	theImportanceSamplingAnalysis 
			= new ImportanceSamplingAnalysis(theReliabilityDomain, 
							 theProbabilityTransformation, 
							 theGFunEvaluator, 
							 theRandomNumberGenerator, 
							 startAtOrigin,
							 interp,
							 numberOfSimulations, 
							 targetCOV,
							 samplingVariance,
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
				theGFunEvaluator,
				theGradGEvaluator,
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
	theReliabilityDomain->getStartPoint(*theDesignPoint);

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
									theGFunEvaluator, 
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
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before a GFunVisualizationAnalysis can be created" << endln;
		return TCL_ERROR;
	}

	if (theProbabilityTransformation == 0 ) {
		opserr << "Need theProbabilityTransformation before a GFunVisualizationAnalysis can be created" << endln;
		//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
		opserr << "Assume all RV's are independent" << endln;
		theProbabilityTransformation = 
		new AllIndependentTransformation(theReliabilityDomain,0);
//		return TCL_ERROR;
//////////////////////////////////////////////////////////////////////////////////
///////////// Modified by K Fujimura /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

		
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
					   theGFunEvaluator, 
					   theProbabilityTransformation, 
					   startAtOrigin,
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
											theGFunEvaluator, 
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
TclReliabilityModelBuilder_rvReduction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

	// Check that this is not done after transformation etc has been created
	if (theProbabilityTransformation != 0 ) {
		opserr << "ERROR: R.v. reduction cannot be performed after the " << endln
			<< " probability transformation has been created." << endln;
		return TCL_ERROR;
	}


	// Initial declarations
	int i, j;
	int rvNum;
	bool result;
	bool isInList, rv1isInList, rv2isInList;
	Vector keepRvs;






	// Read user-given numbers of random variables to keep
	int numInList;
	if (strcmp(argv[1],"-list") == 0) {
		numInList = argc-2;
		Vector tempKeepRvs(numInList);
		for (i=1; i<=numInList; i++) {
			if (Tcl_GetInt(interp, argv[i+1], &rvNum) != TCL_OK) {
				opserr << "ERROR: Invalid input to r.v. reduction command." << endln;
				return TCL_ERROR;
			}
			tempKeepRvs(i-1) = rvNum;
		}
		keepRvs = tempKeepRvs;
	}
	else if (strcmp(argv[1],"-file") == 0) {

		// Check how many entries we should read
		if (Tcl_GetInt(interp, argv[3], &numInList) != TCL_OK) {
			opserr << "ERROR: Invalid input to r.v. reduction command." << endln;
			return TCL_ERROR;
		}

		// Read entries
		Vector tempKeepRvs(numInList);
		ifstream inputFile( argv[2], ios::in );
		if (inputFile.fail()) {
			opserr << "File " << argv[2] << " could not be opened. " << endln;
			return TCL_ERROR;
		}
		for (int i=0; i<numInList; i++) {
				inputFile >> tempKeepRvs(i);
		}
		inputFile.close();
		keepRvs = tempKeepRvs;
	}
	else {
		opserr << "ERROR: Invalid argument to r.v. reduction. " << endln;
		return TCL_ERROR;
	}



	// Handle re-creation of random variables
	ArrayOfTaggedObjects aListOfRandomVariables(256);
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	RandomVariable *aRandomVariable;
	int newTag = 1;
	for (i=1; i<=nrv; i++) {
		isInList = false;
		for (j=1; j<=numInList; j++) {
			if (keepRvs(j-1)==i) {
				isInList = true;
				newTag = j;
				break;
			}
		}
		aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i);
		if (!isInList) {
			delete aRandomVariable;
		}
		else {
			aRandomVariable->setNewTag(newTag);
			result = aListOfRandomVariables.addComponent(aRandomVariable);
		}
		theReliabilityDomain->removeRandomVariable(i);
	}
	// Add the new random variables to the reliability domain
	for (i=1; i<=numInList; i++) {
		theReliabilityDomain->addRandomVariable((RandomVariable*)aListOfRandomVariables.getComponentPtr(i));
	}



	// Handle the re-creation of correlation coefficients
	int nCorr = theReliabilityDomain->getNumberOfCorrelationCoefficients();
	int count = 1;
	ArrayOfTaggedObjects aListOfCorrelationCoefficients(256);
	CorrelationCoefficient *aCorrelationCoefficient;
	CorrelationCoefficient *aNewCorrelationCoefficient;
	int rv1, rv2, newRv1 = 1, newRv2 = 1;
	double rho;
	for (i=1; i<=nCorr; i++) {

		// Get correlation data
		aCorrelationCoefficient = theReliabilityDomain->getCorrelationCoefficientPtr(i);
		rv1 = aCorrelationCoefficient->getRv1();
		rv2 = aCorrelationCoefficient->getRv2();
		rho = aCorrelationCoefficient->getCorrelation();

		// Get random variable numbers and possibly create a new correlation coefficient
		rv1isInList = false;
		rv2isInList = false;
		for (j=1; j<=numInList; j++) {
			if (keepRvs(j-1)==rv1) {
				rv1isInList = true;
				newRv1 = j;
			}
			if (keepRvs(j-1)==rv2) {
				rv2isInList = true;
				newRv2 = j;
			}
		}
		if (rv1isInList && rv2isInList) {
			aNewCorrelationCoefficient = new CorrelationCoefficient(count,newRv1,newRv2,rho);
			aListOfCorrelationCoefficients.addComponent(aNewCorrelationCoefficient);
			count++;
		}
		delete aCorrelationCoefficient;
		theReliabilityDomain->removeCorrelationCoefficient(i);
	}
	// Add the new correlation coefficients to the reliability domain
	for (i=1; i<=(count-1); i++) {
		theReliabilityDomain->addCorrelationCoefficient((CorrelationCoefficient*)aListOfCorrelationCoefficients.getComponentPtr(i));
	}






	// Handle the re-creation of random variable positioners
	int nPos = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	RandomVariablePositioner *aRvPositioner;
	ArrayOfTaggedObjects aListOfPositioners(256);
	count =1;
	int newRvIndex;
	for (i=1; i<=nPos; i++) {
		aRvPositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		rvNum = aRvPositioner->getRvIndex();
		isInList = false;
		for (j=1; j<=numInList; j++) {
			if (keepRvs(j-1)==rvNum) {
				isInList = true;
				newRvIndex = j;
				break;
			}
		}
		if (!isInList) {
			delete aRvPositioner;
		}
		else {
			aRvPositioner->setNewTag(count);
			aRvPositioner->setRvIndex(newRvIndex);
			aListOfPositioners.addComponent(aRvPositioner);
			count++;
		}
		theReliabilityDomain->removeRandomVariablePositioner(i);
	}
	// Add the new random variable positioners to the reliability domain
	for (i=1; i<=numInList; i++) {
		theReliabilityDomain->addRandomVariablePositioner((RandomVariablePositioner*)aListOfPositioners.getComponentPtr(i));
	}





/*
nrv = theReliabilityDomain->getNumberOfRandomVariables();
RandomVariable *aRV;
for (i=1; i<=nrv; i++) {
	aRV = theReliabilityDomain->getRandomVariablePtr(i);
	opserr << "Mean of rv# " << i << ": " << aRV->getMean() << endln;

}
nrv = theReliabilityDomain->getNumberOfCorrelationCoefficients();
CorrelationCoefficient *aaRV;
for (i=1; i<=nrv; i++) {
	aaRV = theReliabilityDomain->getCorrelationCoefficientPtr(i);
	opserr << "Rv1: " << aaRV->getRv1() << endln;
	opserr << "Rv2: " << aaRV->getRv2() << endln;
	opserr << "Corr: " << aaRV->getCorrelation() << endln;
}
nrv = theReliabilityDomain->getNumberOfRandomVariablePositioners();
RandomVariablePositioner *aaaRV;
for (i=1; i<=nrv; i++) {
	aaaRV = theReliabilityDomain->getRandomVariablePositionerPtr(i);
	opserr << "Rv in pos: " << aaaRV->getRvNumber() << endln;
}

while(true) {
}
*/


	return TCL_OK;
}





//////////////////////////////////////////////////////////////////
int 
TclReliabilityModelBuilder_inputCheck(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// Check that tagged objects are consequtive
	int i, num;
	ReliabilityDomainComponent *component;

	/*
	num = theReliabilityDomain->getNumberOfRandomVariables();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getRandomVariablePtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive random variable list." << endln;
			return TCL_ERROR;
		}
	}
	*/

	// Clear out old parameter positioners so we don't produce a memory leak
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

	/*
	num = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive random variable positioner list." << endln;
			return TCL_ERROR;
		}
	}
	*/
	/*
	num = theReliabilityDomain->getNumberOfCorrelationCoefficients();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getCorrelationCoefficientPtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive correlation coefficient list." << endln;
			return TCL_ERROR;
		}
	}
	*/

	num = theReliabilityDomain->getNumberOfFilters();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getFilter(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive filter list." << endln;
			return TCL_ERROR;
		}
	}
	
	/*
	num = theReliabilityDomain->getNumberOfLimitStateFunctions();
	for (i=1; i<=num; i++) {
		component = theReliabilityDomain->getLimitStateFunctionPtr(i);
		if (component == 0) {
			opserr << "ERROR: Non-consequtive limit-state (performance) function list." << endln;
			return TCL_ERROR;
		}
	}
	*/

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
TclReliabilityModelBuilder_printReliability(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc > 1)
    theReliabilityDomain->Print(opserr, 1);
  else
    theReliabilityDomain->Print(opserr);

  return TCL_OK;
}


// ---------- Quan Gu ------------------------
// command: designVariable  1 -name E  < -startPt $E   -lowerBound [expr $E*0.8] -upperBound [expr $E*1.2] >
int 
TclReliabilityModelBuilder_addDesignVariable(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{
  DesignVariable *theDesignVariable = 0;
  int tag;
  char name[20]="";
  char valueString[20]="";
  double value=0;
  double lowerBound=0;
  double upperBound=0;
  int numberOfArguments = argc;

  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 3) {
		opserr << "ERROR: invalid number of arguments to designVariable command \n";
		return TCL_ERROR;
  }


  // GET TAG NUMBER
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-name") == 0)||(strcmp(argv[argvCounter],"-Name") == 0)) {
			argvCounter++;
			strcpy(name,argv[argvCounter]);
			argvCounter++;
		}// if

		else if (strcmp(argv[argvCounter],"-startPt") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &value) != TCL_OK) {
			opserr << "ERROR: invalid input: startPt \n";
			return TCL_ERROR;
			}
			strcpy(valueString,argv[argvCounter]);
			lowerBound=value;
			upperBound=value;
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-lowerBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &lowerBound) != TCL_OK) {
			opserr << "ERROR: invalid input: startPt \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-upperBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &upperBound) != TCL_OK) {
			opserr << "ERROR: invalid input: startPt \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if
		
		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// in tcl, run 'set name value'
  char tclAssignment[50];
  strcpy (tclAssignment, "set ");
  strcat (tclAssignment, name);
  strcat (tclAssignment, " ");

   
 // _gcvt( value, 7, buffer );
   strcat (tclAssignment, valueString);
  
   
   if (Tcl_GetVar(interp, name, TCL_GLOBAL_ONLY) !=NULL)
   {
		opserr<<"Fatal::the variable with name: "<<name <<" is already in system, please use another name!"<<endln;   
		exit(-1);
   }

   if (Tcl_Eval(interp,tclAssignment) !=TCL_OK ){
   		opserr<<"Fatal::can not set varuable with name: "<<name <<"in tcl command!"<<endln;   
		exit(-1);
   }

// here tag ;
  theDesignVariable = new DesignVariable(tag, 
			 name,
			 value,
			 upperBound,
			 lowerBound,
			 interp,
			 theReliabilityDomain,
			 0);

  if (theDesignVariable == 0) {
		opserr << "ERROR: could not create random theDesignVariable number " << tag << endln;
		return TCL_ERROR;
	  }
  

  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addDesignVariable(theDesignVariable) == false) {
	opserr << "ERROR: failed to add theDesignVariable variable to the domain (wrong number of arguments?)\n";
	opserr << "theDesignVariable variable: " << tag << endln;
	delete theDesignVariable; // otherwise memory leak
	return TCL_ERROR;
  }

  return TCL_OK;
}
					   

// command designVariablePositioner 1   -dvNum 1 -element 1     -material E  
int 
TclReliabilityModelBuilder_addDesignVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	DesignVariablePositioner *theDesignVariablePositioner = 0;
	int tag;
	int dvNumber;
	int tagOfObject;
	DomainComponent *theObject;
	int argvCounter = 1;


	// READ THE TAG NUMBER
	if (Tcl_GetInt(interp, argv[argvCounter++], &tag) != TCL_OK) {
		opserr << "ERROR: Invalid input tag to design variable positioner." << endln;
		return TCL_ERROR;
	}


	if (strcmp(argv[argvCounter],"-dvNum") == 0) {
		argvCounter++;
		
		// READ THE RANDOM VARIABLE NUMBER
		if (Tcl_GetInt(interp, argv[argvCounter++], &dvNumber) != TCL_OK) {
			opserr << "ERROR: invalid input: dvNumber \n";
			return TCL_ERROR;
		}

		// CHECK THAT THE RANDOM VARIABLE ACTUALLY EXISTS
		DesignVariable *theDesignVariable = 0;
		theDesignVariable = theReliabilityDomain->getDesignVariablePtr(dvNumber);
		if (theDesignVariable == 0){
			opserr << "ERROR:: A non-existing design variable number " << dvNumber << " is being positioned in the model " << endln;
			return TCL_ERROR;
		}
	}
	else {
		opserr << "ERROR: Illegal design variable specification in  " << endln
			<< " design variable positioner command. " << endln;
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

		theDesignVariablePositioner = new DesignVariablePositioner(tag,
									   theReliabilityDomain,
									   dvNumber,
									   theObject,
									   data,
									   argc-argvCounter);


		int dvnumber = theDesignVariablePositioner->getDVNumber();
	}

	// IF UNCERTAIN *LOAD*
	else if (strcmp(argv[argvCounter],"-loadPattern") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getLoadPattern(tagOfObject);


		theDesignVariablePositioner = new DesignVariablePositioner(tag,
									   theReliabilityDomain,
									   dvNumber,
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

		theDesignVariablePositioner = new DesignVariablePositioner(tag,
							           theReliabilityDomain,
									   dvNumber,
									   theObject,
									   data,
									   argc-argvCounter);
	}
	else {
		opserr << "ERROR: Unknown parameter in designVariablePositioner" << endln;
		return TCL_ERROR;
	}

	delete [] data;

	// ADD THE RANDOMVARIABLEPOSITIONER TO THE DOMAIN
	if (theReliabilityDomain->addDesignVariablePositioner(theDesignVariablePositioner) == false) {
		opserr << "ERROR: failed to add random variable positioner number " << tag << " to the domain." << endln;
		delete theDesignVariablePositioner; // otherwise memory leak
		return TCL_ERROR;
	}

	return TCL_OK;

}


// command: objectiveFunction 1  -name F  -GradientName G -tclFile objective.tcl  ; # -lowerBound -1.e20 -upperBound 1.e20  - multiplier $E -state 0 -linearAdd A ;
int 
TclReliabilityModelBuilder_addObjectiveFunction(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{
  ObjectiveFunction *theObjectiveFunction = 0;
  int tag;
  char name[25] = "" ;
  char * gradientName = 0;
  char tclFileName[35] = "";
  double lowerBound= -1.e20;;
  double upperBound=1.e20;
  double multiplier = 0; 
  double state = 0; 
  int numberOfArguments = argc;

  Vector *linearAdd=0;

  bool isGradProvided=false;

  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 7) {
		opserr << "ERROR: invalid number of arguments to designVariable command "<<endln;
	    opserr <<"command: objectiveFunction 1  -name F  -GradientName G  -tclFile objective.tcl  ; # -lowerBound -1.e20 -upperBound 1.e20 - multiplier $E -state 0 -linearAdd A"<<endln;
		return TCL_ERROR;
  }


  // GET TAG NUMBER
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-name") == 0)||(strcmp(argv[argvCounter],"-Name") == 0)) {

			argvCounter++;
			strcpy(name,argv[argvCounter]);

			if ((Tcl_GetVar(interp, name, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: objectiveFunction Name "<< name << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			
			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-GradientName") == 0)||(strcmp(argv[argvCounter],"-gradientName") == 0)) {
			
			argvCounter++;
			
			gradientName = new char[30];

			strcpy(gradientName,argv[argvCounter]);
			
			if ((Tcl_GetVar(interp, gradientName, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: objectiveFunction gradient Name "<< gradientName << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			isGradProvided = true;
			argvCounter++;
		}// else if
		
		else if ((strcmp(argv[argvCounter],"-tclFile") == 0)||(strcmp(argv[argvCounter],"-TclFile") == 0)) {
			argvCounter++;
			strcpy(tclFileName,argv[argvCounter]);
			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-lowerBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &lowerBound) != TCL_OK) {
			opserr << "ERROR: invalid input: lowerBound \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-upperBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &upperBound) != TCL_OK) {
			opserr << "ERROR: invalid input: upperBound \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-multiplier") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &multiplier) != TCL_OK) {
			opserr << "ERROR: invalid input: multiplier \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if		

		else if (strcmp(argv[argvCounter],"-state") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &state) != TCL_OK) {
			opserr << "ERROR: invalid input: state \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if		

		else if (strcmp(argv[argvCounter],"-linearAdd") == 0) {
			char linearAddName[20];

			argvCounter++;
			strcpy(linearAddName,argv[argvCounter]);


         if( Tcl_GetVar2(interp, linearAddName,"1",TCL_GLOBAL_ONLY ) != NULL)  {

			int numOfDV = theReliabilityDomain ->getNumberOfDesignVariables();
			linearAdd = new Vector(numOfDV);

			const char *  theValue;	
			char index[5];
			for(int i=0; i<numOfDV; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, linearAddName,index,TCL_GLOBAL_ONLY );
				(*linearAdd)(i) = atof(theValue);
			};
			
		 }  //if 
		  else {
			opserr<<"warning: the linearAdd with name "<<linearAddName<< " does not exit"<<endln;
		 } // else

				argvCounter++;
		}// else if	"-linearAdd"

		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// here tag ;
  theObjectiveFunction = new ObjectiveFunction(   tag,
												  theReliabilityDomain,
												  interp,
												  isGradProvided,
											      linearAdd,
												  tclFileName,
												  name, 
												  gradientName, 
												  lowerBound, 
												  upperBound, 
												  multiplier, 									   
												  state		   
												  );




  if (theObjectiveFunction == 0) {
		opserr << "ERROR: could not create random theObjectiveFunction "<< endln;
		return TCL_ERROR;
	  }
 // release memory 
	if (gradientName !=0) delete gradientName;
	if (linearAdd !=0)	 delete linearAdd;


  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addObjectiveFunction(theObjectiveFunction) == false) {
	opserr << "ERROR: failed to add theObjectiveFunction  to the domain (wrong number of arguments?)\n";

	delete theObjectiveFunction; // otherwise memory leak

	return TCL_ERROR;
  }

  return TCL_OK;
}



/*  command 

array set lBound1{1  -2.0  2  -3.0 }
array set uBound1{1   2.0  2   3.0 }

array set M1 {1 0.0  2  0.0 }
array set S1{1  0.0  2    0.0 }

for {set i 1} {$i <=10} {incr i} {
   for {set j 1} {$j < 10} {incr j} {
     set A1($i,$j) [expr $i*10+$j]
   }
 }


constraintFunction 1  -name F1  -GradientName G1 -lowerBound  lBound1 -upperBound uBound1 -tclFile constraint1.tcl -multiplier M1 -state S1 -linearAdd A1 ;



  */
int 
TclReliabilityModelBuilder_addConstraintFunction(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv){

  ConstraintFunction *theConstraintFunction = 0;
  int tag;
  int numberOfConstraints;

  char name[25] = "" ;
  char * gradientName = 0;
  char tclFileName[35] = "";

  Vector * lowerBound = 0;
  Vector * upperBound = 0;
  Vector * multiplier = 0; 
  Vector * state = 0; 
  Matrix *linearAdd=0;

  int numberOfArguments = argc;
  bool isGradProvided=false;



  /* refer
  ConstraintFunction(int passedTag, int passedNumberOfConstraint,
									   ReliabilityDomain * passedReliabilityDomain, 
									   Tcl_Interp *passedTclInterp, 

									   bool passedIsGradProvided, 
									   Matrix * passedLinearAdd,
									   char * passedTclFileName,
									   char * passedName, 
									   char * passedGradientName, 
									   
									   Vector * passedLowerBound, 
									   Vector * passedUpperBound, 
									   Vector * passedMultiplier, 									   
									   Vector * passedState
									   );  */

  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 6) {
		opserr << "ERROR: invalid number of arguments to designVariable command "<<endln;
	    opserr <<"command: constraintFunction 1  -name F1  -GradientName G1 -lowerBound  lBound1 -upperBound uBound1 -tclFile constraint1.tcl -multiplier M1 -state S1 -linearAdd A1 "<<endln;
		return TCL_ERROR;
  }


  // GET TAG NUMBER
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-name") == 0)||(strcmp(argv[argvCounter],"-Name") == 0)) {

			argvCounter++;
			strcpy(name,argv[argvCounter]);

			if ((Tcl_GetVar(interp, name, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: ConstraintFunction Name "<< name << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			
			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-GradientName") == 0)||(strcmp(argv[argvCounter],"-gradientName") == 0)) {
			
			argvCounter++;
			
			gradientName = new char[30];

			strcpy(gradientName,argv[argvCounter]);
			
			if ((Tcl_GetVar(interp, gradientName, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: ConstraintFunction gradient Name "<< gradientName << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			isGradProvided = true;
			argvCounter++;
		}// else if
		
		else if ((strcmp(argv[argvCounter],"-tclFile") == 0)||(strcmp(argv[argvCounter],"-TclFile") == 0)) {
			argvCounter++;
			strcpy(tclFileName,argv[argvCounter]);
			argvCounter++;
		}// else if
		

		else if (strcmp(argv[argvCounter],"-lowerBound") == 0) {
			char lowerBoundName[20];
			argvCounter++;
			strcpy(lowerBoundName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, lowerBoundName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: lowerBound error!"<<endln; exit(-1);}
//			int numOfDV = theReliabilityDomain ->getNumberOfDesignVariables();
			
			lowerBound = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, lowerBoundName,index,TCL_GLOBAL_ONLY );
				(*lowerBound)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-lowerBound"


		else if (strcmp(argv[argvCounter],"-upperBound") == 0) {
			char upperBoundName[20];
			argvCounter++;
			strcpy(upperBoundName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, upperBoundName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: upperBound error!"<<endln; exit(-1);}
//			int numOfDV = theReliabilityDomain ->getNumberOfDesignVariables();
			
			upperBound = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, upperBoundName,index,TCL_GLOBAL_ONLY );
				(*upperBound)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-upperBound"


		else if (strcmp(argv[argvCounter],"-state") == 0) {
			char stateName[20];
			argvCounter++;
			strcpy(stateName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, stateName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: state error!"<<endln; exit(-1);}
//			int numOfDV = theReliabilityDomain ->getNumberOfDesignVariables();
			
			state = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, stateName,index,TCL_GLOBAL_ONLY );
				(*state)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-state"



		else if (strcmp(argv[argvCounter],"-multiplier") == 0) {
			char multiplierName[20];
			argvCounter++;
			strcpy(multiplierName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, multiplierName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: lmultiplier error!"<<endln; exit(-1);}
//			int numOfDV = theReliabilityDomain ->getNumberOfDesignVariables();
			
			multiplier = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, multiplierName,index,TCL_GLOBAL_ONLY );
				(*multiplier)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-multiplier"

		
		
		else if (strcmp(argv[argvCounter],"-linearAdd") == 0) {
			char linearAddName[20];

			argvCounter++;
			strcpy(linearAddName,argv[argvCounter]);

// get Number of constraintFunction
		
			int ii=1;
			char index[10];
			sprintf(index,"%d",ii);
			strcat(index,",1");

			while ( Tcl_GetVar2(interp, linearAddName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
				strcat(index,",1");
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: linearAdd error!"<<endln; exit(-1);}

			int numOfDVs = theReliabilityDomain ->getNumberOfDesignVariables();

 ////////////////////////////////////////////////
         if( Tcl_GetVar2(interp, linearAddName,"1,1",TCL_GLOBAL_ONLY ) != NULL)  {

			linearAdd = new Matrix(numberOfConstraints,numOfDVs);

			const char *  theValue;	
			char temp[5];

			for(int i=0;i<numberOfConstraints; i++){
				for(int j=0; j<numOfDVs; j++){

			       sprintf(temp,"%d",i+1);   // begin with 1
				   strcpy(index,temp);
				   sprintf(temp,"%d",j+1);   // begin with 1
				   strcat(index,",");
				   strcat(index,temp);
					
				   theValue = Tcl_GetVar2(interp, linearAddName,index,TCL_GLOBAL_ONLY );
				   (*linearAdd)(i,j) = atof(theValue);
				}; //for
			} //for


		 }  //if 
		  else {
			opserr<<"warning: the linearAdd with name "<<linearAddName<< " does not exit"<<endln;
		 } // else

				argvCounter++;
		}// else if	"-linearAdd"



		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// here tag ;
/*
  ConstraintFunction *theConstraintFunction = 0;
  int tag;
  int numberOfConstraints;

  char name[25] = "" ;
  char * gradientName = 0;
  char tclFileName[35] = "";

  Vector * lowerBound = 0;
  Vector * upperBound = 0;
  Vector * multiplier = 0; 
  Vector * state = 0; 
  Matrix *linearAdd=0;

  int numberOfArguments = argc;
  bool isGradProvided=false;



   refer
  ConstraintFunction(int passedTag, int passedNumberOfConstraint,
									   ReliabilityDomain * passedReliabilityDomain, 
									   Tcl_Interp *passedTclInterp, 

									   bool passedIsGradProvided, 
									   Matrix * passedLinearAdd,
									   char * passedTclFileName,
									   char * passedName, 
									   char * passedGradientName, 
									   
									   Vector * passedLowerBound, 
									   Vector * passedUpperBound, 
									   Vector * passedMultiplier, 									   
									   Vector * passedState
									   );  */
  theConstraintFunction = new ConstraintFunction(tag, 
									   numberOfConstraints,
									   theReliabilityDomain, 
									   interp, 

									   isGradProvided, 
									   linearAdd,
									   tclFileName,
									   name, 
									   gradientName, 
									   
									   lowerBound, 
									   upperBound, 
									   multiplier, 									   
									   state
									   ); 




  if (theConstraintFunction == 0) {
		opserr << "ERROR: could not create random theConstraintFunction "<< endln;
		return TCL_ERROR;
	  }
  
// delete 

if (gradientName !=0) delete gradientName;
if (linearAdd !=0)	 delete linearAdd;
if (lowerBound !=0) delete lowerBound;
if (upperBound !=0) delete upperBound;
if (multiplier !=0) delete multiplier;
if (state !=0) delete state;



  // ADD THE OBJECT TO THE DOMAIN
  if (theReliabilityDomain->addConstraintFunction(theConstraintFunction) == false) {
	opserr << "ERROR: failed to add theConstraintFunction  to the domain (wrong number of arguments?)\n";

	delete theConstraintFunction; // otherwise memory leak
 


	return TCL_ERROR;
  }




  return TCL_OK;

}






// command: runSNOPTAnalysis -maxNumIter 100 -printOptPointX OptX.out -tclFileToRun tclFileToRun.tcl -printFlag 1
int 
TclReliabilityModelBuilder_runSNOPTAnalysis(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{

#ifndef _SNOPT
  opserr << "TclReliabilityModelBuilder_runSNOPTAnalysis() - NO SNOPT LINKED\n";
  return TCL_ERROR;
#else

  int maxNumberOfIterations=100;
  int printFlag=0;
  char fileNamePrint[25]="";
  char probType[25]="SNOPTAnalysis";
  char * tclFileName=0;

  int numberOfArguments = argc;
  
  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 3) {
		opserr << "ERROR: invalid number of arguments to designVariable command "<<endln;
	    opserr <<"command: runSNOPTAnalysis -maxNumIter 100 -printOptPointX OptX.out -tclFileToRun tclFileToRun.tcl"<<endln;
		return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 1;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-maxNumIter") == 0)||(strcmp(argv[argvCounter],"-maxnumiter") == 0)) {

			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &maxNumberOfIterations) != TCL_OK) {
			opserr << "ERROR: invalid input: maxNumberOfIterations \n";
			return TCL_ERROR;
			}

			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-printFlag") == 0)||(strcmp(argv[argvCounter],"-printflag") == 0)) {

			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &printFlag) != TCL_OK) {
			opserr << "ERROR: invalid input: printFlag \n";
			return TCL_ERROR;
			}

			argvCounter++;
		}// if
		else if ((strcmp(argv[argvCounter],"-printOptPointX") == 0)||(strcmp(argv[argvCounter],"-printoptpointx") == 0)) {
			
			argvCounter++;
			
//  		gradientName = new char[30];

			strcpy(fileNamePrint,argv[argvCounter]);
			
			argvCounter++;
		}// else if
		
		else if ((strcmp(argv[argvCounter],"-tclFileToRun") == 0)||(strcmp(argv[argvCounter],"-TclFileToRun") == 0)) {
			argvCounter++;
			tclFileName = new char[30];
			strcpy(tclFileName,argv[argvCounter]);
			argvCounter++;
		}// else if
		
	
		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// here tag ;
  theSNOPTAnalysis = new SNOPTAnalysis(maxNumberOfIterations, 
					printFlag,
					fileNamePrint,
					probType, 
					theReliabilityDomain,
					interp,
					tclFileName
					);




  if (theSNOPTAnalysis == 0) {
		opserr << "ERROR: could not create random theSNOPTAnalysis "<< endln;
		return TCL_ERROR;
	  }


	theSNOPTAnalysis->runOptAnalysis(theReliabilityDomain);

	if (tclFileName !=0) delete tclFileName;
  return TCL_OK;
#endif 
}



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
	char theCommand[15] = "inputCheck";
	Tcl_Eval( interp, theCommand );


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

/* refer 
			MonteCarloResponseAnalysis(
						ReliabilityDomain *passedReliabilityDomain,
						ProbabilityTransformation *passedProbabilityTransformation,
						RandomNumberGenerator *passedRandomNumberGenerator,
						int passedNumberOfSimulations,
						int printFlag,
						TCL_Char *outputFileName,
						TCL_Char *tclFileToRunFileName)
{
  
  */
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

	/*** FMK MHS QUAN
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


/*
Hessian::Hessian(int pSize,ReliabilityDomain *passedReliabilityDomain, 
				 ProbabilityTransformation *passedProbabilityTransformation,
				 GFunEvaluator *passedGFunEvaluator,
				 GradGEvaluator *passedGradGEvaluator)
				 */

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
	

	Hessian * theHessian = new Hessian(size,theReliabilityDomain,theProbabilityTransformation,theGFunEvaluator,theGradGEvaluator,pTol);
 
	if (FDM) { 
		theHessian->formReducedHessian(designPoint);
		
		
		Matrix hessian=	theHessian->getHessianApproximation();

//		resultsOutputFile <<"Hessian in U space: \n";
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				resultsOutputFile <<hessian(i,j)<<"   ";
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
					theGFunEvaluator,
					theProbabilityTransformation,
					outputFile,
					theGradGEvaluator,
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



int 
TclReliabilityModelBuilder_transformXtoU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{


    if (theProbabilityTransformation ==0){
		opserr<<"Fatal: theProbabilityTransformation does not exist!"<<endln;
		exit(-1);
	}
	
	
	
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	
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

	/*
	theProbabilityTransformation->set_x(pointX);
	theProbabilityTransformation->transform_x_to_u();
	Vector pointU = theProbabilityTransformation->get_u();
	*/
	Vector pointU;
	theProbabilityTransformation->transform_x_to_u(pointX, pointU);

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


int 
TclReliabilityModelBuilder_transformUtoX(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{


    if (theProbabilityTransformation ==0){
		opserr<<"Fatal: theProbabilityTransformation does not exist!"<<endln;
		exit(-1);
	}
	
	
	
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	
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
	/*
	theProbabilityTransformation->set_u(pointU);
	theProbabilityTransformation->transform_u_to_x();
	Vector pointX = theProbabilityTransformation->get_x();
	*/
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


/*ReliabilityDomain *passedReliabilityDomain,
					GFunEvaluator *passedGFunEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					char *passedOutputFileName,
					GradGEvaluator * passedGradGEvaluator, Vector * pDesignPt, int numAxis, 
					char * typeExpPtRule,Tcl_Interp *passedTclInterp, 
					Matrix * passedHessian, char * passedHessianFile, char * typeSurfaceDesign, 
					char * typeRespSurfaceSimulation, Vector * gridInfo,
					RandomNumberGenerator * pRandomNumberGenerator,
					double pTargetCOV,
					int pNumberOfSimulations*/
int 
TclReliabilityModelBuilder_runDP_RSM_SimTimeInvariantAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
 

	if ((theRandomNumberGenerator ==0)||(theGFunEvaluator==0)||(theProbabilityTransformation==0) ||(theGradGEvaluator==0)) {
	
		opserr<<"theRandomNumberGenerator ==0)||(theGFunEvaluator)||(theProbabilityTransformation==0) ||(theGradGEvaluator==0)"<<endln;
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
					theGFunEvaluator,
					theProbabilityTransformation,
					outputFile,
					theGradGEvaluator, 
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
 

	if ((theRandomNumberGenerator ==0)||(theGFunEvaluator==0)||(theProbabilityTransformation==0) ||(theGradGEvaluator==0)) {
	
		opserr<<"theRandomNumberGenerator ==0)||(theGFunEvaluator)||(theProbabilityTransformation==0) ||(theGradGEvaluator==0)"<<endln;
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
					theGFunEvaluator,
					theProbabilityTransformation,
					outputFile,
					theGradGEvaluator, 
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



				/*  ReliabilityDomain *passedReliabilityDomain,
					GFunEvaluator *passedGFunEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					char *passedOutputFileName,
					GradGEvaluator * passedGradGEvaluator, Vector * pDesignPt, int numAxis, 
					char * typeExpPtRule,Tcl_Interp *passedTclInterp, 
					Matrix * passedHessian, char * passedHessianFile, char * typeSurfaceDesign, 
					char * typeRespSurfaceSimulation, Vector * gridInfo,
					RandomNumberGenerator * pRandomNumberGenerator,
					double pTargetCOV,
					int pNumberOfSimulations, double pLittleDt,
					double ImpulseInterval*/

	if (passedHessian !=0) delete passedHessian;
	if (hessianFileName !=0) delete hessianFileName;
	theDP_RSM_Sim_TimeVariant->analyze();


    delete designPoint;
	if (gridInfo !=0) delete gridInfo;
	return 0;

}


///getBetaFORM is not in K.F.
int 
TclReliabilityModelBuilder_getBetaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "ERROR: Invalid number of arguments to getBetaFORM command." << endln;
    return TCL_ERROR;
  }

  int lsfTag;
  if (Tcl_GetInt(interp, argv[1], &lsfTag) != TCL_OK) {
    opserr << "WARNING betaFORM lsfTag? - could not read lsfTag\n";
    return TCL_ERROR;	        
  }   

  LimitStateFunction *theLSF =
    theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);

  if (theLSF == 0) {
    opserr << "WARNING betaFORM LSF with tag " << lsfTag << " not found\n";
    return TCL_ERROR;	        
  }

  double beta = theLSF->getFORM_beta();

  char buffer[40];
  sprintf(buffer,"%35.20f",beta);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

///getGammaFORM is not in K.F.
int 
TclReliabilityModelBuilder_getGammaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "ERROR: Invalid number of arguments to getGammaFORM command." << endln;
    return TCL_ERROR;
  }

  int lsfTag;
  if (Tcl_GetInt(interp, argv[1], &lsfTag) != TCL_OK) {
    opserr << "WARNING gammaFORM lsfTag? rvTag? - could not read lsfTag\n";
    return TCL_ERROR;	        
  }

  LimitStateFunction *theLSF =
    theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);
  
  if (theLSF == 0) {
    opserr << "WARNING gammaFORM LSF with tag " << lsfTag << " not found\n";
    return TCL_ERROR;	        
  }

  const Vector &gammaVec = theLSF->getFORM_gamma();
  
  char buffer[40];

  if (argc > 2) {
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
    double gamma = gammaVec(index);
    
    sprintf(buffer,"%35.20f",gamma);
    
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  }
  else {
    int nrv = gammaVec.Size();
    for (int i = 0; i < nrv; i++) {
      sprintf(buffer, "%35.20f ", gammaVec(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

///getGammaFORM is not in K.F.
int 
TclReliabilityModelBuilder_getAlphaFORM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "ERROR: Invalid number of arguments to getAlphaFORM command." << endln;
    return TCL_ERROR;
  }

  int lsfTag;
  if (Tcl_GetInt(interp, argv[1], &lsfTag) != TCL_OK) {
    opserr << "WARNING alphaFORM lsfTag? rvTag? - could not read lsfTag\n";
    return TCL_ERROR;	        
  }

  LimitStateFunction *theLSF =
    theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);
  
  if (theLSF == 0) {
    opserr << "WARNING alphaFORM LSF with tag " << lsfTag << " not found\n";
    return TCL_ERROR;	        
  }

  const Vector &gammaVec = theLSF->getFORM_alpha();
  
  char buffer[40];

  if (argc > 2) {
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
    double gamma = gammaVec(index);
    
    sprintf(buffer,"%35.20f",gamma);
    
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  }
  else {
    int nrv = gammaVec.Size();
    for (int i = 0; i < nrv; i++) {
      sprintf(buffer, "%35.20f ", gammaVec(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}


///invNormalCDF is not in K.F.
int
TclReliabilityModelBuilder_invNormalCDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  static NormalRV aStdNormal(0, 0.0, 1.0);

  double x;
  if (Tcl_GetDouble(interp, argv[1], &x) != TCL_OK) {
    opserr << "WARNING invNormalCDF x? <mean? stdev?>- could not read x\n";
    return TCL_ERROR;	        
  }

  char buffer[40];

  if (argc < 4) {
    sprintf(buffer,"%35.20f", aStdNormal.getInverseCDFvalue(x));
  }
  else {
    double mean, stdev;
    if (Tcl_GetDouble(interp, argv[2], &mean) != TCL_OK) {
      opserr << "WARNING invNormalCDF x? mean? stdev? - could not read mean\n";
      return TCL_ERROR;	        
    }
    if (Tcl_GetDouble(interp, argv[3], &stdev) != TCL_OK) {
      opserr << "WARNING invNormalCDF x? mean? stdev? - could not read stdev\n";
      return TCL_ERROR;	        
    }
    sprintf(buffer,"%35.20f", mean + stdev*aStdNormal.getInverseCDFvalue(x));
  }

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
	int nfound=0;
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
	if(theGFunEvaluator==0){
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
//								 theGFunEvaluator,
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
								 theGFunEvaluator,
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
		if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an CrossingRateAnalyzer can be created" << endln;
		return TCL_ERROR;
		}
		if (theGradGEvaluator == 0 ) {
		opserr << "Need theGradGEvaluator before an CrossingRateAnalyzer can be created" << endln;
		return TCL_ERROR;
		}
	}
	theCrossingRateAnalyzer = new CrossingRateAnalyzer
								(theReliabilityDomain,
								 theFindDesignPointAlgorithm,
								 theGFunEvaluator,
								 theGradGEvaluator,
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
	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an FirstPassageAnalyzer can be created" << endln;
		return TCL_ERROR;
	}

	if(stationary==1){
		theFirstPassageAnalyzer=new StatFirstPassageAnalyzer
									(theReliabilityDomain,
									 theFindDesignPointAlgorithm,
									 theGFunEvaluator,
									 theFOSeriesSimulation,
									 analysis,
									 twoside,
									 print);
	}else{
		theFirstPassageAnalyzer=new NonStatFirstPassageAnalyzer
								   (theReliabilityDomain,
									theFindDesignPointAlgorithm,
									theGFunEvaluator,
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
	if (theGFunEvaluator == 0 ) {
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

	if (theGFunEvaluator == 0 ) {
		opserr << "Need theGFunEvaluator before an OutCrossingAnalysis can be created" << endln;
		return TCL_ERROR;
	}

	if(stationary==1){
	 	theRandomVibrationSimulation=new StatRandomVibrationSimulation
								(theReliabilityDomain,
								 theStructuralDomain,
						         theGFunEvaluator,
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
						         theGFunEvaluator,
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
	if (theGFunEvaluator == 0 ) {
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
	bool stationary = true;  // ture - stationary problem
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
		if (theGradGEvaluator == 0 ) {
		opserr << "Need theGradGEvaluator before an OutCrossingAnalysis can be created" << endln;
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
//									  theGFunEvaluator);
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
								 theGFunEvaluator,
								 theGradGEvaluator,
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
