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

// $Revision: 1.2 $
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/StaticAnalyzer.cpp,v $

#include <StaticAnalyzer.h>
#include<Integrator.h>
//#include <RandomVariablePositionerIter.h>

StaticAnalyzer::StaticAnalyzer
				(ReliabilityDomain* passedReliabilityDomain,
				 Domain* passedStructuralDomain,
				 InitialStaticAnalysis* passedInitialStaticAnalysis,
				 ReliabilityStaticAnalysis* passedStaticAnalysis,
				// SensitivityAlgorithm* passedSensitivityAlgorithm,
			         Integrator* passedSensitivityAlgorithm,
				 //SensitivityIntegrator *passedSensitivityIntegrator,
				   Integrator *passedSensitivityIntegrator,//Abbas 
				 int passednstep,
				 double passeddelta,
				 int passednumLoadPatterns,
				 int* passedStaticLoadPatterns,
				 bool passedprint)
:Analyzer(passedReliabilityDomain,passedStructuralDomain,
		  passedInitialStaticAnalysis, passedSensitivityAlgorithm,
		  passedSensitivityIntegrator, passednstep,
		  passeddelta,passednumLoadPatterns,passedStaticLoadPatterns,passedprint)
{
    theStaticAnalysis = passedStaticAnalysis;
//	SensitivityAlgorithm* currentSensitivityAlgorithm
	Integrator* currentSensitivityAlgorithm

   	   =theStaticAnalysis->getSensitivityAlgorithm();
	if(currentSensitivityAlgorithm!=0){
		activeSensitivity=false;
	}else{
		activeSensitivity=true;
	}
}
StaticAnalyzer::~StaticAnalyzer()
{
}
int StaticAnalyzer::analyze(Vector x)
{
	if(print) {
		output<<"\n";
		output<<"function StaticAnalyzer::analyzeX\n";
		output<<"\n";
		output.flush();
	}
	if(theInitialStaticAnalysis!=0){
		// initial static analysis including initialization  
		theInitialStaticAnalysis->analyze(x);
		theInitialStaticAnalysis->constandrecoverLoads(0.0);
//		theInitialStaticAnalysis->constLoads(0.0);
//		theInitialStaticAnalysis->recoverLoads();
		if(print) {
		output<<"\n";
		output<<" InitialStaticAnalysisComplete\n";
		output<<"\n";
		output.flush();
		}
	}
	else{
		if(print) {
		output<<"\n";
		output<<"function StaticAnalyzer::analyzeX\n";
		output<<"\n";
		output<<" without initialstaticanalysis\n";
		output<<"\n";
		output.flush();
		}
		theDomain->revertToStart();
	}
	if (activeSensitivity){
//		theSensitivityAlgorithm->revertToStart();
//		int NumGrads = theReliabilityDomain->getNumberOfRandomVariables();
//		theSensitivityIntegrator->revertToStart(NumGrads);
//		if(print) {
//			output<<"\n";
//			output<<" theSensitivityAlgorithm->revertToStart()\n";
//			output<<" NumGrads" << NumGrads <<"\n";
//			output<<" theSensitivityIntegrator->revertToStart()\n";
//			output.flush();
//		}
	}
// setRandomVariables
	/*
	RandomVariablePositioner *theRandomVariablePositioner;
	int rvNumber;
	int i;
	int numberOfRandomVariablePositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	if(print){
		output << "\n";
		output << " Number of Random Variable Positioner ";
		output << numberOfRandomVariablePositioners << "\n";
		output << "\n";
		output.flush();
		output << "\n";
		output << " Set Randomvariable to Vector X \n";
		output << "\n";
		output.flush();
	}
	RandomVariablePositionerIter rvpIter = theReliabilityDomain->getRandomVariablePositioners();
	RandomVariablePositioner *theRVP;
	while ((theRVP = rvpIter()) != 0) {
	//for ( i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
		//theRVP  = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		rvNumber				= theRVP->getRvIndex();
		theRVP->update(x(rvNumber));
		if(print){
			output << " RandomVariable ";
			output << rvNumber << "value=";
			output << x(rvNumber) <<"\n";
		}
		output.flush();
	}
	*/

	if(print){
		output << "\n";
		output << " Static Response Analysis \n";
		output.flush();
	}
	double result=theStaticAnalysis->analyze(Numstep);
//	if( (int)result != 0 ){
//		opserr << "Error in Static Analyzer::analyzen";
//		opserr << "return code" << (int)result << "\n";
//		exit(-1);
//	}
	if(print)this->printresults();
	return (int)result;
}
int StaticAnalyzer::analyzeMean()
{
	if(print) {
		output<<"\n";
		output<<"function StaticAnalyzer::analyzeMean\n";
		output<<"\n";
		output.flush();
	}
	if(theInitialStaticAnalysis!=0){
		// initial static analysis including initialization  
		theInitialStaticAnalysis->analyzeMean();
		theInitialStaticAnalysis->constandrecoverLoads(0.0);
//		theInitialStaticAnalysis->constLoads(0.0);
//		theInitialStaticAnalysis->recoverLoads();
		if(print) {
		output<<"\n";
		output<<" InitialStaticAnalysisComplete\n";
		output<<"\n";
		output.flush();
		}
	}
	else{
		if(print) {
		output<<"\n";
		output<<"function StaticAnalyzer::analyzeMean\n";
		output<<"\n";
		output<<" without initialstaticanalysis\n";
		output<<"\n";
		output.flush();
		}
		theDomain->revertToStart();
	}
	if (activeSensitivity){
//		if(print) {
//		output<<"\n";
//		output<<" Setting For Sensitivity Algorithm";
//		output<<"\n";
//		output.flush();
//		}
//		theSensitivityAlgorithm->revertToStart();
//		int NumGrads = theReliabilityDomain->getNumberOfRandomVariables();
//		theSensitivityIntegrator->revertToStart(NumGrads);
//		if(print) {
//			output<<"\n";
//			output<<" theSensitivityAlgorithm->revertToStart()\n";
//			output<<" NumGrads" << NumGrads <<"\n";
//			output<<" theSensitivityIntegrator->revertToStart()\n";
//			output.flush();
//		}
	}
	/*
	RandomVariablePositioner *theRandomVariablePositioner;
	int rvNumber;
	int i;
	int numberOfRandomVariablePositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	if(print){
		output << "\n";
		output << " Number of Random Variable Positioner ";
		output << numberOfRandomVariablePositioners << "\n";
		output << "\n";
		output.flush();
		output << "\n";
		output << " Set Randomvariable to Vector X \n";
		output << "\n";
		output.flush();
	}
	RandomVariablePositionerIter rvpIter = theReliabilityDomain->getRandomVariablePositioners();
	//for ( i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
	while ((theRandomVariablePositioner = rvpIter()) != 0) {
		//theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		rvNumber				= theRandomVariablePositioner->getRvIndex();
		RandomVariable *theRV	= theRandomVariablePositioner->getRandomVariable();
		theRandomVariablePositioner->update(theRV->getMean());
		if(print){
			double mean=theRV->getMean();
			output << " Pos" <<rvNumber;
			output << " Rv" << rvNumber;
			output << " mean" << mean;
			output.flush();
		}
	}
	*/
	if(print){
		output << "\n";
		output << " Static Response Analysis \n";
		output.flush();
	}
    theDomain->setCurrentTime(0.0);
    theDomain->setCommittedTime(0.0);

	double result=theStaticAnalysis->analyze(Numstep);
//	if( (int)result != 0 ){
//		opserr << "Error in Static Analyzer::analyzen";
//		opserr << "return code" << (int)result << "\n";
//		exit(-1);
//	}
	if(print)this->printresults();
	return (int)result;
}

void StaticAnalyzer::activateSensitivty()
{
	if(theSensitivityAlgorithm==0){
		opserr << "FatatlError \n";
		opserr << "StaticAnalyzer::activateSensitivity\n";
		opserr << "theSensitivityAlgorithm is not defined\n";
		exit(-1);
	}

	if(theInitialStaticAnalysis!=0){
		theInitialStaticAnalysis->activateSensitivity();
	}
	if(print){
		output << "\n";
		output << " StaticAnalyzer::activateSensitivity \n";
		output << "\n";
	}
	//SensitivityAlgorithm* currentSensitivityAlgorithm

	Integrator* currentSensitivityAlgorithm

	   =theStaticAnalysis->getSensitivityAlgorithm();
	if(!activeSensitivity){
//		currently inactive	
		if(currentSensitivityAlgorithm != 0){
			opserr << " SelectLoadStaticAnalysis::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			exit(-1);
		}
		theStaticAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
		activeSensitivity=true;
	}else{
		if(currentSensitivityAlgorithm == 0){
			opserr << " StaticAnalyzer::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			opserr << " must not be the zero pointer \n";
			exit(-1);
		}
	}
}
void StaticAnalyzer::inactivateSensitivty()
{
	if(theInitialStaticAnalysis!=0){
		theInitialStaticAnalysis->inactivateSensitivity();
	}
	if(print){
		output << "\n";
		output << " StaticAnalyzer::inactivateSensitivity \n";
		output << "\n";
	}
//	SensitivityAlgorithm* currentSensitivityAlgorithm
       	Integrator* currentSensitivityAlgorithm

	   =theStaticAnalysis->getSensitivityAlgorithm();
	if(activeSensitivity){
		if(currentSensitivityAlgorithm == 0){
			opserr << " SelectLoadStaticAnalysis::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			opserr << " must not be the zero pointer \n";
			exit(-1);
		}
		activeSensitivity=false;
	//	SensitivityAlgorithm* zeroSensitivityAlgrithm=0;
		Integrator* zeroSensitivityAlgrithm=0;

		theStaticAnalysis->setSensitivityAlgorithm(zeroSensitivityAlgrithm);
	}else{
		if(currentSensitivityAlgorithm != 0){
			opserr << " StaticAnalyzer::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			opserr << " must be the zero pointer \n";
			exit(-1);
		}
	}
}
void StaticAnalyzer::setGFunEachStepEvaluator(GFunEachStepEvaluator *pGFunEachStepEvaluator)
{
	opserr << "ERROR StaticAnalyzer::setGFunEachStepEvaluator\n";
	opserr << "This is not implemented yet\n";
	exit(-1);
}
void StaticAnalyzer::inactivateGFunEachStepEvaluator()
{
	opserr << "ERROR StaticAnalyzer::inactivatesetGFunEachStepEvaluator\n";
	opserr << "This is not implemented yet\n";
	exit(-1);
}
Matrix* StaticAnalyzer::getEachStepResult()
{
	opserr << "ERROR StaticAnalyzer::getEachStepResult\n";
	opserr << "This is not implemented yet\n";
	exit(-1);
}
Matrix* StaticAnalyzer::getEachStepConvFlag()
{
	opserr << "ERROR StaticAnalyzer::getEachStepConvFlag\n";
	opserr << "This is not implemented yet\n";
	exit(-1);
}



