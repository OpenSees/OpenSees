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
// $Date: 2008-08-26 17:34:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/SelectLoadInitialStaticAnalysis.cpp,v $


#include <SelectLoadInitialStaticAnalysis.h>
//#include <RandomVariablePositionerIter.h>

SelectLoadInitialStaticAnalysis::SelectLoadInitialStaticAnalysis
						  (ReliabilityDomain* passedReliabilityDomain,
						   Domain* passedStructuralDomain,
						   int passednstep,
						   int passednumLoadPatterns,
						   int* passedStaticLoadPatterns,
						   bool passedprint)
:InitialStaticAnalysis(passedReliabilityDomain,passedStructuralDomain,
passedprint)
{
	StaticLoadPatterns=NULL;
	theOrgPatternIter=NULL;
	theOrgPatterns=NULL;

	theReliabilityStaticAnalysis=NULL;
	theAnalysisModel=NULL;
	theAlgorithm=NULL;
	theHandler=NULL;
	theNumberer=NULL;
	theSOE=NULL;
	theStaticIntegrator=NULL;
	theTest=NULL;
	theSensitivityAlgorithm=NULL;

	Nstep = passednstep;
	NumLoadPatterns = passednumLoadPatterns;
	if( NumLoadPatterns == 0 ){
		StaticLoadPatterns = NULL;
	}else{
		StaticLoadPatterns = new int[NumLoadPatterns];
		for( int i=0; i<NumLoadPatterns; i++){
			StaticLoadPatterns[i] = passedStaticLoadPatterns[i];
		}
	}
	saveLoads();
	createStaticAnalysis();
	modified=false;

	if(print){
		output << "=====SelectLoadStaticAnalysis=====\n";
		output << "\n";
		output << " Number of Selected Loads" << NumLoadPatterns << "\n";
		output << "\n";
		output << " Load Pattern IDs \n";
		for( int i=0; i<NumLoadPatterns; i++)
			output << StaticLoadPatterns[i] << "\n";
	}
	output.flush();

}
SelectLoadInitialStaticAnalysis::~SelectLoadInitialStaticAnalysis()
{
	if(theReliabilityStaticAnalysis!=NULL)
	{delete theReliabilityStaticAnalysis;theReliabilityStaticAnalysis=0;}
	if(theAlgorithm!=NULL)
	{delete theAlgorithm;theAlgorithm=0;}
	if(theAnalysisModel!=NULL)
	{delete theAnalysisModel;theAnalysisModel=0;}
	if(theTest!=NULL)
	{delete theTest;theTest=0;}
	if(theHandler!=NULL)
	{delete theHandler;theHandler=0;}
	if(theNumberer!=NULL)
	{delete theNumberer;theNumberer=0;}
	if(theStaticIntegrator!=NULL)
	{delete theStaticIntegrator;theStaticIntegrator=0;}
	if(theSOE!=NULL)
	{delete theSOE;theSOE=0;}
	if(theSensitivityAlgorithm!=NULL) 
	{delete theSensitivityAlgorithm;theSensitivityAlgorithm=0;}
	if( theOrgPatterns != NULL ){
		delete theOrgPatterns;
		theOrgPatterns = NULL;
	}
	if( StaticLoadPatterns != NULL ){
		delete [] StaticLoadPatterns;
		StaticLoadPatterns=NULL;
	}
}
void SelectLoadInitialStaticAnalysis::createStaticAnalysis(void)
{
	if(theAnalysisModel!=NULL){delete theAnalysisModel;theAnalysisModel=NULL;}
	theAnalysisModel = new AnalysisModel();
	if(theAnalysisModel == NULL){
		opserr << "Fail to generate theAnalysisModel\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
	if(theTest!=0){delete theTest; theTest=NULL;}
//	theTest = new CTestNormUnbalance(1.0e-12,25,1);      
    theTest = new CTestNormDispIncr(1.0e-8,25,0);             
	if(theTest == NULL) {
		opserr << "Fail to generate theTest\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
    if(theAlgorithm!=NULL){delete theAlgorithm; theAlgorithm=NULL;}
    theAlgorithm = new NewtonRaphson(*theTest); 
	if(theAlgorithm == NULL) {
		opserr << "Fail to generate theAlgorithm\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
	if(theHandler!=NULL){delete theHandler; theHandler=NULL;}
    theHandler = new PlainHandler();       
	if(theHandler == NULL) {
		opserr << "Fail to generate theHandler\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
    RCM *theRCM = new RCM();	
	if(theRCM == NULL) {
		opserr << "Fail to generate theRCM\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
	if(theNumberer!=NULL){delete theNumberer; theNumberer=NULL;}
    theNumberer = new DOF_Numberer(*theRCM);    	
	if(theNumberer == NULL) {
		opserr << "Fail to generate theNumberer\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
	if(theStaticIntegrator!=NULL){delete theStaticIntegrator; theStaticIntegrator=NULL;}
	theStaticIntegrator = new LoadControl(1, 1, 1, 1);     
	if(theStaticIntegrator == NULL) {
  		opserr << "Fail to generate theStaticIntegrator\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}

    BandGenLinSolver *theSolver;
	theSolver = new BandGenLinLapackSolver();
//	ProfileSPDLinSolver *theSolver;
//  theSolver = new ProfileSPDLinDirectSolver(); 	
	if(theSolver == NULL) {
		opserr << "Fail to generate theSolver\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
	if(theSOE!=NULL){delete theSOE; theSOE=NULL;}
    theSOE = new BandGenLinSOE(*theSolver);      
//    theSOE = new ProfileSPDLinSOE(*theSolver);      
	if(theSOE == NULL) {
		opserr << "Fail to generate theSOE\n";
		opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
		exit(-1);
	}
	theReliabilityStaticAnalysis = new ReliabilityStaticAnalysis(*theDomain,
					       *theHandler,
					       *theNumberer,
					       *theAnalysisModel,
					       *theAlgorithm,
					       *theSOE,
					       *theStaticIntegrator);
	if(theReliabilityStaticAnalysis == NULL) {
	opserr << "Fail to generate theStaticAnalysis\n";
	opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
	exit(-1);
	}
	int analysisTypeTag=1;	
	if(theSensitivityAlgorithm!=NULL){delete theSensitivityAlgorithm; theSensitivityAlgorithm=NULL;}

//	theSensitivityAlgorithm = new   // Discuss with Prof Scott
//	    SensitivityAlgorithm(theDomain,
//				 theAlgorithm,
//				 theStaticIntegrator,
//				 analysisTypeTag);
	if(theSensitivityAlgorithm == NULL) {
	opserr << "Fail to generate theSensitivityAlgorithm\n";
	opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
	exit(-1);
	}
	int res = theReliabilityStaticAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
	if(res != 0) {
	opserr << "Fail to set the sensitivity algorithm to static analysis\n";
	opserr << "in SelectLoadStaticAnalysis::createStaticAnalysis \n";
	exit(-1);
	}
	activeSensitivity=true;
}
void SelectLoadInitialStaticAnalysis::analyze(Vector x)
{
	if(NumLoadPatterns !=0 )modifyLoads();
	this->reset();
	if(print){
		output <<"\n";
		output <<" function SelectLoadStaticAnalysis::analyze\n";
		output <<"\n";
		output.flush();
	}
//
	/*
	RandomVariablePositioner *theRandomVariablePositioner;
	int rvNumber;
	int i;
	int numberOfRandomVariablePositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	if(print){
		output << "\n";
		output << "SelectLoadStaticAnalysis::analyze\n";
		output << "\n";
		output << " Number of Random Variable Positioner ";
		output << numberOfRandomVariablePositioners << "\n";
		output << "\n";
		output.flush();
	}
	if(print){
		output << "\n";
		output << " Set Randomvariable to Vector X \n";
		output << "\n";
		output.flush();
	}
	for ( i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
		theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		rvNumber = theRandomVariablePositioner->getRvIndex();
		theRandomVariablePositioner->update(x(rvNumber));
	}
	if(print){
		RandomVariablePositionerIter rvpIter = theReliabilityDomain->getRandomVariablePositioners();
		RandomVariablePositioner *theRVP;
		while ((theRVP = rvpIter()) != 0) {
		//for ( i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
			//theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
			//rvNumber				= theRandomVariablePositioner->getRvIndex();
			RandomVariable *theRV	= theRVP->getRandomVariable();
			theRVP->update(theRV->getMean());
			output << " RandomVariable ";
			output << i << "value=";
			output << x(theRVP->getRvIndex()) <<"\n";
		}
		output.flush();
	}
	*/

	int result;
	double dresult = 0;
	result=theReliabilityStaticAnalysis->analyze(Nstep);
	if( result != 0 ){
		opserr << "Error in SelectLoadStaticAnalysis::analyzen";
		opserr << "return code" << (int)result << "\n";
		exit(-1);
	}
	if(print)this->printResult();
}
void SelectLoadInitialStaticAnalysis::analyzeMean()
{
	if(NumLoadPatterns !=0 )modifyLoads();
	this->reset();
	/*
	RandomVariablePositioner *theRandomVariablePositioner;
	int rvNumber;
	int i;
	int numberOfRandomVariablePositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	if(print){
		output << "\n";
		output << "SelectLoadStaticAnalysis::analyze\n";
		output << "\n";
		output << " Number of Random Variable Positioner ";
		output << numberOfRandomVariablePositioners << "\n";
		output << "\n";
		output.flush();
	}
	  
	RandomVariablePositionerIter rvpIter = theReliabilityDomain->getRandomVariablePositioners();
	RandomVariablePositioner *theRVP;
	while ((theRVP = rvpIter()) != 0) {
	//for ( i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
		//theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		//rvNumber				= theRandomVariablePositioner->getRvIndex();
		RandomVariable *theRV	= theRVP->getRandomVariable();
		theRandomVariablePositioner->update(theRV->getMean());
	}
	
	if(print){
		output << "\n";
		output << " Set Randomvariable to their means \n";
		output << "\n";
		output.flush();
		RandomVariablePositionerIter rvpIter = theReliabilityDomain->getRandomVariablePositioners();
		RandomVariablePositioner *theRVP;
		while ((theRVP = rvpIter()) != 0) {
		//for ( i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
		//theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		//rvNumber				= theRandomVariablePositioner->getRvIndex();
			RandomVariable *theRV	= theRVP->getRandomVariable();
			theRandomVariablePositioner->update(theRV->getMean());
			double mean=theRV->getMean();
			output << " Pos" <<theRVP->getTag();
			output << " Rv" << theRV->getTag();
			output << " mean" << mean;
		}
	}
	*/

	int result;
	double dresult = 0;
	result=theReliabilityStaticAnalysis->analyze(Nstep);
	if( result != 0 ){
		opserr << "Error in SelectLoadStaticAnalysis::analyzen";
		opserr << "return code" << (int)result << "\n";
		exit(-1);
	}
//	theDomain->setLoadConstant();
//	if(NumLoadPatterns !=0 )recoverLoads();
//    theDomain->setCurrentTime(0.0);
//   theDomain->setCommittedTime(0.0);
	if(print)this->printResult();
}
void SelectLoadInitialStaticAnalysis::reset()
{
	if(print) {
		output<<"\n";
		output<<"function SelectLoadStaticAnalysis::reset \n";
		output<<"\n";
		output<<" theDomain->revertToStart\n";
		output<<"\n";
		output.flush();
	}
	theDomain->revertToStart();
	theDomain->unsetLoadConstant();

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
}
void SelectLoadInitialStaticAnalysis::saveLoads()
{
	int numload=0;
	LoadPattern* thePattern;
	if(theOrgPatterns!=NULL){ delete theOrgPatterns; theOrgPatterns=NULL;} 
	theOrgPatterns = new ArrayOfTaggedObjects(32);
	if( theOrgPatterns == NULL ){
		opserr << "OutCrossingAnalysis::saveLoads - out of memory\n";
		opserr << "for theLoadPatternsOrg\n";
		exit(-1);
	}
	LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
	while((thePattern = thePatterns()) != NULL){
		numload++;
		bool result = theOrgPatterns->addComponent(thePattern);
		if(!result) {
			opserr << "OutCrossingAnalysis::saveLoads - out of memory\n";
			opserr << "for copy\n";
			exit(-1);
		}
	}
	if( numload == 0){
		opserr << " FATAL error !! \n";
		opserr << " No loadpattern is defined before InitialStaticAnalysis\n";
		exit(-1);
	}
	if(theOrgPatternIter!=NULL){delete theOrgPatternIter; theOrgPatternIter=NULL;}
	theOrgPatternIter = new LoadPatternIter(theOrgPatterns);
}
void SelectLoadInitialStaticAnalysis::modifyLoads()
{
	modified=true;
	LoadPattern* thePattern;
	LoadPattern* thePat;
	if(NumLoadPatterns != 0){
		LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
		while((thePattern = thePatterns()) != NULL){
			int tag=thePattern->getTag();
			thePat=theDomain->removeLoadPattern(tag);
		}
		theOrgPatternIter->reset();
		while((thePattern = (*theOrgPatternIter)()) != NULL){
			int tag=thePattern->getTag();
			bool found=false;
			for(int i=0; i<NumLoadPatterns; i++){
				if(tag == StaticLoadPatterns[i]){
					found=true;
					break;
				}
			}
			if(found) theDomain->addLoadPattern(thePattern);
		}
	}
	else{
		opserr<< "!!!!! WARNING !!!!\n";
		opserr<< "No Load is selected in the initialstatic analysis\n";
		opserr<< "all Loads are applied to the model\n";
	}

	if(print){
		output << "\n";
		output << " after modify load  in theDomain \n";
		output << "\n";
		LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
		while((thePattern = thePatterns()) != NULL){
			int tag=thePattern->getTag();
			output << " load pattern " << tag <<"\n";
		}
		output << "\n";
		output << " after modify load in theOrgPatterns \n";
		output << "\n";
		theOrgPatternIter->reset();
		while((thePattern = (*theOrgPatternIter)()) != NULL){
			int tag=thePattern->getTag();
			output << " load pattern " << tag <<"\n";
		}
	}
}	
void SelectLoadInitialStaticAnalysis::recoverLoads()
{
	if(!modified){
		opserr << "===== FATAL Error =====\n"; 
		opserr << "It is attempted to recover loads\n"; 
		opserr << "although they are not modified \n"; 
		exit(-1);
	}
	modified=false;
	LoadPattern* thePattern;
	LoadPattern* thePat;
	LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
	while((thePattern = thePatterns()) != NULL){
		int tag=thePattern->getTag();
		thePat=theDomain->removeLoadPattern(tag);
	}
	theOrgPatternIter->reset();
	while((thePattern = (*theOrgPatternIter)()) != NULL){
		theDomain->addLoadPattern(thePattern);
	}
	if(print){
		output << "\n";
		output << " after recover load  in theDomain \n";
		output << "\n";
		LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
		while((thePattern = thePatterns()) != 0){
			int tag=thePattern->getTag();
			output << " load pattern " << tag <<"\n";
		}
	}
}
void SelectLoadInitialStaticAnalysis::activateSensitivity(void)
{
	if(print){
		output << "\n";
		output << " SelectLoadStaticAnalysis::activateSensitivity \n";
		output << "\n";
	}
//	SensitivityAlgorithm* currentSensitivityAlgorithm
	Integrator* currentSensitivityAlgorithm

		=theReliabilityStaticAnalysis->getSensitivityAlgorithm();
	if(!activeSensitivity){
//		currently inactive	
		if(currentSensitivityAlgorithm != NULL){
			opserr << " SelectLoadStaticAnalysis::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			exit(-1);
		}
		theReliabilityStaticAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
		activeSensitivity=true;
	}else{
		if(currentSensitivityAlgorithm == NULL){
			opserr << " SelectLoadStaticAnalysis::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			opserr << " must not be the zero pointer \n";
			exit(-1);
		}
	}
}
void SelectLoadInitialStaticAnalysis::inactivateSensitivity(void)
{
	if(print){
		output << "\n";
		output << " SelectLoadStaticAnalysis::inactivateSensitivity \n";
		output << "\n";
	}
//	SensitivityAlgorithm* currentSensitivityAlgorithm
	Integrator* currentSensitivityAlgorithm

	   =theReliabilityStaticAnalysis->getSensitivityAlgorithm();
	if(activeSensitivity){
		if(currentSensitivityAlgorithm == NULL){
			opserr << " SelectLoadStaticAnalysis::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			opserr << " continue analysis\n";
//			exit(-1);
		}
		activeSensitivity=false;
	//	SensitivityAlgorithm* zeroSensitivityAlgrithm=0;
		Integrator* zeroSensitivityAlgrithm=0;

		theReliabilityStaticAnalysis->setSensitivityAlgorithm(zeroSensitivityAlgrithm);
	}else{
		if(currentSensitivityAlgorithm != NULL){
			opserr << " SelectLoadStaticAnalysis::activateSensitivity \n";
			opserr << " Inconsistency of activeSensitivity \n";
			opserr << " continue analysis \n";
//			exit(-1);
		}
	}
}
void SelectLoadInitialStaticAnalysis::constLoads(double time=0.0)
{
	theDomain->setLoadConstant();
    theDomain->setCurrentTime(time);
    theDomain->setCommittedTime(time);
}
void SelectLoadInitialStaticAnalysis::resetconstLoads(double time=0.0)
{
	theDomain->unsetLoadConstant();
    theDomain->setCurrentTime(time);
    theDomain->setCommittedTime(time);
}
void SelectLoadInitialStaticAnalysis::constandrecoverLoads(double time=0.0)
{
	theDomain->setLoadConstant();
    theDomain->setCurrentTime(time);
    theDomain->setCommittedTime(time);
	if(!modified){
		opserr << "===== FATAL Error =====\n"; 
		opserr << "It is attempted to recover loads\n"; 
		opserr << "although they are not modified \n"; 
		exit(-1);
	}
	modified=false;
	LoadPattern* thePatternOrg;
	LoadPattern* thePattern;
	theOrgPatternIter->reset();
	LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
	while((thePatternOrg = (*theOrgPatternIter)()) != NULL){
		bool found=false;
		int id=thePatternOrg->getTag();
		thePatterns.reset();
		while((thePattern = thePatterns()) != NULL){
			int tag=thePattern->getTag();
			if(id==tag){
				found=true;
				break;
			}
		}
		if(!found){
			theDomain->addLoadPattern(thePatternOrg);
		}
	}
	if(print){
		output << "\n";
		output << " after recover load  in theDomain \n";
		output << "\n";
		LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
		while((thePattern = thePatterns()) != 0){
			int tag=thePattern->getTag();
			output << " load pattern " << tag <<"\n";
		}
	}
}
