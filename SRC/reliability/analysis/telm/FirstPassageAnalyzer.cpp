
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/FirstPassageAnalyzer.cpp,v $



#include <FirstPassageAnalyzer.h>
FirstPassageAnalyzer::FirstPassageAnalyzer
					 (ReliabilityDomain* passedReliabilityDomain,
					  FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
					  FunctionEvaluator* passedGFunEvaluator,
					  FOSeriesSimulation* passedFOSeriesSimulation,
					  int passedanalysisType,
					  bool passedtwoside)
{
	theReliabilityDomain = passedReliabilityDomain;
	theFindDesignPointAlgorithm = passedFindDesignPointAlgorithm;
	theGFunEvaluator = passedGFunEvaluator;
	theFOSeriesSimulation = passedFOSeriesSimulation;

	theRandomProcess = 0;
//	thecomponentResults=0;
//	theOutCrossingResult=0;
	theTimePoints=0;

	analysisType = passedanalysisType;
	numRV	 = theReliabilityDomain->getNumberOfRandomVariables();
	//numRVPos = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	numRVPos = 0;
	delta=theGFunEvaluator->getDt();
	numTimePoints=0;
	twoside=passedtwoside;

	timepoints=0;

	if(analysisType==2){
		if(passedFOSeriesSimulation==0){
			opserr<<"Warning in FirstPassageAnalyzer\n";
			opserr<<"FOSeriesSimulation object is required\n";
			opserr<<"when analysisType==2\n";
			opserr<<"default object is constructed\n";
			theFOSeriesSimulation=new FOSeriesSimulation();
			if(theFOSeriesSimulation==0){
				opserr<<" insufficient memory \n";
				opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				opserr<<" allocation of theFOSeriesSimulation\n";
			}
			theFOSeriesSimulation->setTwoSide(twoside);
		}
	}

}
FirstPassageAnalyzer::~FirstPassageAnalyzer()
{
	if(timepoints==0){delete timepoints; timepoints=0;}
//	if(thecomponentResults !=0){
//		for(int lsf=1;lsf<=numLsf;lsf++){
//			delete thecomponentResults[lsf-1]; thecomponentResults[lsf-1]=0;
//		}
//		delete [] thecomponentResults; thecomponentResults=0;
//	}
}
void
FirstPassageAnalyzer::setRandomProcess(RandomProcess* passedRandomProcess)
{
	theRandomProcess = passedRandomProcess;
	numTotalPulses = theRandomProcess->getNumTotalPulse();
	numOtherRVs=numRV-numTotalPulses;
	delta_Pulse = theRandomProcess->getDeltaPulse();

}
void FirstPassageAnalyzer::setComponentResult
						  (Vector* pudes, Vector* pxdes, Vector* palpha, 
						   Vector* phfunc, double pbeta,double ppf,
						   double panalysistime, int panalysisstep,
						   TimePoints* passedTimePoints,
						   int passedlsf, int passedidfragility,
						   bool stationary)
{
//	theOutCrossingResult=passedOutCrossingResult;

	udesres=pudes;
	xdesres=pxdes;
	alphadesres=palpha;
	hfuncdesres=phfunc;
	betadesres=pbeta;
	pfdesres=ppf;
	analysistime=panalysistime;
	analysisstep=panalysisstep;
	theTimePoints=passedTimePoints;

//	int IdLsf=theOutCrossingResult->getIdLsf();
//	int IdFragility=theOutCrossingResult->getIdFragility();

//	if(IdLsf!=passedlsf||IdFragility!=passedidfragility){
//		opserr<<"Fatal Error in FirstPassageAnalyzer::setComponentResult \n";
//		opserr<<"lsf and idfragility is not consistent with the current outcrossing result \n";
//		opserr<<" passedlsf "<<passedlsf<<" IdLsf" <<IdLsf<<"\n"; 
//		opserr<<" passedidfragility" <<passedidfragility<<" IdFragility " <<IdFragility<<"\n"; 
//		exit(-1);
//	}

//	int numPoints=theOutCrossingResult->getnumPoints();
	numTimePoints= theTimePoints->getNumTimePoints();
//	if(numPoints!=numTimePoints){
//		opserr<<"Fatal Error in FirstPassageAnalyzer::setComponentResult \n";
//		opserr<<" number of timepoints is not consistent with the current outcrossing result \n";
//		opserr<<" numPoints "<<numPoints<<" numTimePoints" <<numTimePoints<<"\n"; 
//		exit(-1);
//	}

	if(timepoints!=0){delete timepoints; timepoints=0;}
	timepoints = new Vector(numTimePoints);
	if(timepoints==0){
		opserr<<" insufficient memory \n";
		opserr<<" FirstPassageAnalyzer::setComponentResult\n";
		opserr<<" allocation of timepoints\n";
	}

//	if(!stationary){
//		(*timepoints)(0)=theOutCrossingResult->getTime(0);
//		for( int i=1;i<=numTimePoints;i++){
//			(*timepoints)(i)=theOutCrossingResult->getTime(i);
//			double time=theTimePoints->getAnalysisTime(i-1);
//			if((*timepoints)(i)!=time){
//				opserr<<"Fatal Error in FirstPassageAnalyzer::setComponentResult \n";
//				opserr<<"analysis time are not consistent\n";
//				opserr<<"(*timepoints)(i)"<<(*timepoints)(i);
//				opserr<<"time" <<time<<"\n"; 
//				exit(-1);
//			}
//		}
//	}else{
//		(*timepoints)(0)=0.0;
		for( int i=0;i<=numTimePoints;i++)
			(*timepoints)(i)=theTimePoints->getAnalysisTime(i);
//	}
}
//void FirstPassageAnalyzer::setSystemResults
//						  (OutCrossingResults* passedOutCrossingResult,
//						   TimePoints* passedTimePoints,
//						   int passednumlsf, int passedidfrag,
//						   bool stationary)
//{
//	if(thecomponentResults!=0){
//		for(int i=0; i<numLsf; i++){
//			delete thecomponentResults[i]; thecomponentResults[i]=0;
//		}
//		delete [] thecomponentResults;
//		thecomponentResults=0;
//	}
//	theOutCrossingResult=passedOutCrossingResult;
//	theTimePoints=passedTimePoints;
//	numLsf=passednumlsf;
//
//	thecomponentResults = new OutCrossingResults*[numLsf];
//	if(thecomponentResults==0){
//		opserr<<" insufficient memory \n";
//		opserr<<" FirstPassageAnalyzer::setSystemResults\n";
//		opserr<<" allocation of thecomponentResults\n";
//		exit(-1);
//	}
//	for(int lsf=1;lsf<=numLsf;lsf++){
//		theOutCrossingResult->readfromFile(lsf,passedidfrag);
//		thecomponentResults[lsf-1]= new OutCrossingResults(*theOutCrossingResult);
//		if(thecomponentResults[lsf-1]==0){
//		opserr<<" insufficient memory \n";
//		opserr<<" FirstPassageAnalyzer::setSystemResults\n";
//		opserr<<" allocation of thecomponentResults\n";
//		exit(-1);
//		}
//	}

//	int numPoints=theOutCrossingResult->getnumPoints();
//	numTimePoints= theTimePoints->getNumTimePoints();
//	if(numPoints!=numTimePoints){
//		opserr<<"Fatal Error in FirstPassageAnalyzer::setComponentResult \n";
//		opserr<<" number of timepoints is not consistent with the current outcrossing result \n";
//		opserr<<" numPoints "<<numPoints<<" numTimePoints" <<numTimePoints<<"\n"; 
//		exit(-1);
//	}
//
//	if(timepoints!=0){delete timepoints; timepoints=0;}
//	timepoints = new Vector(numTimePoints+1);
//	if(timepoints==0){
//		opserr<<" insufficient memory \n";
//		opserr<<" FirstPassageAnalyzer::setComponentResult\n";
//		opserr<<" allocation of timepoints\n";
//	}
//	(*timepoints)(0)=0.0;
//	for( int i=1;i<=numPoints;i++)
//		(*timepoints)(i)=theTimePoints->getAnalysisTime(i-1);
//}
