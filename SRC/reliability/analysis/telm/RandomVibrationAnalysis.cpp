
// $Revision: 1.2 $
// $Date: 2008-05-13 16:30:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/RandomVibrationAnalysis.cpp,v $

#include <RandomVibrationAnalysis.h>
RandomVibrationAnalysis::RandomVibrationAnalysis
						(ReliabilityDomain* passedReliabilityDomain,
						 FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
						 //NewSearchWithStepSizeAndStepDirection* passedFindDesignPointAlgorithm,
						 Domain* passedDomain, 
						 InitialPointBuilder* passedInitialPointBuilder,
						 CrossingRateAnalyzer* passedCrossingRateAnalyzer,
						 FirstPassageAnalyzer* passedFirstPassageAnalyzer,
						 FunctionEvaluator* passedGFunEvaluator,
   						 GradientEvaluator* passedGradGEvaluator,
						 ReliabilityConvergenceCheck* passedReliabilityConvergenceCheck,
						 double passedStartTime,
						 double passedEndTime,
						 double passedTimeInterval,
						 double passedStartAnalysis,
						 double passedFragMin,
						 double passedFragInt,
						 int passednFrag,
						 int passeddesignpoint,
						 bool passedstationary,
						 bool passedmirrorimage,
						 bool passedinitialpoint,
						 bool passedfirstpassage,
						 TCL_Char *passedFileName,
						 char* passedFileBinary,
						 Tcl_Interp* passedTclInterp,
	 					 bool passedprint)
:designpoint(passeddesignpoint), stationary(passedstationary),
 mirrorimage(passedmirrorimage), initialpoint(passedinitialpoint),
 firstpassage(passedfirstpassage), print(passedprint)
{
	theReliabilityDomain = passedReliabilityDomain;
	theFindDesignPointAlgorithm = passedFindDesignPointAlgorithm;
	theDomain = passedDomain;
	theInitialPointBuilder = passedInitialPointBuilder;
	theCrossingRateAnalyzer = passedCrossingRateAnalyzer;
	theFirstPassageAnalyzer = passedFirstPassageAnalyzer;
	theGFunEvaluator = passedGFunEvaluator;
	theGradGEvaluator = passedGradGEvaluator;
	theReliabilityConvergenceCheck = passedReliabilityConvergenceCheck;
	theTclInterp = passedTclInterp;
 
	delta= theGFunEvaluator->getDt();  // time step for the dynamic response analysis
	fileName = new char[256];		   // output file name	
	strcpy(fileName,passedFileName);
	if(passedFileBinary!=0){
		fileBinary = new char[256];		   // output file name	
		strcpy(fileBinary,passedFileBinary);
	}else{
		fileBinary=0;
	}

	numRV	 = theReliabilityDomain->getNumberOfRandomVariables();
	//numRVPos = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	numRVPos = 0;
	numLsf	 = theReliabilityDomain->getNumberOfLimitStateFunctions();

	Scorg=theReliabilityConvergenceCheck->getScaleValue();
	Scfixorg=theReliabilityConvergenceCheck->getScfix();
	////////////////////////////////
	///// find random process /////
	////////////////////////////////
	theRandomProcess = new RandomProcess(theDomain,theReliabilityDomain,print);
	NumTotalPulse = theRandomProcess->getNumTotalPulse();
	delta_Pulse = theRandomProcess->getDeltaPulse();

	if(theCrossingRateAnalyzer!=0)
		theCrossingRateAnalyzer->setRandomProcess(theRandomProcess);
		theCrossingRateAnalyzer->setdelta(delta);
	if(theFirstPassageAnalyzer!=0)
		theFirstPassageAnalyzer->setRandomProcess(theRandomProcess);
	if(initialpoint){
		theInitialPointBuilder->setRandomProcess(theRandomProcess);
		theInitialPointBuilder->setDt(delta);
//		theInitialPointBuilder->setDpulse(delta_Pulse);
	}

	if(passedStartAnalysis<0.0){
		StartSteps = 0;
	}else{
		StartSteps = (int)((passedStartAnalysis+0.01*delta)/delta);
	}
	////////////////////////////////////////////////////
	///////	TimePoints for analysis results ////////////
    ////////////////////////////////////////////////////
	theTimePoints = new TimePoints(passedStartTime,
								   passedEndTime,
								   passedTimeInterval,
								   delta,print);

	NumTimePoints = theTimePoints->getNumTimePoints();

	numFragility=passednFrag;
	Fragility=new Vector(numFragility);
	for(int i=0; i<numFragility; i++){
		(*Fragility)(i)=passedFragMin+passedFragInt*(double)i;
	}
	
	//// save results /////
	numDesPoints=numFragility*numLsf;
	//if(!stationary)numDesPoints=NumTimePoints;
	theOutCrossingResults = new	OutCrossingResults(numLsf,numFragility,numRV, 
												   fileBinary,print);
	
	xDesTmp= new Vector(numRV);   
	uDesTmp= new Vector(numRV);
	alphaTmp= new Vector(numRV);
	xinitial= new Vector(numRV);
	hfunc= new Vector(numRV);

	baseExcitation = 0;
	baseExcitation0 = 0;
	mirExcitation = 0;
	thePerformFuncCoeffs=0;
	thePfCoeffIter=0;

	FirstPassageMatrix=0;
	if(firstpassage){
		FirstPassageMatrix = new Matrix*[numLsf+1];
		if( FirstPassageMatrix == 0){ opserr << "FirstPassageMatrix\n";
									  opserr << "RandomVibrationAnalysis\n";
									  exit(1);}
		for( int i=0; i<=numLsf; i++){
			FirstPassageMatrix[i]=new Matrix(NumTimePoints, numFragility);
			if( FirstPassageMatrix[i] == 0)	{ opserr << "FirstPassageMatrix[i]\n";
											  opserr << "RandomVibrationAnalysis\n";
											  exit(1);}
		}
	}

	if(print){
		output.open("OutCrossingBak.txt", ios::out);
		output << "\n";
		output << "OutCrossing::OutCrossing\n";
		output << "\n";
		output << "delta......................" << delta << "\n"; 
		output << "designpoint................" << designpoint << "\n"; 
		output << "mirrorimage................" << mirrorimage << "\n"; 
		output << "initialpoint..............." << initialpoint << "\n"; 
		output << "stationary................." << stationary << "\n"; 
		output << "firstpassage..............." << firstpassage << "\n"; 
		output << "Scorg......................" << Scorg << "\n"; 
		output << "numRV " << numRV << "\n"; 
		output << "numRVPos " << numRVPos << "\n"; 
		output << "numLsf " << numLsf << "\n"; 
		output << "NumTotalPulse " << NumTotalPulse<< "\n"; 
		output << "delta_Pulse" << delta_Pulse<< "\n"; 
		output << "StartSteps" << StartSteps<< "\n"; 
		output << "NumTimePoints" << NumTimePoints<< "\n"; 
		output << "numFragility" << numFragility<< "\n"; 
		output.flush();
	}
	outputTELs.open("TELs.txt", ios::out);
}
RandomVibrationAnalysis::~RandomVibrationAnalysis()
{
	if(fileName != 0){
		delete [] fileName;
		fileName = 0;
	}

	if(theRandomProcess != 0){
		delete theRandomProcess;
		theRandomProcess = 0;
	}

	if(theTimePoints !=0){
		delete theTimePoints;
		theTimePoints = 0;
	}

	if(Fragility !=0){
		delete Fragility;
		Fragility= 0;
	}

	if(theOutCrossingResults !=0){
		delete theOutCrossingResults;
		theOutCrossingResults = 0;
	}

	if(xDesTmp != 0 ){
		delete xDesTmp;
		xDesTmp = 0;
	}

	if(alphaTmp != 0 ){
		delete alphaTmp;
		alphaTmp = 0;
	}
	if(uDesTmp != 0 ){
		delete uDesTmp;
		uDesTmp= 0;
	}
	if(xinitial != 0 ){
		delete xinitial;
		xinitial= 0;
	}
	if(baseExcitation != 0 ){
		delete baseExcitation;
		baseExcitation= 0;
	}
	if(baseExcitation0 != 0 ){
		delete baseExcitation0;
		baseExcitation0= 0;
	}
	if(mirExcitation!= 0 ){
		delete mirExcitation;
		mirExcitation= 0;
	}
	if(thePerformFuncCoeffs!= 0 ){
		delete thePerformFuncCoeffs;
		thePerformFuncCoeffs= 0;
	}
	if(thePfCoeffIter!= 0 ){
		delete thePfCoeffIter;
		thePfCoeffIter= 0;
	}
}
void RandomVibrationAnalysis::analyze()
{
	outputFile.open(fileName, ios::out);
	// desingpoint = 0 : No design-point analysis
	//             = 1 : design-point & out crossing rate analysis 
	//	           = 2 : read previous results 
	//			   = 3 : read previous results and compute crossing rate
	this->PrintTitle();

	if(designpoint!=0){
		if(designpoint==1){
			this->DesignPoints();
		}else if(designpoint==2 || designpoint==3){
			this->RdDesignPoints();
//			this->readPreviousResults();
			if(designpoint==3){
			}
		}else{
		}
	}
	outputFile.flush();
}	
void 
RandomVibrationAnalysis::DesignPoints(void)
{
 	outputFile<<"\n";
	outputFile<<"###########################################################\n";
	outputFile<<"#####            Design Point Analysis                #####\n";
	outputFile<<"###########################################################\n";

	hfunc= new Vector(numRV);  

	Matrix* hfuncs = new Matrix(numRV,numFragility);
	Vector* threval= new Vector(numFragility);

	double completedThreshold,completedTime;
	int completedSteps;
	
	int ilocres=0;

	/// clear results ////
	theOutCrossingResults->clear(numLsf, numFragility, numRV);

 	for (int lsf=1; lsf<=numLsf; lsf++ ) {
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);
		theLSF = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
		double Performthreshold0=this->setLimitState(theLSF);
		
		if(theInitialPointBuilder!=0){
			theInitialPointBuilder->setLsf(lsf);
			if(theInitialPointBuilder->getType()==1){
				theInitialPointBuilder->setOutCrossingResults
				(outputFile,theOutCrossingResults);
			}
		}

        theOutCrossingResults->printSingleTitle(outputFile,lsf);

		for(int iFrag=0;iFrag<numFragility;iFrag++){
			ilocres++;
//			outputFile<<"\n";
//			outputFile<<"Fragility......................."<<iFrag<<"\n";
//			outputFile<<"Fragility coefficient..........."<<(*Fragility)(iFrag)<<"\n";
//			outputFile<<"\n";
			// ANALYZE LIMIT STATE FUNCTION & SET SCALE VALUE //
			if(iFrag!=0) {
				completedThreshold=Performthreshold;
				completedTime=AnalysisTime;
				completedSteps=AnalysisStep;
			}
			Performthreshold=Performthreshold0*(*Fragility)(iFrag);
			(*threval)(iFrag)=Performthreshold;

			this->setScale(fabs(Performthreshold));
			theGFunEvaluator->setThreshold(Performthreshold);
			theGFunEvaluator->setPerformFuncCoeffs(thePerformFuncCoeffs);
			theGFunEvaluator->setPerformFuncCoeffIter(thePfCoeffIter);
			theGradGEvaluator->setPerformFuncCoeffs(thePerformFuncCoeffs);
			theGradGEvaluator->setPerformFuncCoeffIter(thePfCoeffIter);

			this->printPerformanceFunction(lsf);

			if(!mirrorimage&&StartSteps==0){
				opserr << " Either of MirrorImageBuilder or StartTime\n";
				opserr << " is required for Stationary Analysis\n";
				exit(-1);
			}

			if(mirrorimage)  {
				//// Mirror Image Analysis /////
			 	//theMirrorImageBuilder->setPerformFuncCoeffIter(thePfCoeffIter);
				//theMirrorImageBuilder->setpfthreshold(Performthreshold);
				//opserr << " Mirror Image Excitation Analysis \n";
				//mirStep = theMirrorImageBuilder->buildMirrorImage();
				//if(mirExcitation!=0)
				//{ delete mirExcitation; mirExcitation=0;}
				//mirExcitation = new Vector(mirStep);
				//(*mirExcitation)= theMirrorImageBuilder->getMirrorExcitation();
				//mirTime=delta*(float)mirStep;
				//opserr<< "\n";
				//opserr<< "===== Mirror Image Excitation =====\n";
				//opserr<< "Number of Steps...............";
				//opserr.width(5);
				//opserr<< mirStep << "\n";
				//opserr<< "Duration of Mirror Image......";
				//opserr.width(15);
				//opserr << mirTime << "\n";
				//opserr<< "\n";
				//opserr<< "\n";
				opserr<< " No Mirror Image Excitation Option\n";
			}

			if(StartSteps!=0){
				AnalysisStep=StartSteps;
				AnalysisTime=delta*(float)AnalysisStep;
			}else{
				AnalysisStep=mirStep;
				AnalysisTime=mirTime;
			}

			if(theInitialPointBuilder !=0){
				theInitialPointBuilder->setthreshold(fabs(Performthreshold));
				if(theInitialPointBuilder->getType()==0){
				// mirrorimageinitialpoint//
					if(!mirrorimage){
						opserr << " MirrorImageBuilder is required\n";
						opserr << " For MirrorImageInitialPointBuilder \n";
						exit(-1);
					}
					theInitialPointBuilder->setMirrorImageExcitation(mirStep,(*mirExcitation));
					(*xinitial)=theInitialPointBuilder->buildInitialPoint(AnalysisStep);
					numAnaLineSearch=theInitialPointBuilder->getnumAna();
				}else if(theInitialPointBuilder->getType()==1){
					if(iFrag==0){
					/// starting point //
						if(theInitialPointBuilder->startMirror()){
							theInitialPointBuilder->setMirrorImageExcitation(mirStep,(*mirExcitation));
							(*xinitial)=theInitialPointBuilder->buildInitialPoint(-AnalysisStep);
							numAnaLineSearch=theInitialPointBuilder->getnumAna();
						}else{
							xinitial->Zero();
							numAnaLineSearch=0;
						}
					}else{
					/// 2nd or later point starting point //
						theInitialPointBuilder->setPrevResults
						(completedSteps, completedTime, completedThreshold,
						 theFindDesignPointAlgorithm->get_x());
						theInitialPointBuilder->setOutCrossingResults(outputFile,theOutCrossingResults);
				  		(*xinitial)=theInitialPointBuilder->buildInitialPoint(AnalysisStep);
						numAnaLineSearch=theInitialPointBuilder->getnumAna();
					}
				}
				opserr<< "\n";
				opserr<< "===== Initial Point Analysis =====\n";
		  		opserr<< "analysis step.................";
				opserr.width(5);
				opserr<< AnalysisStep << "\n";
 				opserr<< "Number of Line Search.........";
				opserr.width(5);
 				opserr<< numAnaLineSearch << "\n";
				opserr<< "\n";
			}else{
				// Starting From Origin //
				xinitial->Zero();
				numAnaLineSearch=0;
				opserr<< "\n";
				opserr<< "===== Start Point Analysis =====\n";
				opserr<< "Number of Steps...............";
				opserr.width(5);
				opserr<< AnalysisStep << "\n";
				opserr<< "Time..........................";
				opserr.width(15);
				opserr<< AnalysisTime << "\n";
				opserr<< "\n";
			}
			theGFunEvaluator->setNsteps(AnalysisStep);
			theFindDesignPointAlgorithm->set_x(*xinitial);
			iresult=theFindDesignPointAlgorithm->findDesignPoint();
			opserr<< "check--1\n";
	 		if(iresult<0){ 
				opserr << "Warning!! FORMAnalysis::analyze() - failed while finding the" << endln
					   << " design point for limit-state function number " 
					   << lsf << "." << endln;
				for(int i=0;i<numRV;i++)(*hfunc)(i)=0.0;
			}else{
				this->findTELS();
				//// save results /////	
//				crossRate=this->crossingRate();
			}
			for (int i=0;i<numRV; i++){
				(*hfuncs)(i,iFrag)=(*hfunc)(i);
			}
			opserr<< "check--2\n";
			opserr<< "check--3\n";
			//// find TELS under defined random process /////
			this->saveDesResults(ilocres, iresult, lsf, Performthreshold);
			opserr<< "check--4\n";
			theOutCrossingResults->printSinglePoint(outputFile,ilocres);
			opserr<< "check--5\n";
			theOutCrossingResults->outtoFile();
			////////// component first passage analysis /////
			////////// component first passage analysis /////
//			crossRate=0.0;
//			if(firstpassage) analyzeComponentFirstPasssage(lsf,iFrag, iresult);
		}// Fragility loop end;
	///// print hfunc for the Lsf /////
		outputTELs<<"\n";
		outputTELs<<"###################################################\n";
		outputTELs<<"#####     TELS for Lsf"<<lsf<< "\n";
		outputTELs<<"###################################################\n";
		outputTELs<<"\n";
		outputTELs<<"\n";
		outputTELs<<"\n";
		outputTELs<<"\n";
		outputTELs<<" step"<<"      Time";
		outputTELs.setf(ios::right);
		outputTELs.setf(ios::scientific, ios::floatfield);
		for(int iFrag=0;iFrag<numFragility;iFrag++){
			outputTELs << setw(15)<<setprecision(5)<<(*threval)(iFrag);
		}
		outputTELs<<"\n";
		for(int istep=0;istep<numRV;istep++){
			outputTELs << setw(5) << istep;
			outputTELs << setw(10)<<setprecision(2)<<delta*float(istep);
			for(int iFrag=0;iFrag<numFragility;iFrag++){
			outputTELs << setw(15)<<setprecision(5)<<(*hfuncs)(istep,iFrag);
			}
			outputTELs<<"\n";
		}
		outputTELs.flush();
	}// Limit state Function loop end;
	delete hfunc; hfunc=0;
    delete hfuncs; hfuncs=0;
    delete threval; threval=0;

	//this->systemFirstPassage();
}
void 
RandomVibrationAnalysis::RdDesignPoints(void)
{
 	outputFile<<"\n";
	outputFile<<"###########################################################\n";
	outputFile<<"#####            Read Design Point & Analyze          #####\n";
	outputFile<<"###########################################################\n";

	theOutCrossingResults->readfromFile();

	int pnpoints=theOutCrossingResults->getnumPoints();
	int pnlsf=theOutCrossingResults->getnumLsf();
	int pnfrag=theOutCrossingResults->getnumFragility();
	int pnrv=theOutCrossingResults->getnumRV();

	if(pnrv!=numRV){
		opserr<< "!!!!!Fatal!!!!!\n";
		opserr<< "Inconsistent numbers for randomvariables \n";
		opserr<< "RandomVibrationAnalysis::RdDesignPoints(void)\n";
		opserr<< "pnrv "<<pnrv;
		opserr<< "numRV "<<numRV<<"\n";
		exit(-1);
	}
	if(numLsf!=pnlsf||numFragility!=pnfrag||pnpoints!=numDesPoints){
		opserr<< "!!!!!Warning!!!!!\n";
		opserr<< "Inconsistent numbers for theOutCrossingResults \n";
		opserr<< "passednumLsf "<<pnlsf;
		opserr<< "numLsf "<<numLsf<<"\n";
		opserr<< "passednumFrag "<<pnfrag;
		opserr<< "numFragility "<<numFragility<<"\n";
		opserr<< "passednumPoints "<<pnpoints;
		opserr<< "numDesPoints "<<numDesPoints<<"\n";
		numDesPoints=pnpoints;
		numLsf=pnrv;
		numFragility=pnfrag;
	}

	hfunc= new Vector(numRV);  
	Matrix* hfuncs = new Matrix(numRV,numFragility);
	Vector* threval= new Vector(numFragility);

	double completedThreshold,completedTime;
	int completedSteps;
	
	int ilocres=0;

	/// clear results ////
	//theOutCrossingResults->clear(numLsf*numFragility, numRV);

 	for (int lsf=1; lsf<=numLsf; lsf++ ) {
		for(int iFrag=0;iFrag<numFragility;iFrag++){
			ilocres++;
			int idLsf=theOutCrossingResults->getLsf(ilocres);
			double thretmp=theOutCrossingResults->getthresholdVal(ilocres);
			(*threval)(iFrag)=thretmp;
			AnalysisStep=theOutCrossingResults->getnumSteps(ilocres);
			AnalysisTime=theOutCrossingResults->getTime(ilocres);
			(*xDesTmp)=theOutCrossingResults->getxDesPoints(ilocres);
			(*uDesTmp)=theOutCrossingResults->getuDesPoints(ilocres);
			(*alphaTmp)=theOutCrossingResults->getDesAlpha(ilocres);
			(*hfunc)=theOutCrossingResults->getHfunc(ilocres);
			numEvalinFORM=theOutCrossingResults->getnumAna(ilocres);
			numAnaIncSens=theOutCrossingResults->getnumAnaIncSens(ilocres);
			betaTmp =theOutCrossingResults->getbeta(ilocres);
			pfTmp =theOutCrossingResults->getpf(ilocres);
			double check1=theOutCrossingResults->getcheck1(ilocres);
			double check2=theOutCrossingResults->getcheck2(ilocres);
			double check1_init=theOutCrossingResults->getcheck1_init(ilocres);
			double check2_init=theOutCrossingResults->getcheck2_init(ilocres);
			iresult=theOutCrossingResults->getiresult(ilocres);
			for (int i=0;i<numRV; i++){
				(*hfuncs)(i,iFrag)=(*hfunc)(i);
			}
			////////// component first passage analysis /////
			crossRate=0.0;
			if(firstpassage) analyzeComponentFirstPasssage(lsf,iFrag, iresult);
		}// Fragility loop end;
	///// print hfunc for the Lsf /////
		outputTELs<<"\n";
		outputTELs<<"###################################################\n";
		outputTELs<<"#####     TELS for Lsf"<<lsf<< "\n";
		outputTELs<<"###################################################\n";
		outputTELs<<"\n";
		outputTELs<<"\n";
		outputTELs<<"\n";
		outputTELs<<"\n";
		outputTELs<<" step"<<"      Time";
		outputTELs.setf(ios::right);
		outputTELs.setf(ios::scientific, ios::floatfield);
		for(int iFrag=0;iFrag<numFragility;iFrag++){
			outputTELs << setw(15)<<setprecision(5)<<(*threval)(iFrag);
		}
		outputTELs<<"\n";
		for(int istep=0;istep<numRV;istep++){
			outputTELs << setw(5) << istep;
			outputTELs << setw(10)<<setprecision(2)<<delta*float(istep);
			for(int iFrag=0;iFrag<numFragility;iFrag++){
			outputTELs << setw(15)<<setprecision(5)<<(*hfuncs)(istep,iFrag);
			}
			outputTELs<<"\n";
		}
		outputTELs.flush();
	}// Limit state Function loop end;
	delete hfunc; hfunc=0;
    delete hfuncs; hfuncs=0;
    delete threval; threval=0;
	//this->systemFirstPassage();
}
void RandomVibrationAnalysis::analyzeComponentFirstPasssage
							  (int passedlsf, int passedFrag, int pres)
{

	if(pres<0){
		for(int i=0; i<=NumTimePoints; i++){
			(*FirstPassageMatrix[passedlsf-1])(i,passedFrag)=0.0;
		}
	}else{
		Vector FpProb;
		FpProb =theFirstPassageAnalyzer->componentFisrtPassage
					 (uDesTmp, xDesTmp, alphaTmp, hfunc, 
					  betaTmp, pfTmp, AnalysisTime, AnalysisStep,
					  theTimePoints,
					  passedlsf, passedFrag,
					  outputFile); 
		for(int i=0; i<=NumTimePoints; i++){
			(*FirstPassageMatrix[passedlsf-1])(i,passedFrag)=FpProb(i+1);
		}
	}
}
void RandomVibrationAnalysis::PrintTitle()
{
	outputFile<<"###########################################################\n";
	outputFile<<"#####                                                 #####\n";
	outputFile<<"#####          Random Vibration Analysis              #####\n";
	outputFile<<"#####                                                 #####\n";
	outputFile<<"###########################################################\n";
	outputFile<< "\n";
	outputFile<< "delta......................" << delta << "\n"; 
	outputFile<< "delta_Pulse................" << delta_Pulse<< "\n"; 
	outputFile<< "designpoint................" << designpoint << "\n"; 
	outputFile<< "mirrorimage................" << mirrorimage << "\n"; 
	outputFile<< "initialpoint..............." << initialpoint << "\n"; 
	outputFile<< "stationary................." << stationary << "\n"; 
	outputFile<< "\n";
	outputFile<< "numRV......................" << numRV << "\n"; 
	outputFile<< "numRVPos..................." << numRVPos << "\n"; 
	outputFile<< "numLsf....................." << numLsf << "\n"; 
	outputFile<< "NumTotalPulse.............." << NumTotalPulse << "\n"; 
	outputFile<< "StartSteps................." << StartSteps << "\n"; 
	outputFile<< "NumTimePoints.............." << NumTimePoints << "\n"; 
	outputFile<< "numFragility..............." << numFragility<< "\n"; 
	outputFile.flush();
	theTimePoints->printTimePts(outputFile);
	outputFile.flush();
}	
void RandomVibrationAnalysis::printPerformanceFunction(int idlsf)
{
/* 	outputFile.setf( ios::scientific );
	outputFile<< "\n";
	outputFile<< "===== Limit-State Function "; 
	outputFile<< setw(5) << idlsf; 
	outputFile<< " =====\n"; 
	outputFile<< "\n";
	outputFile<< "Limit-State Function Constant....";
	outputFile<< setw(15) << setprecision(5) << Performthreshold << "\n";
	outputFile << "\n";
	outputFile << "      Node";
	outputFile << " Direction";
	outputFile << "     Coefficient\n";
	thePfCoeffIter->reset();
	PerformanceFunctionCoeff* thePfCoeff;
	while((thePfCoeff = (*thePfCoeffIter)()) != 0){
	    int N=thePfCoeff->getNodeID();
	    int D=thePfCoeff->getDirection();
	    double c=thePfCoeff->getCoefficient();
		outputFile<< setw(10) << N;
		outputFile<< setw(10) << D;
		outputFile<< setw(15) << c;
		outputFile<<"\n";
	}
	outputFile.flush();
*/
}
double RandomVibrationAnalysis::setLimitState(LimitStateFunction* theLSF)
{
	if(thePerformFuncCoeffs !=0 )
	{delete thePerformFuncCoeffs;thePerformFuncCoeffs = 0;}

	thePerformFuncCoeffs = new ArrayOfTaggedObjects(32);
	if(thePerformFuncCoeffs==0){
		opserr << "RandomVibrationAnalysis::setLimitState - out of memory "
			   << "for PerfomanceFunctionCoeff \n";
		exit(-1);}

	if(print){
		output << " limit state function analysis \n";
		output << "\n";
		output << "----- Analyze Limit State Function -----\n";
		output << "\n";
		output.flush();
	}

    PerformanceFunctionCoeff* thePFCoeff;

	const char *theExpression = theLSF->getExpression();
	// This parsing should go away -- MHS 10/7/2011
	//char *theTokExpression = theLSF->getTokenizedExpression();
	char *theTokExpression = "";

	char separators[5] = "}{";
	char *dollarSign = "$";
	char *underscore = "_";
	char tempchar[100]="";

	char lsf_forTokenizing[500];
	strcpy(lsf_forTokenizing,theExpression);
	char tclAssignment[100];

//  1st set all disp to zero and evaluate to have the constant term //
	char *tokenPtr = strtok( lsf_forTokenizing, separators);
	while ( tokenPtr != NULL ) {
		strcpy(tempchar,tokenPtr);
		if ( strncmp(tokenPtr, "ud", 2) == 0) {
			opserr << "Performance Function need to be \n";
			opserr << "a linear function of disp to apply \n";
			opserr << "initial shape analysys \n";
			opserr << "skip initial shape for this lsf \n";
			exit(-1);
		}
		// If a nodal displacement is detected
		else if ( strncmp(tokenPtr, "u", 1) == 0) {
			int nodeNumber, direction;
			sscanf(tempchar,"u_%i_%i", &nodeNumber, &direction);
			sprintf(tclAssignment,"set u_%d_%d 0.0",nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		else if ( strncmp(tokenPtr, "rec",3) == 0) {
			opserr << "Performance Function need to be \n";
			opserr << "a linear function of disp to apply \n";
			opserr << "initial shape analysys \n";
			opserr << "skip initial shape for this lsf \n";
			exit(-1);
		}
		tokenPtr = strtok( NULL, separators);
	}
	double pfthreshold = 0.0;
	double dthreshold = 0.0;
	Tcl_ExprDouble( theTclInterp, theTokExpression, &pfthreshold );

	// Re-create possible recorders for subsequent analyses
	strcpy(lsf_forTokenizing,theExpression);
	tokenPtr = strtok( lsf_forTokenizing, separators);
	while ( tokenPtr != NULL ) {
		strcpy(tempchar,tokenPtr);

        if ( strncmp(tokenPtr, "u", 1) == 0) {
			int nodeNumber, direction;
			sscanf(tempchar,"u_%i_%i", &nodeNumber, &direction);
			sprintf(tclAssignment,"set u_%d_%d 1.0",nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
			Tcl_ExprDouble( theTclInterp, theTokExpression, &dthreshold );
			double coeff=dthreshold-pfthreshold;
			int Tag=thePerformFuncCoeffs->getNumComponents()+1;
            thePFCoeff=new PerformanceFunctionCoeff
		     (Tag, nodeNumber, direction, coeff);
			bool result = 
				thePerformFuncCoeffs->addComponent(thePFCoeff);
			if(!result){
		      opserr << "OutCrossingAnalysis::AnalyzeGfun - out of memory3\n";
			exit(-1);}
			sprintf(tclAssignment,"set u_%d_%d 0.0",nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		
		tokenPtr = strtok( NULL, separators);
	}

	if(thePfCoeffIter !=0 )
	{delete thePfCoeffIter;thePfCoeffIter = 0;}
	thePfCoeffIter=new PerformanceFunctionCoefficientIter(thePerformFuncCoeffs);
	if(thePfCoeffIter==0){
       opserr << "OutCrossingAnalysis::AnalyzeGfun - out of memory3\n";
	   exit(1);
	}

	if(print){
		output.setf( ios::scientific );
		output << "Performance Function Constant............";
		output << setw(15) << setprecision(5) << pfthreshold << "\n";
		output << "\n";
		output << "      Node";
		output << " Direction";
		output << "     Coefficient\n";
		thePfCoeffIter->reset();
		PerformanceFunctionCoeff* thePfCoeff;
		while((thePfCoeff = (*thePfCoeffIter)()) != 0){
		    int N=thePfCoeff->getNodeID();
		    int D=thePfCoeff->getDirection();
		    double c=thePfCoeff->getCoefficient();
			output<< setw(10) << N;
			output<< setw(10) << D;
			output<< setw(15) << c;
		}
		output.flush();
	}
	return pfthreshold;
}
void RandomVibrationAnalysis::setScale(double scale)
{
	theReliabilityConvergenceCheck->Scalefix(false);
	theReliabilityConvergenceCheck->setScaleValue(fabs(Performthreshold));
	theReliabilityConvergenceCheck->Scalefix(true);
}
void  RandomVibrationAnalysis::saveDesResults(int loc, int presult,
											  int plsf, double pthre)
{
  static NormalRV aStdNormRV(1,0.0,1.0);

	(*xDesTmp)=theFindDesignPointAlgorithm->get_x();
	(*uDesTmp)=theFindDesignPointAlgorithm->get_u();
	(*alphaTmp)=theFindDesignPointAlgorithm->get_alpha();
	numEvalinFORM=theFindDesignPointAlgorithm->getNumberOfEvaluations();
	numAnaIncSens=theFindDesignPointAlgorithm->getNumberOfSensAna();
	numIterinFORM=theFindDesignPointAlgorithm->getNumberOfSteps();
	betaTmp =(*alphaTmp)^(theFindDesignPointAlgorithm->get_u());
	pfTmp = 1.0 - aStdNormRV.getCDFvalue(betaTmp);
	double check1=theReliabilityConvergenceCheck->getCheck1();
	double check2=theReliabilityConvergenceCheck->getCheck2();
	double check1_init=theFindDesignPointAlgorithm->get_check1_init();
	double check2_init=theFindDesignPointAlgorithm->get_check2_init();
	double check1_conv=theFindDesignPointAlgorithm->get_check1_conv();
	double check2_conv=theFindDesignPointAlgorithm->get_check2_conv();
//          store the results
	theOutCrossingResults->setLsf(loc,plsf);
	theOutCrossingResults->setnumAna(loc,numEvalinFORM);
	theOutCrossingResults->setnumLinSearch(loc,numAnaLineSearch);
	theOutCrossingResults->setnumAnaIncSens(loc,numAnaIncSens);
	theOutCrossingResults->setxDesPoints(loc,(*xDesTmp));
	theOutCrossingResults->setuDesPoints(loc,(*uDesTmp));
	theOutCrossingResults->setDesAlpha(loc,(*alphaTmp));
	theOutCrossingResults->setHfunc(loc,(*hfunc));
	theOutCrossingResults->setthresholdVal(loc,pthre);
	theOutCrossingResults->setbeta(loc,betaTmp);
	theOutCrossingResults->setpf(loc,pfTmp);
	theOutCrossingResults->setnu(loc,crossRate);
	theOutCrossingResults->settime(loc,AnalysisTime);
	theOutCrossingResults->setnumSteps(loc,AnalysisStep);
	theOutCrossingResults->setcheck1(loc,check1_conv);
	theOutCrossingResults->setcheck2(loc,check2_conv);
	theOutCrossingResults->setcheck1_init(loc,check1_init);
	theOutCrossingResults->setcheck2_init(loc,check2_init);
	theOutCrossingResults->setiresult(loc,presult);
}
double RandomVibrationAnalysis::crossingRate(void)
{
//	Vector udes;
//	Vector alpha;
	opserr<< "check--11\n";

	///// extract design point /////
	(*uDesTmp)=theFindDesignPointAlgorithm->get_u();
	opserr<< "check--22\n";
	(*alphaTmp)=theFindDesignPointAlgorithm->get_alpha();
	opserr<< "check--33\n";
	betaTmp=(*uDesTmp)^(*alphaTmp);
	//// setting CrossingRateAnalyzer ///// 
	opserr<< "check--44\n";
	theCrossingRateAnalyzer->setAnalysis(AnalysisStep,(*uDesTmp),(*alphaTmp), betaTmp);
	//// compute CrossingRateAnalyzer ///// 
	opserr<< "check--55\n";
	RateTmp=theCrossingRateAnalyzer->computeRate();
//	outputFile<< "\n";
//	outputFile<< "Corssing Rate Analyhsis\n";
//	outputFile<< "crossing rate...................."<<RateTmp<<"\n";
	opserr<< "\n";
	opserr<< "Corssing Rate Analyhsis\n";
	opserr<< "crossing rate...................."<<RateTmp<<"\n";
	return RateTmp;
}
void RandomVibrationAnalysis::findTELS(void)
{
//	
//  Find TELS for the obtained design point
//
//
//  step1 make A matrix for the randomprocess
//
//	numRV ! num of RVs
//	numRVPos = num of RVpositioner
//	theRandomProcess ! random process
//	NumTotalPulse
//	delta_Pulse 
//	AnalysisStep=StartSteps;
//	AnalysisTime=delta*(float)AnalysisStep;
//
//  ----- built a-vector -----
//
	(*uDesTmp)=theFindDesignPointAlgorithm->get_u();
	(*alphaTmp)=theFindDesignPointAlgorithm->get_alpha();
	double length=(*uDesTmp).Norm();
	Vector* avector= new Vector(numRV);  
	(*avector)=Performthreshold/length/length*(*uDesTmp);
//
//	----- step1 build the Transpose of Jacobian Matrix
//         J_ij=S_j(t_i) i: time for the excitation j : kickin time for pulse j
//        J^T_ij=S_i(t_j) i : Kick in time for the pulse i 
//
	int nactive_Pulse=theRandomProcess->getNumOfActivePulses(AnalysisTime);
	if(nactive_Pulse!=AnalysisStep){
		opserr<< "\n";
		opserr<< "Fatal\n";
		opserr<< "nactive_Pulse!=AnalysisStep\n";
		opserr<< "nactive_Pulse "<<nactive_Pulse;
		opserr<< "AnalysisStep "<<AnalysisStep;
	}

	Vector* aavector= new Vector(nactive_Pulse);  
	for(int ipulse=0;ipulse<nactive_Pulse;ipulse++){
		int iii=theRandomProcess->getRVseqid(ipulse);
		(*aavector)(ipulse)=(*avector)(iii-1);
	}

	double aa1=(*aavector)(0);
	double aa2=(*aavector)(1);
	double aa3=(*aavector)(2);
	double aan=(*aavector)(nactive_Pulse-1);
	double aan1=(*aavector)(nactive_Pulse-2);
	double aan2=(*aavector)(nactive_Pulse-3);

	Matrix* Jacobian_x_u= new Matrix(numRV,numRV);
	(*Jacobian_x_u)=theFindDesignPointAlgorithm->getJacobian_x_u();
	Matrix* JJacobian_x_u= new Matrix(nactive_Pulse,nactive_Pulse);
	for(int ipulse=0;ipulse<nactive_Pulse;ipulse++){
		int iii=theRandomProcess->getRVseqid(ipulse);
		for(int jpulse=0;jpulse<nactive_Pulse;jpulse++){
			int jjj=theRandomProcess->getRVseqid(jpulse);
			(*JJacobian_x_u)(ipulse,jpulse)=(*Jacobian_x_u)(iii-1,jjj-1);
		}
	}
	double jj1=(*JJacobian_x_u)(0,0);
	double jj2=(*JJacobian_x_u)(1,1);
	double jj3=(*JJacobian_x_u)(2,2);
	double jjn=(*JJacobian_x_u)(nactive_Pulse-1,nactive_Pulse-1);
	double jjn1=(*JJacobian_x_u)(nactive_Pulse-2,nactive_Pulse-2);
	double jjn2=(*JJacobian_x_u)(nactive_Pulse-3,nactive_Pulse-3);
///
/// ------ compute matrix -----
///
	   	 
	Matrix* TJacobian = new Matrix(nactive_Pulse,AnalysisStep);

	double time_pls;
	double time_ext;
	for(int jstep=0;jstep<AnalysisStep;jstep++){
		time_pls=delta_Pulse*(float)(jstep+1);
		for(int istep=0;istep<nactive_Pulse;istep++){
			time_ext=delta*(float)(istep+1);
			double aaa=theRandomProcess->getFactorSensitivity(time_ext,time_pls);
			(*TJacobian)(jstep,istep)=aaa;
		}
	}

	double tjn=(*TJacobian)(AnalysisStep-1,AnalysisStep-1);
	double tjn1=(*TJacobian)(AnalysisStep-2,AnalysisStep-1);
	double tjn2=(*TJacobian)(AnalysisStep-2,AnalysisStep-2);

///
///  ============ for filtered white noise ========
/// 
	for(int jstep=0;jstep<nactive_Pulse;jstep++){
		double tmp=(*aavector)(jstep);
		(*aavector)(jstep)=tmp/delta/(*JJacobian_x_u)(jstep,jstep);
	}
	double aaa1=(*aavector)(0);
	double aaa2=(*aavector)(1);
	double aaa3=(*aavector)(2);
	double aaan=(*aavector)(nactive_Pulse-1);
	double aaan1=(*aavector)(nactive_Pulse-2);
	double aaan2=(*aavector)(nactive_Pulse-3);
	Vector* tmp= new Vector(AnalysisStep);  
	(*tmp)(AnalysisStep-1)=(*aavector)(AnalysisStep-1)/(*TJacobian)(AnalysisStep-1,AnalysisStep-1);
	for(int jstep=AnalysisStep-2;jstep>=0;jstep--){
		double sum=(*aavector)(jstep);
		for(int kstep=jstep+1; kstep<AnalysisStep; kstep++){
			sum=sum-(*TJacobian)(jstep,kstep)*(*tmp)(kstep);
		}
		(*tmp)(jstep)=sum/(*TJacobian)(jstep,jstep);
	}
	double ttt1=(*tmp)(0);
	double ttt2=(*tmp)(1);
	double ttt3=(*tmp)(2);
	double tttn=(*tmp)(AnalysisStep-1);
	double tttn1=(*tmp)(AnalysisStep-2);
	double tttn21=(*tmp)(AnalysisStep-3);
///


			
//	MatrixOperations* theMatrixOperations = 0;
//	theMatrixOperations = new MatrixOperations(*TJacobian);
//	if (theMatrixOperations == 0) {
//		opserr << "RandomVibrationAnalysis::findTELS - could " << endln
//			<< " not create the object to perform matrix operations." << endln;
//	}

	// Cholesky decomposition of correlation matrix
//	int result = theMatrixOperations->computeInverse();
//	if (result < 0) {
//		opserr << "RandomVibrationAnalysis::findTELS() - could not" << endln
//			<< " compute the inverse of Jacobian " << endln;
//	}
//	Matrix* invJacob = new Matrix(nactive_Pulse,nactive_Pulse);
//	(*invJacob)= theMatrixOperations->getInverse();
	
	///// compute h-vector /////
//	Vector* tmp= new Vector(nactive_Pulse);  
//	(*tmp)=(*invJacob)^(*aavector);

	for(int i=0;i<nactive_Pulse;i++){
		(*hfunc)(i)=(*tmp)(nactive_Pulse-i-1);
	}

	if(nactive_Pulse<numRV){
		for(int i=nactive_Pulse;i<numRV;i++){
			(*hfunc)(i)=0.0;
		}
	}
	double hhh1=(*hfunc)(0);
	double hhh2=(*hfunc)(1);
	double hhh3=(*hfunc)(2);
	double hhhn=(*hfunc)(nactive_Pulse-1);
	double hhhn1=(*hfunc)(nactive_Pulse-2);
	double hhhn2=(*hfunc)(nactive_Pulse-3);

	delete avector; avector=0;
	delete aavector; aavector=0;
	delete TJacobian;  TJacobian=0;
//	delete invJacob;  invJacob=0;
	delete tmp; tmp=0;
	delete Jacobian_x_u;
//	delete theMatrixOperations;  theMatrixOperations=0;

}
/////////////////////////////////////////////////////////////////////////////////
////////////////////// previous version /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*
void 
RandomVibrationAnalysis::nonStatDesignPoints(void)
{
	outputFile<<"\n";
	outputFile<<"###########################################################\n";
	outputFile<<"#####    Design Point Analysis (non-stationary)       #####\n";
	outputFile<<"###########################################################\n";

	for(int iFrag=0;iFrag<numFragility;iFrag++){
		//// loop for Fragility steps /////
		outputFile<<"\n";
		outputFile<<"Fragility......................."<<iFrag<<"\n";
		outputFile<<"Fragility coefficient..........."<<(*Fragility)(iFrag)<<"\n";
		outputFile<<"\n";

		for (int lsf=1; lsf<=numLsf; lsf++ ) {
			//// loop for limit state functions /////
			theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);
			theLSF = theReliabilityDomain->getLimitStateFunctionPtr(lsf);

			// ANALYZE LIMIT STATE FUNCTION & SET SCALE VALUE //
			Performthreshold=this->setLimitState(theLSF);
			double fragfactor=(*Fragility)(iFrag);
			Performthreshold*=fragfactor;
			this->setScale(fabs(Performthreshold));
			theGFunEvaluator->setThreshold(Performthreshold);
			theGFunEvaluator->setPerformFuncCoeffs(thePerformFuncCoeffs);
			theGFunEvaluator->setPerformFuncCoeffIter(thePfCoeffIter);
			theGradGEvaluator->setPerformFuncCoeffs(thePerformFuncCoeffs);
			theGradGEvaluator->setPerformFuncCoeffIter(thePfCoeffIter);
			this->printPerformanceFunction(lsf);
			if(theInitialPointBuilder !=0)
			theInitialPointBuilder ->setthreshold(fabs(Performthreshold));

			/// initialize results object ///
			theOutCrossingResults->clear(lsf,iFrag, Performthreshold, fragfactor,
										 NumTimePoints, numRV);
			/// initialize analysis order /////
			theTimePoints->initializeOrder();
			int ipos=0;

			//////////////////////////////////////////////////////////////
			/////
			///// starting point analysis if required 
			///// mirrorimage = true  or  StartSteps !=0
			/////
			///// if mirrorimage = true ( then initialpoint must be true )
			/////    1. construct mirror image exictaion
			/////    2. construct initial point by initial point builder
			/////    if required (StartSteps!=0)  shift the mirror image 
			///// if mirrorimage = false ( then StartSteps must be non zero )
			/////    analyze design point for StartPoint starting from origin
			/////			
			//////////////////////////////////////////////////////////////
			if(mirrorimage||StartSteps!=0){
				if(mirrorimage)  {
					//// Mirror Image Analysis /////
//					theMirrorImageBuilder->setPerformFuncCoeffIter(thePfCoeffIter);
//					theMirrorImageBuilder->setpfthreshold(Performthreshold);
//					opserr << " Mirror Image Excitation Analysis \n";
//					mirStep = theMirrorImageBuilder->buildMirrorImage();
//					if(mirExcitation!=0){ delete mirExcitation; mirExcitation=0;} 
//					mirExcitation = new Vector(mirStep);
//					(*mirExcitation)= theMirrorImageBuilder->getMirrorExcitation();
//					mirTime=delta*(float)mirStep;
//					outputFile<< "\n";
//					outputFile<< "===== Mirror Image Excitation =====\n";
//					outputFile<< "Number of Steps..............." << setw(5) << mirStep << "\n";
//					outputFile<< "Duration of Mirror Image......" << setw(15) << mirTime << "\n";
//					outputFile<< "\n";
//					opserr<< "\n";
//					opserr<< "===== Mirror Image Excitation =====\n";
//					opserr<< "Number of Steps...............";
//					opserr.width(5);
//					opserr<< mirStep << "\n";
//					opserr<< "Duration of Mirror Image......";
//					opserr.width(15); opserr.setFloatField(FIXEDD); opserr.precision(2);
//					opserr<< mirTime << "\n";
//					opserr<< "\n";
//					/// set initialpointbuilder ///
//					theInitialPointBuilder->setMirrorImageExcitation(mirStep,(*mirExcitation));
//					/// set current analysis step and time ///
//					if(StartSteps==0){
//						AnalysisStep=mirStep;
//						AnalysisTime=mirTime;
//					}else{
//						AnalysisStep=StartSteps;
//						AnalysisTime=delta*(float)AnalysisStep;
//					}
//					/// getting initial point ///
//					(*xinitial)=theInitialPointBuilder->buildInitialPoint(AnalysisStep);
//					/// getting number of dynamig response analyses ///
//					numAnaLineSearch=theInitialPointBuilder->getnumAna();
//					outputFile<< "\n";
//					outputFile<< "===== Initial Point Analysis =====\n";
//					outputFile<< "analysis step................." << setw(5) << AnalysisStep << "\n";
//					outputFile<< "Number of Line Search........." << setw(5) << numAnaLineSearch << "\n";
//					outputFile<< "\n";
//					opserr<< "\n";
//					opserr<< "===== Initial Point Analysis =====\n";
//					opserr<< "analysis step.................";
//					opserr.width(5); opserr << AnalysisStep << "\n";
//					opserr<< "Number of Line Search.........";
//					opserr.width(5); opserr << numAnaLineSearch << "\n";
//					opserr<< "\n";
					opserr<< " No mirror image option \n";
				}else {
					//// withoug Mirror Image Analysis    /////
					//// But starting from a middle point /////
					//// setting current analysis step and time /////  
					AnalysisStep=StartSteps;
					AnalysisTime=delta*(float)AnalysisStep;
					numAnaLineSearch=0; /// no initial point builder ///
					xinitial->Zero();   /// starting from the origin 
					outputFile<< "\n";
					outputFile<< "===== Start Point Analysis =====\n";
					outputFile<< "Number of Steps..............." << setw(5) << AnalysisStep << "\n";
					outputFile<< "Time.........................." << setw(15) << AnalysisTime << "\n";
					outputFile<< "\n";
					opserr<< "\n";
					opserr<< "===== Start Point Analysis =====\n";
					opserr<< "Number of Steps...............";
					opserr.width(5); opserr << AnalysisStep << "\n";
					opserr<< "Time..........................";
					opserr.width(15); opserr.setFloatField(FIXEDD); opserr.setPrecision(2);
					opserr<< AnalysisTime << "\n";
					opserr<< "\n";
				}
				//// setting number of steps 
				theGFunEvaluator->setNsteps(AnalysisStep);
				//// design point search 
				theFindDesignPointAlgorithm->set_x(*xinitial);
				iresult=theFindDesignPointAlgorithm->findDesignPoint();
				if(iresult<0){ 
					opserr << "Warning!! FORMAnalysis::analyze() - failed while finding the" << endln
						   << " design point for limit-state function number " << lsf << "." << endln;
				}
				outputFile<< "\n";
				outputFile<< "Finding Design Point Complete\n";
				outputFile<< "iresult..........................." << setw(5) << iresult << "\n";
				outputFile<< "\n";
				//// crossing rate analysis /////	
				crossRate=this->crossingRate();

				//// save results /////	
				ipos=theTimePoints->MakeOrder(AnalysisTime);
				this->saveDesResults(ipos, iresult);

				//// set base excitation for shifting  /////
				if(baseExcitation!=0){ delete baseExcitation; baseExcitation=0;}
				baseExcitation=new Vector(AnalysisStep);
				(*baseExcitation)=theRandomProcess->getExcitation(AnalysisStep,delta);
				baseStep=AnalysisStep;
				baseStep0=baseStep;
				if(baseExcitation0!=0){ delete baseExcitation0; baseExcitation0=0;}
				baseExcitation0=new Vector(AnalysisStep);
				(*baseExcitation0)=(*baseExcitation);
			}

			///// number of Analysis points /////
			int numAnalysis=NumTimePoints;
			if(ipos!=0) numAnalysis--;   /// already analyzed one of them

			for(int ipts=0;ipts<numAnalysis;ipts++){
				///// curremt analysis step & time  /////
				int iord=theTimePoints->getOrder(ipts);
				AnalysisStep=theTimePoints->getAnalysisStep(abs(iord));
				AnalysisTime=delta*(float)AnalysisStep;
				///// print /////
				outputFile<< "\n";
				outputFile<< "Next Time point\n";
				outputFile<< "Order of Time point..............." << setw(5) <<iord<< "\n";
				outputFile<< "AnalysisStep......................" << setw(5) <<AnalysisStep<< "\n";
				outputFile<< "\n";
				opserr<< "\n";
				opserr<< "AnalysisStep......................";
				opserr.width(5); opserr <<AnalysisStep<< "\n";
				opserr<< "\n";
				if(iord<0){
					///// change base excitation /////
					baseStep=baseStep0;
					if(baseExcitation!=0){ delete baseExcitation; baseExcitation=0;}
					baseExcitation=new Vector(baseStep);
					(*baseExcitation)=(*baseExcitation0);
					iord=-iord;
					outputFile<< "\n";
					outputFile<< "Base Point Change\n";
					opserr<< "\n";
					opserr<< "Base Point Change\n";
				}
				if(initialpoint){
					///// build initial point /////
					theInitialPointBuilder->
						setMirrorImageExcitation(baseStep,(*baseExcitation));
					(*xinitial)=theInitialPointBuilder->
						buildInitialPoint(AnalysisStep);
					numAnaLineSearch=theInitialPointBuilder->getnumAna();
					outputFile<< "\n";
					outputFile<< "===== Initial Point Analysis =====\n";
					outputFile<< "Number of Line Search........." << setw(5) << numAnaLineSearch << "\n";
					outputFile<< "\n";
					opserr<< "\n";
					opserr<< "===== Initial Point Analysis =====\n";
					opserr<< "Number of Line Search.........";
					opserr.width(5); opserr << numAnaLineSearch << "\n";
					opserr<< "\n";
				}else{
					xinitial->Zero();
					numAnaLineSearch=0;
				}
				theGFunEvaluator->setNsteps(AnalysisStep);
				theFindDesignPointAlgorithm->set_x(*xinitial);
				iresult=theFindDesignPointAlgorithm->findDesignPoint();
				if(iresult<0){ 
					opserr << "Warning FORMAnalysis::analyze() - failed while finding the" << endln
							<< " design point for limit-state function number " 
							<< lsf << "." << endln;
				}
				outputFile<< "\n";
				outputFile<< "Finding Design Point Complete\n";
				outputFile<< "iresult..........................." << setw(5) << iresult << "\n";
				outputFile<< "\n";
				// crossing rate analysis
				crossRate=this->crossingRate();
				// save results //
				this->saveDesResults(iord, iresult);
				if(baseExcitation!=0){ delete baseExcitation; baseExcitation=0;}
				baseExcitation=new Vector(AnalysisStep);
				(*baseExcitation)=theRandomProcess->getExcitation(AnalysisStep,delta);
				baseStep=AnalysisStep;
			} // Time points loop end 
			theOutCrossingResults->printResults(outputFile);
			theOutCrossingResults->outtoFile();
			////////// component first passage analysis /////
			if(firstpassage) analyzeComponentFirstPasssage(lsf,iFrag);
		}// Limit state Function loop end;
//		this->nonStatSystemFirstPassage(iFrag);
	}// Fragility loop end;
}
*/



