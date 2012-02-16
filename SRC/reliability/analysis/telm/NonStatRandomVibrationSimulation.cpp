// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NonStatRandomVibrationSimulation.cpp,v $

#include <NonStatRandomVibrationSimulation.h>
#include <TimePoints.h>
#include <RandomVibrationSimulatorResult.h>
NonStatRandomVibrationSimulation::NonStatRandomVibrationSimulation
								(ReliabilityDomain* passedReliabilityDomain,
								 Domain* passedDomain,
						         FunctionEvaluator* passedGFunEvaluator,
							     ProbabilityTransformation* passedTransformation,
						         double passedStartTime,
						         double passedEndTime,
						         double passedTimeInterval,
						         double passedFragMin,
						         double passedFragInt,
						         int passednFrag,
								 bool passedtwoside,
								 bool passedsystem,
						         int passedmaxSim,
						         int passedcheckinterval,
						         double passedeps,
						         int passedinstantaneous,
						         int passedfirstpassage,
	  				             TCL_Char *passedFileName,
						         char* passedFileBinary,
 	                             Tcl_Interp *passedTclInterp,
								 bool passedprint)
:RandomVibrationSimulation
(passedReliabilityDomain, passedDomain, passedGFunEvaluator,passedTransformation,
 passedFragMin, passedFragInt, passednFrag, passedtwoside, passedsystem, passedmaxSim,
 passedcheckinterval, passedeps, passedinstantaneous,
 passedfirstpassage, passedFileName, passedFileBinary,passedTclInterp)
{
	print=passedprint;
	if(print){
		output.open("NonStatRandomVibrationSimulation.txt", ios::out);
		output << "\n";
		output << "MaxSim....................." << MaxSim << "\n"; 
		output << "eps........................" << eps << "\n"; 
		output << "checkinterval.............." << checkinterval << "\n"; 
		output << "instantaneous.............." << instantaneous <<"\n"; 
		output << "firstpassage..............." << firstpassage  << "\n"; 
		output << "numRV......................" << numRV  << "\n"; 
		output << "numRVPos..................." << numRVPos << "\n"; 
		output << "numLsf....................." << numLsf << "\n"; 
		output << "delta......................" << delta << "\n"; 
		output << "delta_Pulse................" << delta_Pulse << "\n"; 
		output << "numFragility..............." << numFragility  << "\n"; 
		for(int i=0; i<numFragility; i++){
		output<< setw(5) << i << "   " << setw(15) << (*Fragility)(i) <<"\n";
		}
	}
	TimePoints* theTimePoints = new TimePoints(passedStartTime,
								    passedEndTime,
								    passedTimeInterval,
								    delta,false);
	numTimePoints= theTimePoints->getNumTimePoints();
//	if(timepoints!=0){ delete timepoints; timepoints=0;}
	timepoints=new Vector(numTimePoints);
//	if(anaSteps!=0){ delete anaSteps; anaSteps=0;}
	anaSteps=new int[numTimePoints];
	for(int i=0;i<numTimePoints;i++){
		(*timepoints)(i)=theTimePoints->getAnalysisTime(i);
		anaSteps[i]=theTimePoints->getAnalysisStep(i);
	}
	if(print){
		output << "\n";
		output << "numTimePoints.............." << numTimePoints << "\n"; 
		output << "\n";
		for(int i=0;i<numTimePoints;i++){
			output << setw(5) << i <<"     "<< setw(12) << (*timepoints)(i);
			output <<"     "<< setw(5) << anaSteps[i] <<"\n";
			output.flush();
		}
	}
	delete theTimePoints;
	theTimePoints=0;
}

NonStatRandomVibrationSimulation::~NonStatRandomVibrationSimulation()
{
	if(timepoints!=0){ delete timepoints; timepoints=0;} 
	if(anaSteps!=0){ delete anaSteps; anaSteps=0;} 
}
void NonStatRandomVibrationSimulation::crudeInstantaneousSimulation()
{
	outputFile.open(fileName, ios::out);
	outputFile<<"\n";
	outputFile<<"###########################################################\n";
	outputFile<<"#####    RandomVibraionSimulation(non-stationary)     #####\n";
	outputFile<<"#####    Instantaneous failure probability            #####\n";
	outputFile<<"###########################################################\n";
	outputFile<<"\n";
	outputFile<< "numRV......................" << numRV  << "\n"; 
	outputFile<< "numRVPos..................." << numRVPos << "\n"; 
	outputFile<< "numLSF....................." << numLsf << "\n"; 
	outputFile<< "delta......................" << delta << "\n"; 
	outputFile<< "delta_Pulse................" << delta_Pulse << "\n"; 
	outputFile<< "\n";
	outputFile<< "MaxSim....................." << MaxSim << "\n"; 
	outputFile<< "eps........................" << eps << "\n"; 
	outputFile<< "checkinterval.............." << checkinterval << "\n"; 
	outputFile<< "\n";
	outputFile<< "numFragility..............." << numFragility  << "\n"; 
	outputFile.setf(ios::right);
	outputFile.setf(ios::fixed, ios::floatfield);
	for(int i=0; i<numFragility; i++){
		outputFile << setw(5) << i << "   ";
	 outputFile << setw(15)<< setprecision(3) << (*Fragility)(i) <<"\n";
	}
	outputFile << "\n";
	outputFile << "numTimePoints.............." << numTimePoints << "\n"; 
	outputFile << "\n";
	for(int i=0;i<numTimePoints;i++){
		outputFile << setw(5) << i <<"     ";
		outputFile << setw(15) << setprecision(3) <<(*timepoints)(i)<<"\n";
		outputFile.flush();
	}
	///// check time /////
	int numSteps=checkTimePoints();

	///// allocate Results
	RandomVibrationSimulatorResult** theResults;
	theResults=new RandomVibrationSimulatorResult*[numLsf];
	for(int i=0; i<numLsf; i++) theResults[i]=
		new RandomVibrationSimulatorResult(i,numTimePoints,numFragility,eps);
				
	//// construct GFunEachStepEvaluator /////
	GFunEachStepEvaluator* theGFunEachStepEvaluator;
	theGFunEachStepEvaluator= new GFunEachStepEvaluator
						 (theTclInterp, theReliabilityDomain, theDomain,
						  numSteps,print);
	Vector* PFthreshold=theGFunEachStepEvaluator->getPFthershold();
	theGFunEvaluator->setGFunEachStepEvaluator(theGFunEachStepEvaluator);
	theGFunEvaluator->setNsteps(numSteps);

	int allconv,iconv,lsf,jpt,ifrag;
	Matrix* theLSFValues;
	Matrix* theLSFConv;
	double threshold0,response,frag,qvalue;

	int checkstep=checkinterval;
	int ibreak=0;
	int result;


	for(int isim=1;isim<=MaxSim;){

		///// generate random vector /////
		this->generateRV();

		///// dynamic response analysis /////
		bool FEconvergence = true;
		result = theGFunEvaluator->runAnalysis();
		if (result < 0) {
			// In this case a failure happened during the analysis
			// Hence, register this as failure
			FEconvergence = false;
		}

		///// extract and check 
		theLSFValues=theGFunEvaluator->getEachStepResult();
		theLSFConv=theGFunEvaluator->getEachStepConvFlag();
		allconv=1;
		for (lsf=0;lsf<numLsf; lsf++){
			threshold0=(*PFthreshold)(lsf);
			for(jpt=0; jpt<numTimePoints; jpt++){
				if((*theLSFConv)(lsf,anaSteps[jpt]-1)!=0.0){
					response=(*theLSFValues)(lsf,anaSteps[jpt]-1);
					for(ifrag=0; ifrag<numFragility; ifrag++){
						frag=(*Fragility)(ifrag);
						qvalue=0.0;
						if(response>(frag*threshold0)) qvalue=1.0;
						if(twoside&&(response<((frag*threshold0)*(-1.0)))) qvalue=1.0;
//						iconv=theResults[lsf]->updateq(jpt, ifrag, isim, qvalue);
						iconv=theResults[lsf]->updateq(jpt, ifrag, qvalue);
						allconv*=iconv;
					}
				}else{
					for(ifrag=0; ifrag<numFragility; ifrag++){
						qvalue=1.0;
//						iconv=theResults[lsf]->updateq(jpt, ifrag, isim, qvalue);
						iconv=theResults[lsf]->updateq(jpt, ifrag, qvalue);
						allconv*=iconv;
					}
				}
			}
		}

		if(isim==checkstep){
			if(allconv==1) ibreak=1;
			outputFile<<"\n";
			outputFile<<"----- result at simulation"<<isim<<" -----\n"; 
			for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print1(outputFile);
		}
		if(ibreak==1) break;
	}

	///// simulation ends /////
	outputFile<<"\n";
	outputFile<<"==== simulation results =====\n"; 
	outputFile<<"\n";
	for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print2(outputFile);

	for(int i=0; i<numLsf; i++){ delete theResults[i]; theResults[i]=0; }
	delete [] theResults;
	theResults=0;
	theGFunEvaluator->inactivateGFunEachStepEvaluator();
	delete theGFunEachStepEvaluator;
	theGFunEachStepEvaluator=0;
}
void NonStatRandomVibrationSimulation::crudeFisrtpassageSimulation()
{
	outputFile.open(fileName, ios::out);
	outputFile<<"\n";
	outputFile<<"###########################################################\n";
	outputFile<<"#####    RandomVibraionSimulation(non-stationary)     #####\n";
	outputFile<<"#####    First Passage failure probability            #####\n";
	outputFile<<"###########################################################\n";
	outputFile<<"\n";
	outputFile<< "numRV......................" << numRV  << "\n"; 
	outputFile<< "numRVPos..................." << numRVPos << "\n"; 
	outputFile<< "numLSF....................." << numLsf << "\n"; 
	outputFile<< "delta......................" << delta << "\n"; 
	outputFile<< "delta_Pulse................" << delta_Pulse << "\n"; 
	outputFile<< "\n";
	outputFile<< "MaxSim....................." << MaxSim << "\n"; 
	outputFile<< "eps........................" << eps << "\n"; 
	outputFile<< "checkinterval.............." << checkinterval << "\n"; 
	outputFile<< "\n";
	outputFile<< "numFragility..............." << numFragility  << "\n"; 
	outputFile.setf(ios::right);
	outputFile.setf(ios::fixed, ios::floatfield);
	for(int i=0; i<numFragility; i++){
		outputFile << setw(5) << i << "   ";
	 outputFile << setw(15)<< setprecision(3) << (*Fragility)(i) <<"\n";
	}
	outputFile << "\n";
	outputFile << "numTimePoints.............." << numTimePoints << "\n"; 
	outputFile << "\n";
	for(int i=0;i<numTimePoints;i++){
		outputFile << setw(5) << i <<"     ";
		outputFile << setw(15) << setprecision(3) <<(*timepoints)(i)<<"\n";
		outputFile.flush();
	}
	///// check time /////
	int numSteps=checkTimePoints();

	///// allocate Results
	RandomVibrationSimulatorResult** theResults;
	theResults=new RandomVibrationSimulatorResult*[numLsf];
	for(int i=0; i<numLsf; i++) theResults[i]=
		new RandomVibrationSimulatorResult(i,numTimePoints,numFragility,eps);
				
	//// construct GFunEachStepEvaluator /////
	GFunEachStepEvaluator* theGFunEachStepEvaluator;
	theGFunEachStepEvaluator= new GFunEachStepEvaluator
						 (theTclInterp, theReliabilityDomain, theDomain,
						  numSteps,print);
	Vector* PFthreshold=theGFunEachStepEvaluator->getPFthershold();
	theGFunEvaluator->setGFunEachStepEvaluator(theGFunEachStepEvaluator);
	theGFunEvaluator->setNsteps(numSteps);

	int allconv,iconv,lsf,jpt,kpt,ifrag,iii,istep,kconv;
	Matrix* theLSFValues;
	Matrix* theLSFConv;
	double threshold0,response,frag,conv,threshold;

	int checkstep=checkinterval;
	int ibreak=0;
	int result;

	int** minstep;
	minstep = new int*[numLsf];
	for( lsf=0; lsf<numLsf; lsf++)	minstep[lsf]=new int[numFragility];
	for( lsf=0; lsf<numLsf; lsf++){
		for(ifrag=0; ifrag<numFragility; ifrag++)minstep[lsf][ifrag]=0;
	}
	int iisim=0;
	for(int isim=1;isim<=MaxSim;isim++){
		opserr << "trial " << isim <<"\n";
		///// generate random vector /////
		this->generateRV();
		///// dynamic response analysis /////
		bool FEconvergence = true;
		result = theGFunEvaluator->runAnalysis();
//		Vector* samplExc=new Vector(numSteps);
//		(*samplExc)=theRandomProcess->getExcitation(numSteps, delta);
//		ofstream outputSample2;
///		outputSample2.open("samplexc.txt", ios::out);
//		outputSample2.setf(ios::right);
//		outputSample2.setf(ios::scientific, ios::floatfield);
//		for( int k=0; k<numSteps; k++){
//			outputSample2<<setw(15)<<setprecision(5)<<(*samplExc)(k);
//			outputSample2<<"\n";
//		}
//		outputSample2.flush();
		if (result < 0) {
			// In this case a failure happened during the analysis
			// Hence, register this as failure
			FEconvergence = false;
		}
		if(FEconvergence){
			iisim++;
			///// extract and check 
			theLSFValues=theGFunEvaluator->getEachStepResult();
			theLSFConv=theGFunEvaluator->getEachStepConvFlag();
			allconv=1;
			for (lsf=0;lsf<numLsf; lsf++){
				threshold0=(*PFthreshold)(lsf);
				for(ifrag=0; ifrag<numFragility; ifrag++){
					frag=(*Fragility)(ifrag);
					threshold=threshold0*frag;
					iii=0;
					iconv=0;
					for(istep=0; istep<numSteps; istep++){
						response=(*theLSFValues)(lsf,istep);
						conv=(*theLSFConv)(lsf,istep);
			 			if(!twoside){
							if(response>threshold)iii=istep+1;
//							if(response>threshold0||conv!=0.0)iii=istep+1;
						}else{
//							if(response>threshold0||response<-threshold||conv!=0.0){
							if(response>threshold||response<-threshold){
								iii=istep+1;
							}
						}
						if(iii!=0){
							kpt=0;
							for(jpt=0; jpt<numTimePoints; jpt++){
								if(anaSteps[jpt]>=iii) kpt=jpt;
								if(kpt!=0) break;
							}
							if(kpt!=0){
								iconv=1;
								for(jpt=kpt; jpt<numTimePoints; jpt++){
//									kconv=theResults[lsf]->updateq(jpt, ifrag, isim, 1.0);
									kconv=theResults[lsf]->updateq(jpt, ifrag, 1.0);
									iconv*=kconv;
								}
							}
						}
						if(iii!=0) break;
					}
					allconv*=iconv;
				}
			}
			if(iisim==checkstep){
 				allconv=1;
				iconv=1;
				for (lsf=0;lsf<numLsf; lsf++){
					for(ifrag=0; ifrag<numFragility; ifrag++){
						for(jpt=0; jpt<numTimePoints; jpt++){
					 		kconv=theResults[lsf]->checkconvergence(jpt, ifrag, iisim);
			  				iconv*=kconv;
						}
					}
				}
				allconv*=iconv;
				if(allconv==1) ibreak=1;
				outputFile<<"\n";
				outputFile<<"----- result at simulation"<<iisim<<" -----\n"; 
				for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print1(outputFile);
				outputFile.flush();
				checkstep+=checkinterval;
			}
			if(ibreak==1) break;
		}
	}

	///// simulation ends /////
    outputFile<<"\n";
	outputFile<<"==== simulation results =====\n"; 
	outputFile<<"\n";
	for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print2(outputFile);

	for( lsf=0; lsf<numLsf; lsf++){delete [] minstep[lsf]; minstep[lsf]=0;}
	delete [] minstep; 
	minstep=0;
	for(int i=0; i<numLsf; i++){ delete theResults[i]; theResults[i]=0; }
	delete [] theResults;
	theResults=0;
	theGFunEvaluator->inactivateGFunEachStepEvaluator();
	delete theGFunEachStepEvaluator;
	theGFunEachStepEvaluator=0;
}
int NonStatRandomVibrationSimulation::checkTimePoints()
{
	////// preparation for simulation //////////
	double Tmax=(*timepoints)(numTimePoints-1);
	int maxStep=(int)((Tmax+0.1*delta)/delta);
	
	double pTmax=0.0;
	int numTotalPulses=theRandomProcess->getNumTotalPulse();
	for(int i=0; i<numTotalPulses; i++){
		double tmp=theRandomProcess->getPulseTime(i+1);
		if(tmp>pTmax) pTmax=tmp;
	}
	if(pTmax<Tmax){
		opserr<<"ERROR NonStatRandomVibrationSimulation\n";
		opserr<<"number of random pulses is not large enough\n";
		opserr<<"to describe the random process up to the maximum time points\n";
		opserr<<"Tmax(timpoints)  "<<Tmax<<"\n";
		opserr<<"pTmax(randomProcess)  "<<pTmax<<"\n";
		exit(-1);
	}
	return maxStep;
}
void NonStatRandomVibrationSimulation::samplingInstantaneousSimulation()
{
	opserr<<"ERROR\n";
	opserr<<"NonStatRandomVibrationSimulation::samplingInstantaneousSimulation\n";
	opserr<<"This function is not implemented yet\n";
	exit(-1);
}
void NonStatRandomVibrationSimulation::samplingFisrtpassageSimulation()
{
	opserr<<"ERROR\n";
	opserr<<"NonStatRandomVibrationSimulation::samplingFisrtpassageSimulation\n";
	opserr<<"This function is not implemented yet\n";
	exit(-1);
}
