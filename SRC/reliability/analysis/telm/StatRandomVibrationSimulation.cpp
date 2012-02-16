// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/StatRandomVibrationSimulation.cpp,v $

#include <StatRandomVibrationSimulation.h>
#include <TimePoints.h>
#include <RandomVibrationSimulatorResult.h>
StatRandomVibrationSimulation::StatRandomVibrationSimulation
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
                                 double passedstationarytime,
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
								 bool passedprint,
								 double passedsampleAmp,
								 double passedsampleTime)
//								 Vector *pStartPoint)
:RandomVibrationSimulation
(passedReliabilityDomain, passedDomain, passedGFunEvaluator,passedTransformation,
 passedFragMin, passedFragInt, passednFrag, passedtwoside, passedsystem, passedmaxSim,
 passedcheckinterval, passedeps, passedinstantaneous,
 passedfirstpassage, passedFileName, passedFileBinary,passedTclInterp)
{
	stationaryTime=passedstationarytime;
	print=passedprint;
	sampleAmp=passedsampleAmp;
	sample=false;
	if(sampleAmp!=0.0) {
		sample=true;
		sampleStep=(int)((passedsampleTime+delta*0.1)/delta);
	}
//	startPoint = pStartPoint;
	if(print){
		output.open("StatRandomVibrationSimulation.txt", ios::out);
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
		output << "stationaryTime............." << stationaryTime << "\n"; 
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
	numTimePoints++;
	timepoints=new Vector(numTimePoints);
	anaSteps=new int[numTimePoints];
	(*timepoints)(0)=0.0;
	anaSteps[0]=0;
	
//	int addstep=(int)((stationaryTime+delta*0.1)/delta);
	for(int i=1;i<numTimePoints;i++){
		(*timepoints)(i)=theTimePoints->getAnalysisTime(i-1);
		anaSteps[i]=theTimePoints->getAnalysisStep(i-1);
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
	samplExc=0;
	samplResp=0;
	respMat=0;
	excMat=0;

}

StatRandomVibrationSimulation::~StatRandomVibrationSimulation()
{
	if(timepoints!=0){ delete timepoints; timepoints=0;} 
	if(anaSteps!=0){ delete anaSteps; anaSteps=0;} 
}
void StatRandomVibrationSimulation::crudeInstantaneousSimulation()
{
/*	outputFile.open(fileName, ios::out);
	outputFile<<"\n";
	outputFile<<"###########################################################\n";
	outputFile<<"#####    RandomVibraionSimulation(stationary)         #####\n";
	outputFile<<"#####    Instantaneous failure probability            #####\n";
	outputFile<<"###########################################################\n";
	outputFile<<"\n";
	outputFile<< "numRV......................" << numRV  << "\n"; 
	outputFile<< "numRVPos..................." << numRVPos << "\n"; 
	outputFile<< "numLSF....................." << numLsf << "\n"; 
	outputFile<< "delta......................" << delta << "\n"; 
	outputFile<< "delta_Pulse................" << delta_Pulse << "\n"; 
	outputFile<< "stationaryTime............." << stationaryTime << "\n"; 
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
	for(i=0;i<numTimePoints;i++){
		outputFile << setw(5) << i <<"     ";
		outputFile << setw(15) << setprecision(3) <<(*timepoints)(i)<<"\n";
		outputFile.flush();
	}
	///// check time /////
	int numTimePoints0=numTimePoints;
	numTimePoints=1;
	int numSteps=checkTimePoints();
	
	bool systemAna=false;
	if(system&&numLsf>1)systemAna=true;

	///// allocate Results
	RandomVibrationSimulatorResult** theResults;
	theResults=new RandomVibrationSimulatorResult*[numLsf];
	for(i=0; i<numLsf; i++) theResults[i]=
		new RandomVibrationSimulatorResult(i,numTimePoints,numFragility,eps);

	if(systemAna){
		RandomVibrationSimulatorResult* theSysResults;
		theSysResults=
			new RandomVibrationSimulatorResult(numLsf+1,numTimePoints,numFragility,eps);
	}

				
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

	int nsim=0;
	for(int isim=1;isim<=MaxSim;){
		nsim++;
		opserr << "trial " << isim <<"\n";
		///// generate random vector /////
		this->generateRV();

		///// dynamic response analysis /////
		bool FEconvergence = true;
		result = theGFunEvaluator->runGFunAnalysis(*xRV);
		if (result < 0) {
			// In this case a failure happened during the analysis
			// Hence, register this as failure
			FEconvergence = false;
		}
//		} else {
			///// extract and check 
		theLSFValues=theGFunEvaluator->getEachStepResult();
		theLSFConv=theGFunEvaluator->getEachStepConvFlag();
		allconv=1;
		for (lsf=0;lsf<numLsf; lsf++){
			threshold0=(*PFthreshold)(lsf);
			for(jpt=0; jpt<numTimePoints; jpt++){
				if((*theLSFConv)(lsf,anaSteps[jpt]-1)!=0.0){
					response=(*theLSFValues)(lsf,anaSteps[jpt]-1);
		 			for(ifrag=0; i<numFragility; i++){
						frag=(*Fragility)(ifrag);
						qvalue=0.0;
						if(response>(frag*threshold0)) qvalue=1.0;
						if(twoside&&(response<((frag*threshold0)*(-1.0)))) qvalue=1.0;
//							iconv=theResults[lsf]->updateq(jpt, ifrag, isim, qvalue);
						iconv=theResults[lsf]->updateq(jpt, ifrag, nsim, qvalue);
						allconv*=iconv;
					}
				}else{
					for(ifrag=0; i<numFragility; i++){
						qvalue=1.0;
//							iconv=theResults[lsf]->updateq(jpt, ifrag, isim, qvalue);
						iconv=theResults[lsf]->updateq(jpt, ifrag, nsim, qvalue);
						allconv*=iconv;
					}
				}
			}
		}
//		if(isim==checkstep){
		if(nsim==checkstep){
			if(allconv==1) ibreak=1;
			outputFile<<"\n";
			outputFile<<"----- result at simulation"<<isim<<" -----\n"; 
			for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print1(outputFile);
			outputFile.flush();
			checkstep+=checkinterval;
		}
		if(ibreak==1) break;
		numTimePoints=numTimePoints0;
	}

	///// simulation ends /////
/*	outputFile<<"\n";
	outputFile<<"==== simulation results =====\n"; 
	outputFile<<"\n";
	for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print2(outputFile);

	for(i=0; i<numLsf; i++){ delete theResults[i]; theResults[i]=0; }
	delete [] theResults;
	theResults=0;
	theGFunEvaluator->inactivateGFunEachStepEvaluator();
	delete theGFunEachStepEvaluator;
	theGFunEachStepEvaluator=0;
*/
}

void StatRandomVibrationSimulation::crudeFisrtpassageSimulation()
{
	outputFile.open(fileName, ios::out);
	outputFile<<"\n";
	outputFile<<"###########################################################\n";
	outputFile<<"#####    RandomVibraionSimulation(stationary)         #####\n";
	outputFile<<"#####    First Passage failure probability            #####\n";
	outputFile<<"###########################################################\n";
	outputFile<<"\n";
	outputFile<< "numRV......................" << numRV  << "\n"; 
	outputFile<< "numRVPos..................." << numRVPos << "\n"; 
	outputFile<< "numLSF....................." << numLsf << "\n"; 
	outputFile<< "delta......................" << delta << "\n"; 
	outputFile<< "delta_Pulse................" << delta_Pulse << "\n"; 
	outputFile<< "stationaryTime............." << stationaryTime << "\n"; 
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
		outputFile << setw(5) << anaSteps[i] <<"     ";
		outputFile << setw(15) << setprecision(3) <<(*timepoints)(i)<<"\n";
		outputFile.flush();
	}
	///// check time /////
	int numSteps=checkTimePoints();
	// number of steps to be analyzed //
	ifsample=0;
    if(sample){
		if(samplExc!=0){ delete samplExc; samplExc=0; } 
		if(samplResp!=0){ delete samplResp; samplResp=0; } 
		samplExc=new Vector(numSteps);
		samplResp=new Vector(numSteps);
		samplRec.open("work_outcross.bin" , ios::binary);
		nsample=0;
		double difmin=1.0e10;
		for(int ifrag=0; ifrag<numFragility; ifrag++){
			double frag=(*Fragility)(ifrag);
			double dif=fabs(frag-sampleAmp);
			if(dif<difmin){
				difmin=dif;
				ifsample=ifrag;
			}
		}

	}

	///// allocate Results
	RandomVibrationSimulatorResult** theResults=0;
	RandomVibrationSimulatorResult* theSysResults=0;
	theResults=new RandomVibrationSimulatorResult*[numLsf];
	for(int i=0; i<numLsf; i++) theResults[i]=
		new RandomVibrationSimulatorResult(i,numTimePoints,numFragility,eps);

	bool systemAna=false;
	if(system&&numLsf>1)systemAna=true;
	if(systemAna){
		theSysResults=
			new RandomVibrationSimulatorResult(numLsf+1,numTimePoints,numFragility,eps);
	}
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
	for( lsf=0; lsf<numLsf; lsf++){
		minstep[lsf]=new int[numFragility];
	}
	for( lsf=0; lsf<numLsf; lsf++){
		for(ifrag=0; ifrag<numFragility; ifrag++)minstep[lsf][ifrag]=0;
	}
	int convmax;
	int nsim=0;
	bool count;
	for(int isim=1;isim<=MaxSim;isim++){
		opserr << "trial " << isim <<"\n";
			///// generate random vector /////
		this->generateRV();
		///// dynamic response analysis /////
		bool FEconvergence = true;
//		xRV=startPoint;
		result = theGFunEvaluator->runAnalysis();
//		(*samplExc)=theRandomProcess->getExcitation(numSteps, delta);
//		ofstream outputSample2;
//		outputSample2.open("samplexc.txt", ios::out);
//		outputSample2.setf(ios::right);
//		outputSample2.setf(ios::scientific, ios::floatfield);
//		for( int k=0; k<numSteps; k++){
//			outputSample2<<setw(15)<<setprecision(5)<<(*samplExc)(k);
//			outputSample2<<"\n";
//		}

		convmax=numSteps;
		theLSFConv=theGFunEvaluator->getEachStepConvFlag();
		theLSFValues=theGFunEvaluator->getEachStepResult();
		count=true;
	   	if (result < 0) {
			// In this case a failure happened during the analysis
			// Hence, register this as failure
			FEconvergence = false;
			convmax=0;
			bool tbreak=false;
			count=true;
			for(istep=0; istep<numSteps; istep++){
				conv=(*theLSFConv)(0,istep);
				if(conv<1.0) {
					convmax=istep;
					for (lsf=0;lsf<numLsf; lsf++){
	 					double respmax=(*theLSFValues)(lsf,istep);
						double thresholdmax=(*PFthreshold)(lsf)*(*Fragility)(numFragility-1);
						if(!twoside){
		  					if(respmax<=thresholdmax) count=false;
						}else{
							if(respmax<=thresholdmax||respmax>=-thresholdmax) count=false;
						}
					}
					tbreak=true;
				}
				if(tbreak) break;
			}

		}
	 	if(count&&convmax>=stationaryStep){
			nsim++;
			///// extract and check 
			theLSFValues=theGFunEvaluator->getEachStepResult();
//			ofstream output;
//			output.open("randomvibrationsimulation.txt", ios::out);
//			for(int j=0;j<numSteps;j++){
//				output <<  setw(30) << (*theLSFValues)(0,j) << "\n";
//			}
//			output.close();
			for(ifrag=0; ifrag<numFragility; ifrag++){
				frag=(*Fragility)(ifrag);
				int iexceed=0;
				int jptmin=-999;
				int jptstart;
				for (lsf=0;lsf<numLsf; lsf++){
					threshold0=(*PFthreshold)(lsf);
					threshold=threshold0*frag;
					iii=-999;
//					for(istep=stationaryStep; istep<numSteps; istep++){
					for(istep=stationaryStep; istep<convmax; istep++){
		   				response=(*theLSFValues)(lsf,istep);
						conv=(*theLSFConv)(lsf,istep);
						if(!twoside){
							if(response>threshold){
			  		 			iii=istep-stationaryStep;
								iexceed++;
							}
						}else{
							if(response>threshold||response<-threshold){
								iii=istep-stationaryStep;
								iexceed++;
							}
						}
				 		if(iii>=0){
							/// for check print ///
					  		if(sample && ifrag==ifsample && lsf==0 && iii<=sampleStep) {
								(*samplExc)=theRandomProcess->getExcitation(numSteps, delta);
								for(int j=0;j<numSteps;j++){
									(*samplResp)(j)=(*theLSFValues)(0,j);
								}
								nsample++;
								sampletofile(numSteps, samplRec);
							}
 							kpt=0;
							for(jpt=0; jpt<numTimePoints; jpt++){
								if(anaSteps[jpt]>=iii) kpt=jpt+1;
								if(kpt!=0) break;
							}
							if(kpt!=0){
								jptstart=kpt-1;
								if(jptmin<0) jptmin=jptstart;
								if(jptmin>=0&&jptstart<jptmin) jptmin=jptstart;
								for(jpt=kpt-1; jpt<numTimePoints; jpt++){
							 		kconv=theResults[lsf]->updateq(jpt, ifrag, 1.0);
								}
							}
						}
						if(iii>=0) break;
						if(conv<1.0) break;
					}
					allconv*=iconv;
				}
			  	if(systemAna){
					if(jptmin>=0){
						for(jpt=jptmin; jpt<numTimePoints; jpt++){
							kconv=theSysResults->updateq(jpt, ifrag, 1.0);
						}
					}
				}
//				if(iexceed==0) break;
			}
		}
		if(nsim==checkstep){
 			allconv=1;
			iconv=1;
			for (lsf=0;lsf<numLsf; lsf++){
				for(ifrag=0; ifrag<numFragility; ifrag++){
					for(jpt=0; jpt<numTimePoints; jpt++){
					 	kconv=theResults[lsf]->checkconvergence(jpt, ifrag, nsim);
						iconv*=kconv;
					}
				}
			}
			allconv*=iconv;
			if(systemAna){
				iconv=1;
				for(jpt=0; jpt<numTimePoints; jpt++){
					kconv=theSysResults->checkconvergence(jpt, ifrag, nsim);
					iconv*=kconv;
				}
				allconv*=iconv;
			}
			if(allconv==1) ibreak=1;
			outputFile<<"\n";
			outputFile<<"----- result at simulation"<<nsim<<" -----\n"; 
			for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print1(outputFile);
			outputFile.flush();
			checkstep+=checkinterval;
		}
		if(ibreak==1) break;
	}

	///// simulation ends /////
	outputFile<<"\n";
	outputFile<<"==== simulation results =====\n"; 
	outputFile<<"\n";
	outputFile.flush();
	for (lsf=0;lsf<numLsf; lsf++)theResults[lsf]->print2(outputFile);
	if(systemAna){
		outputFile<<"\n";
		outputFile<<"==== simulation results (system) =====\n"; 
		outputFile<<"\n";
		theSysResults->print2(outputFile);
	}
	outputFile.flush();

	if(sample){
		if(nsample>=50) nsample=50;
		respMat=new Matrix(nsample,numSteps);
		excMat=new Matrix(nsample,numSteps);
		samplRead.open("work_outcross.bin" , ios::binary);
		samplRead.seekg(0);
		for (int j=0; j<nsample; j++){
			samplefromfile(numSteps, samplRead);
			for( int k=0; k<numSteps; k++){
				(*respMat)(j,k)=(*samplResp)(k);
				(*excMat)(j,k)=(*samplExc)(k);
			}
		}
		ofstream outputSample1;
		ofstream outputSample2;
		outputSample1.open("samplresp.txt", ios::out);
		outputSample2.open("samplexc.txt", ios::out);
		outputSample1.setf(ios::right);
		outputSample2.setf(ios::right);
		outputSample1.setf(ios::scientific, ios::floatfield);
		outputSample2.setf(ios::scientific, ios::floatfield);
		for( int k=0; k<numSteps; k++){
			for (int j=0; j<nsample; j++){
				outputSample1<<setw(15)<<setprecision(5)<<(*respMat)(j,k);
				outputSample2<<setw(15)<<setprecision(5)<<(*excMat)(j,k);
			}
			outputSample1<<"\n";
			outputSample2<<"\n";
		}
		outputSample2.flush();
		outputSample1.flush();
 		delete respMat; respMat=0;
		delete excMat; excMat=0;
	}

	for( lsf=0; lsf<numLsf; lsf++){delete [] minstep[lsf]; minstep[lsf]=0;}
	delete [] minstep; 
	for(int i=0; i<numLsf; i++){ delete theResults[i]; theResults[i]=0; }
	delete [] theResults;
	theResults=0;
	theGFunEvaluator->inactivateGFunEachStepEvaluator();
	delete theGFunEachStepEvaluator;	
	theGFunEachStepEvaluator=0;
}
int StatRandomVibrationSimulation::checkTimePoints()
{
	////// preparation for simulation //////////
	double Tmax=(*timepoints)(numTimePoints-1)+stationaryTime;
	stationaryStep=(int)((stationaryTime+0.1*delta)/delta)-1;
	int maxStep=(int)((Tmax+0.1*delta)/delta);
	
	double pTmax=0.0;
	int numTotalPulses=theRandomProcess->getNumTotalPulse();
	for(int i=0; i<numTotalPulses; i++){
		double tmp=theRandomProcess->getPulseTime(i);
		if(tmp>pTmax) pTmax=tmp;
	}
//	if(pTmax<Tmax){
//		opserr<<"ERROR NonStatRandomVibrationSimulation\n";
//		opserr<<"number of random pulses is not large enough\n";
//		opserr<<"to describe the random process up to the maximum time points\n";
//		opserr<<"Tmax(timpoints)  "<<Tmax<<"\n";
//		opserr<<"pTmax(randomProcess)  "<<pTmax<<"\n";
//		exit(-1);
//	}
	return maxStep;
}
void StatRandomVibrationSimulation::samplingInstantaneousSimulation()
{
	opserr<<"ERROR\n";
	opserr<<"NonStatRandomVibrationSimulation::samplingInstantaneousSimulation\n";
	opserr<<"This function is not implemented yet\n";
	exit(-1);
}
void StatRandomVibrationSimulation::samplingFisrtpassageSimulation()
{
	opserr<<"ERROR\n";
	opserr<<"NonStatRandomVibrationSimulation::samplingFisrtpassageSimulation\n";
	opserr<<"This function is not implemented yet\n";
	exit(-1);
}


int
StatRandomVibrationSimulation::sampletofile(int nstep, ofstream& record)
{
	long streampos=record.tellp();
//nt sz=U->Size();
//	record.write(reinterpret_cast<const char*>(&time), sizeof(time));
//
	double val;
	for( int i=0 ; i<nstep; i++){
		val=(*samplExc)(i);
		record.write(
		reinterpret_cast<const char*>(&val), sizeof( val));
	}
	for(int  i=0 ; i<nstep; i++){
		val=(*samplResp)(i);
		record.write(
		reinterpret_cast<const char*>(&val), sizeof( val));
	}
	record.flush();
	return streampos;
}

long
StatRandomVibrationSimulation::samplefromfile(int nstep,ifstream& record)
{
	long streampos=record.tellg();
//	int sz=U->Size();
//	record.read(reinterpret_cast<char*>(&time), sizeof(time));

	double val;
	for( int i=0 ; i<nstep; i++){
		record.read(
		reinterpret_cast<char*>(&val), sizeof( val));
		(*samplExc)(i)=val;
	}
	for( int i=0 ; i<nstep; i++){
		record.read(
		reinterpret_cast<char*>(&val), sizeof( val));
		(*samplResp)(i)=val;
	}
	return streampos;
}














