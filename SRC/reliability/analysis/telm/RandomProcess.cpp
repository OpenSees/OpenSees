// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/RandomProcess.cpp,v $

#include <RandomProcess.h>
RandomProcess::RandomProcess(Domain* passedDomain
							,ReliabilityDomain* passedReliabilityDomain,
							 bool passedprint)
{
	theReliabilityDomain=0;
	theDiscreSeriese=0;
	iRVPulse=0; 
	timePulse=0;
	theReliabilityDomain = passedReliabilityDomain;
	print=passedprint;
	if(print){
		output.open("RandomProcess.txt", ios::out);
		output << "\n";
		output << "RandomProcess::RandomProcess\n";
		output << "\n";
	}

	FindRandomProcess(passedDomain);
	AnalyzeProcess();
} 
RandomProcess::~RandomProcess()
{
	if( iRVPulse !=0 ) {
		delete [] iRVPulse;
		iRVPulse = 0;
	}
	if( timePulse !=0 ){
		delete timePulse;
		timePulse=0;
	}

} 
void RandomProcess::AnalyzeProcess()
{
	numTotalPulse = theDiscreSeriese->getNumPulses();

	if(print){
		output << "\n";
		output << "RandomProcess::AnalyzeProcess\n";
		output << "\n";
		output << "  numTotalPulse " << numTotalPulse << "\n";
		output.flush();
	}

	if( iRVPulse !=0 ) {delete [] iRVPulse;iRVPulse = 0;}
	if( timePulse !=0 ){delete timePulse;timePulse=0;}
	iRVPulse = new int [numTotalPulse];
	if(iRVPulse == 0) {
		opserr << " memory allocation error 1 in AnalyzeProcess \n";
		exit(-1);
	}
	timePulse = new Vector(numTotalPulse);
	if(timePulse == 0) {
		opserr << " memory allocation error 2 in AnalyzeProcess \n";
		exit(-1);
	}

	for(int i=0; i<numTotalPulse; i++) iRVPulse[i]=0;

	int numRV=theReliabilityDomain->getNumberOfRandomVariables();

	for(int i=1; i<=numRV; i++){
		int ii=theDiscreSeriese->getPulseSequentialID(i);
		if(ii>=0){
			iRVPulse[ii]=i;
			double time=theDiscreSeriese->getkickInTimes(i);
			(*timePulse)(ii) = time;
		}
	}


	delta_Pulse=(*timePulse)(1)-(*timePulse)(0);

	for(int i=0; i<	numTotalPulse; i++){
		for(int j=i+1; j<numTotalPulse; j++){
			if((*timePulse)(i)>(*timePulse)(j)){
				int iwork = iRVPulse[i];
				double twork = (*timePulse)(i);
				iRVPulse[i]=iRVPulse[j];
				(*timePulse)(i)=(*timePulse)(j);
				iRVPulse[j]=iwork;
				(*timePulse)(j)=twork;
			}
		}
	}

	if(print){
		output << "\n";
		output << " link pulse to RandomVariable \n";
		output << "\n";
		output << "  numTotalPulse " << numTotalPulse << "\n";
		output << "\n";
		for(int i=0; i<	numTotalPulse; i++){
			output << " pulse " << i ;
			output << " randomVariable " << iRVPulse[i]; 
			output << " time " << (*timePulse)(i); 
		}
		output.flush();
	}

	RandomVariable* theRV;
	int irv;
	for (int i=0; i< numTotalPulse; i++ ) {
		irv=iRVPulse[i];
		theRV = theReliabilityDomain->getRandomVariablePtr(irv);
		theRV -> setStartValue(0.0);
		theDiscreSeriese->updateRV(irv,0.0);
	}
	theDiscreSeriese->activateParameter(0);
}
void RandomProcess::FindRandomProcess(Domain* theStrDom)
{
	int Direction;
	NewDiscretizedRandomProcessSeries* theSeriese;
	LoadPatternIter& thePatterns = theStrDom->getLoadPatterns();
    LoadPattern *thePattern;
	UniformExcitation *theUniform;
	int ifound=0;
    while((thePattern = thePatterns()) != 0){
		int loadpatterntag=thePattern->getTag();
		int classtag=thePattern->getClassTag();
		if( classtag==PATTERN_TAG_UniformExcitation){
			theUniform=(UniformExcitation *)thePattern;
			const GroundMotion *theMotion=theUniform->getGroundMotion();
			Direction=theUniform->getDirection();
			const TimeSeries *theAccel=theMotion->getAccelSeries();
			if( theAccel != 0) classtag=theAccel->getClassTag();
			if( classtag == TSERIES_TAG_DiscretizedRandomProcessSeries){
				ifound++;
				theSeriese=(NewDiscretizedRandomProcessSeries*)theAccel;
				opserr<< " Random Excitation Problem \n";
			}
		}
	}
	if(ifound==0 || ifound > 1){
		opserr<< " For Current implementation of outcrossing rate problem \n";
		opserr<< " the problem must include only one uniform excitataion \n"; 
		opserr<< " with the acceleration time series defined by \n";
		opserr<< " discretized random process \n";
		opserr<< " For problems of other types, some adjustment might be required \n";
		opserr<< " Get contact with the developer !!\n";
		exit(-1);
	}
	theDiscreSeriese = theSeriese;
	Direction++;

	if(print){
		output << "\n";
		output << " Direction " << Direction;
		output << "\n";
		output.flush();
	}

	DirExcitation=Direction;
}
int RandomProcess::getNumOfActivePulses(double time)
{
	int i;
	for(i=0; i<numTotalPulse; i++){
 		if(time-(*timePulse)(i) < -1.0e-7) break;;
	}
	return i;
}
double RandomProcess::getPulseTime(int ipulse){
	return (*timePulse)(ipulse-1);
}
double RandomProcess::getFactorSensitivity(double ctime, double ktime){
	double ddd=theDiscreSeriese->getFactorSensitivity(ctime,ktime);
	return ddd;
}
Vector RandomProcess::getExcitation(int Nstep, double delta){
	Vector Excitation(Nstep);
	for(int i=0;i<Nstep;i++){
		double time=delta*(double)(i+1);
		Excitation(i)=theDiscreSeriese->getFactor(time);
	}
//	ofstream outputSample2;
//	outputSample2.open("samplexc.txt", ios::out);
//	outputSample2.setf(ios::right);
//	outputSample2.setf(ios::scientific, ios::floatfield);
//	for( int k=0; k<Nstep; k++){
//		outputSample2<<setw(15)<<setprecision(5)<<Excitation(k);
//		outputSample2<<"\n";
//	}
	return Excitation;
}
int RandomProcess::getRVseqid(int ipulse){

	int irv=iRVPulse[ipulse];
	return irv;
}
