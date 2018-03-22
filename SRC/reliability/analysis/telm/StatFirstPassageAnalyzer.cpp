// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/StatFirstPassageAnalyzer.cpp,v $

#include <StatFirstPassageAnalyzer.h>
StatFirstPassageAnalyzer::StatFirstPassageAnalyzer
						  (ReliabilityDomain* passedReliabilityDomain,
						   FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
						   FunctionEvaluator* passedGFunEvaluator,
						   FOSeriesSimulation* passedFOSeriesSimulation,
						   int passedanalysisType,
						   bool passedtwoside,
						   bool passedprint)
:FirstPassageAnalyzer(passedReliabilityDomain, passedFindDesignPointAlgorithm,
					  passedGFunEvaluator, passedFOSeriesSimulation, 
					  passedanalysisType,passedtwoside)
{
	print=passedprint;
	if(print){
		output.open("StatFirstPassageAnalyzer.txt", ios::out);
		output << "\n";
		output << "NonStatFirstPassageAnalyzer::NonStatFirstPassageAnalyzer\n";
		output << "\n";
		output << "analysisType"<<analysisType<<"\n";
		output << "numRV"<<numRV<<"\n";
		output << "numrvpos"<<numRVPos<<"\n";
		output << "detla"<<delta<<"\n";
		output << "twoside"<<twoside<<"\n";
		output.flush();
	}
	FPprob=0;
	covres=0;
	betares=0;
	numSim=0;
}
StatFirstPassageAnalyzer::~StatFirstPassageAnalyzer()
{
	if(FPprob!=0){delete FPprob; FPprob=0;}
	if(covres!=0){delete covres; covres=0;}
	if(betares!=0){delete betares; betares=0;}
	if(numSim!=0){delete [] numSim; numSim=0;}
}
Vector
StatFirstPassageAnalyzer::componentFisrtPassage
					 (Vector* pudes, Vector* pxdes, Vector* palpha, Vector* phfunc, 
					  double pbeta, double ppf, double panalysistime, int panalysisstep,
					  TimePoints* passedTimePoints,
					  int passedlsf, int passedidfragility,
					  ofstream& outputFile) 
{
	this->setComponentResult(pudes, pxdes,palpha,phfunc,
		                     pbeta,ppf,panalysistime,panalysisstep,
							 passedTimePoints,
						     passedlsf,passedidfragility,true); 
	if(print){
		output << "\n";
		output << "StatFirstPassageAnalyzer::componentFisrtPassage\n";
		output << "\n";
		output << "after setComponentResult\n";
		output << "\n";
		output << "numTimePoints "<<numTimePoints<<"\n";
		output << "\n";
		output << "   seq. ID"<<"        Time\n";
//		output.setf( ios::scientific, ios::floatfield );
		output.setf(ios::fixed, ios::floatfield);
		for(int i=0;i<=numTimePoints;i++){
			output << setw(10) << i;
			output << setw(12) << setprecision(2) << (*timepoints)(i);
			output <<"\n";
		}
		output.flush();
	}

	outputFile <<"\n";
	outputFile <<"----------------------------------------------------------------\n";
	outputFile <<"--- Component First-passage Probability Analysis(stationary) ---\n";
	outputFile <<"----------------------------------------------------------------\n";
	outputFile <<"\n";
	outputFile <<"Lsf......................................"<<passedlsf<<"\n";
	outputFile <<"Fagility................................."<<twoside<<"\n";
	outputFile <<"\n";
	outputFile <<"analysisType............................."<<analysisType<<"\n";
	outputFile <<"twoside.................................."<<passedidfragility<<"\n";
	outputFile <<"\n";
	outputFile.flush();

	if(analysisType==0){
		this->componentFisrtPassage1();
		outputFile <<"--- Integration of up-crossing rates ---\n";
		outputFile <<"\n";
		outputFile <<"        Time"<<"    Probability\n";
		for(int i=0; i<=numTimePoints; i++){
			outputFile.setf(ios::fixed, ios::floatfield);
			outputFile << setw(12) << setprecision(2) << (*timepoints)(i);
			outputFile.setf( ios::scientific, ios::floatfield);
			outputFile << setw(15) << setprecision(5) << (*FPprob)(i);
			outputFile<<"\n";
		}
		outputFile.flush();
	}else{
		this->componentFisrtPassage2();
		outputFile <<"--- First Order Series System Simulation ---\n";
		outputFile <<"\n";
		outputFile <<"        Time"<<"    Probability";
		outputFile <<"           beta"<<"  numsimulation";
		outputFile <<"       c.o.v"<<"\n";
		for(int i=0; i<=numTimePoints; i++){
			outputFile.setf(ios::fixed, ios::floatfield);
			outputFile << setw(12) << setprecision(2) << (*timepoints)(i);
			outputFile.setf(ios::scientific, ios::floatfield );
			outputFile << setw(15) << setprecision(5) << (*FPprob)(i);
			outputFile << setw(15) << setprecision(5) << (*betares)(i);
			outputFile << setw(15) << numSim[i];
			outputFile << setw(15) << setprecision(5) << (*covres)(i);
			outputFile<<"\n";
		}
		outputFile.flush();
	}
	return (*FPprob);
}
void
StatFirstPassageAnalyzer::componentFisrtPassage2()
{
 	int numSteps0,iRV,ntemp,i,j;
	int NactivePulse0,lastStep,NactiveMax,iflag,sizeU,numComp;
	double Time0,lastPulseTime,lastTime,compTime,beta0,pf0;
	double beta02,curtime,pfres,cvar,betapt,pulse2,xxx,uOthers2;
    int curNrv,curComp,numSimulation;

	if(theRandomProcess==0){
		opserr<< "fatal in StatFirstPassageAnalyzer::componentFisrtPassage2\n";
		opserr<< "theRandomProcess is not set yet\n";
		exit(-1);
	}

	if(print){
		output.setf( ios::scientific, ios::floatfield );
		output << "\n";
		output << "StatFirstPassageAnalyzer::componentFisrtPassage2\n";
		output << "\n";
		output << " numTotalPulses "<<numTotalPulses<<"\n";
		output << " numOtherRVs "<<numOtherRVs<<"\n";
		output << " delta_Pulse "<<delta_Pulse<<"\n";
		output << "\n";
		output.flush();
	}


	Vector* uOthers=0;
	uOthers=new Vector(numOtherRVs);
	if(uOthers==0){opserr<<" insufficient memory \n";
				   opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				   opserr<<" allocation of uOthers\n";exit(-1);}

	Vector* Pulse0=0;
	Vector* Tpulse0=0;
	Vector* shiftedPulse=0;
	Vector** alphaComp=0;
	Vector** uDesComp=0;
	Vector* betaVec=0;

	numSteps0=analysisstep;
	Time0=analysistime;
	beta0=betadesres;
	pf0=pfdesres;
	beta02=beta0*beta0;
	if(print){
		output << "\n";
		output << "numSteps0 "<<numSteps0<<"\n";
		output << "Time0 "<<Time0<<" beta0 "<<beta0<<"\n";
		output.flush();
	}

     ///////// extract random variable for structures ////////
	Vector uDesing0=(*udesres);
	int* itemp=0;
	itemp = new int[numRV];
	if(itemp==0){opserr<<" insufficient memory \n";
				 opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				 opserr<<" allocation of itemp\n";exit(-1);}

	for(i=0;i<numRV;i++) itemp[i]=1;
	for(i=0; i<numTotalPulses; i++){iRV=theRandomProcess->getRVID(i);
									itemp[iRV-1]=0;}
	ntemp=0;
	for(i=0;i<numRV;i++){
		if(itemp[i]==1){
			(*uOthers)(ntemp)=uDesing0(i);
			ntemp=ntemp+1;
		}
	}
	delete [] itemp;
	itemp=0;

	uOthers2=uOthers->Norm();
	uOthers2*=uOthers2;

	if(ntemp!=numOtherRVs){
		opserr<< "fatal in componentFisrtPassage2\n";
		opserr<< "ntemp != numOtherRVs\n";
		opserr<< "ntemp "<<ntemp<<" numOtherRVs "<<numOtherRVs<<"\n";
		exit(-1);
	}
	if(print && numOtherRVs!=0){
		output <<"\n";
		output <<"other radom variables \n";
		output <<"\n";
		output <<"   seq. ID"<<"     coordinate"<<"\n";
		for(int i=0; i<numOtherRVs; i++){
			output<< setw(10) <<i;
			output<< setw(15) << setprecision(5) << (*uOthers)(i);
			output<<"\n";
		}
		output.flush();
	}

	///////// extract random variable for excitation ////////
    NactivePulse0=theRandomProcess->getNumOfActivePulses(Time0);

	Pulse0=new Vector(NactivePulse0);
	if(Pulse0==0){opserr<<" insufficient memory \n";
			      opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				  opserr<<" allocation of Pulse0\n";exit(-1);}
	Tpulse0=new Vector(NactivePulse0);
	if(Tpulse0==0){opserr<<" insufficient memory \n";
			       opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				   opserr<<" allocation of Tpulse0\n";exit(-1);}

	for(i=0; i<NactivePulse0; i++){
		iRV=theRandomProcess->getRVID(i);
		(*Pulse0)(i)=uDesing0(iRV-1);
		(*Tpulse0)(i)=theRandomProcess->getPulseTime(i+1);
	}
	lastPulseTime=(*Tpulse0)(NactivePulse0-1);

	if(print){
		output <<"\n";
		output <<"input Pulses\n";
		output <<"\n";
		output <<"   seq. ID"<<"     pulse time"<<"     coordinate"<<"\n";
		for(int i=0; i<numOtherRVs; i++){
			output<< setw(10) <<i;
			output<< setw(15) << setprecision(5) << (*Tpulse0)(i);
			output<< setw(15) << setprecision(5) << (*Pulse0)(i);
			output<<"\n";
		}
		output <<"\n";
		output <<" lastPulseTime "<<lastPulseTime<<"\n";
		output.flush();
	}
	////////// problem size ////////////////
	lastStep=theTimePoints->getAnalysisStep(numTimePoints-1);
	lastTime=theTimePoints->getAnalysisTime(numTimePoints-1)+Time0;

	NactiveMax=NactivePulse0;
	iflag=0;
	if(lastPulseTime>=lastTime) iflag=1; 
	while(iflag==0){
		NactiveMax++;
		lastPulseTime+=delta_Pulse;
		if(lastPulseTime>=lastTime) iflag=1; 
	}
	shiftedPulse=new Vector(NactiveMax);
	if(shiftedPulse==0){opserr<<" insufficient memory \n";
					    opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				        opserr<<" allocation of shiftedPulse\n";exit(-1);}
	sizeU=numOtherRVs+NactiveMax;
	numComp=lastStep+1;

	if(print){
		output <<"\n";
		output <<" lastStep "<<lastStep<<"\n";
		output <<" lastTime "<<lastTime<<"\n";
		output <<" sizeU "<<sizeU<<"\n";
		output <<" numComp "<<numComp<<"\n";
		output.flush();
	}

	alphaComp=new Vector*[numComp];
	if(alphaComp==0){opserr<<" insufficient memory \n";
					 opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				     opserr<<" allocation of alphaComp\n";exit(-1);}
	int ialloc=0;
	for(i=0;i<numComp;i++){alphaComp[i]=new Vector(sizeU);if(alphaComp[i]==0) ialloc=1;}
	if(ialloc!=0){opserr<<" insufficient memory \n";
			      opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				  opserr<<" allocation of alphaComp2\n";exit(-1);}

	uDesComp=new Vector*[numComp];
	if(uDesComp==0){opserr<<" insufficient memory \n";
			        opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				    opserr<<" allocation of uDesComp\n";exit(-1);}
	ialloc=0;
	for(i=0;i<numComp;i++){uDesComp[i]=new Vector(sizeU);if(uDesComp[i]==0) ialloc=1;}
	if(ialloc!=0){opserr<<" insufficient memory \n";
			      opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				  opserr<<" allocation of uDesComp2\n";exit(-1);}

	betaVec=new Vector(numComp);
	if(betaVec==0){opserr<<" insufficient memory \n";
				   opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				   opserr<<" allocation of betaVec\n";exit(-1);}

	/////////// construct component //////////
	for( i=0;i<numComp;i++){
		compTime=Time0+delta*(float)i;
		for(j=0;j<numOtherRVs;j++) {
			(*alphaComp[i])(j)=(*uOthers)(j);
			(*uDesComp[i])(j)=(*uOthers)(j);
		}
		shiftPulse(NactivePulse0, Time0, *Pulse0, *Tpulse0,
				   NactiveMax, compTime, *shiftedPulse);
		pulse2=shiftedPulse->Norm();
		pulse2*=pulse2;
		xxx=(beta02-uOthers2)/pulse2;
		xxx=sqrt(xxx);
		for(j=numOtherRVs;j<sizeU;j++){
			(*alphaComp[i])(j)=(*shiftedPulse)(j-numOtherRVs);
			(*uDesComp[i])(j)=xxx*(*shiftedPulse)(j-numOtherRVs);
		}
		for(j=0;j<sizeU;j++)(*alphaComp[i])(j)=(*uDesComp[i])(j);
		alphaComp[i]->Normalize();
		(*betaVec)(i)=beta0;
	}

	if(print){ 
		int* iprint=new int[100];
		int interval=numComp/10;
		int nprint=1;
		iprint[0]=0;
		int istep=0;
		while(nprint>0){
			istep+=interval;
			if(istep<=numComp){
				iprint[nprint]=istep;
				nprint++;
			}else break;
		}
		output<<"\n";
		output<<"--- constructed points -----\n";
		output<<"               ";
		for(i=0;i<nprint;i++)output<<setw(15)<<iprint[i];
		output<<"\n";
		output<<"               ";
		for(i=0;i<nprint;i++) output<<setw(15)<<Time0+delta*(float)iprint[i];
		output<<"\n";
		for(j=0;j<NactiveMax;j++){
			for(i=0;i<nprint;i++) output<<setw(15)<<(*uDesComp[iprint[i]])(j);
			output<<"\n";
		}
		output.flush();
		output<<"\n";
		output<<"--- constructed alpha -----\n";
		output<<"               ";
		for(i=0;i<nprint;i++) output<<setw(15)<<iprint[i];
		output<<"\n";
		output<<"               ";
		for(i=0;i<nprint;i++) output<<setw(15)<<Time0+delta*(float)iprint[i];
		output<<"\n";
		for(j=0;j<NactiveMax;j++){
		  for(i=0;i<nprint;i++) output<<setw(15)<<(*alphaComp[iprint[i]])(j);
		  output<<"\n";
		}
		output.flush();
		delete [] iprint;
		iprint=0;
	}
	if(print){
		output <<"\n";
		output <<"FO series simulation \n";
		output <<"\n";
	}
	theFOSeriesSimulation->setAlphaVec(alphaComp);
	theFOSeriesSimulation->setBetaVec(betaVec);
	theFOSeriesSimulation->setuDesVec(uDesComp);

	if(FPprob!=0) {delete FPprob;FPprob=0;}
	if(covres!=0) {delete covres;covres=0;}
	if(betares!=0) {delete betares;betares=0;}
	if(numSim!=0) {delete[] numSim; numSim=0;}
	FPprob=new Vector(numTimePoints+1);
	covres=new Vector(numTimePoints+1);
	betares=new Vector(numTimePoints+1);
	numSim=new int[numTimePoints+1];
	if(FPprob==0){opserr<<" insufficient memory \n";
			      opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				  opserr<<" allocation of FPprob\n";exit(-1);}
	if(covres==0){opserr<<" insufficient memory \n";
			      opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				  opserr<<" allocation of covres\n";exit(-1);}
	if(betares==0){opserr<<" insufficient memory \n";
			      opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				  opserr<<" allocation of betares\n";exit(-1);}
	if(numSim==0){opserr<<" insufficient memory \n";
			      opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage2\n";
				  opserr<<" allocation of numSim\n";exit(-1);}

//	double pftime0;
//	double betatime0;
//	NormalRV* aStdNormRV = new NormalRV(1,0.0,1.0,0.0);
//	if(!twoside){
//		pftime0=pf0;
//	}else{
//		pftime0=pf0*2.0;
//	}
//	(*FPprob)(0)=pftime0;
//	numSim[0]=0;
//	(*covres)(0)=0.0;
//	betatime0=-aStdNormRV->getInverseCDFvalue(pftime0);
//	(*betares)(0)=betatime0;

	for(i=0;i<=numTimePoints;i++){
		curtime=(*timepoints)(i);
	    curNrv=theRandomProcess->getNumOfActivePulses(curtime);
	    curNrv+=NactivePulse0;
		curNrv+=numOtherRVs;
		theFOSeriesSimulation->setNrv(curNrv);
		curComp=(int)((curtime+0.1*delta)/delta)+1;
		theFOSeriesSimulation->setNcomp(curComp);
		theFOSeriesSimulation->analyze();
		pfres=theFOSeriesSimulation->getpfres();
		cvar=theFOSeriesSimulation->getcvar();
		betapt=theFOSeriesSimulation->getbetares();
		numSimulation=theFOSeriesSimulation->getnumSimulation();
		(*FPprob)(i)=pfres;
		(*covres)(i)=cvar;
	 	(*betares)(i)=betapt;
		numSim[i]=numSimulation;
		if(print){
			output<<"point id "<<i;
			output<<" time "<<curtime;
			output<<" curNrv "<<curNrv<<"\n";
			output<<"results\n";
			output<<" pfres "<<pfres;
			output<<" cvar "<<cvar;
			output<<" numSimulation "<<numSimulation;
			output<<" betares "<<betares<<"\n";
		}
	}

	if(uOthers!=0) {delete uOthers;uOthers=0;}
	if(Pulse0!=0) {delete Pulse0;Pulse0=0;}
	if(Tpulse0!=0) {delete Tpulse0;Tpulse0=0;}
	if(shiftedPulse!=0) {delete shiftedPulse;shiftedPulse=0;}
	if(alphaComp!=0) {
		for(i=0;i<numComp;i++){delete alphaComp[i];alphaComp[i]=0;}
		delete [] alphaComp;
		alphaComp=0;
	}
	if(uDesComp!=0) {
		for(i=0;i<numComp;i++){delete uDesComp[i];uDesComp[i]=0;}
		delete [] uDesComp;
		uDesComp=0;
	}
	if(betaVec!=0) {delete betaVec;betaVec=0;}
}
void StatFirstPassageAnalyzer::shiftPulse
(int Npulse0, double time0, Vector& Pulse0, Vector& PulseTime0,
 int NpulseNew, double timeNew, Vector& shiftedPulse)
{
	int i,ileft,iright,iileft,iiright,icurrent,ibreak;
	double ptime,pleft,pright,ppp,shiftTime,tleft,tright;

	shiftTime=timeNew-time0;
	if(shiftTime<0) {
		opserr<<"Fatal in StatFirstPassageAnalyzer::shiftPulse\n";
		opserr<<"shiftTime must be positive\n";
		opserr<<"timeNew "<<timeNew;
		opserr<<"time0 "<<time0;
		opserr<<"shiftTime "<<shiftTime;
		exit(-1);
	}

	if(shiftTime==0.0){
		for(i=0; i<Npulse0; i++) shiftedPulse(i)=Pulse0(i);
		for(i=Npulse0; i<NpulseNew; i++) shiftedPulse(i)=0.0;
	}else{
		for(i=0; i<NpulseNew; i++){
			ileft=i+1;
			iright=i;
			icurrent=i;
			ptime=delta_Pulse*double(i+1)+1.0e-9;
			if(ptime>=shiftTime) break;
			shiftedPulse(i)=0.0;
		}
		ibreak=0;
		for(i=icurrent; i<NpulseNew; i++){
			iileft=i-ileft;
			iiright=i-iright;
			double Ptime=delta_Pulse*double(i+1)+1.0e-9;
			if(iileft<0){ 
				pleft=0.0;
				tleft=shiftTime;
			}else{
				pleft=Pulse0(iileft);
				tleft=PulseTime0(iileft)+shiftTime;
			}
			if(iiright<Npulse0){ 
				tright=PulseTime0(iiright)+shiftTime;
				pright=Pulse0(iiright);
			}else{
				tright=PulseTime0(Npulse0-1)+shiftTime+delta_Pulse;
				pright=0.0;
				ibreak=1;
			}
			if(ibreak==1) {
				icurrent=i+1;
				break;
			}
			ppp=pleft+(pright-pleft)/(tright-tleft)*(Ptime-tleft);
			shiftedPulse(i)=ppp;
		}
		for(i=icurrent; i<NpulseNew; i++){
			shiftedPulse(i)=0.0;
		}
	}
}
void StatFirstPassageAnalyzer::systemFisrtPassage(void){
	opserr<<"StatFirstPassageAnalyzer::systemFisrtPassage";
	opserr<<"Is not implemented yet\n";
}
void
StatFirstPassageAnalyzer::componentFisrtPassage1()
{
	opserr<<"StatFirstPassageAnalyzer::componentFisrtPassage1";
	opserr<<"Is not implemented yet\n";
//	double statnu=theOutCrossingResult->getnu(0);
//	double statpf=theOutCrossingResult->getpf(0);
//	double aaa=1.0-statpf;
/*	if(print){
		output.setf( ios::scientific, ios::floatfield );
		output << "\n";
		output << "StatFirstPassageAnalyzer::componentFisrtPassage1\n";
		output << "\n";
		output << "statnu"<<statnu<<"\n";
		output << "statpf"<<statpf<<"\n";
		output << "aaa"<<aaa<<"\n";
		output.flush();
	}

	///// allocation of FPprob /////
	if(FPprob!=0){ delete FPprob; FPprob=0; }
	FPprob=new Vector(numTimePoints+1);
	if(FPprob==0){
		opserr<<" insufficient memory \n";
		opserr<<" StatFirstPassageAnalyzer::componentFisrtPassage1\n";
		opserr<<" allocation of FPprob\n";
	}

	double amp=1.0;
	if(twoside) amp=2.0;

	for(int i=0;i<=numTimePoints; i++){
		double curtime=(*timepoints)(i);
		double nuint=statnu*curtime;
		double fp=1.0-aaa*exp(-nuint);
		if(print){
			output.setf( ios::scientific, ios::floatfield );
			output << "\n";
			output << "point "<<i<<" time "<<curtime<<"\n";
			output << "nuint "<<nuint<<" fp "<<fp<<"\n";
			output.flush();
		}
		double ppp=fp*amp;
		if( ppp >= 1.0) ppp=1.0;
		(*FPprob)(i)=ppp;
	}
*/
}	
  
