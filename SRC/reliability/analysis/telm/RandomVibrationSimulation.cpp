#include <RandomVibrationSimulation.h>
RandomVibrationSimulation::RandomVibrationSimulation
						 (ReliabilityDomain* passedReliabilityDomain,
						  Domain* passedDomain,
						  FunctionEvaluator* passedGFunEvaluator,
						  ProbabilityTransformation* passedTransformation,
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
                          Tcl_Interp *passedTclInterp)
{
	theReliabilityDomain = passedReliabilityDomain;
	theDomain=passedDomain;
	numRV	 = theReliabilityDomain->getNumberOfRandomVariables();
	//numRVPos = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	numRVPos = 0;
	numLsf	 = theReliabilityDomain->getNumberOfLimitStateFunctions();
	theRandomProcess = new RandomProcess(theDomain,theReliabilityDomain,false);
	NumTotalPulse = theRandomProcess->getNumTotalPulse();
	delta_Pulse = theRandomProcess->getDeltaPulse();

	theGFunEvaluator = passedGFunEvaluator;
	delta= theGFunEvaluator->getDt(); 
	theTransformation =passedTransformation;
    theRandomNumberGenerator = new GeneralRandGenerator();
	if (theRandomNumberGenerator==0) {
		opserr << "RandomVibrationSimulationn() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}
	uRV=0;
	xRV=0;
	uRV=new Vector(numRV);
	xRV=new Vector(numRV);
	if (uRV==0||xRV==0) {
		opserr << "RandomVibrationSimulation() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}

	twoside=passedtwoside;
	system=passedsystem;

	MaxSim=passedmaxSim;
	eps=passedeps;
	checkinterval=passedcheckinterval;
	instantaneous=passedinstantaneous;
	firstpassage=passedfirstpassage;

	numFragility=passednFrag;
	Fragility=new Vector(numFragility);
	if (Fragility==0) {
		opserr << "RandomVibrationSimulation() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}
	for(int i=0; i<numFragility; i++){
		(*Fragility)(i)=passedFragMin+passedFragInt*(double)i;
	}

	fileName = new char[256];		   // output file name	
	strcpy(fileName,passedFileName);
	if(passedFileBinary!=0){
		fileBinary = new char[256];		   // output file name	
		strcpy(fileBinary,passedFileBinary);
	}else{
		fileBinary=0;
	}
	theTclInterp=passedTclInterp;


}
RandomVibrationSimulation::~RandomVibrationSimulation()
{
	if(Fragility!=0){ delete Fragility; Fragility=0;}
	if(theRandomProcess!=0){
		delete theRandomProcess;
		theRandomProcess=0;
	}
	if(theRandomNumberGenerator!=0){
		delete theRandomNumberGenerator;
		theRandomNumberGenerator=0;
	}
	if(fileName!=0){ delete [] fileName; fileName=0;}
	if(fileBinary!=0){ delete [] fileBinary; fileBinary=0;}
	if(uRV!=0){delete uRV; uRV=0;}
	if(xRV!=0){delete xRV; xRV=0;}
}
void RandomVibrationSimulation::analyze()
{

	theGFunEvaluator->inactivateSensitivty();

	if(instantaneous!=0){
  		if(instantaneous==1){
			this->crudeInstantaneousSimulation();
		}else{
			this->samplingInstantaneousSimulation();
		}
	}

	if(firstpassage!=0){
		if(firstpassage==1){
			this->crudeFisrtpassageSimulation();
		}else{
			this->samplingFisrtpassageSimulation();
		}
	}
}
void RandomVibrationSimulation::generateRV()
{
	int result;
	for(int j=0;j<numRV;j++) 
		(*uRV)(j)=theRandomNumberGenerator->generate_singleStdNormalNumber();

//	ofstream output;
//	output.open("randomvibrationsimulation.txt", ios::out);
//	for(j=0;j<numRV;j++){
//		output <<  setw(30) << (*uRV)(j) << "\n";
//	}
//	output.close();
    /* FMK
	result = theTransformation->set_u(*uRV);
	if (result < 0) {
		opserr << "NonStatRandomVibrationSimulation::crudeinstantaneousSimulation()";
		opserr << "- could not " << endln
			<< " set the u-vector for xu-transformation. " << endln;
		exit(-1);
	}
	result = theTransformation->transform_u_to_x();
	if (result < 0) {
		opserr << "SamplingAnalysis::analyze() - could not " << endln
			<< " transform u to x. " << endln;
		exit(-1);
	}
	(*xRV) = theTransformation->get_x();
	FMK */
	result = theTransformation->transform_u_to_x(*uRV, *xRV);
	if (result < 0) {
		opserr << "SamplingAnalysis::analyze() - could not " << endln
			<< " transform u to x. " << endln;
		exit(-1);
	}
}

