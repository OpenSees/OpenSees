

// $Revision: 1.3 $
// $Date: 2008-08-27 17:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/CrossingRateAnalyzer.cpp,v $
#include <math.h>
#include <CrossingRateAnalyzer.h>
CrossingRateAnalyzer::CrossingRateAnalyzer(ReliabilityDomain* passedReliabilityDomain,
					   FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
					   FunctionEvaluator* passedGFunEvaluator,
					   GradientEvaluator* passedGradGEvaluator,
					   int passedanalysisType,
					   double passedlittleDt,
  			         bool passedprint)
{
	///// INITIALIZE /////
	numSteps=0;
	numRV=0;
	beta0=0.0;
	betaShifted=0.0;
	beta1=0.0;


	theReliabilityDomain=passedReliabilityDomain;
	theFindDesignPointAlgorithm=passedFindDesignPointAlgorithm;
	theGFunEvaluator=passedGFunEvaluator;
	theGradGEvaluator = passedGradGEvaluator;
	theRandomProcess=0;
	analysisType=passedanalysisType;
	littleDt=passedlittleDt;
	delta=0.0;
	delta_pulse=0.0;
    print=passedprint;
	uDesign0=0;
	uShifted=0;

	if(analysisType==1){
		if(theReliabilityDomain==0){
			opserr << "Need theReliabilityDomain before an "
				   << "crossingrateanalyzer can be created" << endln;
			exit(-1);
		}
		if(theFindDesignPointAlgorithm==0){
			opserr << "Need theFindDesignPointAlgorithm before an "
				   << "crossingrateanalyzer can be created" << endln;
			exit(-1);
		}
		if(theGFunEvaluator==0){
			opserr << "Need GFunEvaluator before an "
				   << "crossingrateanalyzer can be created" << endln;
			exit(-1);
		}
	}

	if(print){
		output.open("CrossingRateAnalyzer.txt", ios::out);
	}
}
CrossingRateAnalyzer::~CrossingRateAnalyzer()
{
	if(uDesign0!=0){ delete uDesign0; uDesign0=0;}
	if(uShifted!=0){ delete uShifted; uShifted=0;}
}
void CrossingRateAnalyzer::clear()
{
	if(uDesign0!=0){ delete uDesign0; uDesign0=0;}
	if(uShifted!=0){ delete uShifted; uShifted=0;}
}
void CrossingRateAnalyzer::setAnalysis
                           (int passedNstep, Vector& udes, 
						    Vector& alphades, double betades)
{
	opserr<< "check--111\n";
	if(theRandomProcess==0){
		opserr<<" theRandomProcess is not set\n";
		opserr<<" CrossingRateAnalyzer::setAnalysis\n";
		opserr<<"\n";
		exit(-1);
	}
	if(delta==0.0){
		opserr<<" delta is not set\n";
		opserr<<" CrossingRateAnalyzer::setAnalysis\n";
		opserr<<"\n";
		exit(-1);
	}
	opserr<< "check--222\n";
	numSteps=passedNstep;
	numRV=udes.Size();
	opserr<< "check--333\n";
	if(uDesign0!=0){ delete uDesign0; uDesign0=0;}
	opserr<< "check--337\n";
	if(uShifted!=0){ delete uShifted; uShifted=0;}
	opserr<< "check--444\n";
	int num=udes.Size();
	uDesign0=new Vector(num);
	opserr<< "check--555\n";
	if(uDesign0==0){opserr << " Insufficient memory /n";
					opserr << " CrossingRateAnalyzer uDesing0 /n";
					exit(-1);}
	opserr<< "check--666\n";
	uShifted=new Vector(num);
	opserr<< "check--777\n";
	if(uShifted==0){opserr << " Insufficient memory /n";
					opserr << " CrossingRateAnalyzer uShifted /n";
					exit(-1);}
	*uDesign0=udes;
	*uShifted=udes;


	opserr<< "check--888\n";
	beta0=betades;
	betaShifted=beta0;
	opserr<< "check--999\n";
    int NactivePulses=
		theRandomProcess->getNumOfActivePulses(delta*(float)numSteps);
	int NumTotalPulse=
		theRandomProcess->getNumTotalPulse();

	Vector* Pulse0=new Vector(NumTotalPulse);
	if(Pulse0==0){opserr << " Insufficient memory /n";
				  opserr << " CrossingRateAnalyzer Pulse0 /n";
				  exit(-1);}

	for(int i=0; i<NactivePulses; i++){
		int iRV=theRandomProcess->getRVID(i);
		(*Pulse0)(i)=udes(iRV-1);
	}

	// Construct Shifted Pulses //
	Vector* ShiftedPulse=new Vector(NumTotalPulse);
	if(ShiftedPulse==0){opserr << " Insufficient memory /n";
						opserr << " CrossingRateAnalyzer ShiftedPulse /n";
					    exit(-1);}

	double dt_small=delta*littleDt;
	(*ShiftedPulse)(0)=(*Pulse0)(0)*(1.0-dt_small/delta_pulse);
	for (int j=1; j<NactivePulses; j++) 
	(*ShiftedPulse)(j)=
		(*Pulse0)(j)-((*Pulse0)(j)-(*Pulse0)(j-1))*dt_small/delta_pulse;
	for (int j=NactivePulses; j<NumTotalPulse; j++) 
		(*ShiftedPulse)(j)=0.0;

	// Construct Shifted design point //
	for(int i=0; i<NumTotalPulse; i++){
		int iRV=theRandomProcess->getRVID(i);
		if(i<NactivePulses) (*uShifted)(iRV-1)=(*ShiftedPulse)(i);
		else (*uShifted)(iRV-1)=0.0;
	}

	betaShifted=(*uShifted).Norm();
	if(beta0<0){
		betaShifted=-betaShifted;
	}

	if(print){
		output<<" CrossingRateAnalyzer::setAnalysis \n";
		output<<"\n";
		output<<" numSteps " << numSteps <<"\n";
		output<<" numRV " << numRV <<"\n";
		output<<" NactivePulses " << NactivePulses <<"\n";
		output<<" NumTotalPulse " << NumTotalPulse <<"\n";
		output<<"\n";
		output<<" Original and Shifted beta \n";
		output<<"\n";
		output<<" Origianl Beta " << setw(15) << setprecision(5); 
		output<<beta0<<"\n";
		output<<" Shifted Beta " << setw(15) << setprecision(5);
		output<<betaShifted<<"\n";
		output<<" InnerProduct " << setw(15) << setprecision(5); 
		double inner=(*uDesign0)^(*uShifted);
		inner/=((*uDesign0).Norm()*(*uShifted).Norm());
		output<< inner  <<"\n";
		output<<"\n";
		output<<" Original and Shifted Pulse \n";
		for(int i=0; i<NumTotalPulse; i++){
			output << setw(10) << i+1;
			output << setw(15) << setprecision(5) << (*Pulse0)(i);
			output << setw(15) << setprecision(5) << (*ShiftedPulse)(i);
			output << "\n";
		}
		output<<"\n";
		output<<" Original and Shifted design point \n";
		for(int i=0; i<numRV; i++){
			output << setw(10) << i+1;
			output << setw(15) << setprecision(5) << (*uDesign0)(i);
			output << setw(15) << setprecision(5) << (*uShifted)(i);
			output << "\n";
		}
		output.flush();
	}
	delete ShiftedPulse; ShiftedPulse=0;
	delete Pulse0;Pulse0=0;
	
}
double CrossingRateAnalyzer::computeRate()
{
	if(uDesign0==0||uShifted==0){
		opserr<< "ERROR in ICrossingRateAnalyzer::computeRate\n";
		opserr<< "uDesign0 and uShifted are not set yet \n";
		opserr<< "Consult with developer \n";
		exit(-1);
	}

	opserr<< "check--1\n";
	double result;
	if(analysisType==1){
		result=computeRate1();
		result/=(delta_pulse*littleDt);
	}else{
		result=computeRate2();
		result/=(delta_pulse*littleDt);
	}
	return result;
}
double CrossingRateAnalyzer::computeRate1()
{
	char string[500];
	static NormalRV aStdNormRV(1,0.0,1.0);

	// Get the 'dgdu' vector from the sensitivity evaluator
	// (The returned matrix contains 'node#' 'dir#' 'dgdu' in rows)

	// This will go away -- MHS 10/7/2011
	//Matrix DgDdispl = theGradGEvaluator->getDgDdispl();
	Matrix DgDdispl(5,5); // So it will compile


	// Add extra term to limit-state function
	// (should add an alternative option where user give additional limit-state functions)
	// (... because this works only when g is linear in u)
    int lsf=theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction* theLSF = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
	int numVel = DgDdispl.noRows();
	double accuSum = 0.0;
	double dt_small=delta*littleDt;

	char *expressionPtr = new char[100];
	for (int j=0; j<numVel; j++) {
		int nodeNumber = (int)DgDdispl(j,0);
		int dofNumber = (int)DgDdispl(j,1);
		double dgduValue = DgDdispl(j,2);
		char expression[100];
		sprintf(expression,"+(%10.8f)*(%8.5f)*\\$ud(%d,%d)",dt_small, dgduValue, nodeNumber, dofNumber);
		strcpy(expressionPtr,expression);
		// Add it to the limit-state function
		// Need a better way to do this -- MHS 10/7/2011
		//theLSF->addExpression(expressionPtr);
	}
	delete [] expressionPtr;
	expressionPtr=0;
	// Inform the user
	strcpy(string,theLSF->getExpression());
	opserr << " ...evaluating G2=" << string << endln;
	bool DSPTfailed = false;

	theGFunEvaluator->setNsteps(numSteps);
	theFindDesignPointAlgorithm->set_u(*uShifted);
	int iresult=theFindDesignPointAlgorithm->findDesignPoint();

	if (iresult < 0){
		opserr << "OutCrossingAnalysis::analyze() - failed while finding the" << endln
	           << " design point for limit-state function number " << 
			   theLSF->getTag() << "." << endln;
		DSPTfailed = true;
	}
	// Zero out the added expression in the limit-state function
	// Need a better way to do this -- MHS 10/7/2011
	//theLSF->removeAddedExpression();

	double pf2;
	if (!DSPTfailed) {
		// Get results from the "find design point algorithm"
		(*uShifted) = theFindDesignPointAlgorithm->get_u();
		// Postprocessing (remember; here is an assumption that the mean point is in the safe domain)
		beta1 = theFindDesignPointAlgorithm->get_beta();
		pf2 = 1.0 - aStdNormRV.getCDFvalue(beta1);
				
//   Post-processing to find parallel system probability
	// should use CorrelatedStandardNormal class
		double	a = -((*uDesign0)^(*uShifted));	// Interval start
		a/=((*uDesign0).Norm()*(*uShifted).Norm());
		double 	b = 0.0;				// Interval end
		double  n_2 = 100;				// Half the number of intervals
		double 	h = b-a;
		double 	fa = functionToIntegrate(a,beta0,beta1);
		double 	fb = functionToIntegrate(b,beta0,beta1);
		double 	sum_fx2j = 0.0;
		double 	sum_fx2j_1 = 0.0;
		for (int j=1;  j<=n_2;  j++) {
			sum_fx2j = sum_fx2j + functionToIntegrate(   (double) (a+(j*2)*h/(2*n_2)) ,beta0, beta1  );
			sum_fx2j_1 = sum_fx2j_1 + functionToIntegrate(   (double)(a+(j*2-1)*h/(2*n_2)) , beta0, beta1  );
		}
		sum_fx2j = sum_fx2j - functionToIntegrate((double)(b),beta0,beta1);
		double integral = h/(2*n_2)/3.0*(fa + 2.0*sum_fx2j + 4.0*sum_fx2j_1 + fb);
		double 	Pmn1 = aStdNormRV.getCDFvalue(beta0)*aStdNormRV.getCDFvalue(-beta1) - integral;
		return Pmn1;
	}else {
//		outputFile << "#  Second limit-state function did not converge.                      #" << endln;
		return -1.0;
	}
	delete [] expressionPtr;
	expressionPtr=0;
}
double CrossingRateAnalyzer::computeRate2()
{
	opserr<< "check--1\n";
	beta1 = -beta0;
	opserr<< "check--2\n";
	// Post-processing to find parallel system probability
//	double	a = (*alpha1)^(*alpha0);
	double	a = (*uDesign0)^(*uShifted);
	opserr<< "check--3\n";
	a/=((*uDesign0).Norm()*(*uShifted).Norm());
	opserr<< "check--4\n";
	a*=-1.0;
	opserr<< "check--5\n";
    double pi = 4.0*atan(1.0);
//		3.14159265358979;
	opserr<< "check--6\n";
	double bbb = exp(-beta1*beta1*0.5);
 	double ccc = asin(a);
	opserr<< "check--7\n";
	double Pmn1 = 1.0/(2.0*pi) * exp(-beta1*beta1*0.5) * (asin(a)+0.5*pi);
	return Pmn1;
}

double
CrossingRateAnalyzer::functionToIntegrate(double rho, double beta1, double beta2)
{
	double result;

	if (fabs(rho-1.0) < 1.0e-8) {
		result = 0.0;
	}
	else {
		double pi = 3.14159265358979;
		result = 1.0/(2.0*pi*sqrt(1.0-rho*rho)) 
			* exp(-(beta1*beta1+beta2*beta2+2.0*rho*beta1*beta2)
			/(2.0*(1.0-rho*rho)));
	}
	return result;
}
// analyze performance function //
void CrossingRateAnalyzer::setRandomProcess(RandomProcess* passedRandomProcess){
	theRandomProcess=passedRandomProcess;
	delta_pulse=theRandomProcess->getDeltaPulse();
}
