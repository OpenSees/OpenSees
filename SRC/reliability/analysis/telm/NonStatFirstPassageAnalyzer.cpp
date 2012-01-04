// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NonStatFirstPassageAnalyzer.cpp,v $

#include <NonStatFirstPassageAnalyzer.h>
NonStatFirstPassageAnalyzer::NonStatFirstPassageAnalyzer(ReliabilityDomain* passedReliabilityDomain,
							 FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
							 FunctionEvaluator* passedGFunEvaluator,
							 FOSeriesSimulation* passedFOSeriesSimulation,
							 int passedanalysisType,
							 int passedinterpolationType,
                             bool passedtwoside,
							 bool passedprint)
:FirstPassageAnalyzer(passedReliabilityDomain, passedFindDesignPointAlgorithm,
		      passedGFunEvaluator, passedFOSeriesSimulation, 
		      passedanalysisType,passedtwoside)
{
  interpolationType=passedinterpolationType;
  print=passedprint;
  
  if(print){
    output.open("NonStatFirstPassageAnalyzer.txt", ios::out);
    output << "\n";
    output << "NonStatFirstPassageAnalyzer::NonStatFirstPassageAnalyzer\n";
    output << "\n";
    output << "analysisType"<<analysisType<<"\n";
    output << "numRV"<<numRV<<"\n";
    output << "numrvpos"<<numRVPos<<"\n";
    output << "detla"<<delta<<"\n";
    output << "interpolationType"<<interpolationType<<"\n";
    output.flush();
  }
  FPprob=0;
  covres=0;
  betares=0;
  numSim=0;
}
NonStatFirstPassageAnalyzer::~NonStatFirstPassageAnalyzer()
{
	if(FPprob!=0){delete FPprob; FPprob=0;}
	if(covres!=0){delete covres; covres=0;}
	if(betares!=0){delete betares; betares=0;}
	if(numSim!=0){delete [] numSim; numSim=0;}
}
Vector 
NonStatFirstPassageAnalyzer::componentFisrtPassage
                            (Vector* pudes, Vector* pxdes, Vector* palpha, Vector* phfunc, 
							 double pbeta,double ppf, double panalysistime, int panalysisstep,
							 TimePoints* passedTimePoints,
							 int passedlsf, int passedidfragility,
							 ofstream& outputFile)  
{
	this->setComponentResult(pudes, pxdes,palpha,phfunc,
		                     pbeta,ppf,panalysistime,panalysisstep,
							 passedTimePoints,
						     passedlsf,passedidfragility,false); 

	if(print){
		output << "\n";
		output << "NonStatFirstPassageAnalyzer::componentFisrtPassage\n";
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
	outputFile <<"-------------------------------------------------------------------\n";
	outputFile <<"--- Component First-passage Probability Analysis(nonstationary) ---\n";
	outputFile <<"-------------------------------------------------------------------\n";
	outputFile <<"\n";
	outputFile <<"Lsf......................................"<<passedlsf;
	outputFile <<"Fagility................................."<<passedidfragility;
	outputFile <<"\n";
	outputFile <<"analysisType............................."<<analysisType;
	outputFile <<"twoside.................................."<<twoside;
	outputFile <<"\n";

	if(analysisType==0){
		this->componentFisrtPassage1();
		outputFile <<"--- Integration of up-crossing rates ---\n";
		outputFile <<"\n";
		outputFile <<"        Time"<<"    Probability\n";
		for(int i=0; i<numTimePoints; i++){
			output.setf(ios::fixed, ios::floatfield);
			outputFile << setw(12) << setprecision(2) << (*timepoints)(i);
			output.setf( ios::scientific, ios::floatfield );
			outputFile << setw(15) << setprecision(15) << (*FPprob)(i);
			outputFile<<"\n";
		}
	}else{
		this->componentFisrtPassage2();
		outputFile <<"--- First Order Series System Simulation ---\n";
		outputFile <<"\n";
		outputFile <<"        Time"<<"    Probability";
		outputFile <<"        beta"<<"  numsimulation";
		outputFile <<"       c.o.v"<<"\n";
		for(int i=0; i<numTimePoints; i++){
			output.setf(ios::fixed, ios::floatfield);
			outputFile << setw(12) << setprecision(2) << (*timepoints)(i);
			output.setf( ios::scientific, ios::floatfield );
			outputFile << setw(15) << setprecision(15) << (*FPprob)(i);
			outputFile << setw(15) << setprecision(15) << (*betares)(i);
			outputFile << setw(15) << numSim[i];
			outputFile << setw(15) << setprecision(15) << (*covres)(i);
			outputFile<<"\n";
		}
	}
	return (*FPprob);
}
void
NonStatFirstPassageAnalyzer::componentFisrtPassage2()
{
	opserr<< " NonStatFirstPassageAnalyzer::componentFisrtPassage2 \n";
	opserr<< " This function is not implemented yet\n";
	opserr<< " return zeros\n";

	if(FPprob!=0){ delete FPprob; FPprob=0; }
	FPprob=new Vector(numTimePoints+1);
	if(FPprob==0){
		opserr<<" insufficient memory \n";
	  	opserr<<" NonStatFirstPassageAnalyzer::componentFisrtPassage1\n";
		opserr<<" allocation of FPprob\n";
	}
}
void
NonStatFirstPassageAnalyzer::systemFisrtPassage()
{
	opserr<< " NonStatFirstPassageAnalyzer::systemFisrtPassage \n";
	opserr<< " This function is not implemented yet\n";
	opserr<< " return zeros\n";
}
void
NonStatFirstPassageAnalyzer::componentFisrtPassage1()
{
/*	if(FPprob!=0){ delete FPprob; FPprob=0; }
	FPprob=new Vector(numTimePoints+1);
	if(FPprob==0){
		opserr<<" insufficient memory \n";
		opserr<<" NonStatFirstPassageAnalyzer::componentFisrtPassage1\n";
		opserr<<" allocation of FPprob\n";
	}

	Vector* nupt = new Vector(numTimePoints+1);
	if(nupt==0){
		opserr<<" insufficient memory \n";
		opserr<<" NonStatFirstPassageAnalyzer::componentFisrtPassage1\n";
		opserr<<" allocation of nupt\n";
	}

	for( int i=0;i<=numTimePoints;i++)
		(*nupt)(i)=theOutCrossingResult->getnu(i);

	if(print){
		output << "\n";
		output << "NonStatFirstPassageAnalyzer::componentFisrtPassage1\n";
		output << "\n";
		output << "   seq. ID"<<"           time"<<"             nu\n";
		output.setf( ios::scientific, ios::floatfield );
//		output.setf(ios::fixed, ios::floatfield);
		for(int i=0;i<=numTimePoints;i++){
			output << setw(10) << i;
			output << setw(15) << setprecision(5) << (*timepoints)(i);
			output << setw(15) << setprecision(5) << (*nupt)(i);
			output <<"\n";
		}
		output.flush();
	}

	double timeleft=0.0;
	double nuleft=0.0;
	double timeright;
	double nuright;
	double time0=(*timepoints)(0);
	double nu0=(*nupt)(0);
	double fp;
	double nuintegral=0.0;
	double amp=1.0;
	if(twoside) amp=2.0;

	for(int i=1;i<=numTimePoints;i++){
		if(time0>timeleft && time0< (*timepoints)(i)){
			timeright=time0;
			nuright=nu0;
			nuintegral+=(nuleft+nuright)*(timeright-timeleft)/2.0;
			fp=1.0-exp(-nuintegral);
			(*FPprob)(0)=fp*amp;
			if(print){
				output<<"\n";
				output<<"time0>timeleft && time0< (*timepoints)(i)-1-\n";
				output<<"timeleft "<<timeleft<<" timeright "<<timeright<<"\n";
				output<<"nuleft "<<nuleft<<" nuright "<<nuright<<"\n";
				output<<"nuintegral "<<nuintegral<<" fp "<<fp<<"\n";
				output.flush();
			}
			timeleft=time0;
			nuleft=nu0;
			timeright=(*timepoints)(i);
			nuright=(*nupt)(i);
			nuintegral+=(nuleft+nuright)*(timeright-timeleft)/2.0;
			fp=1.0-exp(-nuintegral);
			(*FPprob)(i)=fp*amp;
			if(print){
				output<<"time0>timeleft && time0< (*timepoints)(i)-2-\n";
				output<<"\n";
				output<<"timeleft "<<timeleft<<" timeright "<<timeright<<"\n";
				output<<"nuleft "<<nuleft<<" nuright "<<nuright<<"\n";
				output<<"nuintegral "<<nuintegral<<" fp "<<fp<<"\n";
				output.flush();
			}
		}else{
			timeright=(*timepoints)(i);
			nuright=(*nupt)(i);
			nuintegral+=(nuleft+nuright)*(timeright-timeleft)/2.0;
			fp=1.0-exp(-nuintegral);
			(*FPprob)(i)=fp*amp;
			if(print){
				output<<"timeleft "<<timeleft<<" timeright "<<timeright<<"\n";
				output<<"nuleft "<<nuleft<<" nuright "<<nuright<<"\n";
				output<<"nuintegral "<<nuintegral<<" fp "<<fp<<"\n";
				output.flush();
			}
		}
		timeleft=timeright;
		nuleft=nuright;
	}
	if(time0 > (*timepoints)(numTimePoints)){
		timeright=time0;
		nuright=nu0;
		nuintegral+=(nuleft+nuright)*(timeright-timeleft)/2.0;
		fp=1.0-exp(-nuintegral);
		(*FPprob)(0)=fp*amp;
	}
	delete nupt;
	nupt=0;
*/
}	

