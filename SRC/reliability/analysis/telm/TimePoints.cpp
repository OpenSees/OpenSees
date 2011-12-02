/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/TimePoints.cpp,v $

#include <TimePoints.h>
TimePoints::TimePoints(double passedStartTime,
					   double passedEndTime,
					   double passedInterval,
					   double passedDelta,
					   bool passedprint)
{
	delta=passedDelta;
	StepsToStart=(int)((passedStartTime+0.01*delta)/delta);
	StepsToEnd=(int)((passedEndTime+0.01*delta)/delta);
	SamplingFrequency=(int)((passedInterval+0.01*delta)/delta);
	TmPts=0;
	anaStps=0;
	OrdAna=0;
	print=passedprint;
	if(print){
		output.open("TimePoints.txt", ios::out);
		output << "\n";
		output << "TimePoints::TimePoints\n";
		output << "\n";
		output << "StepsToStart " << StepsToStart << "\n"; 
		output << "StepsToEnd " << StepsToEnd << "\n"; 
		output << "SamplingFrequency " << SamplingFrequency << "\n"; 
		output << "delta " << delta << "\n"; 
		output.flush();
	}
	MakeTimePoints();
} 
TimePoints::~TimePoints()
{
	if(TmPts !=0 ) {delete TmPts; TmPts=0;}
	if(anaStps !=0 ){ delete [] anaStps; anaStps=0;}
	if(OrdAna !=0 ){ delete [] OrdAna; OrdAna = 0; }
} 
void TimePoints::MakeTimePoints()
{
  	NumPoints=1;
	int tmpSteps=StepsToStart;
	int iflag=0;
	while(iflag==0){
		NumPoints++;
		tmpSteps+=SamplingFrequency;
		if(tmpSteps>=StepsToEnd) iflag=1;
	}

	if(TmPts !=0 ) {delete TmPts; TmPts=0;}
	if(anaStps !=0 ){ delete [] anaStps; anaStps=0;}
	if(OrdAna !=0 ){ delete [] OrdAna; OrdAna = 0; }
  	TmPts = new Vector(NumPoints);
	anaStps = new int[NumPoints];
	OrdAna = new int[NumPoints];
	if(TmPts == 0 || anaStps == 0 || OrdAna == 0 ) {
		opserr << " memory allocation error MakeTimePoints \n";
		exit(-1);
	}
		// construct time points //
	for(int i=0; i<NumPoints; i++){ 
		tmpSteps=StepsToStart+i*SamplingFrequency;
		if(tmpSteps>StepsToEnd) tmpSteps=StepsToEnd; 
		(*TmPts)(i)=delta*(double)tmpSteps;
		anaStps[i]=tmpSteps;
		OrdAna[i]=i;
	}
	if(print){
	output << "\n\n";
	output << " numPoints............................." <<NumPoints<< "\n";
	output << " dt...................................." <<delta<< "\n";
	output << "\n\n";
	output << "   ID        Time     Steps     Order\n";
	output.setf( ios::fixed, ios::floatfield );
	for(int i=0; i<NumPoints; i++){
		output<<setw(5)<<i;
		output<<setw(12)<<setprecision(3)<<(*TmPts)(i);	
		output<<setw(10)<<anaStps[i];
		output<<setw(10)<<OrdAna[i];
		output<<"\n";
	}
	output.flush();
	}
}
int TimePoints::MakeOrder(double MirTime)
{
	Vector* Tafter = 0;
	Tafter =new Vector(NumPoints);
	Vector* Tbefore = 0;
	Tbefore = new Vector(NumPoints);
	int* Ibefore = 0;
	Ibefore = new int[NumPoints];
	int* Iafter = 0;
	Iafter =new int[NumPoints];

	int Nbefore=0;
	int Nafter=0;
	int Nequal=0;
	int Ieqaul=0;

	for(int i=0;i<NumPoints;i++) OrdAna[i]=0;
	for(int i=0;i<NumPoints;i++){
		if((*TmPts)(i)==MirTime) {
			Nequal=Nequal+1;
			Ieqaul=i;
		}else if((*TmPts)(i)>MirTime){
		    (*Tafter)(Nafter)=(*TmPts)(i); 
		    Iafter[Nafter]=i;
		    Nafter=Nafter+1;
		}else{
		    Nbefore=Nbefore+1;
		    (*Tbefore)(Nbefore)=(*TmPts)(i);
		    Ibefore[Nbefore]=i;
		}
	}
	for(int i=0; i<Nbefore-1; i++){
		for( int j=i+1; j<Nbefore; j++){
			if((*Tbefore)(i)<(*Tbefore)(j)) {
				double twork=(*Tbefore)(i);
				int iwork=Ibefore[i];
				(*Tbefore)(i)=(*Tbefore)(j);
				Ibefore[i]=Ibefore[j];
				(*Tbefore)(j)=twork;
				Ibefore[j]=iwork;
			}
		}
	}
	for(int i=0; i<Nafter-1; i++){
		for( int j=i+1; j<Nafter; j++){
			if((*Tafter)(i)>(*Tafter)(j)) {
				double twork=(*Tafter)(i);
				int iwork=Iafter[i];
				(*Tafter)(i)=(*Tafter)(j);
				Iafter[i]=Iafter[j];
				(*Tafter)(j)=twork;
				Iafter[j]=iwork;
			}
		}
	}

	for(int i=0; i<Nbefore; i++){
	  int ii=Ibefore[i];
	  OrdAna[i]=ii;
	}	

	for(int i=0; i<Nafter; i++){
	  int ii=Iafter[i];
	  if(i==0) OrdAna[i+Nbefore]=-ii;
	  else OrdAna[i+Nbefore]=ii;
	}

	if(print){
		opserr << "Number of analysis points " <<NumPoints<<"\n";
		opserr << "Mirror time " <<MirTime<<"\n";
		for(int i=0; i<NumPoints; i++){
		opserr << "order"; 
		opserr.width(5);
		opserr << i << "  ";
		opserr.width(5);
		opserr << OrdAna[i] << "  ";
		opserr <<" time ";
		opserr.width(12);
		opserr << (*TmPts)(abs(OrdAna[i])) ;
		opserr << " step ";
		opserr.width(10);
		opserr << anaStps[abs(OrdAna[i])] << "\n";
		}
	}
	if(Tafter!=0){delete Tafter;Tafter=0;}
	if(Tbefore!=0){delete Tbefore;Tbefore = 0;}
	if(Ibefore != 0){delete [] Ibefore;Ibefore = 0;}
	if(Iafter != 0){delete [] Iafter;Iafter = 0;}
	return Ieqaul;
}
int TimePoints::getAnalysisStep(int jj)
{
	int istep=anaStps[jj];
	return istep;
}
int TimePoints::getOrder(int ii)
{
	int jj=OrdAna[ii];
	return jj;
}
int TimePoints::getNumTimePoints(void)
{
	return NumPoints;
}
void TimePoints::printTimePts(ofstream& outputFile)
{
	outputFile<< "\n";
	outputFile<< "number of time points......" << NumPoints << "\n"; 
	outputFile<< "\n";
	outputFile<< "      seq.      step        time\n"; 
	for(int i=0; i<NumPoints; i++){
		outputFile << setw(10) << i;
		outputFile << setw(10) << anaStps[i];
		outputFile << setw(12) << setprecision(5) << (*TmPts)(i);
		outputFile <<"\n";
	}
	outputFile.flush();
}
void TimePoints::initializeOrder(){
	for(int i=0; i<NumPoints; i++)OrdAna[i]=i;
}
