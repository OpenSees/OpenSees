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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/RandomVibrationSimulatorResult.cpp,v $
#include <math.h>
#include <RandomVibrationSimulatorResult.h>
RandomVibrationSimulatorResult::RandomVibrationSimulatorResult
								(int plsf, int pnumtimepoints,
								 int pnumfragility, double passedeps)
{
	lsf=plsf;
	numTimePoints=pnumtimepoints;
	numFragility=pnumfragility;

	sum_of_q	=0;
	sum_of_qsquared	=0;
	probability	=0;
	cov			=0;
	probability_at_convergence=0;

	sum_of_q	=new Matrix(numTimePoints,numFragility);
	sum_of_qsquared	=new Matrix(numTimePoints,numFragility);
	probability	=new Matrix(numTimePoints,numFragility);
	cov	=new Matrix(numTimePoints,numFragility);
	probability_at_convergence	=new Matrix(numTimePoints,numFragility);

	convergenceFlag=new int*[numTimePoints];
	num_of_simulation_at_convergence=new int*[numTimePoints];
	for(int i=0; i<numTimePoints; i++){
		convergenceFlag[i]=new int[numFragility];
		num_of_simulation_at_convergence[i]=new int[numFragility];
	}
	
	for(int i=0; i<numTimePoints; i++){
		for(int j=0; j<numFragility; j++){
			(*sum_of_q)(i,j)=0.0;
			(*sum_of_qsquared)(i,j)=0.0;
			(*probability)(i,j)=0.0;
			(*cov)(i,j)=0.0;
			(*probability_at_convergence)(i,j)=0;
			convergenceFlag[i][j]=0;
			num_of_simulation_at_convergence[i][j]=0;
		}
	}
	eps=passedeps;
}

RandomVibrationSimulatorResult::~RandomVibrationSimulatorResult()
{
	if(sum_of_q!=0){ delete sum_of_q; sum_of_q=0; }
	if(sum_of_qsquared!=0){ delete sum_of_qsquared; sum_of_qsquared=0; }
	if(probability!=0){ delete probability; probability=0; }
	if(cov!=0){ delete cov; cov=0; }
	if(probability_at_convergence!=0){ delete probability_at_convergence; 
										probability_at_convergence=0; }
	if(convergenceFlag!=0){
		for(int i=0;i<numTimePoints;i++){
			delete [] convergenceFlag[i];
			convergenceFlag[i]=0;
		}
		delete [] convergenceFlag;
		convergenceFlag=0;
	}
	if(num_of_simulation_at_convergence!=0){
		for(int i=0;i<numTimePoints;i++){
			delete [] num_of_simulation_at_convergence[i];
			num_of_simulation_at_convergence[i]=0;
		}
		delete [] num_of_simulation_at_convergence;
		num_of_simulation_at_convergence=0;
	}
}

int 
RandomVibrationSimulatorResult::updateq(int itime, int ifrag, double qvalue)
{
	(*sum_of_q)(itime,ifrag)+=qvalue;
	(*sum_of_qsquared)(itime,ifrag)+=qvalue*qvalue;
	return 0;
}
int 
RandomVibrationSimulatorResult::checkconvergence
(int itime, int ifrag,int numsimulation)
{
	double sumq=(*sum_of_q)(itime,ifrag);	
	double sumq2=(*sum_of_qsquared)(itime,ifrag);	
	double amp=1.0/(double)numsimulation;

	double qbar=0.0;
	double varqbar=0.0;
	double covqbar=0.0;

	if (sumq > 0.0) {
		qbar = amp*sumq;
		varqbar = amp*(amp*sumq2-(sumq*amp)*(sumq*amp));
		if (varqbar< 0.0) varqbar= 0.0;
		covqbar = sqrt(varqbar)/qbar;
	}
	(*probability)(itime,ifrag)=qbar;
 	(*cov)(itime,ifrag)=covqbar;

	if(convergenceFlag[itime][ifrag]==0){
		(*probability_at_convergence)(itime,ifrag)=qbar;
		num_of_simulation_at_convergence[itime][ifrag]=numsimulation;
	}

	int iconv=0;
	if(covqbar>0.0&&covqbar<eps) {
		convergenceFlag[itime][ifrag]=1;
		iconv=1;
	}
	return iconv;
}
void RandomVibrationSimulatorResult::print1(ofstream& output)
{
	output.setf(ios::right);
	output.setf(ios::scientific, ios::floatfield);
	output << "\n";
	output << " ===== LSF "<<lsf<<" =====\n";
	output << "\n";
	output << " time";
	for(int i=0; i< numFragility; i++) {
		output << " ---Fragility " <<setw(10)<<i <<"  ----"; 
	}
	output << "\n";
	output << "     ";
	for(int i=0; i< numFragility; i++) {
		output << "    probability";
		output << "          c.o.v";
	}
	output << "\n";
	for(int i=0; i<numTimePoints; i++){
		output<<setw(5)<< i;
		for(int j=0; j<numFragility; j++){
			output << setw(15) << setprecision(5) << (*probability)(i,j);
			output << setw(15) << setprecision(5) << (* cov)(i,j);
		}
		output << "\n";
	}
}
void RandomVibrationSimulatorResult::print2(ofstream& output)
{
	output.setf(ios::right);
	output.setf(ios::scientific, ios::floatfield);
	output << "\n";
	output << " ===== LSF "<<lsf<<" =====\n";
	output << "\n";
	output << " time";
	for(int i=0; i< numFragility; i++) {
		output << "  --Fragility " <<setw(10)<<i <<"  --  "; 
	}
	output << "\n";
	output << "     ";
	for(int i=0; i< numFragility; i++) {
		output << "    probability";
		output << "      nim/c.o.v";
	}
	output << "\n";
	for(int  i=0; i<numTimePoints; i++){
		output<<setw(5)<< i;
		for(int j=0; j<numFragility; j++){
			if(convergenceFlag[i][j]==1){
				output << setw(15) << setprecision(5) << (*probability_at_convergence)(i,j);
				output << setw(15) << num_of_simulation_at_convergence[i][j];
			}else{
				output << setw(15) << setprecision(5) << (*probability)(i,j);
				output << setw(15) << setprecision(5) << (*cov)(i,j);
			}
		}
		output << "\n";
	}
	output.flush();
}
