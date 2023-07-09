/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//
#include <stdio.h> 
#include <math.h>
#include "SafeGuardedZeroFindingAlgorithm.h"

#include <OPS_Globals.h>

SafeGuardedZeroFindingAlgorithm::SafeGuardedZeroFindingAlgorithm()
{

}

SafeGuardedZeroFindingAlgorithm::~SafeGuardedZeroFindingAlgorithm()
{
	if (x_1 !=0) delete [] x_1;
	if (x_2 !=0) delete [] x_2;
	if (G_1 !=0) delete [] G_1;
	if (G_2 !=0) delete [] G_2;

}
int SafeGuardedZeroFindingAlgorithm::findZeroPoint(double x0)
{
	if (! theSamplingAnalysis->getContribution()) return 0;  //Gu 2007 Dec 23
 

	// note:  1, the point x0 is not used for this algorithm
	//        2, a & Fa must be set before calling. Fa must be positive (safe domain).
	//        3, either b& Fb set,  or FEconvergence is set to be false (unsafe domain, either G<0 or diverge)
	//        4, isFbValid decide whether Fb is converged value;
    
	//create new array to save history
	double * Xs;
	double * Fs;
	int ii;

	if (functionType ==1) {
		Xs = x_1;
	    Fs = G_1;
		ii=ii_1;
	}

	else if (functionType ==2) {
		Xs = x_2;
	    Fs = G_2;
		ii=ii_2;
	}

/////////////////////////////////////////////////////////

	double TOL=variableTol/1000.0;
	

	double x1,x_new, last_x_new;
	double Fx0, Fx1,Fx_new;
	
	if (Fa < -functionTol) {
		opserr<<"error in SafeGuardedZeroFindingAlgorithm::findZeroPoint. Fa<0"<<endln; 
		return -1;
	}

	

	if (FEconvergence) {  //save data
		isFbValid = true; 

		x_new = -Fb*(b-a)/(Fb-Fa)+b;
		if ((x_new<a-TOL)||(x_new<b+TOL)){  
			x_new=(a+b)/2.0;
		}
	}
	else {  // !FEconvergence, bisection
		x_new = (a+b)/2.0;
		// improve  ??????????????????????????????????????????????????????????????????????
//		x_new = a+0.15;
		isFbValid = false; 
	}

	if (this->functionType == 1){
		iterationNum =2; // here is function evaluation number
	}
	else if (this->functionType == 2){
		iterationNum =0; // here is function evaluation number
	}

	last_x_new = x_new+999;
 

	Fx_new=myFunction(x_new);

	iterationNum++;

	if (FEconvergence) {

		if (fabs(Fx_new)<functionTol) {
			zeroPoint = x_new; 
			opserr<<"Return: x_new is:"<<x_new<<" Fx_new is: "<<Fx_new<<endln;
			return 1;
		}

		if (Fx_new * Fa >functionTol*functionTol){ //good
			
			a = x_new;
			Fa = Fx_new;
			
		}
		else if (Fx_new * Fa < -functionTol*functionTol){
			b = x_new;
			Fb = Fx_new;
			isFbValid = true;

		}
		else { 
			opserr<<"wrong bound!"<<endln; return -1;
		}
	}
	else {//diverge
		b=x_new; Fb=-1.0;
		isFbValid = false;
	
	}
	
	
	
	if (maxIterNum>100) {
		opserr<<"warning: SafeGuardedZeroFindingAlgorithm::findZeroPoint may cause member ";
		opserr<<" error since maxIterNum:"<<maxIterNum<<" exceed the limit 100..."<<endln;
	}

	while (theSamplingAnalysis->getContribution() && /*fabs(Fx_new)>functionTol &&*/ iterationNum < maxIterNum && (b-a)>variableTol) {
	
		iterationNum++;

		if (isFbValid ) { // select x0 and x1..
			if(fabs(Fa)>fabs(Fb)+functionTol) {
				x1=b; Fx1=Fb;

			}
			else {
				x1=a; Fx1=Fa;

			}
		}
		else {// ! isFbValid , begin from a
			x1=a; Fx1=Fa;
		}
		//	--------- select best x0, which is closest to x1 ---------
		double  interval =999;
		for(int i=0;i<ii;i++){
			if ((fabs(Xs[i]-x1)<interval)&&(fabs(Xs[i]-x1)>TOL/1000.0)) 
			{x0= Xs[i];Fx0=Fs[i]; interval=fabs(Xs[i]-x1); }
		}
			

		last_x_new = x_new;  // save old x_new
		// --- newton secant method 
		if (fabs(Fx1-Fx0)> 1.0e-10){
			x_new = -Fx1*(x1-x0)/(Fx1-Fx0)+x1;
		}
		else {
			x_new = (a+b)/2.0;	
		}




		// --- make sure it is in trust region	
		if ((x_new>b+TOL)||(x_new<a-TOL)) {
			x_new=(a+b)/2.0;
		}
		//correction if starked! ----

		if ((fabs(x_new-a)< TOL/1000.)||(fabs(x_new-b)< TOL/1000.))  x_new=(a+b)/2.0;
		if (fabs(last_x_new-x_new)<TOL/1000.) x_new=(a+b)/2.0;
		
		

		// check x_new, is not solution 
		Fx_new=myFunction(x_new);

		if (FEconvergence){
		


		
			if(fabs(Fx_new)<functionTol){
					zeroPoint = x_new;
					opserr<<"Return: x_new is:"<<x_new<<" Fx_new is: "<<Fx_new<<endln;				
					return 1;
			} 

			else if (Fx_new*Fa > functionTol*functionTol){ //good
				
				a = x_new;
				Fa = Fx_new;
				
			}
			else if (Fx_new * Fa < -functionTol*functionTol){
				b = x_new;
				Fb = Fx_new;
				isFbValid = true;

			}
		
			else { opserr<<"bound wrong"<<endln; return -1; }//opserr<<"wrong bound!"<<endln; return -1;}
			printf("a: %f, b: %f, Fa: %f, Fb: %f,x_new: %9.4f, Fx_new: %14.7f \n",a,b,Fa,Fb,x_new,Fx_new );

			}
		else { //FEdivergence
		
			b=x_new;
			Fb=-1.0;
			isFbValid = false;
			printf("diverge-- a: %f, b: %f, Fa: %f, Fb: %f,x_new: %9.4f, Fx_new: %14.7f \n",a,b,Fa,Fb,x_new,Fx_new );	
		}			
		
		
	};//while


	if (iterationNum == maxIterNum) {//not converge
		return -1;
	}
	
	zeroPoint = x_new;


	return 0;
} ;

double SafeGuardedZeroFindingAlgorithm::myFunction(double x){
	return theSamplingAnalysis->getSampledValue(x);
};


SafeGuardedZeroFindingAlgorithm::SafeGuardedZeroFindingAlgorithm(SamplingAnalysis * thePassedAnalysis)
{
	theSamplingAnalysis = thePassedAnalysis;
	functionType=1;
	FEconvergence = true;
	
	if (functionType ==1){  // failure prob.
		x_1 = new double [100];
		G_1 = new double [100];
	}
	else if (functionType ==2){  // upcrossing
		x_1 = new double [100];
		G_1 = new double [100];
		x_2 = new double [100];
		G_2 = new double [100];
	}	
	
	

}

