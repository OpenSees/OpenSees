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
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2008-02-29 19:47:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/convergenceCheck/OptimalityConditionReliabilityConvergenceCheck.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <OptimalityConditionReliabilityConvergenceCheck.h>
#include <ReliabilityConvergenceCheck.h>
#include <math.h>
#include <Vector.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


OptimalityConditionReliabilityConvergenceCheck::OptimalityConditionReliabilityConvergenceCheck(double passedE1, double passedE2, double pscaleValue, int print)
:ReliabilityConvergenceCheck()
{
	e1 = passedE1;
	e2 = passedE2;
	criterium1 = 0.0;
	criterium2 = 0.0;
	scaleValue = pscaleValue;
	printFlag = print;

 	//S  Modified by K Fujimura 10/10/2004
 	fixscale = false;
 	if(scaleValue != 0.0) 	fixscale = true;
 	//E  Modified by K Fujimura 10/10/2004
 	
	logfile.open("ConvergenceCheckLog.txt", ios::out);//not in K.F.
}

OptimalityConditionReliabilityConvergenceCheck::~OptimalityConditionReliabilityConvergenceCheck()
{
	logfile.close();//not in K.F.
}



int	
OptimalityConditionReliabilityConvergenceCheck::setScaleValue(double passedScaleValue)
{
	/////S modified by K Fujimura /////
	if(!fixscale){
		scaleValue = passedScaleValue;
	}
	//if (scaleValue == 0.0) {
	//	scaleValue = passedScaleValue;
	//}
	/////E modified by K Fujimura /////

	return 0;
}


int	
OptimalityConditionReliabilityConvergenceCheck::check(const Vector &u, double g, const Vector &gradG)
{
	if (scaleValue == 0.0) {
		opserr << "OptimalityConditionReliabilityConvergenceCheck::check() --" << endln
			<< " scale value has not been set!" << endln;
	}

	// Convergence criteria
	criterium1 = fabs(g / scaleValue);
	criterium2 = 1.0 - fabs(1.0/(gradG.Norm()*u.Norm()) * (gradG^u));
	
 	/////S Modified by K Fujimura///////
 	if(u.Norm()!=0.0){
 		criterium2 = 1.0+1.0/(gradG.Norm()*u.Norm()) * (gradG^u);
 	}else criterium2 = 1.0;
	/////E Modified by K Fujimura///////


	// Inform user about convergence status 
	static ofstream logfile( "SearchLog.out", ios::out );///// Modified by K Fujimura
	char outputString[100];
	sprintf(outputString,"check1=(%11.3e), check2=(%10.3e), dist=%16.14f",criterium1,criterium2,u.Norm());
	if (printFlag != 0) {
		opserr << outputString << endln;
	}
	logfile << outputString << endln;
	logfile.flush();//not in K.F.

	// Return '1' if the analysis converged ('-1' otherwise)
	if ( ( criterium1 < e1 ) && ( criterium2 < e2 ) ) {
		return 1;
	}
	else {
/////S Modified by K Fujimura///////
		if(criterium2 < e2 ){
 			return -1;
 		}else if(criterium1 < e1 ){
 			return -2;
 		}else{
 			return -3;
		}
/////E Modified by K Fujimura///////
	}
}

int 
OptimalityConditionReliabilityConvergenceCheck::getNumberOfCriteria()
{
	return 2;
}

double 
OptimalityConditionReliabilityConvergenceCheck::getCriteriaValue(int whichCriteria)
{
	if (whichCriteria == 1) {
		return criterium1; 
	}
	else if (whichCriteria == 2) {
		return criterium2;
	}
	else {
		opserr << "OptimalityConditionReliabilityConvergenceCheck::getCriteriaValue() -- " << endln
			<< " criterium number " << whichCriteria << " does not exist!" << endln;
		return 0.0;
	}

}

/////S added by K Fujimura /////

double 
OptimalityConditionReliabilityConvergenceCheck::getCheck1()
{ return criterium1;}
double 
OptimalityConditionReliabilityConvergenceCheck::getCheck2()
{ return criterium2;}
void
OptimalityConditionReliabilityConvergenceCheck::Scalefix(bool fix)
{
	fixscale=fix;
}
int	
OptimalityConditionReliabilityConvergenceCheck::checkG(double g)
{
	if (scaleValue == 0.0) {
		opserr << "OptimalityConditionReliabilityConvergenceCheck::check() --" << endln
			<< " scale value has not been set!" << endln;
	}
	// Convergence criteria
	criterium1 = fabs(g / scaleValue);
	// Return '1' if the analysis converged ('-1' otherwise)
	if (  criterium1 < e1 ) {
		return 1;
	}else {
		return -1;
	}
}

/////E added by K Fujimura /////