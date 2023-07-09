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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-05-13 18:37:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/convergenceCheck/StandardReliabilityConvergenceCheck.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <StandardReliabilityConvergenceCheck.h>
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


StandardReliabilityConvergenceCheck::StandardReliabilityConvergenceCheck(double passedE1, 
																		 double passedE2, 
																		 double pscaleValue,
																		 int print)
:ReliabilityConvergenceCheck()
{
	e1 = passedE1;
	e2 = passedE2;
	criterium1 = 0.0;
	criterium2 = 0.0;
	scaleValue = pscaleValue;
	printFlag = print;
 	/////////////////////////////////////
 	// S Modified by K Fujimura 10/10/2004
 	/////////////////////////////////////
 	fixscale = false;
 	if (scaleValue != 0.0)
		fixscale = true;
 	/////////////////////////////////////
 	//E  Modified by K Fujimura 10/10/2004
 	/////////////////////////////////////
	logfile.open("ConvergenceCheckLog.txt", ios::out);//not in K.F.

}

StandardReliabilityConvergenceCheck::~StandardReliabilityConvergenceCheck()
{
	logfile.close();//not in K.F.
}



int	
StandardReliabilityConvergenceCheck::setScaleValue(double passedScaleValue)
{
 	/////////////////////////////////////
 	//S  Modified by K Fujimura 10/10/2004
 	/////////////////////////////////////
 	if ( fixscale == false ) {
		scaleValue = passedScaleValue;
	}

	//if (scaleValue == 0.0) {
	//	scaleValue = passedScaleValue;
	//}
	//E  Modified by K Fujimura 10/10/2004

	return 0;
}


int	
StandardReliabilityConvergenceCheck::check(const Vector &u, double g, const Vector &gradG)
{
	// Alpha vector
	//Vector alpha = gradG *  ( (-1.0) / gradG.Norm() );
	double gradGnorm = gradG.Norm();
	double oneOver_gradGnorm = 1.0/gradGnorm;

	// The scaling factor (rather new development)
	double temp0=u.Norm();
	if (temp0 < 1.0) {
		temp0 = 1.0; 
	}
	double oneOver_temp0 = 1.0/temp0;

	// Scaled u-vector
	//Vector u_scaled = (1.0/temp0) * u;
	//Vector temp1 = u_scaled - (alpha^u_scaled)*alpha;
	double alpha_dot_uScaled = -(gradG^u)*oneOver_gradGnorm*oneOver_temp0;

	//criterium2 = temp1.Norm();
	//opserr << criterium2 << endln;

	criterium2 = 0.0;
	double tmp1;
	double tmp2 = -alpha_dot_uScaled*oneOver_gradGnorm;
	// This loop replaces 'Vector temp1' above -- MHS
	for (int i = 0; i < u.Size(); i++) {
	  tmp1 = oneOver_temp0*u(i) - tmp2*gradG(i);
	  criterium2 += tmp1*tmp1;
	}
	criterium2 = sqrt(criterium2);
	//opserr << criterium2 << ' ';

	// Convergence criteria
	if (scaleValue == 0.0) {
		opserr << "StandardReliabilityConvergenceCheck::check() -- scale value is zero or has " << 
			"not been set!" << endln;
		criterium1 = fabs(g);
	}
	else
		criterium1 = fabs(g / scaleValue);
	
	// Inform user about convergence status 
	char outputString[100];
	sprintf(outputString,"check1=(%11.3e), check2=(%10.3e), dist=%16.14f",criterium1,criterium2,u.Norm());
	if (printFlag != 0) {
		opserr << outputString << endln;
	}
	logfile << outputString << endln;
	logfile.flush();

	// Return '1' if the analysis converged ('-1' otherwise)
	if ( ( criterium1 < e1 ) && ( criterium2 < e2 ) ) {
		return 1;
	}
	else {
		/////S modified by K Fujimura /////
 		if(criterium2 < e2 ){
 			return -1;
 		}else if(criterium1 < e1 ){
 			return -2;
 		}else{
 			return -3;
 		}
/////		return -1;
		/////E modified by K Fujimura /////
	}
}

int 
StandardReliabilityConvergenceCheck::getNumberOfCriteria()
{
	return 2;
}

double 
StandardReliabilityConvergenceCheck::getCriteriaValue(int whichCriteria)
{
	if (whichCriteria == 1) {
		return criterium1; 
	}
	else if (whichCriteria == 2) {
		return criterium2;
	}
	else {
		opserr << "StandardReliabilityConvergenceCheck::getCriteriaValue() -- " << endln
			<< " criterium number " << whichCriteria << " does not exist!" << endln;
		return 0.0;
	}

}
////////////////////////////////////
//S Modified by K Fujimura 10/10/2004
////////////////////////////////////
void
StandardReliabilityConvergenceCheck::Scalefix(bool fix)
{
	fixscale = fix;
}

double 
StandardReliabilityConvergenceCheck::getCheck1()
{
	return criterium1;
}

double 
StandardReliabilityConvergenceCheck::getCheck2()
{
	return criterium2;
}

int	
StandardReliabilityConvergenceCheck::checkG(double g)
{
	// Convergence criterium 1 only
	if (scaleValue == 0.0) {
		opserr << "StandardReliabilityConvergenceCheck::checkG() -- scale value is zero or has " << 
			"not been set!" << endln;
		criterium1 = fabs(g);
	}
	else
		criterium1 = fabs(g / scaleValue);
	
	// Return '1' if the analysis converged ('-1' otherwise)
	if ( criterium1 < e1 ) {
		return 1;
	} else {
		return -1;
	}
}
////////////////////////////////////
//E Modified by K Fujimura 10/10/2004
////////////////////////////////////
