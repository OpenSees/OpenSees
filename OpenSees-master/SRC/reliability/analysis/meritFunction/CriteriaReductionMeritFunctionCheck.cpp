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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-02-29 19:47:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/meritFunction/CriteriaReductionMeritFunctionCheck.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <CriteriaReductionMeritFunctionCheck.h>
#include <MeritFunctionCheck.h>
#include <ReliabilityConvergenceCheck.h>
#include <math.h>
#include <Vector.h>


CriteriaReductionMeritFunctionCheck::CriteriaReductionMeritFunctionCheck(ReliabilityConvergenceCheck *passedReliabilityConvergenceCheck)
:MeritFunctionCheck()
{
	theReliabilityConvergenceCheck = passedReliabilityConvergenceCheck;
}

CriteriaReductionMeritFunctionCheck::~CriteriaReductionMeritFunctionCheck()
{
}







int
CriteriaReductionMeritFunctionCheck::check(const Vector &u_old, 
					   double g_old, 
					   const Vector &grad_G_old, 
					   double stepSize,
					   const Vector &stepDirection,
					   double g_new, 
					   const Vector &grad_G_new,
					   /// added by K Fujimura /////
								  int reschk)
{
	// New point in standard normal space    (// not in K.F. version)
	//Vector u_new = u_old + stepSize*stepDirection;
  Vector u_new(u_old);
  u_new.addVector(1.0, stepDirection, stepSize);


	// Number of convergence criteria
	int numCrit = theReliabilityConvergenceCheck->getNumberOfCriteria();


	// Initial declarations
	int i;
	Vector oldCriteriaValues(numCrit);
	Vector newCriteriaValues(numCrit);


	// Convergence checks for the old point
	theReliabilityConvergenceCheck->check(u_old,g_old,grad_G_old);
	for (i=1; i<=numCrit; i++) {
		oldCriteriaValues(i-1) = theReliabilityConvergenceCheck->getCriteriaValue(i);
	}
	

	// Convergence checks for the new point
	theReliabilityConvergenceCheck->check(u_new,g_new,grad_G_new);
	for (i=1; i<=numCrit; i++) {
		newCriteriaValues(i-1) = theReliabilityConvergenceCheck->getCriteriaValue(i);
	}


	// All criteria must have improved for the step to be OK
	bool OK = true;
	for (i=1; i<=numCrit; i++) {
		if (newCriteriaValues(i-1)>oldCriteriaValues(i-1)) {
			OK = false;
		}
	}


	// Return
	if (OK) {
		return 1;
	}
	else {
		return -1;
	}
}





double
CriteriaReductionMeritFunctionCheck::getMeritFunctionValue(const Vector &u, 
							   double g, 
							   const Vector &grad_G)
{
	opserr << "CriteriaReductionMeritFunctionCheck::getMeritFunctionValue() -- this method" << endln
		<< " is not implemented in this specific class." << endln;

	return 0.0;
}

int
CriteriaReductionMeritFunctionCheck::updateMeritParameters(const Vector &u, 
							   double g, 
							   const Vector &grad_G, int reschk)
{
	opserr << "CriteriaReductionMeritFunctionCheck::updateMeritParameters() -- this method" << endln
		<< " is not implemented in this specific class." << endln;

	return 0;
}
