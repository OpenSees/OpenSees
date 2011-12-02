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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-03-04 00:39:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/FiniteDifferenceGradGEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <FiniteDifferenceGradGEvaluator.h>
#include <Vector.h>
#include <GradGEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <GFunEvaluator.h>
#include <RandomVariable.h>
#include <tcl.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


FiniteDifferenceGradGEvaluator::FiniteDifferenceGradGEvaluator(
					GFunEvaluator *passedGFunEvaluator,
					ReliabilityDomain *passedReliabilityDomain,
					Tcl_Interp *passedTclInterp,
					double passedPerturbationFactor,
					bool PdoGradientCheck)
:GradGEvaluator(passedReliabilityDomain, passedTclInterp)
{
	theGFunEvaluator = passedGFunEvaluator;
	perturbationFactor = passedPerturbationFactor;
	doGradientCheck = PdoGradientCheck;

	int nrv = passedReliabilityDomain->getNumberOfRandomVariables();
	grad_g = new Vector(nrv);

	DgDdispl = 0;
	DgDpar = 0;
}

FiniteDifferenceGradGEvaluator::~FiniteDifferenceGradGEvaluator()
{
	delete grad_g;

	if (DgDdispl != 0)
		delete DgDdispl;
	if (DgDpar != 0)
		delete DgDpar;
}



Vector
FiniteDifferenceGradGEvaluator::getGradG()
{
	return (*grad_g);
}



int
FiniteDifferenceGradGEvaluator::evaluateGradG(double gFunValue, Vector passed_x)
{
	// Call base class method
	computeParameterDerivatives(gFunValue);

	
	// Initial declarations
	int numberOfRandomVariables = passed_x.Size();
	Vector perturbed_x(numberOfRandomVariables);
	RandomVariable *theRandomVariable;
	int i;
	double h;
	double gFunValueAStepAhead;
	double stdv;


	// For each random variable: perturb and run analysis again
	for ( i=0 ; i<numberOfRandomVariables ; i++ )
	{
		// Get random variable from domain
		theRandomVariable = theReliabilityDomain->getRandomVariablePtr(i+1);


		// Get the standard deviation
		stdv = theRandomVariable->getStdv();


		// Compute perturbation
		h = stdv/perturbationFactor;


		// Compute perturbed vector of random variables realization
		perturbed_x = passed_x;
		perturbed_x(i) = perturbed_x(i) + h;


		// Evaluate limit-state function
		double result = theGFunEvaluator->runGFunAnalysis(perturbed_x);
		if (result < 0) {
			opserr << "FiniteDifferenceGradGEvaluator::evaluate_grad_g() - " << endln
				<< " could not run analysis to evaluate limit-state function. " << endln;
			return -1;
		}
		result = theGFunEvaluator->evaluateG(perturbed_x);
		if (result < 0) {
			opserr << "FiniteDifferenceGradGEvaluator::evaluate_grad_g() - " << endln
				<< " could not tokenize limit-state function. " << endln;
			return -1;
		}
		gFunValueAStepAhead = theGFunEvaluator->getG();


		// Compute the derivative by finite difference
		(*grad_g)(i) = (gFunValueAStepAhead - gFunValue) / h;
	}





	if (doGradientCheck) {
		char myString[100];
		ofstream outputFile( "FFDgradients.out", ios::out );
		opserr << endln;
		for (int ffd=0; ffd<grad_g->Size(); ffd++) {
			opserr << "FFD("<< (ffd+1) << ") = " << (*grad_g)(ffd) << endln;
			sprintf(myString,"%20.16e ",(*grad_g)(ffd));
			outputFile << myString << endln;
		}
		outputFile.close();
		opserr << "PRESS Ctrl+C TO TERMINATE APPLICATION!" << endln;
		while(true) {
		}
	}






	return 0;

}





Matrix
FiniteDifferenceGradGEvaluator::getDgDdispl()
{
	opserr << "ERROR: FiniteDifferenceGradGEvaluator::getDgDdispl() is " << endln
		<< " currently not implemented." << endln;
	Matrix dummy(1,1);
	return dummy;
}

