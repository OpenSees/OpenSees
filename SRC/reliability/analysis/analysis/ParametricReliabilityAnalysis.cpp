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
                                                                        
// $Revision: 1.2 $
// $Date: 2004-08-27 17:55:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/ParametricReliabilityAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ParametricReliabilityAnalysis.h>
#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <FindDesignPointAlgorithm.h>
#include <GradGEvaluator.h>
#include <ParameterPositioner.h>
#include <NormalRV.h>
#include <tcl.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

ParametricReliabilityAnalysis::ParametricReliabilityAnalysis(ReliabilityDomain *passedReliabilityDomain,
									 FindDesignPointAlgorithm *passedFindDesignPointAlgorithm,
									 GradGEvaluator *passedGradGEvaluator,
									 int pParameterNumber,
									 double pFirst,
									 double pLast,
									 int pNumIntervals,
									 TCL_Char *passedFileName,
									 Tcl_Interp *passedTclInterp)
:ReliabilityAnalysis()
{
	parameterNumber = pParameterNumber;
	first = pFirst;
	last = pLast;
	numIntervals = pNumIntervals;

	theReliabilityDomain = passedReliabilityDomain;
	theFindDesignPointAlgorithm = passedFindDesignPointAlgorithm;
	theGradGEvaluator = passedGradGEvaluator;

	fileName = new char[256];
	strcpy(fileName,passedFileName);

	theTclInterp = passedTclInterp;
}


ParametricReliabilityAnalysis::~ParametricReliabilityAnalysis()
{
	if (fileName != 0)
	delete [] fileName;
}



int 
ParametricReliabilityAnalysis::analyze(void)
{
	
	// Alert the user that the FORM analysis has started
	opserr << "Fragility Analysis is running ... " << endln;


	// The relevant commands are now, for instance, given like this
	// performanceFunction 1 "{par_1}-{u_7_1}"
	// runParametricReliabilityAnalysis output.out -par 1 -range 14.0 16.0 -num 10


	// Open output file
	ofstream outputFile( fileName, ios::out );


	// Get number of limit-state functions
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();

	// Initial declarations
	Vector pf(numIntervals+1);
	Vector pdf(numIntervals+1);
	Vector uStar, alpha;
	double beta;
	NormalRV aStdNormRV(1,0.0,1.0,0.0);
	Matrix dGdPar;
	int numPars;
	double thedGdPar;
	Vector gradient;
	double currentValue;
	Vector currentValues(numIntervals+1);
	int numPos;
	ParameterPositioner *theParameterPositioner = 0;


	// Loop over number of limit-state functions and perform FORM analysis
	for (int lsf=1; lsf<=numLsf; lsf++ ) {


		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsf << endln;


		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);


		// "Download" limit-state function from reliability domain
		// fmk to Terje: you just set it so why do you need the tag again
		//     also you can't do a redef inside a loop with the same def as loop variable!!!!
		// => changing int lst to int newLsf in line below
		int newLsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
		LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
		if (theLimitStateFunction == 0) {
			opserr << "ParametricReliabilityAnalysis::analyze() - could not find" << endln
				<< " limit-state function with tag #" << lsf << "." << endln;
			return -1;
		}


		// Print results to output file
		outputFile << "#######################################################################" << endln;
		outputFile << "#  FORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"            #" << endln;
		outputFile << "#                                                                     #" << endln;
		outputFile << "#                                                                     #" << endln;
		outputFile << "#                    Failure probability     Estimated probability    #" << endln;
		outputFile << "#    Parameter       estimate (fragility)     densitity function      #" << endln;
		outputFile << "#     value               (CDF)                    (PDF)              #" << endln;
		outputFile.setf(ios::scientific, ios::floatfield);
		outputFile.flush();


		// Range over parameter values
		currentValue = first;
		for (int counter=0; counter<(numIntervals+1); counter++) {

			currentValues(counter) = currentValue;

			// Detect parameter and set value, first the 
			// case where the parameter is represented by a 
			// parameterPositioner
			numPos = theReliabilityDomain->getNumberOfParameterPositioners();
			if (numPos != 0) {
				theParameterPositioner = theReliabilityDomain->getParameterPositionerPtr(parameterNumber);
				if (theParameterPositioner == 0) {
					opserr << "ParametricReliabilityAnalysis::analyze() -- The parameter number in the " << endln
						<< " fragility analysis object does not match the parameter number " << endln
						<< " being set by the parameter positioner." << endln;
				}
			}
			else {
				char separators[5] = "_}{";
				char *theExpression = theLimitStateFunction->getExpression();
				char *lsf_forTokenizing = new char[1000];
				strcpy(lsf_forTokenizing,theExpression);
				char *tokenPtr = strtok( lsf_forTokenizing, separators);
				while ( tokenPtr != NULL ) {
					if ( strcmp(tokenPtr, "par") == 0) {
						tokenPtr = strtok( NULL, separators);
						int par = atoi( tokenPtr );
						if (par==parameterNumber) {
							char tclAssignment[100];
							sprintf(tclAssignment , "set par_%d  %15.5f", parameterNumber, currentValue);
							Tcl_Eval( theTclInterp, tclAssignment);
						}
						else {
							opserr << "ParametricReliabilityAnalysis::analyze() -- The parameter number " << endln
								<< " in the limit-state function does not match the parameter " << endln
								<< " number in the fragility analysis object." << endln;
						}
					}
					tokenPtr = strtok( NULL, separators);
				}
			}


			// Possibly set the parameter value in the FE domain
			if (theParameterPositioner != 0) {
				theParameterPositioner->update(currentValue);
			}


			// Find the design point
			if (theFindDesignPointAlgorithm->findDesignPoint(theReliabilityDomain) < 0){

				// Set detectable 'crazy' values when the design point wasn't found
				pf(counter) = -1.0;
				pdf(counter) = -1.0;
			}
			else {

				// Store probability of failure
				uStar				= theFindDesignPointAlgorithm->get_u();
				alpha				= theFindDesignPointAlgorithm->get_alpha();
				beta = alpha ^ uStar;
				pf(counter) = 1.0 - aStdNormRV.getCDFvalue(beta);

				// Compute PDF, first; derivative of lsf wrt. parameter
				dGdPar = theGradGEvaluator->getDgDpar();

				// Find the right element
				numPars = dGdPar.noRows();
				for (int i=0; i<numPars; i++) {
					if (((int)(dGdPar(i,0))) == parameterNumber) {
						thedGdPar = dGdPar(i,1);
					}
				}

				// Gradient of limit-state function
				gradient = theFindDesignPointAlgorithm->getGradientInStandardNormalSpace();

				// Compute PDF value
				pdf(counter) = fabs( aStdNormRV.getPDFvalue(-beta) / gradient.Norm() * thedGdPar );

			}

			currentValue = currentValue + (last-first)/numIntervals;


			// Print results to output file
			outputFile << "#   " << setprecision(3)<<setw(11)<<currentValues(counter)<<"         ";
			if (pf(counter)==-1.0) {
				outputFile << "--failed--              ";
				outputFile << "--failed--            #" << endln;
			}
			else {
				outputFile <<setprecision(3)<<setw(11)<<pf(counter)<<"             ";
				outputFile <<setprecision(3)<<setw(11)<<pdf(counter)<<"           #" << endln;
			}
			outputFile.flush();

		} // Done looping over parameter range


		outputFile << "#                                                                     #" << endln;
		outputFile << "#######################################################################" << endln << endln << endln;

	} // Done looping over limit-state functions
		

	// Print summary of results to screen (more here!!!)
	opserr << "Fragility Analysis completed." << endln;


	// Clean up
	outputFile.close();

	return 0;
}

