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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-05-13 16:30:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/FragilityAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <FragilityAnalysis.h>
#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <FindDesignPointAlgorithm.h>
#include <GradGEvaluator.h>
#include <ParameterPositioner.h>
#include <NormalRV.h>
#include <tcl.h>

#include <LimitStateFunctionIter.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

FragilityAnalysis::FragilityAnalysis(ReliabilityDomain *passedReliabilityDomain,
									 FindDesignPointAlgorithm *passedFindDesignPointAlgorithm,
									 GradGEvaluator *passedGradGEvaluator,
									 int pParameterNumber,
									 double pFirst,
									 double pLast,
									 int pNumIntervals,
									 const char *passedFileName,
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

	strcpy(fileName,passedFileName);

	theTclInterp = passedTclInterp;
}


FragilityAnalysis::~FragilityAnalysis()
{

}



int 
FragilityAnalysis::analyze(void)
{
	
	// Alert the user that the FORM analysis has started
	opserr << "Fragility Analysis is running ... " << endln;


	// The relevant commands are now, for instance, given like this
	// performanceFunction 1 "{par_1}-{u_7_1}"
	// runFragilityAnalysis output.out -par 1 -range 14.0 16.0 -num 10


	// Open output file
	ofstream outputFile( fileName, ios::out );


	// Initial declarations
	Vector pf(numIntervals+1);
	Vector pdf(numIntervals+1);
	static NormalRV aStdNormRV(1,0.0,1.0);
	int numPars, numPos;
	double thedGdPar = 0, currentValue, beta;
	Vector currentValues(numIntervals+1);
	ParameterPositioner *theParameterPositioner = 0;


	// Loop over number of limit-state functions and perform FORM analysis
	LimitStateFunctionIter &lsfIter = theReliabilityDomain->getLimitStateFunctions();
	LimitStateFunction *theLimitStateFunction;
	//for (lsf = 1 ; lsf <= numLsf; lsf++ ) {
	while ((theLimitStateFunction = lsfIter()) != 0) {
	  int lsf = theLimitStateFunction->getTag();

		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsf << endln;
		Tcl_SetVar2Ex(theTclInterp,"RELIABILITY_lsf",NULL,Tcl_NewIntObj(lsf),TCL_NAMESPACE_ONLY);

		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);

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
		for (int counter=0; counter < numIntervals+1; counter++) {

			currentValues(counter) = currentValue;

			// Detect parameter and set value, first the 
			// case where the parameter is represented by a 
			// parameterPositioner
			numPos = theReliabilityDomain->getNumberOfParameterPositioners();
			if (numPos != 0) {
				theParameterPositioner = theReliabilityDomain->getParameterPositionerPtr(parameterNumber);
				if (theParameterPositioner == 0) {
					opserr << "FragilityAnalysis::analyze() -- The parameter number in the " << endln
						<< " fragility analysis object does not match the parameter number " << endln
						<< " being set by the parameter positioner." << endln;
				}
			}
			else {
				char separators[5] = "_}{";
				const char *theExpression = theLimitStateFunction->getExpression();
				char lsf_forTokenizing[1000];
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
							opserr << "FragilityAnalysis::analyze() -- The parameter number " << endln
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
			if (theFindDesignPointAlgorithm->findDesignPoint() < 0){

				// Set detectable 'crazy' values when the design point wasn't found
				pf(counter) = -1.0;
				pdf(counter) = -1.0;
			}
			else {

				// Store probability of failure
				const Vector &uStar				= theFindDesignPointAlgorithm->get_u();
				const Vector &alpha				= theFindDesignPointAlgorithm->get_alpha();
				beta = alpha ^ uStar;
				pf(counter) = 1.0 - aStdNormRV.getCDFvalue(beta);

				// Compute PDF, first; derivative of lsf wrt. parameter
				const Matrix &dGdPar = theGradGEvaluator->getDgDpar();

				// Find the right element
				numPars = dGdPar.noRows();
				for (int i=0; i<numPars; i++) {
					if (((int)(dGdPar(i,0))) == parameterNumber) {
						thedGdPar = dGdPar(i,1);
					}
				}

				// Gradient of limit-state function
				const Vector &gradient = theFindDesignPointAlgorithm->getGradientInStandardNormalSpace();

				// Compute PDF value
				pdf(counter) = aStdNormRV.getPDFvalue(-beta) / gradient.Norm() * thedGdPar;

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

