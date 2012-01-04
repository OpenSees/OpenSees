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
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/ParametricReliabilityAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ParametricReliabilityAnalysis.h>
#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <FindDesignPointAlgorithm.h>
#include <GradientEvaluator.h>
#include <ParameterPositioner.h>
#include <NormalRV.h>
#include <tcl.h>
#include <LimitStateFunctionIter.h>
#include <math.h>

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
									 GradientEvaluator *passedGradGEvaluator,
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

	strcpy(fileName,passedFileName);

	theTclInterp = passedTclInterp;
}


ParametricReliabilityAnalysis::~ParametricReliabilityAnalysis()
{

}



int 
ParametricReliabilityAnalysis::analyze(void)
{
	
	// Alert the user that the FORM analysis has started
	opserr << "Fragility Analysis is running ... " << endln;


	// The relevant commands are now, for instance, given like this
	// performanceFunction 1 "{par_1}-{u_7_1}"
	// KRM note:	{par_1}  -->  \$par(1)
	//				{u_7_1}  -->  \$u(7,1)
	// runParametricReliabilityAnalysis output.out -par 1 -range 14.0 16.0 -num 10


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

	LimitStateFunctionIter lsfIter = theReliabilityDomain->getLimitStateFunctions();
	LimitStateFunction *theLimitStateFunction;
	// Loop over number of limit-state functions and perform FORM analysis
	//for (lsf=1; lsf <= numLsf; lsf++ ) {
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

		// Detect parameter, first the case where the parameter is represented by a parameterPositioner
		numPos = theReliabilityDomain->getNumberOfParameterPositioners();
		if (numPos != 0) {
			theParameterPositioner = theReliabilityDomain->getParameterPositionerPtr(parameterNumber);
			if (theParameterPositioner == 0) {
				opserr << "ParametricReliabilityAnalysis::analyze() -- The parameter number in the " << endln
					<< " fragility analysis object does not match the parameter number " << endln
					<< " being set by the parameter positioner." << endln;
			}
		}
		
		// Now when the original \$par(i) syntax is used
		else {
			
			// check that par(i) exists in the LSF
			int llength = 0;
			int foundit = 0;
			// This Tcl stuff should disappear -- MHS 10/7/2011
			//Tcl_Obj *passedList = theLimitStateFunction->getParameters();
			Tcl_Obj *passedList = 0; // Added by MHS 10/7/2011

			Tcl_Obj *paramList = Tcl_DuplicateObj(passedList);
			Tcl_Obj *objPtr;
			Tcl_ListObjLength(theTclInterp,paramList,&llength);
				
			for (int jk = 0; jk < llength; jk++) {
				Tcl_ListObjIndex(theTclInterp,paramList,jk,&objPtr);
				char *listStr = Tcl_GetStringFromObj(objPtr,NULL);
				
				int tempNumber;
				int args = sscanf(listStr,"par(%i)",&tempNumber);
				if (args == 1 && tempNumber == parameterNumber)
					foundit = 1;
			}
			
			// reclaim Tcl object space
			Tcl_DecrRefCount(paramList);
			
			if (foundit == 0) {
				opserr << "ParametricReliabilityAnalysis::analyze() -- The parameter number " << endln
									<< " in the limit-state function does not match the parameter " << endln
									<< " number in the fragility analysis object." << endln;
			}
		}

		// Range over parameter values
		currentValue = first;
		for (int counter=0; counter<(numIntervals+1); counter++) {

			currentValues(counter) = currentValue;
			
			if (theParameterPositioner != 0) {
				// Possibly set the parameter value in the FE domain
				theParameterPositioner->update(currentValue);
			}
			else {
				// set parameter value in Tcl domain
				char theIndex[20];
				sprintf(theIndex,"%d",parameterNumber);
							
				Tcl_Obj *outp = Tcl_SetVar2Ex(theTclInterp,"par",theIndex,Tcl_NewDoubleObj(currentValue),TCL_LEAVE_ERR_MSG);
				if (outp == NULL) {
					opserr << "ERROR ParametricReliabilityAnalysis -- error setting parameter with tag " << parameterNumber << endln;
					opserr << theTclInterp->result << endln;
					return -1;
				}
			}

			// Find the design point
			if (theFindDesignPointAlgorithm->findDesignPoint() < 0) {

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
				// Needs to be retooled -- MHS 10/7/2011
				//const Matrix &dGdPar = theGradGEvaluator->getDgDpar();
				Matrix dGdPar(5,5); // So it will compile

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

