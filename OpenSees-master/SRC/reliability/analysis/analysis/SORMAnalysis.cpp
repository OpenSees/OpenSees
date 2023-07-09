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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-04-10 18:10:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SORMAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SORMAnalysis.h>
#include <FORMAnalysis.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>
#include <ReliabilityAnalysis.h>
#include <FindCurvatures.h>
#include <LimitStateFunction.h>
#include <LimitStateFunctionIter.h>
#include <NormalRV.h>
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


SORMAnalysis::SORMAnalysis(ReliabilityDomain *passedReliabilityDomain,
                           FunctionEvaluator *passedEvaluator,
                           FORMAnalysis *passedFORM,
                           FindCurvatures *passedCurvaturesAlgorithm,
                           TCL_Char *passedFileName)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
    theFunctionEvaluator = passedEvaluator;
    theFORManalysis = passedFORM;
	theCurvaturesAlgorithm = passedCurvaturesAlgorithm;
	strcpy(fileName,passedFileName);
}


SORMAnalysis::~SORMAnalysis()
{
  
}



int 
SORMAnalysis::analyze(void)
{
	// Alert the user that the SORM analysis has started
	opserr << "SORM Analysis is running ... " << endln;


	// Declare variables used in this method
	int numberOfCurvatures = 0;
	double beta, pf1;
	double psi_beta, product;
	double pf2Breitung, betaBreitung;
	static NormalRV aStdNormRV(1,0.0,1.0);


	// Number of limit-state functions
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();
    int numRV = theReliabilityDomain->getNumberOfRandomVariables();

	// Open output file
	ofstream outputFile( fileName, ios::out );

	LimitStateFunctionIter &lsfIter = theReliabilityDomain->getLimitStateFunctions();
	LimitStateFunction *theLimitStateFunction;

	// Loop over number of limit-state functions
	//while ((theLimitStateFunction = lsfIter()) != 0) {
    for (int lsf = 0; lsf < numLsf; lsf++ ) {
        theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtrFromIndex(lsf);
		int lsfTag = theLimitStateFunction->getTag();
        
		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsfTag << endln;

		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsfTag);

		// Compute curvature(s)
		if (theCurvaturesAlgorithm->computeCurvatures() < 0) {
			opserr << "SORMAnalysis::analyze() - failed while finding " << endln
				<< " curvatures for limit-state function number " << lsf << "." << endln;
			return -1;
		}

		// Get results
		const Vector &curvatures = theCurvaturesAlgorithm->getCurvatures();
		numberOfCurvatures = curvatures.Size();
        opserr << "SORM curvatures: " << curvatures << endln;
			
		// Get FORM results from the function evaluator
        beta = theFunctionEvaluator->getResponseVariable("betaFORM", lsfTag);
		pf1 = 1.0 - aStdNormRV.getCDFvalue(beta);
        

		// Compute failure probability by "Breitung"
		double denominator = aStdNormRV.getCDFvalue(-beta);
		if (denominator == 0.0) {
			opserr << "SORMAnalysis::analyze() - denominator zero " << endln
				<< " due to too large reliability index value." << endln;
			return -1;
		}
		psi_beta = aStdNormRV.getPDFvalue(beta)/denominator;
		product = 1.0;
		for (int i = 0; i < numberOfCurvatures; i++ )
			product = product / sqrt(1.0+psi_beta*curvatures(i));
		pf2Breitung = pf1 * product;


		// Compute corresponding beta's
		betaBreitung = -aStdNormRV.getInverseCDFvalue(pf2Breitung);


		// store key results using function evaluator
        theFunctionEvaluator->setResponseVariable("betaSORM", lsfTag, betaBreitung);
        theFunctionEvaluator->setResponseVariable("pfSORM", lsfTag, pf2Breitung);


		// Print SORM results to the output file
		outputFile << "#######################################################################" << endln;
		outputFile << "#  SORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsfTag <<"            #" << endln;
		outputFile << "#  (Curvatures found from search algorithm.)                          #" << endln;
		outputFile << "#                                                                     #" << endln;
		outputFile << "#  Number of principal curvatures used: ............... " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<numberOfCurvatures
			<< "  #" << endln;
		outputFile << "#  Reliability index beta (impr. Breitung's formula):.. " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<betaBreitung 
			<< "  #" << endln;
		outputFile << "#  Corresponding estimated probability of failure pf2:.." 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf2Breitung 
			<< "  #" << endln;
		outputFile << "#                                                                     #" << endln;
		outputFile << "#######################################################################" << endln << endln << endln;
	}


	// Inform user on screen
	opserr << "SORM analysis completed. " << endln;

	// Clean up
	outputFile.close();

	return 0;
}

