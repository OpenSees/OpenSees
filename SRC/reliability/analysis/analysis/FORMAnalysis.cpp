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
// $Date: 2003-10-27 23:45:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/FORMAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <FORMAnalysis.h>
#include <FindDesignPointAlgorithm.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixOperations.h>
#include <NormalRV.h>
#include <RandomVariable.h>
#include <math.h>
#include <ProbabilityTransformation.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


FORMAnalysis::FORMAnalysis(ReliabilityDomain *passedReliabilityDomain,
						   FindDesignPointAlgorithm *passedFindDesignPointAlgorithm,
						   ProbabilityTransformation *passedProbabilityTransformation,
						   TCL_Char *passedFileName,
						   int p_relSensTag)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	theFindDesignPointAlgorithm = passedFindDesignPointAlgorithm;
	theProbabilityTransformation = passedProbabilityTransformation;
	fileName = new char[256];
	strcpy(fileName,passedFileName);
	relSensTag = p_relSensTag;
}


FORMAnalysis::~FORMAnalysis()
{
	if (fileName != 0)
		delete [] fileName;
}



int 
FORMAnalysis::analyze(void)
{

	// Alert the user that the FORM analysis has started
	opserr << "FORM Analysis is running ... " << endln;


	// Declare variables used in this method
	Vector xStar;
	Vector uStar;
	Vector alpha;
	Vector gamma;
	double stdv,mean;
	int i;
	double Go, Glast;
	Vector uSecondLast;
	Vector alphaSecondLast;
	Vector lastSearchDirection;
	double beta;
	double pf1;
	int numberOfEvaluations;
	int lsf;
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();
	RandomVariable *aRandomVariable;
	LimitStateFunction *theLimitStateFunction;
	NormalRV *aStdNormRV=0;
	aStdNormRV = new NormalRV(1,0.0,1.0,0.0);
	Vector delta(numRV); 
	Vector eta(numRV);
	Vector kappa(numRV);


	// Check if computer ran out of memory
	if (aStdNormRV==0) {
		opserr << "FORMAnalysis::analyze() - out of memory while instantiating internal objects." << endln;
		return -1;
	}


	// Open output file
	ofstream outputFile( fileName, ios::out );


	// Loop over number of limit-state functions and perform FORM analysis
	for (lsf=1; lsf<=numLsf; lsf++ ) {


		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsf << endln;


		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);


		// Get the limit-state function pointer
		theLimitStateFunction = 0;
		lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
		theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
		if (theLimitStateFunction == 0) {
			opserr << "FORMAnalysis::analyze() - could not find" << endln
				<< " limit-state function with tag #" << lsf << "." << endln;
			return -1;
		}


		// Find the design point
		if (theFindDesignPointAlgorithm->findDesignPoint(theReliabilityDomain) < 0){
			opserr << "FORMAnalysis::analyze() - failed while finding the" << endln
				<< " design point for limit-state function number " << lsf << "." << endln;
			outputFile << "#######################################################################" << endln;
			outputFile << "#  FORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
				<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"            #" << endln;
			outputFile << "#                                                                     #" << endln;
			outputFile << "#  No convergence!                                                    #" << endln;
			outputFile << "#                                                                     #" << endln;
			outputFile << "#######################################################################" << endln << endln << endln;
		}
		else {
		
			// Get results from the "find desingn point algorithm"
			xStar				= theFindDesignPointAlgorithm->get_x();
			uStar				= theFindDesignPointAlgorithm->get_u();
			alpha				= theFindDesignPointAlgorithm->get_alpha();
			gamma				= theFindDesignPointAlgorithm->get_gamma();
			i					= theFindDesignPointAlgorithm->getNumberOfSteps();
			Go					= theFindDesignPointAlgorithm->getFirstGFunValue();
			Glast				= theFindDesignPointAlgorithm->getLastGFunValue();
			uSecondLast			= theFindDesignPointAlgorithm->getSecondLast_u();
			alphaSecondLast		= theFindDesignPointAlgorithm->getSecondLast_alpha();
			lastSearchDirection	= theFindDesignPointAlgorithm->getLastSearchDirection();
			numberOfEvaluations	= theFindDesignPointAlgorithm->getNumberOfEvaluations();


			// Postprocessing
			beta = alpha ^ uStar;
			pf1 = 1.0 - aStdNormRV->getCDFvalue(beta);


			// Reliability sensitivity analysis wrt. mean/stdv
			if (relSensTag == 1) {
				Vector DuStarDmean;
				Vector DuStarDstdv; 
				double dBetaDmean;
				double dBetaDstdv;

				for ( int j=1; j<=numRV; j++ )
				{
					DuStarDmean = theProbabilityTransformation->meanSensitivityOf_x_to_u(xStar,j);
					DuStarDstdv = theProbabilityTransformation->stdvSensitivityOf_x_to_u(xStar,j);
					dBetaDmean = alpha^DuStarDmean;
					dBetaDstdv = alpha^DuStarDstdv;
					aRandomVariable = theReliabilityDomain->getRandomVariablePtr(j);
					stdv = aRandomVariable->getStdv();
					mean = aRandomVariable->getMean();
					delta(j-1) = stdv * dBetaDmean;
					eta(j-1) = stdv * dBetaDstdv;
					// (Kappa is the sensitivity wrt. the coefficient of variation)
					if (mean != 0.0) {
						kappa(j-1) = -dBetaDmean*stdv/((stdv/mean)*(stdv/mean))+dBetaDstdv*mean;
					}
					else {
						kappa(j-1) = 0.0;
					}
				}
				// Don't scale them:
//				delta = delta * (1.0/delta.Norm());
//				eta = eta * (1.0/eta.Norm());
//				kappa = kappa * (1.0/kappa.Norm());
			}


			// Store key results in the limit-state functions
			theLimitStateFunction->FORMReliabilityIndexBeta				= beta;
			theLimitStateFunction->FORMProbabilityOfFailure_pf1			= pf1;
			theLimitStateFunction->designPoint_x_inOriginalSpace		= xStar;
			theLimitStateFunction->designPoint_u_inStdNormalSpace		= uStar;
			theLimitStateFunction->normalizedNegativeGradientVectorAlpha= alpha;
			theLimitStateFunction->importanceVectorGamma				= gamma;
			theLimitStateFunction->numberOfStepsToFindDesignPointAlgorithm		= i;
			theLimitStateFunction->GFunValueAtStartPt					= Go;
			theLimitStateFunction->GFunValueAtEndPt						= Glast;
			theLimitStateFunction->secondLast_u							= uSecondLast;
			theLimitStateFunction->secondLastAlpha						= alphaSecondLast;
			theLimitStateFunction->lastSearchDirection					= lastSearchDirection;


			// Print FORM results to the output file
			outputFile << "#######################################################################" << endln;
			outputFile << "#  FORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
				<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"            #" << endln;
			outputFile << "#                                                                     #" << endln;
			outputFile << "#  Limit-state function value at start point: ......... " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<Go 
				<< "  #" << endln;
			outputFile << "#  Limit-state function value at end point: ........... " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<Glast 
				<< "  #" << endln;
			outputFile << "#  Number of steps: ................................... " 
				<<setiosflags(ios::left)<<setw(12)<<i 
				<< "  #" << endln;
			outputFile << "#  Number of g-function evaluations: .................. " 
				<<setiosflags(ios::left)<<setw(12)<<numberOfEvaluations 
				<< "  #" << endln;
			outputFile << "#  Reliability index beta: ............................ " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<beta 
				<< "  #" << endln;
			outputFile << "#  FO approx. probability of failure, pf1: ............ " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf1 
				<< "  #" << endln;
			outputFile << "#                                                                     #" << endln;
			outputFile << "# rv#     x*          u*         alpha    gamma    delta    eta       #" << endln;
			outputFile.setf( ios::scientific, ios::floatfield );
			for (int i=0;  i<xStar.Size(); i++) {
			outputFile << "#  " <<setw(3)<<(i+1)<<" ";
			outputFile.setf(ios::scientific, ios::floatfield);
			if (xStar(i)<0.0) { outputFile << "-"; }
			else { outputFile << " "; }
			outputFile <<setprecision(3)<<setw(11)<<fabs(xStar(i));
			if (uStar(i)<0.0) { outputFile << "-"; }
			else { outputFile << " "; }
			outputFile <<setprecision(3)<<setw(11)<<fabs(uStar(i));
			outputFile.unsetf( ios::scientific );
			outputFile.setf(ios::fixed, ios::floatfield);

			if (alpha(i)<0.0) { outputFile << "-"; }
			else { outputFile << " "; }
			outputFile<<setprecision(5)<<setw(8)<<fabs(alpha(i));

			if (gamma(i)<0.0) { outputFile << "-"; }
			else { outputFile << " "; }
			outputFile<<setprecision(5)<<setw(8)<<fabs(gamma(i));		
			
			if (relSensTag == 1) {
				if (delta(i)<0.0) { outputFile << "-"; }
				else { outputFile << " "; }
				outputFile<<setprecision(5)<<setw(8)<<fabs(delta(i));		
				
				if (eta(i)<0.0) { outputFile << "-"; }
				else { outputFile << " "; }
				outputFile<<setprecision(5)<<setw(8)<<fabs(eta(i));		

//              Printing reliability sensitivity wrt. coefficient of variation
//				aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i+1);
//				mean = aRandomVariable->getMean();
//				if (mean==0.0) { outputFile << "    -    "; }
//				else {
//					if (kappa(i)<0.0) { outputFile << "-"; }
//					else { outputFile << " "; }
//					outputFile<<setprecision(5)<<setw(8)<<fabs(kappa(i));		
//				}
			}
			else {
				outputFile << "    -        -    ";
			}
			
			outputFile<<"   #" << endln;
			}
			outputFile << "#                                                                     #" << endln;
			outputFile << "#######################################################################" << endln << endln << endln;


			// Inform the user that we're done with this limit-state function
			opserr << "Done analyzing limit-state function " << lsf << ", beta=" << beta << endln;
		}

	}


	// Clean up
	outputFile.close();
	delete aStdNormRV;

	// Print summary of results to screen (more here!!!)
	opserr << "FORMAnalysis completed." << endln;

	return 0;
}

