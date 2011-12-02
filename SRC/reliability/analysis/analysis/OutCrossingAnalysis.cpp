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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-10-27 23:45:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/OutCrossingAnalysis.cpp,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <OutCrossingAnalysis.h>
#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <FindDesignPointAlgorithm.h>
#include <math.h>
#include <tcl.h>
#include <string.h>
#include <NormalRV.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


OutCrossingAnalysis::OutCrossingAnalysis(
				ReliabilityDomain *theRelDom,
				GFunEvaluator *theGFunEval,
				GradGEvaluator *theSensEval,
				FindDesignPointAlgorithm *theFindDesPt,
				int pAnalysisType,
				int p_stepsToStart,
				int p_stepsToEnd,
				int p_sampleFreq,
				double p_littleDeltaT,
				TCL_Char *passedFileName)
:ReliabilityAnalysis()
{
	theFindDesignPointAlgorithm = theFindDesPt;
	theReliabilityDomain = theRelDom;
	theGFunEvaluator = theGFunEval;
	theGradGEvaluator = theSensEval;
	analysisType = pAnalysisType;
	stepsToStart = p_stepsToStart;
	stepsToEnd = p_stepsToEnd;
	sampleFreq = p_sampleFreq;
	littleDeltaT = p_littleDeltaT;
	fileName = new char[256];
	strcpy(fileName,passedFileName);
}

OutCrossingAnalysis::~OutCrossingAnalysis()
{
	if (fileName != 0)
		delete [] fileName;
}

int 
OutCrossingAnalysis::analyze(void)
{

	// Alert the user that the analysis has started
	opserr << "Out-Crossing Analysis is running ... " << endln;

	// Declare variables used in this method
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();
	LimitStateFunction *theLimitStateFunction;
	NormalRV aStdNormRV(1,0.0,1.0,0.0);
	Matrix DgDdispl;
	int numVel, i, j, k, kk, nodeNumber, dofNumber, lsf; 
	double dgduValue, accuSum;
	Vector uStar, uStar2;
	Vector alpha(numRV), alpha2(numRV), alpha_k(numRV), alpha_kk(numRV);
	double beta1, beta2;
	double pf1, pf2;
	int n_2;
	double a, b, integral, h, fa, fb, sum_fx2j, sum_fx2j_1, Pmn1;
	char string[500];


	// Determine number of points
	double nsteps = stepsToEnd-stepsToStart;
	int numPoints = (int)floor(nsteps/sampleFreq);
	numPoints++;
	double dt = theGFunEvaluator->getDt();
	double Dt = dt*sampleFreq;
	double T = dt*sampleFreq*numPoints;


	// Allocate reult vectors
	Vector nu(numPoints);
	Vector ED(numPoints);
	Vector pf(numPoints);
	Vector beta(numPoints);
	Matrix Pmn2(numPoints,numPoints);
	Matrix allAlphas(numRV,numPoints);


	// Open output file and start writing to it
	ofstream outputFile( fileName, ios::out );


	// Loop over number of limit-state functions and perform analysis
	for (lsf=1; lsf<=numLsf; lsf++ ) {


		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsf << endln;


		// Start printing results to the output file
		outputFile << "########################################################R##############" << endln;
		outputFile << "#  OUT-CROSSING RESULTS, LIMIT-STATE FUNCTION NUMBER      "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"        #" << endln;
		outputFile << "#                                                                     #" << endln;



		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);


		// Get the limit-state function pointer
		theLimitStateFunction = 0;
		lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
		theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
		if (theLimitStateFunction == 0) {
			opserr << "OutCrossingAnalysis::analyze() - could not find" << endln
				<< " limit-state function with tag #" << lsf << "." << endln;
			return -1;
		}


		// Loop over the intervals where probabilities are to be computed
		bool oneFailed = false;
		bool DSPTfailed;
		for (i=1; i<=numPoints; i++) {

			// First assume that we find the design point(s)
			DSPTfailed = false;

			// Set 'nsteps' in the GFunEvaluator 
			theGFunEvaluator->setNsteps(stepsToStart+(i-1)*sampleFreq);


			// Inform the user
			char printTime[10];
			sprintf(printTime, "%5.3f", ((stepsToStart+(i-1)*sampleFreq)*dt) );
			if (analysisType == 1) {
				strcpy(string,theLimitStateFunction->getExpression());
				opserr << " ...evaluating -G1=" << string << " at time " << printTime << " ..." << endln;
			}
			else {
				opserr << " ...evaluating performance function at time " << printTime << " ..." << endln;
			}


			// Find the design point for the original limit-state function
			if (theFindDesignPointAlgorithm->findDesignPoint(theReliabilityDomain) < 0){
				opserr << "OutCrossingAnalysis::analyze() - failed while finding the" << endln
					<< " design point for limit-state function number " << lsf << "." << endln;
				DSPTfailed = true;
				oneFailed = true;
			}


			// Start writing to output file
			if (i==1) {
				outputFile << "#         Reliability    Estimated         Mean         Duration      #" << endln;
				outputFile << "#           index         failure       out-crossing    of single     #" << endln;
				outputFile << "#  Time      beta        probability       rate         excursion     #" << endln;
				outputFile << "#                                                                     #" << endln;
				outputFile.flush();
			}

			if (!DSPTfailed) {


				// Get results from the "find design point algorithm"
				uStar = theFindDesignPointAlgorithm->get_u();
				alpha = theFindDesignPointAlgorithm->get_alpha();


				// Put all alphas into rows in a matrix
				for (j=0; j<alpha.Size(); j++) {
					allAlphas(j,i-1) = alpha(j);
				}

		
				// Postprocessing (note that g = -g1)
				beta(i-1) = alpha ^ uStar;
				pf(i-1) = 1.0 - aStdNormRV.getCDFvalue(beta(i-1));
				beta1 = -beta(i-1);
				pf1 = 1.0 - pf(i-1);


				if (analysisType == 1) {
					// Get the 'dgdu' vector from the sensitivity evaluator
					// (The returned matrix containes 'node#' 'dir#' 'dgdu' in rows)
					DgDdispl = theGradGEvaluator->getDgDdispl();


					// Add extra term to limit-state function
					// (should add an alternative option where user give additional limit-state functions)
					// (... because this works only when g is linear in u)
					numVel = DgDdispl.noRows();
					accuSum = 0.0;
					char *expressionPtr = new char[100];
					for (j=0; j<numVel; j++) {

						nodeNumber = (int)DgDdispl(j,0);
						dofNumber = (int)DgDdispl(j,1);
						dgduValue = DgDdispl(j,2);

						char expression[100];
						sprintf(expression,"+(%10.8f)*(%8.5f)*{ud_%d_%d}",littleDeltaT, dgduValue, nodeNumber, dofNumber);
						strcpy(expressionPtr,expression);

						// Add it to the limit-state function
						theLimitStateFunction->addExpression(expressionPtr);
					}
					delete []expressionPtr;


					// Inform the user
					strcpy(string,theLimitStateFunction->getExpression());
					opserr << " ...evaluating G2=" << string << endln;


					// Find the design point for the edited limit-state function
					if (theFindDesignPointAlgorithm->findDesignPoint(theReliabilityDomain) < 0){
						opserr << "OutCrossingAnalysis::analyze() - failed while finding the" << endln
							<< " design point for limit-state function number " << lsf << "." << endln;
						DSPTfailed = true;
						oneFailed = true;
					}


					// Zero out the added expression in the limit-state function
					theLimitStateFunction->removeAddedExpression();


					if (!DSPTfailed) {

						// Get results from the "find design point algorithm"
						uStar2 = theFindDesignPointAlgorithm->get_u();
						alpha2 = theFindDesignPointAlgorithm->get_alpha();


						// Postprocessing (remember; here is an assumption that the mean point is in the safe domain)
						beta2 = (alpha2 ^ uStar2);
						pf2 = 1.0 - aStdNormRV.getCDFvalue(beta2);
					}
					else {
						outputFile << "#  Second limit-state function did not converge.                      #" << endln;
					}
				
					// Post-processing to find parallel system probability
					a = -(alpha ^ alpha2);	// Interval start
					b = 0.0;				// Interval end
					n_2 = 600;				// Half the number of intervals
					h = b-a;
					fa = functionToIntegrate(a,beta1,beta2);
					fb = functionToIntegrate(b,beta1,beta2);
					sum_fx2j = 0.0;
					sum_fx2j_1 = 0.0;
					for (int j=1;  j<=n_2;  j++) {
						sum_fx2j = sum_fx2j + functionToIntegrate(   (double) (a+(j*2)*h/(2*n_2)) ,beta1, beta2  );
						sum_fx2j_1 = sum_fx2j_1 + functionToIntegrate(   (double)(a+(j*2-1)*h/(2*n_2)) , beta1, beta2  );
					}
					sum_fx2j = sum_fx2j - functionToIntegrate((double)(b),beta1,beta2);
					integral = h/(2*n_2)/3.0*(fa + 2.0*sum_fx2j + 4.0*sum_fx2j_1 + fb);
					Pmn1 = aStdNormRV.getCDFvalue(-beta1)*aStdNormRV.getCDFvalue(-beta2) - integral;
				
				}
				else {
					// Use Heonsang's method to find new alpha and beta
					beta2 = -beta1;
					alpha2(0) = alpha(0) - littleDeltaT/Dt * ( alpha(0)-0.0 );
					for (j=1; j<alpha.Size(); j++) {
						alpha2(j) = alpha(j) - littleDeltaT/Dt * ( alpha(j)-alpha(j-1) );
					}

					// Post-processing to find parallel system probability
					a = -1.0*(alpha ^ alpha2);
					double pi = 3.14159265358979;
					Pmn1 = 1.0/(2.0*pi) * exp(-beta2*beta2*0.5) * (asin(a)+1.570796326794897);
				
				
				}


				// POST-PROCESSING

				// Mean out-crossing rate
				if (Pmn1<1.0e-10) {
					opserr << "WARNING: Zero or negative parallel probability: " << endln
						<< " The correlation is probably too high! " << endln;
				}
				nu(i-1) = Pmn1 / littleDeltaT;


				// Duration of single excursion
				if (fabs(nu(i-1))<1.0e-9) {
					ED(i-1) = 999999.0;
				}
				else {
					ED(i-1) = pf(i-1)/nu(i-1);
				}

				
				// Print the results to file right away...
				outputFile.setf( ios::fixed, ios::floatfield );

				outputFile << "#  " <<setprecision(2)<<setw(9)<<((stepsToStart+(i-1)*sampleFreq)*dt);
				
				if (beta(i-1)<0.0) { outputFile << "-"; }
				else { outputFile << " "; }
				outputFile <<setprecision(2)<<setw(11)<<fabs(beta(i-1));
				
				outputFile.setf( ios::scientific, ios::floatfield );
				if (pf(i-1)<0.0) { outputFile << "-"; }
				else { outputFile << " "; }
				outputFile <<setprecision(4)<<setw(16)<<fabs(pf(i-1));
				
				if (nu(i-1)<0.0) { outputFile << "-"; }
				else { outputFile << " "; }
				outputFile <<setprecision(5)<<setw(13)<<fabs(nu(i-1));
				
				if (ED(i-1)<0.0) { outputFile << "-"; }
				else { outputFile << " "; }
				outputFile <<setprecision(3)<<setw(10)<<fabs(ED(i-1));
				
				outputFile.setf( ios::fixed, ios::floatfield );
				outputFile<<"    #" << endln;
				outputFile.flush();


				// Inform the user
				char mystring[200];
				sprintf(mystring, " ... beta(G1)=%5.3f, beta(G2)=%5.3f, rate=%10.3e, rho=%10.4e",beta1,beta2,nu(i-1),a);
				opserr << mystring << endln;

			}
			else {
				outputFile << "#  First limit-state function did not converge.                       #" << endln;
			}


		} // Done looping over points on the time axis


		if (!oneFailed) {

			// Upper bound to probability of excursion during the interval 
			double Upper = 0.0;
			for (j=0; j<nu.Size(); j++) {
				Upper += nu(j)*Dt;
			}
			if (Upper >= 1.0) {
				Upper = 1.0;
			}


			// Approximation to true probability of failure
			double pTrue = 1.0 - exp(-Upper);
			if (pTrue >= 1.0) {
				pTrue = 1.0;
			}


			// Mean occupancy time
			double Eeta = 0.0;
			for (j=0; j<pf.Size(); j++) {
				Eeta += pf(j)*Dt;
			}


			// Matrix of intersection probabilities between 'g's (not g1 or g2) at all times
			for (k=0; k<pf.Size(); k++) {
				for (kk=0; kk<pf.Size(); kk++) {
					// Extract alpha vectors
					for (j=0; j<alpha.Size(); j++) {
						alpha_k(j) = allAlphas(j,k);
						alpha_kk(j) = allAlphas(j,kk);
					}
					a = 0.0;				// Interval start
					b = (alpha_k ^ alpha_kk);	// Interval end
					n_2 = 100;					// Half the number of intervals
					h = b-a;
					fa = functionToIntegrate(a, beta(k), beta(kk));
					fb = functionToIntegrate(b, beta(k), beta(kk));
					sum_fx2j = 0.0;
					sum_fx2j_1 = 0.0;
					for (int j=1;  j<=n_2;  j++) {
						sum_fx2j = sum_fx2j + functionToIntegrate(   (double) (a+(j*2)*h/(2*n_2)) , beta(k), beta(kk) );
						sum_fx2j_1 = sum_fx2j_1 + functionToIntegrate(   (double)(a+(j*2-1)*h/(2*n_2)) , beta(k), beta(kk)  );
					}
					sum_fx2j = sum_fx2j - functionToIntegrate((double)(b), beta(k), beta(kk));
					integral = h/(2*n_2)/3.0*(fa + 2.0*sum_fx2j + 4.0*sum_fx2j_1 + fb);
					Pmn2(k,kk) = aStdNormRV.getCDFvalue(-beta(k))*aStdNormRV.getCDFvalue(-beta(kk)) + integral;
				
				}
			}

			// Mean square of occupancy time
			double EetaSquared = 0.0;
			for (j=0; j<pf.Size(); j++) {
				for (k=0; k<pf.Size(); k++) {
					EetaSquared += Pmn2(j,k)*Dt;
				}
			}


			// Variance of occupancy time
			double VarEta = EetaSquared - Eeta*Eeta;


			// Cumulative area of excursion
			// (Require gfun-parameter sensitivity + evaluation of more lsf's)


			

			outputFile << "#                                                                     #" << endln;
			outputFile << "#                                                                     #" << endln;
			outputFile << "#  ACCUMULATED RESULTS:                                               #" << endln;
			outputFile << "#                                                                     #" << endln;
			outputFile << "#  Total time T: ...................................... " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<< T << "  #" << endln;
			
			outputFile.setf( ios::scientific, ios::floatfield );
			outputFile << "#  Upper bound to probability of excursion during T:... " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<< Upper << "  #" << endln;
			
			outputFile << "#  Approximation to true failure probability:.......... " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<< pTrue << "  #" << endln;
			
			outputFile << "#  Mean occupancy time: ............................... " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<< Eeta << "  #" << endln;

			outputFile << "#  Variance of occupancy time: ........................ " 
				<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<< VarEta << "  #" << endln;

			outputFile << "#                                                                     #" << endln;
			outputFile << "#######################################################################" << endln << endln << endln;
		}

	} // Done looping over limit-state functions


	// Clean up
	outputFile.close();


	// Print summary of results to screen (more here!!!)
	opserr << "Out-Crossing Analysis completed." << endln;

	return 0;
}



double
OutCrossingAnalysis::functionToIntegrate(double rho, double beta1, double beta2)
{
	double result;

	if (fabs(rho-1.0) < 1.0e-8) {
		result = 0.0;
	}
	else {
		double pi = 3.14159265358979;
		result = 1.0/(2.0*pi*sqrt(1.0-rho*rho)) 
			* exp(-(beta1*beta1+beta2*beta2-2.0*rho*beta1*beta2)
			/(2.0*(1.0-rho*rho)));
	}
	return result;
}
