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
// $Date: 2008-05-27 15:42:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/OutCrossingAnalysis.cpp,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//
// Quan and Michele Feb 2006

#include <OutCrossingAnalysis.h>
#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <FindDesignPointAlgorithm.h>
#include <CorrelatedStandardNormal.h>
#include <math.h>
#include <tcl.h>
#include <string.h>
#include <NormalRV.h>
#include <LimitStateFunctionIter.h>

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
				FunctionEvaluator *theGFunEval,
				GradientEvaluator *theSensEval,
				FindDesignPointAlgorithm *theFindDesPt,
				int pAnalysisType,
				int p_stepsToStart,
				int p_stepsToEnd,
				int p_sampleFreq,
				int p_impulseFreq,
				double p_littleDeltaT,
				TCL_Char *passedFileName,
				double p_integralTolerance,
				bool p_useFirstDesignPt)
:ReliabilityAnalysis()
{
	theFindDesignPointAlgorithm = theFindDesPt;
	theReliabilityDomain = theRelDom;
	theGFunEvaluator = theGFunEval;
	theGradGEvaluator = theSensEval;
	analysisType = pAnalysisType;
	stepsToStart = p_stepsToStart;  //  integration time steps
	stepsToEnd = p_stepsToEnd;    // integration time steps
	sampleFreq = p_sampleFreq;   // number of integration time steps between two meanOutcrossingAnalysis.
	littleDeltaT = p_littleDeltaT;
	strcpy(fileName,passedFileName);
	impulseFreq = p_impulseFreq;

	integralTolerance = p_integralTolerance; // default value. 
	useFirstDesignPt = p_useFirstDesignPt; // default 
}

OutCrossingAnalysis::~OutCrossingAnalysis()
{

}

int 
OutCrossingAnalysis::analyze(void)
{

	// Alert the user that the analysis has started
	opserr << "Out-Crossing Analysis is running ... " << endln;

	// Declare variables used in this method
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	static NormalRV aStdNormRV(1,0.0,1.0);
	int numVel, i, j, k, kk, nodeNumber, dofNumber; 
	double dgduValue, accuSum;
	Vector uStar2(numRV);
	Vector alpha(numRV), alpha2(numRV), alpha_k(numRV), alpha_kk(numRV);
	double beta1, beta2;
	double pf1, pf2;
//	int n_2;
//	double a, b, integral, h, fa, fb, sum_fx2j, sum_fx2j_1, Pmn1;
//	double a, b, c, integral, fa, fb, fc;
	double a, Pmn1;
	char string[500];


	// Determine number of points
	double nsteps = stepsToEnd-stepsToStart;       //total integration steps
	int numPoints = (int)floor(nsteps/sampleFreq);  //  how many time points of computation of outcrossingRate
	numPoints++;   
	double dt = theGFunEvaluator->getDt();  // integration time step
	double Dt = dt*impulseFreq;              // time between impuses.
	double T = dt*sampleFreq*numPoints;     // total time of last analysis


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
	LimitStateFunctionIter lsfIter = theReliabilityDomain->getLimitStateFunctions();
	LimitStateFunction *theLimitStateFunction;
	//for (lsf=1; lsf<=numLsf; lsf++ ) {
	while ((theLimitStateFunction = lsfIter()) != 0) {

	  int lsf = theLimitStateFunction->getTag();

		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsf << endln;


		// Start printing results to the output file
		outputFile << "########################################################R##############" << endln;
		outputFile << "#  OUT-CROSSING RESULTS, LIMIT-STATE FUNCTION NUMBER      "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"        #" << endln;
		outputFile << "#                                                                     #" << endln;



		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);


		// Loop over the intervals where probabilities are to be computed
		bool oneFailed = false;
		bool DSPTfailed;
		for (i=1; i<=numPoints; i++) {

			// First assume that we find the design point(s)
			DSPTfailed = false;

			// Set 'nsteps' in the GFunEvaluator 
			theGFunEvaluator->setNsteps(stepsToStart+(i-1)*sampleFreq);


			// Inform the user
			char printTime[23];
			sprintf(printTime, "%14.8f", ((stepsToStart+(i-1)*sampleFreq)*dt) );
			if (analysisType == 1) {
				strcpy(string,theLimitStateFunction->getExpression());
				opserr << " ...evaluating -G1=" << string << " at time " << printTime << " ..." << endln;
			}
			else {
				opserr << " ...evaluating performance function at time " << printTime << " ..." << endln;
			}


			// Find the design point for the original limit-state function
			if (theFindDesignPointAlgorithm->findDesignPoint() < 0){
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
				const Vector &uStar = theFindDesignPointAlgorithm->get_u();
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

				//  ====== possible options: using first design point as the beginning of this second search.. 
					if (useFirstDesignPt) {

						Vector xStar = theFindDesignPointAlgorithm->get_x();
                        // KRM -- parameters now handle all values
						//theFindDesignPointAlgorithm->setStartPt(&xStar);
						
					}
					
				 // ============ by default,start from origin of the standard normal space ====================	


					// Get the 'dgdu' vector from the sensitivity evaluator
					// (The returned matrix contains 'node#' 'dir#' 'dgdu' in rows)
					
					// This should be retooled -- MHS 10/7/2011
					//const Matrix &DgDdispl = theGradGEvaluator->getDgDdispl();
					Matrix DgDdispl(5,5); // So that it compiles


					// Add extra term to limit-state function
					// (should add an alternative option where user give additional limit-state functions)
					// (... because this works only when g is linear in u)
					numVel = DgDdispl.noRows();
					accuSum = 0.0;

					char expressionPtr[200];

					for (j=0; j<numVel; j++) {

						nodeNumber = (int)DgDdispl(j,0);
						dofNumber = (int)DgDdispl(j,1);
						dgduValue = DgDdispl(j,2);

						char expression[200];
						sprintf(expression,"+(%20.14e)*(%20.14e)*\\$ud(%d,%d)",littleDeltaT, dgduValue, nodeNumber, dofNumber);
						strcpy(expressionPtr,expression);

						// Add it to the limit-state function
						// There's a better way to do this -- MHS 10/7/2011
						//theLimitStateFunction->addExpression(expressionPtr);
					}

					// Inform the user
					strcpy(string,theLimitStateFunction->getExpression());
					opserr << " ...evaluating G2=" << string << endln;


					// Find the design point for the edited limit-state function
					if (theFindDesignPointAlgorithm->findDesignPoint() < 0) {
						opserr << "OutCrossingAnalysis::analyze() - failed while finding the" << endln
							<< " design point for limit-state function number " << lsf << "." << endln;
						DSPTfailed = true;
						oneFailed = true;
					}


					// Zero out the added expression in the limit-state function
					// There's a better way to do this -- MHS 10/7/2011
					//theLimitStateFunction->removeAddedExpression();


					if (!DSPTfailed) {

						// Get results from the "find design point algorithm"
						uStar2 = theFindDesignPointAlgorithm->get_u();
						alpha2 = theFindDesignPointAlgorithm->get_alpha();


						// Postprocessing (remember; here is an assumption that the mean point is in the safe domain)
						beta2 = (alpha2 ^ uStar2);  /////// ??????????????????????????
						pf2 = 1.0 - aStdNormRV.getCDFvalue(beta2);
					}
					else {
						outputFile << "#  Second limit-state function did not converge.                      #" << endln;
						
						// NOTE - need to provide a value for beta2 or initialize in declaration
					}
				
					// Post-processing to find parallel system probability
					// use CorrelatedStandardNormal class
					CorrelatedStandardNormal phi2(0.0);
					phi2.setCorrelation( alpha ^ alpha2 );
					Pmn1 = phi2.getCDF(-beta1,-beta2);
					
					a = -(alpha ^ alpha2);	// Interval start  ??????????????????????
					/*
					c = 0.0;				// Interval end
					b = (c+a)/2.0;
					fa = functionToIntegrate(a,beta1,beta2);
					fb = functionToIntegrate(b,beta1,beta2);
					fc = functionToIntegrate(c,beta1,beta2);
					integral = this->getAdaptiveIntegralValue(integralTolerance,a,c,fa,fb,fc,beta1,beta2);
					Pmn1 = aStdNormRV.getCDFvalue(-beta1)*aStdNormRV.getCDFvalue(-beta2) - integral;
					*/
				
				}
				else {
					// Use Heonsang's method to find new alpha and beta
					beta2 = -beta1;
					alpha2(0) = alpha(0) - littleDeltaT/Dt * ( alpha(0)-0.0 );
					for (j=1; j<alpha.Size(); j++) {
						alpha2(j) = alpha(j) - littleDeltaT/Dt * ( alpha(j)-alpha(j-1) );
					}
					
					// Post-processing to find parallel system probability
					a = -1.0*(alpha ^ alpha2)/alpha2.Norm();  
				
					double pi = acos(-1.0);
					Pmn1 = 1.0/(2.0*pi) * exp(-beta2*beta2*0.5) * (asin(a) + pi/2.0);

				}


				// POST-PROCESSING

				// Mean out-crossing rate
				if (Pmn1<1.0e-10) {
					opserr << "WARNING: Zero or negative parallel probability: " << endln
						<< " The correlation is probably too high! " << endln;
				}
				nu(i-1) = Pmn1 / littleDeltaT;
				
				opserr<<"nu(i-1) is:" <<nu(i-1)<<endln;

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
				outputFile <<setprecision(7)<<setw(11)<<fabs(beta(i-1));
				
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
				char mystring[400];
				sprintf(mystring, " ... beta(G1)=%16.8f, beta(G2)=%16.8f, rate=%16.8e, rho=%16.8e",beta1,beta2,nu(i-1),a);
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
			CorrelatedStandardNormal phi2(0.0);
			for (k=0; k<pf.Size(); k++) {
				for (kk=0; kk<pf.Size(); kk++) {
					// Extract alpha vectors
					for (j=0; j<alpha.Size(); j++) {
						alpha_k(j) = allAlphas(j,k);
						alpha_kk(j) = allAlphas(j,kk);
					}
					
					// use CorrelatedStandardNormal class
					phi2.setCorrelation( alpha_k ^ alpha_kk );
					Pmn2(k,kk) = phi2.getCDF(-beta(k),-beta(kk));
					
					/*
					a = 0.0;
					c = (alpha_k ^ alpha_kk);				// Interval end
					b = (c+a)/2.0;
					fa = functionToIntegrate(a,beta(k),beta(kk));
					fb = functionToIntegrate(b,beta(k),beta(kk));
					fc = functionToIntegrate(c,beta(k),beta(kk));
					integral = this->getAdaptiveIntegralValue(integralTolerance, a,c,fa,fb,fc,beta(k),beta(kk));

					Pmn2(k,kk) = aStdNormRV.getCDFvalue(-beta(k))*aStdNormRV.getCDFvalue(-beta(kk)) + integral;
					*/
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


/*
//Quan & Michele January '06
double 
OutCrossingAnalysis::getAdaptiveIntegralValue(double tol, double lowerBound, double upperBound, double fa, double fb, double fc, double beta1, double beta2 ){
	double b= (lowerBound+upperBound)/2.0;
	double fd = functionToIntegrate((lowerBound+b)/2.0, beta1, beta2);
	double fe = functionToIntegrate((upperBound+b)/2.0, beta1, beta2);
//	numOfEvaluations +=2;
	double S1 = (upperBound-lowerBound)*(fa+4*fb+fc)/6.0;
	double S2 = (upperBound-lowerBound)*(fa+4*fd+2*fb+4*fe+fc)/12.0;
	if (fabs(S2-S1)/15.0 <= tol) 
		return (16.0*S2-S1)/15.0;
	else {
		
		double int1 = getAdaptiveIntegralValue( tol/2.0,lowerBound, b, fa, fd, fb, beta1, beta2);
		double int2 = getAdaptiveIntegralValue(  tol/2.0, b, upperBound, fb, fe, fc, beta1, beta2);
		return int1+int2;
	
	}
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
*/
