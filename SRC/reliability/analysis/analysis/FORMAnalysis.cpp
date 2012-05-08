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
                                                                        
// $Revision: 1.18 $
// $Date: 2008-08-27 17:12:23 $
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

#include <RandomVariableIter.h>
#include <LimitStateFunctionIter.h>

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
                           FunctionEvaluator *passedFunctionEvaluator,
						   ProbabilityTransformation *passedProbabilityTransformation,
						   Tcl_Interp *passedInterp, TCL_Char *passedFileName, int p_relSensTag)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	theFindDesignPointAlgorithm = passedFindDesignPointAlgorithm;
    theFunctionEvaluator			= passedFunctionEvaluator;
	theProbabilityTransformation = passedProbabilityTransformation;
	strcpy(fileName,passedFileName);
	relSensTag = p_relSensTag;
	interp = passedInterp;
}


FORMAnalysis::~FORMAnalysis()
{

}



int 
FORMAnalysis::analyze()
{

	// Alert the user that the FORM analysis has started
	opserr << "FORM Analysis is running ... " << endln;

	// Declare variables used in this method
	double stdv, mean, Go, Glast, beta, pf1;
	int numberOfSteps, numberOfEvaluations;
	static NormalRV aStdNormRV(1,0.0,1.0);

	// reliability domain info
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();
	LimitStateFunction *theLimitStateFunction;
	
	Vector delta(numRV);
	Vector eta(numRV);
	Vector kappa(numRV);
	
	// Open output file
	ofstream outputFile( fileName, ios::out );

	// Loop over number of limit-state functions and perform FORM analysis
	//LimitStateFunctionIter lsfIter = theReliabilityDomain->getLimitStateFunctions();
	//while ((theLimitStateFunction = lsfIter()) != 0) {
	for (int lsf = 0; lsf < numLsf; lsf++ ) {
		theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtrFromIndex(lsf);
		int lsfTag = theLimitStateFunction->getTag();

		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsfTag << endln;
		Tcl_SetVar2Ex(interp,"RELIABILITY_lsf",NULL,Tcl_NewIntObj(lsfTag),TCL_NAMESPACE_ONLY);

		// Set tag of "active" limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsfTag);
	  
		outputFile << "#######################################################################" << endln;
		outputFile << "#  FORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
				   << setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsfTag <<"            #" << endln;
		outputFile << "#                                                                     #" << endln;

		// Find the design point
		if (theFindDesignPointAlgorithm->findDesignPoint() < 0) {
            opserr << "FORMAnalysis::analyze() - failed while finding the" << endln
                   << " design point for limit-state function number " << lsfTag << "." << endln;
			
			outputFile << "#  No convergence!                                                    #" << endln;
			outputFile << "#  Results shown for last iteration of analysis                       #" << endln;
		}

		// Get results from the "find desingn point algorithm"
		const Vector &xStar = theFindDesignPointAlgorithm->get_x();
		const Vector &uStar = theFindDesignPointAlgorithm->get_u();
		const Vector &alpha = theFindDesignPointAlgorithm->get_alpha();
		const Vector &gamma = theFindDesignPointAlgorithm->get_gamma();
		numberOfSteps		= theFindDesignPointAlgorithm->getNumberOfSteps();
		Go					= theFindDesignPointAlgorithm->getFirstGFunValue();
		Glast				= theFindDesignPointAlgorithm->getLastGFunValue();
		//const Vector &uSecondLast = theFindDesignPointAlgorithm->getSecondLast_u();
		//const Vector &alphaSecondLast = theFindDesignPointAlgorithm->getSecondLast_alpha();
		//const Vector &lastSearchDirection = theFindDesignPointAlgorithm->getLastSearchDirection();
		numberOfEvaluations	= theFindDesignPointAlgorithm->getNumberOfEvaluations();
	  
		// Postprocessing
		beta = alpha ^ uStar;
		pf1 = 1.0 - aStdNormRV.getCDFvalue(beta);
        RandomVariable *theRV;
	  
		// Reliability sensitivity analysis wrt. mean/stdv
		if (relSensTag == 1) {
			double dBetaDmean;
			double dBetaDstdv;
			Vector DuStarDmean(numRV);
			Vector DuStarDstdv(numRV);

			RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
			for (int j = 0; j < numRV; j++) {
                theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
			//while ((theRV = rvIter()) != 0) {
				//int j = theRV->getIndex();
				int rvTag = theRV->getTag();
				//int j = theReliabilityDomain->getRandomVariableIndex(rvTag);
				
				DuStarDmean = theProbabilityTransformation->meanSensitivityOf_x_to_u(xStar,rvTag);
				DuStarDstdv = theProbabilityTransformation->stdvSensitivityOf_x_to_u(xStar,rvTag);
				dBetaDmean = alpha^DuStarDmean;
				dBetaDstdv = alpha^DuStarDstdv;
				stdv = theRV->getStdv();
				mean = theRV->getMean();
				
				delta(j) = stdv * dBetaDmean;
				eta(j) = stdv * dBetaDstdv;
				// (Kappa is the sensitivity wrt. the coefficient of variation)
				if (mean != 0.0)
					kappa(j) = -dBetaDmean*stdv/((stdv/mean)*(stdv/mean))+dBetaDstdv*mean;
				else
					kappa(j) = 0.0;
			}
		}
		
        // store key results using function evaluator
        for (int j = 0; j < numRV; j++) {
            theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
            int rvTag = theRV->getTag();
            theFunctionEvaluator->setResponseVariable("gammaFORM", lsfTag, rvTag, gamma(j));
            theFunctionEvaluator->setResponseVariable("alphaFORM", lsfTag, rvTag, alpha(j));
        }
        theFunctionEvaluator->setResponseVariable("betaFORM", lsfTag, beta);


        // report the rest to output file
	  outputFile << "#  Limit-state function value at start point: ......... " 
		     <<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<Go 
		     << "  #" << endln;
	  outputFile << "#  Limit-state function value at end point: ........... " 
		     <<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<Glast 
		     << "  #" << endln;
	  outputFile << "#  Number of steps: ................................... " 
		     <<setiosflags(ios::left)<<setw(12)<<numberOfSteps
		     << "  #" << endln;
	  outputFile << "#  Number of g-function evaluations: .................. " 
		     <<setiosflags(ios::left)<<setw(12)<<numberOfEvaluations 
		     << "  #" << endln;
	  outputFile << "#  Reliability index beta: ............................ " 
		     <<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<beta 
		     << "  #" << endln;
	  outputFile << "#  FO approx. probability of failure, pf1: ............ " 
		     <<setiosflags(ios::left)<<setiosflags(ios::scientific)<<setprecision(5)<<setw(12)<<pf1 
		     << "  #" << endln;
	  outputFile << "#                                                                     #" << endln;
	  outputFile << "# rvtag   x*          u*         alpha    gamma    delta    eta       #" << endln;
	  outputFile.setf( ios::scientific, ios::floatfield );
	  
	  RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
	  //for (int i=0;  i<xStar.Size(); i++) {
	  while ((theRV = rvIter()) != 0) {
	    //int i = theRV->getIndex();
	    int tag = theRV->getTag();
	    int i = theReliabilityDomain->getRandomVariableIndex(tag);
	    outputFile << "#  " <<setw(3)<<tag<<" ";
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
	  opserr << "Done analyzing limit-state function " << lsfTag << ", beta=" << beta << endln;
	}
	
	// Clean up
	outputFile.close();

	// Print summary of results to screen (more here!!!)
	opserr << "FORMAnalysis completed." << endln;

	return 0;
}

