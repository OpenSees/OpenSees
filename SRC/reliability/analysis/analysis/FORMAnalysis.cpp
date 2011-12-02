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
                                                                        
// $Revision: 1.15 $
// $Date: 2008-04-10 18:10:29 $
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
						   ProbabilityTransformation *passedProbabilityTransformation,
						   Tcl_Interp *passedInterp, TCL_Char *passedFileName, int p_relSensTag)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	theFindDesignPointAlgorithm = passedFindDesignPointAlgorithm;
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
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	static NormalRV aStdNormRV(1,0.0,1.0,0.0);
	Vector delta(numRV); 
	Vector eta(numRV);
	Vector kappa(numRV);


	// Open output file
	ofstream outputFile( fileName, ios::out );

	LimitStateFunctionIter lsfIter = theReliabilityDomain->getLimitStateFunctions();
	LimitStateFunction *theLimitStateFunction;

	// Loop over number of limit-state functions and perform FORM analysis
	//for (lsf=1; lsf<=numLsf; lsf++ ) {
	while ((theLimitStateFunction = lsfIter()) != 0) {

	  int lsf = theLimitStateFunction->getTag();

	  // Inform the user which limit-state function is being evaluated
	  opserr << "Limit-state function number: " << lsf << endln;
	  Tcl_SetVar2Ex(interp,"RELIABILITY_lsf",NULL,Tcl_NewIntObj(lsf),TCL_NAMESPACE_ONLY);

	  // Set tag of "active" limit-state function
	  theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);
	  
	  outputFile << "#######################################################################" << endln;
	  outputFile << "#  FORM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
		     <<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"            #" << endln;
	  outputFile << "#                                                                     #" << endln;
	  
	  // Find the design point
	  if (theFindDesignPointAlgorithm->findDesignPoint(theReliabilityDomain) < 0){
	    opserr << "FORMAnalysis::analyze() - failed while finding the" << endln
		   << " design point for limit-state function number " << lsf << "." << endln;
	    
			
			outputFile << "#  No convergence!                                                    #" << endln;
			outputFile << "#  Results shown for last iteration of analysis                       #" << endln;
	  }

	  // Get results from the "find desingn point algorithm"
	  const Vector &xStar = theFindDesignPointAlgorithm->get_x();
	  const Vector &uStar = theFindDesignPointAlgorithm->get_u();
	  const Vector &alpha = theFindDesignPointAlgorithm->get_alpha();
	  const Vector &gamma = theFindDesignPointAlgorithm->get_gamma();
	  numberOfSteps					= theFindDesignPointAlgorithm->getNumberOfSteps();
	  Go					= theFindDesignPointAlgorithm->getFirstGFunValue();
	  Glast				= theFindDesignPointAlgorithm->getLastGFunValue();
	  const Vector &uSecondLast = theFindDesignPointAlgorithm->getSecondLast_u();
	  const Vector &alphaSecondLast = theFindDesignPointAlgorithm->getSecondLast_alpha();
	  const Vector &lastSearchDirection = theFindDesignPointAlgorithm->getLastSearchDirection();
	  numberOfEvaluations	= theFindDesignPointAlgorithm->getNumberOfEvaluations();
	  
	  // Postprocessing
	  beta = alpha ^ uStar;
	  pf1 = 1.0 - aStdNormRV.getCDFvalue(beta);
	  
	  
	  // Reliability sensitivity analysis wrt. mean/stdv
	  if (relSensTag == 1) {
	    double dBetaDmean;
	    double dBetaDstdv;
	    Vector DuStarDmean(numRV);
	    Vector DuStarDstdv(numRV);
	    
	    RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
	    RandomVariable *theRV;
	    //for ( int j=1; j<=numRV; j++ ) {
	    while ((theRV = rvIter()) != 0) {
	      //int j = theRV->getIndex();
	      int rvTag = theRV->getTag();
	      int j = theReliabilityDomain->getRandomVariableIndex(rvTag);
	      DuStarDmean = theProbabilityTransformation->meanSensitivityOf_x_to_u(xStar,rvTag);
	      DuStarDstdv = theProbabilityTransformation->stdvSensitivityOf_x_to_u(xStar,rvTag);
	      dBetaDmean = alpha^DuStarDmean;
	      dBetaDstdv = alpha^DuStarDstdv;
	      stdv = theRV->getStdv();
	      mean = theRV->getMean();
	      delta(j) = stdv * dBetaDmean;
	      eta(j) = stdv * dBetaDstdv;
	      // (Kappa is the sensitivity wrt. the coefficient of variation)
	      if (mean != 0.0) {
		kappa(j) = -dBetaDmean*stdv/((stdv/mean)*(stdv/mean))+dBetaDstdv*mean;
	      }
	      else {
		kappa(j) = 0.0;
	      }
	    }
	    
	    // Don't scale them: WHY????? they give gibberish otherwise
	    delta /= delta.Norm();
	    eta /= eta.Norm();
	    kappa /= kappa.Norm();
	  }

	  
	  // Store key results in the limit-state functions
	  theLimitStateFunction->setFORM_beta(beta);
	  theLimitStateFunction->setFORM_pf1(pf1);
	  theLimitStateFunction->setFORM_x(xStar);
	  theLimitStateFunction->setFORM_u(uStar);
	  theLimitStateFunction->setFORM_alpha(alpha);
	  theLimitStateFunction->setFORM_gamma(gamma);
	  theLimitStateFunction->setFORM_numSteps(numberOfSteps);
	  theLimitStateFunction->setFORM_startGFun(Go);
	  theLimitStateFunction->setFORM_endGFun(Glast);
	  theLimitStateFunction->setSecondLast_u(uSecondLast);
	  theLimitStateFunction->setSecondLast_alpha(alphaSecondLast);
	  theLimitStateFunction->setLastSearchDirection(lastSearchDirection);


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
		     <<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf1 
		     << "  #" << endln;
	  outputFile << "#                                                                     #" << endln;
	  outputFile << "# rvtag   x*          u*         alpha    gamma    delta    eta       #" << endln;
	  outputFile.setf( ios::scientific, ios::floatfield );
	  
	  RandomVariableIter rvIter = theReliabilityDomain->getRandomVariables();
	  RandomVariable *theRV;
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
	  opserr << "Done analyzing limit-state function " << lsf << ", beta=" << beta << endln;
	}
	
	// Clean up
	outputFile.close();

	// Print summary of results to screen (more here!!!)
	opserr << "FORMAnalysis completed." << endln;

	return 0;
}

