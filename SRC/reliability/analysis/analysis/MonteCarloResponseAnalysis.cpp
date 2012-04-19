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
// $Date: 2008-03-13 22:22:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/MonteCarloResponseAnalysis.cpp,v $


//
// Written by Quan Gu (qgu@ucsd.edu)
//

#include "MonteCarloResponseAnalysis.h"
#include <ProbabilityTransformation.h>
#include <NatafProbabilityTransformation.h>
#include <RandomNumberGenerator.h>
#include <RandomVariable.h>
//#include <RandomVariablePositioner.h>
#include <NormalRV.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>

using std::ifstream;
using std::ios;

//using std::ios;
MonteCarloResponseAnalysis::MonteCarloResponseAnalysis(
						ReliabilityDomain *passedReliabilityDomain,
						Tcl_Interp *passedTclInterp,
						ProbabilityTransformation *passedProbabilityTransformation,
						RandomNumberGenerator *passedRandomNumberGenerator,
						int passedNumberOfSimulations,
						int passedPrintFlag,
						TCL_Char *passedFileName,
						TCL_Char *pTclFileToRunFileName,
						int pSeed
						)
{
	theReliabilityDomain = passedReliabilityDomain;
	theTclInterp = passedTclInterp;
	theProbabilityTransformation = passedProbabilityTransformation;
	theRandomNumberGenerator = passedRandomNumberGenerator;
	numberOfSimulations = passedNumberOfSimulations;
	printFlag = passedPrintFlag;
	strcpy(fileName,passedFileName);
	seed = pSeed;

	if (pTclFileToRunFileName !=0){
		tclFileToRun=new char [30];
		strcpy(tclFileToRun,pTclFileToRunFileName);
	}
	else tclFileToRun = 0;

}



MonteCarloResponseAnalysis::~MonteCarloResponseAnalysis()
{
	if (tclFileToRun !=0) delete [] tclFileToRun;
}

int MonteCarloResponseAnalysis::analyze(){


	opserr << "Monte Carlo Response Analysis is running ... " << endln;
	

	int result;
	int kk = 0;
	bool isFirstSimulation = true;
	//int seed = 1;
	if (printFlag ==2) {
	
	  	 // check whether the restart file '_restart.tmp' exist, (this file is wrote by openSees only, not by user)
		 //                      if yes, read data;{ success reading: set values above; otherwise: do noting} 
    	 //		                 if no, create file '_restart.tmp'
	     //	data format:                 seed           numOfGFunEvaluations
			
	  ifstream inputFile( "_restart.tmp", ios::in );
	  if (!inputFile) {
	    // file doesn't exist; so it's safe to create it and write to it
	    
	  }
	  else { 
	    inputFile >> seed;
	    inputFile >> kk;
	    inputFile.close();
	    isFirstSimulation = false;
	  }
	  
	}
	
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();

	Vector x(numRV);

	Vector u(numRV);
	Vector randomArray(numRV);

	ofstream *outputFile = 0;


	
	// Prepare output file
	ofstream resultsOutputFile( fileName, ios::out );



	while( kk< numberOfSimulations){ // && govCov>targetCOV || k<=2) ) {

		// Keep the user posted
		if (printFlag == 1 || printFlag == 2) {
			opserr << "Sample #" << kk << ":" << endln;
//			resultsOutputFile<< "Sample #" << kk << ":" << endln;

		}

		
		// Create array of standard normal random numbers
		if (isFirstSimulation) {
			result = theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numRV,seed);
		}
		else {
			result = theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numRV);
		}
		seed = theRandomNumberGenerator->getSeed();
		if (result < 0) {
			opserr << "MonteCarloResponseAnalysis::analyze() - could not generate" << endln
				<< " random numbers for simulation." << endln;
			return -1;
		}
		randomArray = theRandomNumberGenerator->getGeneratedNumbers();


		// Compute the point in standard normal space

		
		u = randomArray;   // Quan

		// Transform into original space
		/*
		result = theProbabilityTransformation->set_u(u);
		if (result < 0) {
			opserr << "MonteCarloResponseAnalysis::analyze() - could not " << endln
				<< " set the u-vector for xu-transformation. " << endln;
			return -1;
		}

		
		result = theProbabilityTransformation->transform_u_to_x();
		if (result < 0) {
			opserr << "MonteCarloResponseAnalysis::analyze() - could not " << endln
				<< " transform u to x. " << endln;
			return -1;
		}
		x = theProbabilityTransformation->get_x();
		*/

		int result = theProbabilityTransformation->transform_u_to_x(u, x);
		if (result < 0) {
			opserr << "MonteCarloResponseAnalysis::analyze() - could not " << endln
			       << " transform u to x. " << endln;
			return -1;
		}


      // ------ here recorder x ----
//		opserr << "RV x is: " << x << endln;
//		resultsOutputFile << "RV x is: " <<endln;
		resultsOutputFile.precision(15);
		for (int ii=0;ii<numRV;ii++)
		   resultsOutputFile << x(ii)<<endln ;
//		resultsOutputFile <<endln;


		// --------------- update structure parameter -----------------
		
		// No longer using positioners -- MHS 4/2012
		/*
		int numberOfRandomVariablePositioners = theReliabilityDomain->getNumberOfRandomVariablePositioners();
		RandomVariablePositioner *theRandomVariablePositioner;
		int rvNumber;
		//FMK
		for (int i=1 ; i<=numberOfRandomVariablePositioners ; i++ )  {
			theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
			rvNumber				= theRandomVariablePositioner->getRvIndex();
			theRandomVariablePositioner->update(x(rvNumber));
		}
		*/

		// ---------------------- run tcl file and  recorder ---------------------


		if (tclFileToRun != 0) {     
			char theRevertToStartCommand[10] = "reset";
			Tcl_Eval( theTclInterp, theRevertToStartCommand );
			char theWipeAnalysis[15] = "wipeAnalysis";
			Tcl_Eval( theTclInterp, theWipeAnalysis );

			if(Tcl_EvalFile(theTclInterp, tclFileToRun) !=TCL_OK){
				opserr<<"MonteCarloResponseAnalysis: the file "<<tclFileToRun<<" can not be run!"<<endln;
				exit(-1);
			}  //if

		}  //if

		kk++;
		isFirstSimulation = false;	



		if (printFlag ==2){
 
			// write necessary data into file '_restart.tmp' .... close file
			ofstream resultsOutputFile5( "_restart.tmp");
			resultsOutputFile5<< seed        <<endln;
			resultsOutputFile5<< kk <<endln;
			
			resultsOutputFile5.flush();
			resultsOutputFile5.close();
		}




		
		
	}// while 
		
		
	opserr << endln;


	// Delete possible 'new' objects
	if (outputFile != 0) {
		delete outputFile;
	}


	// Print summary of results to screen 
	opserr << "Simulation Analysis completed." << endln;

	// Clean up
	resultsOutputFile.close();
//	delete theMatrixOperations;

	return 0;


};
