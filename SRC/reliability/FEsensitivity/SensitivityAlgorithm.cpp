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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-03-04 00:46:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/SensitivityAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SensitivityAlgorithm.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <EquiSolnAlgo.h>
#include <ReliabilityDomain.h>
#include <RandomVariablePositioner.h>


SensitivityAlgorithm::SensitivityAlgorithm(ReliabilityDomain *passedReliabilityDomain,
										   EquiSolnAlgo *passedAlgorithm,
										   SensitivityIntegrator *passedSensitivityIntegrator,
										   int passedAnalysisTypeTag)
{
	// The reliability domain is needed to get hold 
	// of the random variable positioners:
	theReliabilityDomain = passedReliabilityDomain;

	// The finite element equation solution algorithm is 
	// needed to get hold of the system of equations (SOE):
	theAlgorithm = passedAlgorithm;

	// The sensitivity integrator is needed to assemble the 
	// new right-hand side of the system of equations:
	theSensitivityIntegrator = passedSensitivityIntegrator;

	// Tag to tell whether grads should be computed at each step
	// and whether they should be computed wrt. random variables
	analysisTypeTag = passedAnalysisTypeTag;

}




SensitivityAlgorithm::~SensitivityAlgorithm()
{
}



int 
SensitivityAlgorithm::computeSensitivities(void)
{
	// Meaning of analysisTypeTag:
	// 1: compute at each step wrt. random variables
	// 2: compute at each step wrt. parameters
	// 3: compute by command wrt. random variables
	// 4: compute by command wrt. parameters

	
	// Initial declarations
	int numGrads, numPos, i, gradNumber, posNumber;
	RandomVariablePositioner *theRandomVariablePositioner;
	ParameterPositioner *theParameterPositioner;


	// Get pointer to the system of equations (SOE)
	LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();


	// Get pointer to incremental integrator
	IncrementalIntegrator *theIncInt = theAlgorithm->getIncrementalIntegratorPtr();


	// Form current tangent at converged state
	// (would be nice with an if-statement here in case
	// the current tangent is already formed)
	if (theIncInt->formTangent(CURRENT_TANGENT) < 0){
		opserr << "WARNING SensitivityAlgorithm::computeGradients() -";
		opserr << "the Integrator failed in formTangent()\n";
		return -1;
	}


	// Get number of random variables and random variable positioners
	if (analysisTypeTag==1 || analysisTypeTag==3) {
		numGrads = theReliabilityDomain->getNumberOfRandomVariables();
		numPos = theReliabilityDomain->getNumberOfRandomVariablePositioners();
	}
	else {
		numPos = theReliabilityDomain->getNumberOfParameterPositioners();
		numGrads = numPos;
	}
	

	// Zero out the old right-hand side of the SOE
	theSOE->zeroB();
		

	// Form the part of the RHS which are indepent of parameter
	theSensitivityIntegrator->formIndependentSensitivityRHS();


	// Loop where one gradient vector is computed per iteration
	for (gradNumber=1; gradNumber<=numGrads; gradNumber++ )  {


		// Clear sensitivity flags
		for (i=1; i<=numPos; i++ ) {
			if (analysisTypeTag==1 || analysisTypeTag==3) {
				theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
				theRandomVariablePositioner->activate(false);
			}
			else {
				theParameterPositioner = theReliabilityDomain->getParameterPositionerPtr(i);
				theParameterPositioner->activate(false);
			}
		}


		// Loop over positioners to check which are "active" this time
		// and set flag for phase 1
		for (posNumber=1; posNumber<=numPos; posNumber++ ) {

			if (analysisTypeTag==1 || analysisTypeTag==3) {
				// Get the random positioner and its rv#
				theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(posNumber);
				int rvNumber = theRandomVariablePositioner->getRvNumber();

				// Check if rv# of this rv positioner is the same as the gradient#
				if ( rvNumber==gradNumber ) {

					// Set sensitivity flag so that this one contributes to the RHS
					theRandomVariablePositioner->activate(true);

				} // End if rv# == gradient#
			}
			else {
				theParameterPositioner = theReliabilityDomain->getParameterPositionerPtr(posNumber);
				theParameterPositioner->activate(true);
			}

		} // End loop over r.v. positioners to check if it contributes


		// Zero out the old right-hand side
		theSOE->zeroB();

		
		// Form new right-hand side
		theSensitivityIntegrator->formSensitivityRHS(gradNumber);


		// Solve the system of equation with the new right-hand side
		theSOE->solve();


		// Save 'v' to the nodes for a "sensNodeDisp node? dof?" command
		theSensitivityIntegrator->saveSensitivity( theSOE->getX(), gradNumber, numGrads );


		// Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
		theSensitivityIntegrator->commitSensitivity(gradNumber, numGrads);

	} // End loop for each gradient

    return 0;
}

bool 
SensitivityAlgorithm::shouldComputeAtEachStep(void)
{
	if (analysisTypeTag==1 || analysisTypeTag==2) {
		return true;
	}
	else {
		return false;
	}
}