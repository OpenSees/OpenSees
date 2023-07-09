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
                                                                        
// $Revision: 1.13 $
// $Date: 2008-08-26 16:15:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/SensitivityAlgorithm.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SensitivityAlgorithm.h>
#include <ReliabilityDomain.h>
#include <Integrator.h>

#include <LinearSOE.h>
#include <EquiSolnAlgo.h>
#include <Vector.h>
#include <Domain.h>
#include <Parameter.h>
#include <ParameterIter.h>


SensitivityAlgorithm::SensitivityAlgorithm(Domain *passedDomain,
					   EquiSolnAlgo *passedAlgorithm,
					   Integrator *passedSensitivityIntegrator,
					   int passedAnalysisTypeTag):
  theDomain(passedDomain), theAlgorithm(passedAlgorithm), 
  analysisTypeTag(passedAnalysisTypeTag)
{

   
}

SensitivityAlgorithm::~SensitivityAlgorithm()
{

}
/*
int 
SensitivityAlgorithm::computeSensitivities(void)
{
   
  if (theAlgorithm == 0) {
    opserr << "ERROR the FE algorithm must be defined before ";
    opserr << "the sensitivity algorithm\n";
    return -1;
  }

	// Get pointer to the system of equations (SOE)
	LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();
	if (theSOE == 0) {
	  opserr << "ERROR the FE linearSOE must be defined before ";
	  opserr << "the sensitivity algorithm\n";
	  return -1;
	}

	// Get pointer to incremental integrator
	IncrementalIntegrator *theIncInt = theAlgorithm->getIncrementalIntegratorPtr();
//	IncrementalIntegrator *theIncIntSens=theAlgorithm->getIncrementalIntegratorPtr();//Abbas
	if (theIncInt == 0 ) {
	  opserr << "ERROR the FE integrator must be defined before ";
	  opserr << "the sensitivity algorithm\n";
	  return -1;
	}

	// Form current tangent at converged state
	// (would be nice with an if-statement here in case
	// the current tangent is already formed)
	if (theIncInt->formIndependentSensitivityLHS(CURRENT_TANGENT) < 0){
		opserr << "WARNING SensitivityAlgorithm::computeGradients() -";
		opserr << "the Integrator failed in formTangent()\n";
		return -1;
	}
	// Zero out the old right-hand side of the SOE
	theSOE->zeroB();
		
	Integrator* theSensitivityIntegrator = theIncInt;

	// Form the part of the RHS which are independent of parameter
	theSensitivityIntegrator->formIndependentSensitivityRHS();

	ParameterIter &paramIter = theDomain->getParameters();
//	opserr<<" get parameters "<<theDomain->getParameters()<<endln;//Abbas.......
	Parameter *theParam;
	// De-activate all parameters
	while ((theParam = paramIter()) != 0)
	  theParam->activate(false);

	// Now, compute sensitivity wrt each parameter
	int numGrads = theDomain->getNumParameters();
	//opserr<<"the numGrads is "<<numGrads<<endln;//Abbas...............................
	paramIter = theDomain->getParameters();
	
	while ((theParam = paramIter()) != 0) {

	  // Activate this parameter
	  theParam->activate(true);

	  // Zero the RHS vector
	  theSOE->zeroB();

	  // Get the grad index for this parameter
	  int gradIndex = theParam->getGradIndex();

	  // Form the RHS
	  theSensitivityIntegrator->formSensitivityRHS(gradIndex);
         
	  // Solve for displacement sensitivity
	 
	  theSOE->solve();

	  // Save sensitivity to nodes
	  theSensitivityIntegrator->saveSensitivity( theSOE->getX(), gradIndex, numGrads );
	 


	  // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
	  theSensitivityIntegrator->commitSensitivity(gradIndex, numGrads);
	  
	  // De-activate this parameter for next sensitivity calc
	  theParam->activate(false);
	}

	return 0;
}
*/

/*
bool 
SensitivityAlgorithm::shouldComputeAtEachStep(void)
{
  return (analysisTypeTag == 1);
}
*/
