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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-27 23:05:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/NewmarkSensitivityIntegrator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SensitivityIntegrator.h>
#include <Newmark.h>
#include <NewmarkSensitivityIntegrator.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <Information.h>
#include <DOF_GrpIter.h>
#include <FE_EleIter.h>
#include <LoadPattern.h>
#include <Domain.h>
#include <LoadPatternIter.h>
#include <Node.h>
#include <NodeIter.h>


NewmarkSensitivityIntegrator::NewmarkSensitivityIntegrator()
:Newmark(),SensitivityIntegrator(),parameterID(0),sensitivityFlag(0),gradNumber(0)
{
	massMatrixMultiplicator = 0;
	dampingMatrixMultiplicator = 0;
	assemblyFlag = 0;
}

NewmarkSensitivityIntegrator::NewmarkSensitivityIntegrator(int passedAssemblyFlag, double theGamma, double theBeta, bool dispFlag)
:Newmark(theGamma, theBeta, dispFlag),SensitivityIntegrator(),
parameterID(0),sensitivityFlag(0),gradNumber(0)
{
	massMatrixMultiplicator = 0;
	dampingMatrixMultiplicator = 0;
	assemblyFlag = passedAssemblyFlag;
}

NewmarkSensitivityIntegrator::NewmarkSensitivityIntegrator(int passedAssemblyFlag, double theGamma, double theBeta, 
		 double alpham, double betak, 
		 double betaki, double betakc,
		 bool dispFlag)
:Newmark(theGamma, theBeta, alpham, betak, betaki, betakc, dispFlag),
SensitivityIntegrator(),
parameterID(0),sensitivityFlag(0),gradNumber(0)
{
	massMatrixMultiplicator = 0;
	dampingMatrixMultiplicator = 0;
	assemblyFlag = passedAssemblyFlag;
}

NewmarkSensitivityIntegrator::~NewmarkSensitivityIntegrator()
{
	if (massMatrixMultiplicator!=0)
		delete massMatrixMultiplicator;

	if (dampingMatrixMultiplicator!=0)
		delete dampingMatrixMultiplicator;
}

int
NewmarkSensitivityIntegrator::formEleResidual(FE_Element *theEle)
{

	if (sensitivityFlag == 0) {  // NO SENSITIVITY ANALYSIS

		this->Newmark::formEleResidual(theEle);

	}
	else {  // (ASSEMBLE ALL TERMS)

		theEle->zeroResidual();

		// Compute the time-stepping parameters on the form
		// udotdot = a1*ui+1 + a2*ui + a3*udoti + a4*udotdoti
		// udot    = a5*ui+1 + a6*ui + a7*udoti + a8*udotdoti
		// (see p. 166 of Chopra)

		// The constants are:
		// a1 = 1.0/(beta*dt*dt)
		// a2 = -1.0/(beta*dt*dt)
		// a3 = -1.0/beta*dt
		// a4 = 1.0 - 1.0/(2.0*beta)
		// a5 = gamma/(beta*dt)
		// a6 = -gamma/(beta*dt)
		// a7 = 1.0 - gamma/beta
		// a8 = 1.0 - gamma/(2.0*beta)

		// We can make use of the data members c2 and c3 of this class. 
		// As long as disp==true, they are defined as:
		// c2 = gamma/(beta*dt)
		// c3 = 1.0/(beta*dt*dt)

		// So, the constants can be computed as follows:
		if (displ==false) {
			opserr << "ERROR: Newmark::formEleResidual() -- the implemented"
				<< " scheme only works if the displ variable is set to true." << endln;
		}
		double a2 = -c3;
		double a3 = -c2/gamma;
		double a4 = 1.0 - 1.0/(2.0*beta);
		double a6 = -c2;
		double a7 = 1.0 - gamma/beta;
		double dt = gamma/(beta*c2);
		double a8 = dt*(1.0 - gamma/(2.0*beta));


		// Obtain sensitivity vectors from previous step
		int vectorSize = U->Size();
		Vector V(vectorSize);
		Vector Vdot(vectorSize);
		Vector Vdotdot(vectorSize);
		int i, loc;

		AnalysisModel *myModel = this->getAnalysisModelPtr();
		DOF_GrpIter &theDOFs = myModel->getDOFs();
		DOF_Group *dofPtr;
		while ((dofPtr = theDOFs()) != 0) {

			const ID &id = dofPtr->getID();
			int idSize = id.Size();
			const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);	
			for (i=0; i < idSize; i++) {
				loc = id(i);
				if (loc >= 0) {
					V(loc) = dispSens(i);		
				}
			}

			const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
			for (i=0; i < idSize; i++) {
				loc = id(i);
				if (loc >= 0) {
					Vdot(loc) = velSens(i);
				}
			}

			const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);	
			for (i=0; i < idSize; i++) {
				loc = id(i);
				if (loc >= 0) {
					Vdotdot(loc) = accelSens(i);
				}
			}
		}


		// Pre-compute the vectors involving a2, a3, etc.
		Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
		Vector tmp2 = V*a6 + Vdot*a7 + Vdotdot*a8;

		if (massMatrixMultiplicator == 0)
			massMatrixMultiplicator = new Vector(tmp1.Size());
		if (dampingMatrixMultiplicator == 0)
			dampingMatrixMultiplicator = new Vector(tmp2.Size());

		(*massMatrixMultiplicator) = tmp1;
		(*dampingMatrixMultiplicator) = tmp2;


		// Now we're ready to make calls to the FE Element:

		// The term -dPint/dh|u fixed
		theEle->addResistingForceSensitivity(gradNumber); 

		// The term -dM/dh*acc
		theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

		// The term -M*(a2*v + a3*vdot + a4*vdotdot)
		theEle->addM_Force(*massMatrixMultiplicator,-1.0);

		// The term -C*(a6*v + a7*vdot + a8*vdotdot)
		theEle->addD_Force(*dampingMatrixMultiplicator,-1.0);

		// The term -dC/dh*vel
		theEle->addD_ForceSensitivity(gradNumber, *Udot,-1.0);
		
	}

	return 0;
}



int
NewmarkSensitivityIntegrator::formNodUnbalance(DOF_Group *theDof)
{

	if (sensitivityFlag == 0) {  // NO SENSITIVITY ANALYSIS

		this->Newmark::formNodUnbalance(theDof);

	}
	else {  // ASSEMBLE ALL TERMS

		theDof->zeroUnbalance();


		// The term -M*(a2*v + a3*vdot + a4*vdotdot)
		theDof->addM_Force(*massMatrixMultiplicator,-1.0);


		// The term -dM/dh*acc
		theDof->addM_ForceSensitivity(*Udotdot, -1.0);


		// The term -C*(a6*v + a7*vdot + a8*vdotdot)
		theDof->addD_Force(*dampingMatrixMultiplicator,-1.0);


		// The term -dC/dh*vel
		theDof->addD_ForceSensitivity(*Udot,-1.0);


		// In case of random loads (have already been formed by 'applyLoadSensitivity')
		theDof->addPtoUnbalance();

	}


	return 0;
}


int 
NewmarkSensitivityIntegrator::formSensitivityRHS(int passedGradNumber)
{
	sensitivityFlag = 1;


	// Set a couple of data members
	gradNumber = passedGradNumber;

	// Get pointer to the SOE
	LinearSOE *theSOE = this->getLinearSOEPtr();


	// Possibly set the independent part of the RHS
	if (assemblyFlag != 0) {
		theSOE->setB(independentRHS);
	}

	// Get the analysis model
	AnalysisModel *theModel = this->getAnalysisModelPtr();



	// Randomness in external load (including randomness in time series)
	// Get domain
	Domain *theDomain = theModel->getDomainPtr();

	// Loop through nodes to zero the unbalaced load
	Node *nodePtr;
	NodeIter &theNodeIter = theDomain->getNodes();
	while ((nodePtr = theNodeIter()) != 0)
	nodePtr->zeroUnbalancedLoad();


	// Loop through load patterns to add external load sensitivity
	LoadPattern *loadPatternPtr;
	LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
	double time;
	while((loadPatternPtr = thePatterns()) != 0) {
		time = theDomain->getCurrentTime();
		loadPatternPtr->applyLoadSensitivity(time);
	}


	// Randomness in element/material contributions
	// Loop through FE elements
	FE_Element *elePtr;
	FE_EleIter &theEles = theModel->getFEs();    
	while((elePtr = theEles()) != 0) {
		theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
	}


	// Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
	DOF_Group *dofPtr;
	DOF_GrpIter &theDOFs = theModel->getDOFs();
	while((dofPtr = theDOFs()) != 0) {
		theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
	}


	// Reset the sensitivity flag
	sensitivityFlag = 0;

	return 0;
}
		





int 
NewmarkSensitivityIntegrator::formIndependentSensitivityRHS()
{
	// For now; don't use this
/*
	sensitivityFlag = 2; // Tell subsequent methods what to be assembled

	// Get pointer to the SOE
	LinearSOE *theSOE = this->getLinearSOEPtr();


	// Get the analysis model
	AnalysisModel *theModel = this->getAnalysisModelPtr();

	
	// Loop through FE elements
	FE_Element *elePtr;
	FE_EleIter &theEles = theModel->getFEs();    
	while((elePtr = theEles()) != 0) {
		theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
	}


	// Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
	DOF_Group *dofPtr;
	DOF_GrpIter &theDOFs = theModel->getDOFs();
	while((dofPtr = theDOFs()) != 0) {
		theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
	}


	// Set the data member of this class
	independentRHS = theSOE->getB();


	// Reset the sensitivity flag
	sensitivityFlag = 0;
*/

	return 0;
}
		




int 
NewmarkSensitivityIntegrator::saveSensitivity(const Vector & vNew,int gradNum,int numGrads)
{

	// Compute Newmark parameters in general notation
	double a1 = c3;
	double a2 = -c3;
	double a3 = -c2/gamma;
	double a4 = 1.0 - 1.0/(2.0*beta);
	double a5 = c2;
	double a6 = -c2;
	double a7 = 1.0 - gamma/beta;
	double dt = gamma/(beta*c2);
	double a8 = dt*(1.0 - gamma/(2.0*beta));


	// Recover sensitivity results from previous step
	int vectorSize = U->Size();
	Vector V(vectorSize);
	Vector Vdot(vectorSize);
	Vector Vdotdot(vectorSize);
	int i, loc;

	AnalysisModel *myModel = this->getAnalysisModelPtr();
	DOF_GrpIter &theDOFs = myModel->getDOFs();
	DOF_Group *dofPtr;
	while ((dofPtr = theDOFs()) != 0) {

		const ID &id = dofPtr->getID();
		int idSize = id.Size();
		const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);	
		for (i=0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				V(loc) = dispSens(i);		
			}
		}

		const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
		for (i=0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				Vdot(loc) = velSens(i);
			}
		}

		const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);	
		for (i=0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				Vdotdot(loc) = accelSens(i);
			}
		}
	}


	// Compute new acceleration and velocity vectors:
	Vector *vNewPtr = new Vector(vectorSize);
	Vector *vdotNewPtr = new Vector(vectorSize);
	Vector *vdotdotNewPtr = new Vector(vectorSize);
	(*vdotdotNewPtr) = vNew*a1 + V*a2 + Vdot*a3 + Vdotdot*a4;
	(*vdotNewPtr) = vNew*a5 + V*a6 + Vdot*a7 + Vdotdot*a8;
	(*vNewPtr) = vNew;


	// Now we can save vNew, vdotNew and vdotdotNew
    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group 	*dofPtr1;
    while ( (dofPtr1 = theDOFGrps() ) != 0)  {
		dofPtr1->saveSensitivity(vNewPtr,vdotNewPtr,vdotdotNewPtr,gradNum,numGrads);
	}

	delete vNewPtr;
	delete vdotNewPtr;
	delete vdotdotNewPtr;

	return 0;
}



int 
NewmarkSensitivityIntegrator::commitSensitivity(int gradNum, int numGrads)
{

	// Loop through the FE_Elements and set unconditional sensitivities
	AnalysisModel *theAnalysisModel = this->getAnalysisModelPtr();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
		elePtr->commitSensitivity(gradNum, numGrads);
	}

	return 0;
}




int
NewmarkSensitivityIntegrator::setParameter(char **argv, int argc, Information &info)
{
    if (strcmp(argv[0],"alphaM") == 0) {
        info.theType = DoubleType;
        return 1;
    }
    if (strcmp(argv[0],"betaK") == 0) {
        info.theType = DoubleType;
        return 2;
    }
    // otherwise parameter is unknown
    else {
		opserr << "ERROR: Unknown random parameter in Newmark::setParameter()" << endln;
      return -1;
	}
}

int
NewmarkSensitivityIntegrator::updateParameter   (int parameterID, Information &info)
{
  switch (parameterID) {
      
    case 1:
        this->alphaM = info.theDouble;
        return 0;

    case 2:
        this->betaK = info.theDouble;
        return 0;

    default:
		return 0;
  }
}

int
NewmarkSensitivityIntegrator::activateParameter (int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}









