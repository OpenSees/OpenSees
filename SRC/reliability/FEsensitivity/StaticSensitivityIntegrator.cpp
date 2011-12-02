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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-03-04 00:46:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/StaticSensitivityIntegrator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SensitivityIntegrator.h>
#include <StaticSensitivityIntegrator.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <classTags.h>


StaticSensitivityIntegrator::StaticSensitivityIntegrator(AnalysisModel *theModel, LinearSOE *theLinSOE)
:SensitivityIntegrator(), StaticIntegrator(INTEGRATOR_TAGS_StaticSensitivity),
theAnalysisModel(theModel),theSOE(theLinSOE)
{
}


StaticSensitivityIntegrator::~StaticSensitivityIntegrator()
{
}





int
StaticSensitivityIntegrator::formEleResidual(FE_Element *theEle)
{
	theEle->zeroResidual();
	theEle->addResistingForceSensitivity(gradNumber); 

	return 0;
}




int 
StaticSensitivityIntegrator::formIndependentSensitivityRHS()
{
	// For now everything is done each time in the static case
	return 0;
}


int 
StaticSensitivityIntegrator::formSensitivityRHS(int passedGradNumber)
{
	// Set a couple of data members
	gradNumber = passedGradNumber;


	// Loop through elements
	FE_Element *elePtr;
	FE_EleIter &theEles = theAnalysisModel->getFEs();    
	while((elePtr = theEles()) != 0) {
		theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
	}

	// Loop through the loadPatterns and add the dPext/dh contributions
	Vector oneDimVectorWithOne(1);
	oneDimVectorWithOne(0) = 1.0;
	ID oneDimID(1);
	Node *aNode;
	DOF_Group *aDofGroup;
	int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;
	LoadPattern *loadPatternPtr;
	Domain *theDomain = theAnalysisModel->getDomainPtr();
	LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    while((loadPatternPtr = thePatterns()) != 0) {
		const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
		sizeRandomLoads = randomLoads.Size();
		if (sizeRandomLoads == 1) {
			// No random loads in this load pattern
		}
		else {
			// Random loads: add contributions to the 'B' vector
			numRandomLoads = (int)(sizeRandomLoads/2);
			for (i=0; i<numRandomLoads*2; i=i+2) {
				nodeNumber = (int)randomLoads(i);
				dofNumber = (int)randomLoads(i+1);
				aNode = theDomain->getNode(nodeNumber);
				aDofGroup = aNode->getDOF_GroupPtr();
				const ID &anID = aDofGroup->getID();
				relevantID = anID(dofNumber-1);
				oneDimID(0) = relevantID;
				theSOE->addB(oneDimVectorWithOne, oneDimID);
			}
		}
	}

	return 0;
}
		



int
StaticSensitivityIntegrator::saveSensitivity(const Vector &v, int gradNum, int numGrads)
{

    DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
    DOF_Group 	*dofPtr;

	Vector *vNewPtr = new Vector(v.Size());
	(*vNewPtr) = v;
    while ( (dofPtr = theDOFGrps() ) != 0)  {
		dofPtr->saveSensitivity(vNewPtr,0,0,gradNum,numGrads);
	}

	delete vNewPtr;
    
    return 0;
}



int 
StaticSensitivityIntegrator::commitSensitivity(int gradNum, int numGrads)
{

	// Loop through the FE_Elements and set unconditional sensitivities
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
		elePtr->commitSensitivity(gradNum, numGrads);
	}

	return 0;
}





int 
StaticSensitivityIntegrator::newStep(void)
{
	return 0;
}
int 
StaticSensitivityIntegrator::update(const Vector &deltaU)
{
	return 0;
}
int 
StaticSensitivityIntegrator::setDeltaLambda(double newDeltaLambda)
{
	return 0;
}
int 
StaticSensitivityIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}
int 
StaticSensitivityIntegrator::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return 0;
}
void 
StaticSensitivityIntegrator::Print(OPS_Stream &s, int flag)  
{
}





