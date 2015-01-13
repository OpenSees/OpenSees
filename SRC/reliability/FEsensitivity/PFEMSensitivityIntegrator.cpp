
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
                               

// $Revision: 1.0 $
// $Date: 2014-01-16 10:03:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/PFEMSensitivityIntegrator.h,v $


//
// Written by Minjie Zhu (Oregon State University)
//                                         


#include <SensitivityIntegrator.h>
#include <PFEMSensitivityIntegrator.h>
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


PFEMSensitivityIntegrator::PFEMSensitivityIntegrator()
    :SensitivityIntegrator(),PFEMIntegrator(),
     parameterID(0),sensitivityFlag(0),
     gradNumber(0),assemblyFlag(0),independentRHS(),dVn()
{
}

PFEMSensitivityIntegrator::PFEMSensitivityIntegrator(int passedAssemblyFlag)
    :SensitivityIntegrator(),PFEMIntegrator(),
     parameterID(0),sensitivityFlag(0),gradNumber(0),
     assemblyFlag(passedAssemblyFlag),independentRHS(),dVn()
{
}


PFEMSensitivityIntegrator::~PFEMSensitivityIntegrator()
{
}

int
PFEMSensitivityIntegrator::formEleResidual(FE_Element *theEle)
{

    if (sensitivityFlag == 0) {  // NO SENSITIVITY ANALYSIS

        this->PFEMIntegrator::formEleResidual(theEle);

    }
    else {  // (ASSEMBLE ALL TERMS)

        theEle->zeroResidual();

        // Compute the time-stepping parameters on the form
        // udotdot = 1/dt*vn+1 - 1/dt*vn
        // u       = un + dt*vn+1


        // Obtain sensitivity vectors from previous step
        dVn.resize(U->Size()); dVn.Zero();
        Vector dUn(U->Size());

        AnalysisModel *myModel = this->getAnalysisModel();
        DOF_GrpIter &theDOFs = myModel->getDOFs();
        DOF_Group *dofPtr = 0;
        while ((dofPtr = theDOFs()) != 0) {

            const ID &id = dofPtr->getID();
            int idSize = id.Size();

            const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
            for (int i=0; i < idSize; i++) {
                int loc = id(i);
                if (loc >= 0) {
                    dUn(loc) = dispSens(i);
                }
            }

            const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
            for (int i=0; i < idSize; i++) {
                int loc = id(i);
                if (loc >= 0) {
                    dVn(loc) = velSens(i);
                }
            }
        }

        // Now we're ready to make calls to the FE Element:

        // The term -dPint/dh|u fixed
        theEle->addResistingForceSensitivity(gradNumber); 

        // The term -dM/dh*acc
        theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

        // The term -M*(-1/dt*dvn)
        theEle->addM_Force(dVn, c3);

        // The term -K*(dun)
        theEle->addK_Force(dUn, -1.0);

        // The term -dC/dh*vel
        theEle->addD_ForceSensitivity(gradNumber, *Udot,-1.0);
		
    }

    return 0;
}



int
PFEMSensitivityIntegrator::formNodUnbalance(DOF_Group *theDof)
{

    if (sensitivityFlag == 0) {  // NO SENSITIVITY ANALYSIS

        this->PFEMIntegrator::formNodUnbalance(theDof);

    }
    else {  // ASSEMBLE ALL TERMS

        theDof->zeroUnbalance();

        // The term -M*(-1/dt*dvn)
        theDof->addM_Force(dVn, c3);

        // The term -dM/dh*acc
        theDof->addM_ForceSensitivity(*Udotdot, -1.0);

        // The term -C*(dvn)
        //theDof->addD_Force(dVn ,-1.0);


        // The term -dC/dh*vel
        theDof->addD_ForceSensitivity(*Udot,-1.0);


        // In case of random loads (have already been formed by 'applyLoadSensitivity')
        theDof->addPtoUnbalance();

    }


    return 0;
}


int 
PFEMSensitivityIntegrator::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;


    // Set a couple of data members
    gradNumber = passedGradNumber;

    // Get pointer to the SOE
    LinearSOE *theSOE = this->getLinearSOE();


    // Possibly set the independent part of the RHS
    if (assemblyFlag != 0) {
        theSOE->setB(independentRHS);
    }

    // Get the analysis model
    AnalysisModel *theModel = this->getAnalysisModel();



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
PFEMSensitivityIntegrator::formIndependentSensitivityRHS()
{
    return 0;
}
		




int 
PFEMSensitivityIntegrator::saveSensitivity(const Vector & dVNew,int gradNum,int numGrads)
{
    // Recover sensitivity results from previous step
    int vectorSize = U->Size();
    Vector dUn(vectorSize);
    dVn.resize(vectorSize); dVn.Zero();

    AnalysisModel *myModel = this->getAnalysisModel();
    DOF_GrpIter &theDOFs = myModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
	  
        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);	
        for (int i=0; i < idSize; i++) {
	    int loc = id(i);
	    if (loc >= 0) {
                dUn(loc) = dispSens(i);		
	    }
        }
	  
        const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
        for (int i=0; i < idSize; i++) {
	    int loc = id(i);
	    if (loc >= 0) {
                dVn(loc) = velSens(i);
	    }
        }

    }


    // Compute new acceleration and velocity vectors:
    Vector dUNew(vectorSize);
    Vector dANew(vectorSize);

    // dudotdot = 1/dt*dv{n+1} - 1/dt*dvn
    dANew.addVector(0.0, dVNew, c3);
    dANew.addVector(1.0, dVn, -c3);

    // du       = dun + dt*dv{n+1}
    dUNew.addVector(0.0, dVNew, c1);
    dUNew.addVector(1.0, dUn, 1.0);

    // Now we can save vNew, vdotNew and vdotdotNew
    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group 	*dofPtr1;
    while ( (dofPtr1 = theDOFGrps() ) != 0)  {
        dofPtr1->saveSensitivity(dUNew,dVNew,dANew,gradNum,numGrads);
    }
	
    return 0;
}



int 
PFEMSensitivityIntegrator::commitSensitivity(int gradNum, int numGrads)
{

    // Loop through the FE_Elements and set unconditional sensitivities
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
        elePtr->commitSensitivity(gradNum, numGrads);
    }

    return 0;
}




int
PFEMSensitivityIntegrator::setParameter(char **argv, int argc, Information &info)
{
    opserr<<"ERROR: Unknown random parameter in ";
    opserr<<"PFEMSensitivityIntegrator::setParameter()\n";
    return -1;
}

int
PFEMSensitivityIntegrator::updateParameter(int parameterID, Information &info)
{
    return 0;
}

int
PFEMSensitivityIntegrator::activateParameter (int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}


int 
PFEMSensitivityIntegrator::updateGradNumber(int passedGradNumber)
{
    opserr << "Fatal PFEMSensitivityIntegrator::updateGradNumber";
    opserr << "This function should not be called from this object\n";
    opserr << "This function should be called from NewPFEMSensitivityIntegrator\n";
    return -1;
}
int 
PFEMSensitivityIntegrator::sensitivityDomainChanged(int NumGrads)
{
    opserr << "Fatal PFEMSensitivityIntegrator::updateGradNumber";
    opserr << "This function should not be called from this object\n";
    opserr << "This function should be called from NewPFEMSensitivityIntegrator\n";
    return -1;
}
bool 
PFEMSensitivityIntegrator::staticSensitivity(void)
{
    return false;
}
bool 
PFEMSensitivityIntegrator::NewSensitivity(void)
{
    return false;
}







