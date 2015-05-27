/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2007-04-02 23:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/TransientIntegrator.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/TransientIntegrator.C
// 
// Written: fmk 
// Created: Tue Sept 17 15:54:47: 1996
// Revision: A
//
// Description: This file contains the class definition for TransientIntegrator.
// TransientIntegrator is an algorithmic class for setting up the finite element
// equations for a static analysis and for Incrementing the nodal displacements
// with the values in the soln vector to the LinearSOE object. 
//
// What: "@(#) TransientIntegrator.C, revA"

#include <TransientIntegrator.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>

TransientIntegrator::TransientIntegrator(int clasTag)
:IncrementalIntegrator(clasTag)
{

}

TransientIntegrator::~TransientIntegrator()
{

}

int 
TransientIntegrator::formTangent(int statFlag)
{
    int result = 0;
    statusFlag = statFlag;

    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0) {
	opserr << "WARNING TransientIntegrator::formTangent() ";
	opserr << "no LinearSOE or AnalysisModel has been set\n";
	return -1;
    }
    
    // the loops to form and add the tangents are broken into two for 
    // efficiency when performing parallel computations
    
    theLinSOE->zeroA();

    // do modal damping
    bool inclModalMatrix=theModel->inclModalDampingMatrix();
    if (inclModalMatrix == true) {
      const Vector *modalValues = theModel->getModalDampingFactors();
      if (modalValues != 0) {
	this->addModalDampingMatrix(modalValues);
      }
    }


    // loop through the DOF_Groups and add the unbalance
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    
    while ((dofPtr = theDOFs()) != 0) {
	if (theLinSOE->addA(dofPtr->getTangent(this),dofPtr->getID()) <0) {
	    opserr << "TransientIntegrator::formTangent() - failed to addA:dof\n";
	    result = -1;
	}
    }    

    // loop through the FE_Elements getting them to add the tangent    
    FE_EleIter &theEles2 = theModel->getFEs();    
    FE_Element *elePtr;    
    while((elePtr = theEles2()) != 0)     {
	if (theLinSOE->addA(elePtr->getTangent(this),elePtr->getID()) < 0) {
	    opserr << "TransientIntegrator::formTangent() - failed to addA:ele\n";
	    result = -2;
	}
    }
    return result;
}


    
int
TransientIntegrator::formUnbalance(void) {
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();

    if (theModel == 0 || theLinSOE == 0) {
 	opserr << "WARNING IncrementalIntegrator::formUnbalance -";
	opserr << " no AnalysisModel or LinearSOE has been set\n";
	return -1;
    }
    
    theLinSOE->zeroB();

    // do modal damping
    const Vector *modalValues = theModel->getModalDampingFactors();
    if (modalValues != 0) {
      this->addModalDampingForce(modalValues);
    }
    
    if (this->formElementResidual() < 0) {
	opserr << "WARNING IncrementalIntegrator::formUnbalance ";
	opserr << " - this->formElementResidual failed\n";
	return -1;
    }
    
    if (this->formNodalUnbalance() < 0) {
	opserr << "WARNING IncrementalIntegrator::formUnbalance ";
	opserr << " - this->formNodalUnbalance failed\n";
	return -2;
    }    

    return 0;
}
    
int
TransientIntegrator::formEleResidual(FE_Element *theEle)
{
  theEle->zeroResidual();
  theEle->addRIncInertiaToResidual();
  return 0;
}    

int
TransientIntegrator::formNodUnbalance(DOF_Group *theDof)
{
  theDof->zeroUnbalance();
  theDof->addPIncInertiaToUnbalance();
  return 0;
}    



