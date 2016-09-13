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
                                                                        
// $Revision: 1.15 $
// $Date: 2009-05-14 22:50:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DirectIntegrationAnalysis.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of the
// DirectIntegrationAnalysis class.
//
// What: "@(#) DirectIntegrationAnalysis.C, revA"


#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <FE_EleIter.h>

#include <DirectIntegrationAnalysis.h>
#include <EquiSolnAlgo.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <TransientIntegrator.h>
#include <Domain.h>

#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>
// AddingSensitivity:BEGIN //////////////////////////////////
#ifdef _RELIABILITY
#include <SensitivityAlgorithm.h>
#include<Integrator.h>//Abbas
#endif
// AddingSensitivity:END ////////////////////////////////////

// Constructor
//    sets theModel and theSysOFEqn to 0 and the Algorithm to the one supplied

DirectIntegrationAnalysis::DirectIntegrationAnalysis(Domain &the_Domain,
						     ConstraintHandler &theHandler,
						     DOF_Numberer &theNumberer,
						     AnalysisModel &theModel,
						     EquiSolnAlgo &theSolnAlgo,		   
						     LinearSOE &theLinSOE,
						     TransientIntegrator &theTransientIntegrator,
						     ConvergenceTest *theConvergenceTest)
:TransientAnalysis(the_Domain), 
 theConstraintHandler(&theHandler),
 theDOF_Numberer(&theNumberer), 
 theAnalysisModel(&theModel), 
 theAlgorithm(&theSolnAlgo), 
 theSOE(&theLinSOE), 
 theEigenSOE(0),
 theIntegrator(&theTransientIntegrator), 
 theTest(theConvergenceTest),
 domainStamp(0)
{
  // first we set up the links needed by the elements in the 
  // aggregation
  theAnalysisModel->setLinks(the_Domain, theHandler);
  theConstraintHandler->setLinks(the_Domain, theModel, theTransientIntegrator);
  theDOF_Numberer->setLinks(theModel);
  theIntegrator->setLinks(theModel, theLinSOE, theTest);
  theAlgorithm->setLinks(theModel, theTransientIntegrator, theLinSOE, theTest);
  theSOE->setLinks(theModel);

  if (theTest != 0)
    theAlgorithm->setConvergenceTest(theTest);
  else
    theTest = theAlgorithm->getConvergenceTest();
  
// AddingSensitivity:BEGIN ////////////////////////////////////
#ifdef _RELIABILITY
//	theSensitivityAlgorithm = 0;
#endif
// AddingSensitivity:END //////////////////////////////////////
}    

DirectIntegrationAnalysis::~DirectIntegrationAnalysis()
{
  // we don't invoke the destructors in case user switching
  // from a static to a direct integration analysis 
  // clearAll() must be invoked if user wishes to invoke destructor
}    

void
DirectIntegrationAnalysis::clearAll(void)
{
  // invoke the destructor on all the objects in the aggregation
  if (theAnalysisModel != 0)     
    delete theAnalysisModel;
  if (theConstraintHandler != 0) 
    delete theConstraintHandler;
  if (theDOF_Numberer != 0)      
    delete theDOF_Numberer;
  if (theIntegrator != 0) 
    delete theIntegrator;
  if (theAlgorithm != 0)  
    delete theAlgorithm;
  if (theSOE != 0)
    delete theSOE;
  if (theEigenSOE != 0)
    delete theEigenSOE;
  if (theTest != 0)
    delete theTest;

// AddingSensitivity:BEGIN ////////////////////////////////////
#ifdef _RELIABILITY
//	delete theSensitivityAlgorithm;
#endif
// AddingSensitivity:END //////////////////////////////////////

    theAnalysisModel =0;
    theConstraintHandler =0;
    theDOF_Numberer =0;
    theIntegrator =0;
    theAlgorithm =0;
    theSOE =0;
    theEigenSOE =0;
    theTest =0;
}    

#include <NodeIter.h>
#include <Node.h>

int 
DirectIntegrationAnalysis::initialize(void)
{
    Domain *the_Domain = this->getDomainPtr();

    // check if domain has undergone change
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
      domainStamp = stamp;	
      if (this->domainChanged() < 0) {
	opserr << "DirectIntegrationAnalysis::initialize() - domainChanged() failed\n";
	return -1;
      }	
    }
    if (theIntegrator->initialize() < 0) {
	opserr << "DirectIntegrationAnalysis::initialize() - integrator initialize() failed\n";
	return -2;
    } else
      theIntegrator->commit();

    return 0;
}

int 
DirectIntegrationAnalysis::analyze(int numSteps, double dT)
{
  int result = 0;
  Domain *the_Domain = this->getDomainPtr();
 // if (theEigenSOE != 0)
 //   theIntegrator->setEigenSOE(theEigenSOE);

  for (int i=0; i<numSteps; i++) {

    if (theAnalysisModel->analysisStep(dT) < 0) {
      opserr << "DirectIntegrationAnalysis::analyze() - the AnalysisModel failed";
      opserr << " at time " << the_Domain->getCurrentTime() << endln;
      the_Domain->revertToLastCommit();
      return -2;
    }
    
    // check if domain has undergone change
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
      domainStamp = stamp;	
      if (this->domainChanged() < 0) {
	opserr << "DirectIntegrationAnalysis::analyze() - domainChanged() failed\n";
	return -1;
      }	
    }

    if (theIntegrator->newStep(dT) < 0) {
      opserr << "DirectIntegrationAnalysis::analyze() - the Integrator failed";
      opserr << " at time " << the_Domain->getCurrentTime() << endln;
      the_Domain->revertToLastCommit();
      theIntegrator->revertToLastStep();
      return -2;
    }
    
    result = theAlgorithm->solveCurrentStep();
    if (result < 0) {
      opserr << "DirectIntegrationAnalysis::analyze() - the Algorithm failed";
      opserr << " at time " << the_Domain->getCurrentTime() << endln;
      the_Domain->revertToLastCommit();	    
      theIntegrator->revertToLastStep();
      return -3;
    }    
    
// AddingSensitivity:BEGIN ////////////////////////////////////
#ifdef _RELIABILITY
   /*
    if (theSensitivityAlgorithm != 0) {
      result = theIntegrator->computeSensitivities();
      if (result < 0) {
	opserr << "StaticAnalysis::analyze() - the SensitivityAlgorithm failed";
	opserr << " at iteration: " << i << " with domain at load factor ";
	opserr << the_Domain->getCurrentTime() << endln;
	the_Domain->revertToLastCommit();	    
	theIntegrator->revertToLastStep();
	return -5;
      }    
    }
    */
#endif
    // AddingSensitivity:END //////////////////////////////////////
    
    result = theIntegrator->commit();
    if (result < 0) {
      opserr << "DirectIntegrationAnalysis::analyze() - ";
      opserr << "the Integrator failed to commit";
      opserr << " at time " << the_Domain->getCurrentTime() << endln;
      the_Domain->revertToLastCommit();	    
      theIntegrator->revertToLastStep();
      return -4;
    } 
  }    
  return result;
}

int 
DirectIntegrationAnalysis::eigen(int numMode, bool generalized, bool findSmallest)
{
    if (theAnalysisModel == 0 || theEigenSOE == 0) {
      opserr << "WARNING DirectIntegrationAnalysis::eigen() - no EigenSOE has been set\n";
      return -1;
    }

    int result = 0;
    Domain *the_Domain = this->getDomainPtr();

    result = theAnalysisModel->eigenAnalysis(numMode, generalized, findSmallest);

    int stamp = the_Domain->hasDomainChanged();

    if (stamp != domainStamp) {
      domainStamp = stamp;
  
      result = this->domainChanged();
      
      if (result < 0) {
	     opserr << "DirectIntegrationAnalysis::eigen() - domainChanged failed";
	return -1;
      }	
    }


    //
    // zero A and M
    //
    theEigenSOE->zeroA();
    theEigenSOE->zeroM();

    //
    // form K
    //


    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    FE_Element *elePtr;

    while((elePtr = theEles()) != 0) {
      elePtr->zeroTangent();
      elePtr->addKtToTang(1.0);
      if (theEigenSOE->addA(elePtr->getTangent(0), elePtr->getID()) < 0) {
	opserr << "WARNING DirectIntegrationAnalysis::eigen() -";
	opserr << " failed in addA for ID " << elePtr->getID();	    
	result = -2;
      }
    }

    //
    // if generalized is true, form M
    //

    if (generalized == true) {
      int result = 0;
      FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
      while((elePtr = theEles2()) != 0) {     
	elePtr->zeroTangent();
	elePtr->addMtoTang(1.0);
	if (theEigenSOE->addM(elePtr->getTangent(0), elePtr->getID()) < 0) {
	  opserr << "WARNING DirectIntegrationAnalysis::eigen() -";
	  opserr << " failed in addA for ID " << elePtr->getID();	    
	  result = -2;
	}
      }
      
      DOF_Group *dofPtr;
      DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();    
      while((dofPtr = theDofs()) != 0) {
	dofPtr->zeroTangent();
	dofPtr->addMtoTang(1.0);
	if (theEigenSOE->addM(dofPtr->getTangent(0),dofPtr->getID()) < 0) {
	  opserr << "WARNING DirectIntegrationAnalysis::eigen() -";
	  opserr << " failed in addM for ID " << dofPtr->getID();	    
	  result = -3;
	}
      }
    }
    
    // 
    // solve for the eigen values & vectors
    //

    if (theEigenSOE->solve(numMode, generalized, findSmallest) < 0) {
	opserr << "WARNING DirectIntegrationAnalysis::eigen() - EigenSOE failed in solve()\n";
	return -4;
    }
	
    //
    // now set the eigenvalues and eigenvectors in the model
    //

    theAnalysisModel->setNumEigenvectors(numMode);
    Vector theEigenvalues(numMode);
    for (int i = 1; i <= numMode; i++) {
      theEigenvalues[i-1] = theEigenSOE->getEigenvalue(i);
      theAnalysisModel->setEigenvector(i, theEigenSOE->getEigenvector(i));
    }    
    theAnalysisModel->setEigenvalues(theEigenvalues);
  
    return 0;
}



int
DirectIntegrationAnalysis::domainChanged(void)
{
    Domain *the_Domain = this->getDomainPtr();
    int stamp = the_Domain->hasDomainChanged();
    domainStamp = stamp;

    theAnalysisModel->clearAll();    
    theConstraintHandler->clearAll();
    
    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.
    theConstraintHandler->handle();

    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.
    theDOF_Numberer->numberDOF();

    theConstraintHandler->doneNumberingDOF();

    // we invoke setGraph() on the LinearSOE which
    // causes that object to determine its size
    Graph &theGraph = theAnalysisModel->getDOFGraph();

    int result = theSOE->setSize(theGraph);
    if (result < 0) {
	opserr << "DirectIntegrationAnalysis::handle() - ";
	opserr << "LinearSOE::setSize() failed";
	return -3;
    }	    

    if (theEigenSOE != 0) {
      result = theEigenSOE->setSize(theGraph);
      if (result < 0) {
	opserr << "DirectIntegrationAnalysis::handle() - ";
	opserr << "EigenSOE::setSize() failed";
	return -3;
      }	    
    }

    theAnalysisModel->clearDOFGraph();

    // we invoke domainChange() on the integrator and algorithm
    theIntegrator->domainChanged();
    theAlgorithm->domainChanged();

    return 0;
}    

// AddingSensitivity:BEGIN //////////////////////////////
#ifdef _RELIABILITY
int 
DirectIntegrationAnalysis::setSensitivityAlgorithm(/*SensitivityAlgorithm*/ Integrator *passedSensitivityAlgorithm)
{
    int result = 0;

    // invoke the destructor on the old one
  //  if (theSensitivityAlgorithm != 0) {
   //   delete theSensitivityAlgorithm;
  //  }

 //   theSensitivityAlgorithm = passedSensitivityAlgorithm;
    
    return 0;
}
#endif
// AddingSensitivity:END ///////////////////////////////


int 
DirectIntegrationAnalysis::setNumberer(DOF_Numberer &theNewNumberer) 
{
    // invoke the destructor on the old one
    if (theDOF_Numberer != 0)
	delete theDOF_Numberer;

    // first set the links needed by the Algorithm
    theDOF_Numberer = &theNewNumberer;
    theDOF_Numberer->setLinks(*theAnalysisModel);

    // invoke domainChanged() either indirectly or directly
    domainStamp = 0;
    return 0;
}



int 
DirectIntegrationAnalysis::setAlgorithm(EquiSolnAlgo &theNewAlgorithm) 
{
  // invoke the destructor on the old one
  if (theAlgorithm != 0)
    delete theAlgorithm;
  
  // first set the links needed by the Algorithm
  theAlgorithm = &theNewAlgorithm;

  if (theAnalysisModel != 0 && theIntegrator != 0 && theSOE != 0)
    theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
  // invoke domainChanged() either indirectly or directly
  // domainStamp = 0;
  if (domainStamp != 0)
    theAlgorithm->domainChanged();

  return 0;
}


int 
DirectIntegrationAnalysis::setIntegrator(TransientIntegrator &theNewIntegrator)
{
  // invoke the destructor on the old one
  if (theIntegrator != 0) {
      delete theIntegrator;
  }
  // set the links needed by the other objects in the aggregation
  Domain *the_Domain = this->getDomainPtr();
  theIntegrator = &theNewIntegrator;
  theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);
  theConstraintHandler->setLinks(*the_Domain, *theAnalysisModel, *theIntegrator);
  theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);

  // cause domainChanged to be invoked on next analyze
  //  domainStamp = 0;
  if (domainStamp != 0)
    theIntegrator->domainChanged();
   
  return 0;
}
int 
DirectIntegrationAnalysis::setLinearSOE(LinearSOE &theNewSOE)
{
  // invoke the destructor on the old one
  if (theSOE != 0)
    delete theSOE;

  // set the links needed by the other objects in the aggregation
  theSOE = &theNewSOE;
  theIntegrator->setLinks(*theAnalysisModel,*theSOE, theTest);
  theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
  theSOE->setLinks(*theAnalysisModel);
  
  if (theEigenSOE != 0) 
    theEigenSOE->setLinearSOE(*theSOE);
  
  // cause domainChanged to be invoked on next analyze
  domainStamp = 0;
  
  return 0;
}

int 
DirectIntegrationAnalysis::setEigenSOE(EigenSOE &theNewSOE)
{
  // invoke the destructor on the old one if not the same!
  if (theEigenSOE != 0) {
    if (theEigenSOE->getClassTag() != theNewSOE.getClassTag()) {
      delete theEigenSOE;
      theEigenSOE = 0;
    }
  }

  if (theEigenSOE == 0) {
    theEigenSOE = &theNewSOE;
    theEigenSOE->setLinks(*theAnalysisModel);
    theEigenSOE->setLinearSOE(*theSOE);

    /*
    if (domainStamp != 0) {
      Graph &theGraph = theAnalysisModel->getDOFGraph();
      theEigenSOE->setSize(theGraph);
    }
    */
    domainStamp = 0;
  }
 
  return 0;
}

int 
DirectIntegrationAnalysis::setConvergenceTest(ConvergenceTest &theNewTest)
{
  // invoke the destructor on the old one
  if (theTest != 0)
    delete theTest;
  
  // set the links needed by the other objects in the aggregation
  theTest = &theNewTest;

  if (theIntegrator != 0)
    theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);

  if (theAlgorithm != 0)
    theAlgorithm->setConvergenceTest(theTest);
  
  return 0;
}


int
DirectIntegrationAnalysis::checkDomainChange(void)
{
  Domain *the_Domain = this->getDomainPtr();

  // check if domain has undergone change
  int stamp = the_Domain->hasDomainChanged();
  if (stamp != domainStamp) {
    domainStamp = stamp;	
    if (this->domainChanged() < 0) {
      opserr << "DirectIntegrationAnalysis::initialize() - domainChanged() failed\n";
      return -1;
    }	
  }

  return 0;
}


EquiSolnAlgo *
DirectIntegrationAnalysis::getAlgorithm(void)
{
  return theAlgorithm;
}

AnalysisModel *
DirectIntegrationAnalysis::getModel(void)
{
  return theAnalysisModel;
}


TransientIntegrator *
DirectIntegrationAnalysis::getIntegrator(void)
{
  return theIntegrator;
}

ConvergenceTest *
DirectIntegrationAnalysis::getConvergenceTest(void)
{
  return theTest;
}




