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

// $Revision: 1.2 $
// $Date: 2009-05-11 22:42:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/ReliabilityDirectIntegrationAnalysis.cpp,v $
                                                                        
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <FE_EleIter.h>

#include <ReliabilityDirectIntegrationAnalysis.h>
#include <EquiSolnAlgo.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
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

// AddingSensitivity:BEGIN ////////////////////////////////////
//#ifdef _RELIABILITY
#include <SensitivityAlgorithm.h>
#include<Integrator.h>
//#endif
// AddingSensitivity:END ////////////////////////////////////

// Constructor
//    sets theModel and theSysOFEqn to 0 and the Algorithm to the one supplied

ReliabilityDirectIntegrationAnalysis::ReliabilityDirectIntegrationAnalysis(Domain &the_Domain,
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
 theIntegrator(&theTransientIntegrator), 
 theTest(theConvergenceTest),
 domainStamp(0)
{
  // first we set up the links needed by the elements in the 
  // aggregation
  theAnalysisModel->setLinks(the_Domain, theHandler);
  theConstraintHandler->setLinks(the_Domain,theModel,theTransientIntegrator);
  theDOF_Numberer->setLinks(theModel);
  theIntegrator->setLinks(theModel,theLinSOE, theTest);
  theAlgorithm->setLinks(theModel,theTransientIntegrator,theLinSOE, theTest);

  if (theTest != 0)
    theAlgorithm->setConvergenceTest(theTest);
  else
    theTest = theAlgorithm->getConvergenceTest();
  
// AddingSensitivity:BEGIN ////////////////////////////////////
	theSensitivityAlgorithm = 0;
	theGFunEachStepEvaluator = 0;
// AddingSensitivity:END //////////////////////////////////////
}    

ReliabilityDirectIntegrationAnalysis::~ReliabilityDirectIntegrationAnalysis()
{
  // we don't invoke the destructors in case user switching
  // from a static to a direct integration analysis 
  // clearAll() must be invoked if user wishes to invoke destructor
}    

void
ReliabilityDirectIntegrationAnalysis::clearAll(void)
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
  if (theTest != 0)
    delete theTest;

// AddingSensitivity:BEGIN ////////////////////////////////////
//#ifdef _RELIABILITY
	delete theSensitivityAlgorithm;
//#endif
// AddingSensitivity:END //////////////////////////////////////

    theAnalysisModel =0;
    theConstraintHandler =0;
    theDOF_Numberer =0;
    theIntegrator =0;
    theAlgorithm =0;
    theSOE =0;
    theTest =0;
}    

#include <NodeIter.h>
#include <Node.h>

int 
ReliabilityDirectIntegrationAnalysis::initialize(void)
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
ReliabilityDirectIntegrationAnalysis::analyze(int numSteps, double dT)
{

  int result = 0;
  Domain *the_Domain = this->getDomainPtr();

    for (int i=0; i<numSteps; i++) {

	  if( i%100 ==0) opserr << "step "<<i << endln;

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
//#ifdef _RELIABILITY
      if (theSensitivityAlgorithm != 0) {
	result = theSensitivityAlgorithm->computeSensitivities();
	if (result < 0) {
	  opserr << "StaticAnalysis::analyze() - the SensitivityAlgorithm failed";
	  opserr << " at iteration: " << i << " with domain at load factor ";
	  opserr << the_Domain->getCurrentTime() << endln;
	  the_Domain->revertToLastCommit();	    
	  theIntegrator->revertToLastStep();
	  return -5;
	}    
      }
//#endif
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
      
      // opserr << "DirectIntegrationAnalysis - time: " << the_Domain->getCurrentTime() << endln;
    }    
    return result;
}


int
ReliabilityDirectIntegrationAnalysis::domainChanged(void)
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

    // opserr << theAnalysisModel->getDOFGroupGraph();

    /*
    DOF_GrpIter &theDOFs = theAnalysisModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) 
      opserr << dofPtr->getID();
    */

    theDOF_Numberer->numberDOF();


    theConstraintHandler->doneNumberingDOF();

    /*
    DOF_GrpIter &theDOFs1 = theAnalysisModel->getDOFs();
    while ((dofPtr = theDOFs1()) != 0) 
      opserr << dofPtr->getID();


    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    FE_Element *elePtr;
    while((elePtr = theEles()) != 0)     
       opserr << elePtr->getID();
    */

    // we invoke setGraph() on the LinearSOE which
    // causes that object to determine its size

    theSOE->setSize(theAnalysisModel->getDOFGraph());

    // we invoke domainChange() on the integrator and algorithm
    theIntegrator->domainChanged();
    theAlgorithm->domainChanged();


	// AddingSensitivity:BEGIN ////////////////////////////////////
	if (theSensitivityAlgorithm != 0){ 
		int result = theSensitivityAlgorithm->sensitivityDomainChanged();
	}
	// AddingSensitivity:END //////////////////////////////////////

    return 0;
}    

// AddingSensitivity:BEGIN //////////////////////////////
//#ifdef _RELIABILITY
int 
ReliabilityDirectIntegrationAnalysis::setSensitivityAlgorithm(/*SensitivityAlgorithm*/ Integrator *passedSensitivityAlgorithm)
{
    int result = 0;

    // invoke the destructor on the old one
//    if (theSensitivityAlgorithm != 0) {
//      delete theSensitivityAlgorithm;
//    }

    theSensitivityAlgorithm = passedSensitivityAlgorithm;
    
    return 0;
}
//#endif
// AddingSensitivity:END ///////////////////////////////


int 
ReliabilityDirectIntegrationAnalysis::setNumberer(DOF_Numberer &theNewNumberer) 
{
    int result = 0;

    // invoke the destructor on the old one
    if (theDOF_Numberer != 0)
	delete theDOF_Numberer;

    // first set the links needed by the Algorithm
    theDOF_Numberer = &theNewNumberer;
    theDOF_Numberer->setLinks(*theAnalysisModel);

    // invoke domainChanged() either indirectly or directly
    Domain *the_Domain = this->getDomainPtr();
    int stamp = the_Domain->hasDomainChanged();
    domainStamp = stamp;
    result = this->domainChanged();    
    if (result < 0) {
      opserr << "StaticAnalysis::setNumberer() - setNumberer() failed";
      return -1;
    }	

    return 0;
}



int 
ReliabilityDirectIntegrationAnalysis::setAlgorithm(EquiSolnAlgo &theNewAlgorithm) 
{
    // invoke the destructor on the old one
    if (theAlgorithm != 0)
	delete theAlgorithm;

    // first set the links needed by the Algorithm
    theAlgorithm = &theNewAlgorithm;
    theAlgorithm->setLinks(*theAnalysisModel,*theIntegrator,*theSOE, theTest);

    // invoke domainChanged() either indirectly or directly
    Domain *the_Domain = this->getDomainPtr();
    // check if domain has undergone change
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;	    
	if (this->domainChanged() < 0) {
	    opserr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	    opserr << "domainChanged failed";
	    return -1;
	}	
    } else {
	if (theAlgorithm->domainChanged() < 0) {
	    opserr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	    opserr << "algorithm::domainChanged() failed";
	    return -2;
	}
    }

    return 0;
}


int 
ReliabilityDirectIntegrationAnalysis::setIntegrator(TransientIntegrator &theNewIntegrator)
{
  // set the links needed by the other objects in the aggregation
  Domain *the_Domain = this->getDomainPtr();
  theIntegrator = &theNewIntegrator;
  theIntegrator->setLinks(*theAnalysisModel,*theSOE, theTest);
  theConstraintHandler->setLinks(*the_Domain,*theAnalysisModel,*theIntegrator);
  theAlgorithm->setLinks(*theAnalysisModel,*theIntegrator,*theSOE, theTest);

  // invoke domainChanged() either indirectly or directly
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;	    
	if (this->domainChanged() < 0) {
	  opserr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  opserr << "domainChanged failed";
	  return -1;
      }	
  }
  else {
      if (theIntegrator->domainChanged() < 0) {
	  opserr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  opserr << "Integrator::domainChanged failed";
	  return -1;
      }	
  }
   
  return 0;
}
int 
ReliabilityDirectIntegrationAnalysis::setLinearSOE(LinearSOE &theNewSOE)
{
  // invoke the destructor on the old one
  if (theSOE != 0)
    delete theSOE;

  // set the links needed by the other objects in the aggregation
  theSOE = &theNewSOE;
  theIntegrator->setLinks(*theAnalysisModel,*theSOE, theTest);
  theAlgorithm->setLinks(*theAnalysisModel,*theIntegrator,*theSOE, theTest);

  // set the size either indirectly or directly
  Domain *the_Domain = this->getDomainPtr();
  int stamp = the_Domain->hasDomainChanged();
  if (stamp != domainStamp) {
      domainStamp = stamp;	    
      if (this->domainChanged() < 0) {
	  opserr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  opserr << "domainChanged failed";
	  return -1;
      }	
  } else {
      Graph &theGraph = theAnalysisModel->getDOFGraph();
      if (theSOE->setSize(theGraph) < 0) {
	  opserr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  opserr << "LinearSOE::setSize() failed";
	  return -2;	
      }
  }
  
  return 0;
}


int 
ReliabilityDirectIntegrationAnalysis::setConvergenceTest(ConvergenceTest &theNewTest)
{
  // invoke the destructor on the old one
  if (theTest != 0)
    delete theTest;
  
  // set the links needed by the other objects in the aggregation
  theTest = &theNewTest;
  theAlgorithm->setConvergenceTest(theTest);
  
  return 0;
}


int
ReliabilityDirectIntegrationAnalysis::checkDomainChange(void)
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
ReliabilityDirectIntegrationAnalysis::getAlgorithm(void)
{
  return theAlgorithm;
}


TransientIntegrator *
ReliabilityDirectIntegrationAnalysis::getIntegrator(void)
{
  return theIntegrator;
}

ConvergenceTest *
ReliabilityDirectIntegrationAnalysis::getConvergenceTest(void)
{
  return theTest;
}



Matrix* ReliabilityDirectIntegrationAnalysis::getEachStepResult()
{
	if(theGFunEachStepEvaluator!=0){
		return theGFunEachStepEvaluator->getLSFValues();
	}else{
		opserr<<"ERROR in DirectIntegrationAnalysis::getEachStepResult\n";
		opserr<<"theGFunEachStepEvaluator is NULL\n";
		opserr<<"when calling getEachStepResult()\n";
		exit(-1);
	}
}
Matrix* ReliabilityDirectIntegrationAnalysis::getEachStepConvFlag()
{
	if(theGFunEachStepEvaluator!=0){
		return theGFunEachStepEvaluator->getConvFlag();
	}else{
		opserr<<"ERROR in DirectIntegrationAnalysis::getEachStepResult\n";
		opserr<<"theGFunEachStepEvaluator is NULL\n";
		opserr<<"when calling getEachStepConvFlag()\n";
		exit(-1);
	}
}





