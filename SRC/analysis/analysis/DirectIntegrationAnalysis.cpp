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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DirectIntegrationAnalysis.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/DirectIntegrationAnalysis.C
// 
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
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <TransientIntegrator.h>
#include <Domain.h>

#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>

// Constructor
//    sets theModel and theSysOFEqn to 0 and the Algorithm to the one supplied

DirectIntegrationAnalysis::DirectIntegrationAnalysis(
			      Domain &the_Domain,
			      ConstraintHandler &theHandler,
			      DOF_Numberer &theNumberer,
			      AnalysisModel &theModel,
			      EquiSolnAlgo &theSolnAlgo,		   
			      LinearSOE &theLinSOE,
			      TransientIntegrator &theTransientIntegrator)
:TransientAnalysis(the_Domain), 
 theConstraintHandler(&theHandler),
 theDOF_Numberer(&theNumberer), 
 theAnalysisModel(&theModel), 
 theAlgorithm(&theSolnAlgo), 
 theSOE(&theLinSOE),
 theIntegrator(&theTransientIntegrator), domainStamp(0)
{
    // first we set up the links needed by the elements in the 
    // aggregation
    theAnalysisModel->setLinks(the_Domain);
    theConstraintHandler->setLinks(the_Domain,theModel,theTransientIntegrator);
    theDOF_Numberer->setLinks(theModel);
    theIntegrator->setLinks(theModel,theLinSOE);
    theAlgorithm->setLinks(theModel,theTransientIntegrator,theLinSOE);
}    

DirectIntegrationAnalysis::~DirectIntegrationAnalysis()
{

}    

void
DirectIntegrationAnalysis::clearAll(void)
{
    // invoke the destructor on all the objects in the aggregation
    delete theAnalysisModel;
    delete theConstraintHandler;
    delete theDOF_Numberer;
    delete theIntegrator;
    delete theAlgorithm;
    delete theSOE;
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
	cerr << "DirectIntegrationAnalysis::initialize() - domainChanged() failed\n";
	return -1;
      }	
    }
    if (theIntegrator->initialize() < 0) {
	cerr << "DirectIntegrationAnalysis::initialize() - integrator initialize() failed\n";
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

    for (int i=0; i<numSteps; i++) {
	// check if domain has undergone change
	int stamp = the_Domain->hasDomainChanged();
	if (stamp != domainStamp) {
	    domainStamp = stamp;	
	    if (this->domainChanged() < 0) {
		cerr << "DirectIntegrationAnalysis::analyze() - domainChanged() failed\n";
		return -1;
	    }	
	}

	if (theIntegrator->newStep(dT) < 0) {
	    cerr << "DirectIntegrationAnalysis::analyze() - the Integrator failed";
	    cerr << " at time " << the_Domain->getCurrentTime() << endl;
	    the_Domain->revertToLastCommit();
	    return -2;
	}
	
	result = theAlgorithm->solveCurrentStep();
	if (result < 0) {
	    cerr << "DirectIntegrationAnalysis::analyze() - the Algorithm failed";
	    cerr << " at time " << the_Domain->getCurrentTime() << endl;
	    the_Domain->revertToLastCommit();	    
	    return -3;
	}    
	
	result = theIntegrator->commit();
	if (result < 0) {
	    cerr << "DirectIntegrationAnalysis::analyze() - ";
	    cerr << "the Integrator failed to commit";
	    cerr << " at time " << the_Domain->getCurrentTime() << endl;
	    the_Domain->revertToLastCommit();	    
	    return -4;
	} 
   
	// cerr << "DirectIntegrationAnalysis - time: " << the_Domain->getCurrentTime() << endl;
	}    
    return result;
}


int
DirectIntegrationAnalysis::domainChanged(void)
{

    theAnalysisModel->clearAll();    
    theConstraintHandler->clearAll();
    
    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.

    theConstraintHandler->handle();
    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.

    // cerr << theAnalysisModel->getDOFGroupGraph();

    /*
    DOF_GrpIter &theDOFs = theAnalysisModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) 
      cerr << dofPtr->getID();
    */

    theDOF_Numberer->numberDOF();

    // cerr << theAnalysisModel->getDOFGraph();


    /*
    DOF_GrpIter &theDOFs1 = theAnalysisModel->getDOFs();
    while ((dofPtr = theDOFs1()) != 0) 
      cerr << dofPtr->getID();


    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    FE_Element *elePtr;
    while((elePtr = theEles()) != 0)     
       cerr << elePtr->getID();
    */

    // we invoke setGraph() on the LinearSOE which
    // causes that object to determine its size

    theSOE->setSize(theAnalysisModel->getDOFGraph());



    // we invoke domainChange() on the integrator and algorithm
    theIntegrator->domainChanged();
    theAlgorithm->domainChanged();

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
    theAlgorithm->setLinks(*theAnalysisModel,*theIntegrator,*theSOE);

    // invoke domainChanged() either indirectly or directly
    Domain *the_Domain = this->getDomainPtr();
    // check if domain has undergone change
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;	    
	if (this->domainChanged() < 0) {
	    cerr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	    cerr << "domainChanged failed";
	    return -1;
	}	
    } else {
	if (theAlgorithm->domainChanged() < 0) {
	    cerr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	    cerr << "algorithm::domainChanged() failed";
	    return -2;
	}
    }

    return 0;
}


int 
DirectIntegrationAnalysis::setIntegrator(TransientIntegrator &theNewIntegrator)
{
  // set the links needed by the other objects in the aggregation
  Domain *the_Domain = this->getDomainPtr();
  theIntegrator = &theNewIntegrator;
  theIntegrator->setLinks(*theAnalysisModel,*theSOE);
  theConstraintHandler->setLinks(*the_Domain,*theAnalysisModel,*theIntegrator);
  theAlgorithm->setLinks(*theAnalysisModel,*theIntegrator,*theSOE);

  // invoke domainChanged() either indirectly or directly
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;	    
	if (this->domainChanged() < 0) {
	  cerr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  cerr << "domainChanged failed";
	  return -1;
      }	
  }
  else {
      if (theIntegrator->domainChanged() < 0) {
	  cerr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  cerr << "Integrator::domainChanged failed";
	  return -1;
      }	
  }
   
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
  theIntegrator->setLinks(*theAnalysisModel,*theSOE);
  theAlgorithm->setLinks(*theAnalysisModel,*theIntegrator,*theSOE);

  // set the size either indirectly or directly
  Domain *the_Domain = this->getDomainPtr();
  int stamp = the_Domain->hasDomainChanged();
  if (stamp != domainStamp) {
      domainStamp = stamp;	    
      if (this->domainChanged() < 0) {
	  cerr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  cerr << "domainChanged failed";
	  return -1;
      }	
  } else {
      Graph &theGraph = theAnalysisModel->getDOFGraph();
      if (theSOE->setSize(theGraph) < 0) {
	  cerr << "DirectIntegrationAnalysis::setAlgorithm() - ";
	  cerr << "LinearSOE::setSize() failed";
	  return -2;	
      }
  }
  
  return 0;
}









