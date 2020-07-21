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
                                                                        
// $Revision: 1.10 $
// $Date: 2009-05-14 22:50:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/TransientDomainDecompositionAnalysis.cpp,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation of TransientDomainDecompositionAnalysis.
//
// What: "@(#) TransientDomainDecompositionAnalysis.C, revA"

#include <TransientDomainDecompositionAnalysis.h>
#include <EquiSolnAlgo.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <LinearSOESolver.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <TransientIntegrator.h>
#include <ConvergenceTest.h>
#include <Subdomain.h>
//#include <Timer.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// AddingSensitivity:BEGIN //////////////////////////////////
#ifdef _RELIABILITY
#include <SensitivityAlgorithm.h>
#endif
// AddingSensitivity:END ////////////////////////////////////


#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>

#include <Vector.h>
#include <Matrix.h>

TransientDomainDecompositionAnalysis::TransientDomainDecompositionAnalysis(Subdomain &the_Domain)
  :DomainDecompositionAnalysis(ANALYSIS_TAGS_TransientDomainDecompositionAnalysis, the_Domain), 
   theConstraintHandler(0),
   theDOF_Numberer(0), 
   theAnalysisModel(0), 
   theAlgorithm(0), 
   theSOE(0),
   theEigenSOE(0),
   theIntegrator(0), 
   theTest(0),
   domainStamp(0)
{

}    

TransientDomainDecompositionAnalysis::TransientDomainDecompositionAnalysis(Subdomain &the_Domain,
									   ConstraintHandler &theHandler,
									   DOF_Numberer &theNumberer,
									   AnalysisModel &theModel,
									   EquiSolnAlgo &theSolnAlgo,		   
									   LinearSOE &theLinSOE,
									   TransientIntegrator &theTransientIntegrator,
									   ConvergenceTest *theConvergenceTest,
									   bool setLinks)
  :DomainDecompositionAnalysis(ANALYSIS_TAGS_TransientDomainDecompositionAnalysis, the_Domain), 
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
  
  if (setLinks == true) {
    // set up the links needed by the elements in the aggregation
    theAnalysisModel->setLinks(the_Domain, theHandler);
    theConstraintHandler->setLinks(the_Domain,theModel,theTransientIntegrator);
    theDOF_Numberer->setLinks(theModel);
    theIntegrator->setLinks(theModel,theLinSOE, theTest);
    theAlgorithm->setLinks(theModel,theTransientIntegrator,theLinSOE, theTest);
    theSOE->setLinks(*theAnalysisModel);
  }

}    


TransientDomainDecompositionAnalysis::~TransientDomainDecompositionAnalysis()
{
  // we don't invoke the destructors in case user switching
  // from a static to a direct integration analysis 
  // clearAll() must be invoked if user wishes to invoke destructor
}    

void
TransientDomainDecompositionAnalysis::clearAll(void)
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
  
  // now set the pointers to NULL
  theAnalysisModel =0;
  theConstraintHandler =0;
  theDOF_Numberer =0;
  theIntegrator =0;
  theAlgorithm =0;
  theSOE =0;
  theEigenSOE =0;
  theTest = 0;
}    


bool 
TransientDomainDecompositionAnalysis::doesIndependentAnalysis(void)
{
  return true;
}


int 
TransientDomainDecompositionAnalysis::analyze(double dT)
{
  int result = 0;
  Domain *the_Domain = this->getDomainPtr();

  // check for change in Domain since last step. As a change can
  // occur in a commit() in a domaindecomp with load balancing
  // this must now be inside the loop

  //opserr << "TransientDomainDecomp - hasDomainChanged\n";
  int stamp = the_Domain->hasDomainChanged();

  if (stamp != domainStamp) {
    domainStamp = stamp;
    result = this->domainChanged();
    if (result < 0) {
      opserr << "TransientDomainDecompositionAnalysis::analyze() - domainChanged failed";
      return -1;
    }	
  }
  // opserr << "TransientDomainDecomp - newStep\n";
  // result = theAnalysisModel->newStepDomain();
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::analyze() - the AnalysisModel failed";
    opserr << " with domain at load factor ";
    opserr << the_Domain->getCurrentTime() << endln;
    the_Domain->revertToLastCommit();
    
    return -2;
  }

  // opserr << "TransientDomainDecomp - integrator newStep\n";
  result = theIntegrator->newStep(dT);

  // barrierCheck in PartitionedDomain
 // result = this->checkAllResult(result);

  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::analyze() - the Integrator failed";
    opserr << " with domain at load factor ";
    opserr << the_Domain->getCurrentTime() << endln;
    the_Domain->revertToLastCommit();
    theIntegrator->revertToLastStep();    
    return -2;
  }
 //  opserr << "TransientDomainDecomp - solveCurrentStep\n";
  result = theAlgorithm->solveCurrentStep();
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::analyze() - the Algorithm failed";
    opserr << " with domain at load factor ";
    opserr << the_Domain->getCurrentTime() << endln;
    the_Domain->revertToLastCommit();	    
    theIntegrator->revertToLastStep();
    
    return -3;
  }    

  result = theIntegrator->commit();
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::analyze() - ";
    opserr << "the Integrator failed to commit";
    opserr << " with domain at load factor ";
    opserr << the_Domain->getCurrentTime() << endln;
    the_Domain->revertToLastCommit();	    
    theIntegrator->revertToLastStep();
    
    return -4;
  }    	

  return 0;
}


int 
TransientDomainDecompositionAnalysis::eigen(int numMode, bool generalized, bool findSmallest)
{
  int result = 0;
  Domain *the_Domain = this->getDomainPtr();


  if (theEigenSOE == 0) {
    opserr << "TransientDomainDecompositionAnalysis::eigen() - no eigen solver has been set\n";
    return -1;
  }

  // check for change in Domain since last step. As a change can


  // occur in a commit() in a domaindecomp with load balancing
  // this must now be inside the loop
  int stamp = the_Domain->hasDomainChanged();

  if (stamp != domainStamp) {
    domainStamp = stamp;
    result = this->domainChanged();
    if (result < 0) {
      opserr << "TransientDomainDecompositionAnalysis::eigen() - domainChanged failed";
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
	opserr << "WARNING TransientDomainDecomposition::eigen() -";
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
	  opserr << "WARNING TransientDomainDecomposition::eigen() -";
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
	  opserr << "WARNING TransientDomainDecomposition::eigen() -";
	  opserr << " failed in addM for ID " << dofPtr->getID();	    
	  result = -3;
	}
      }
    }
    
    // 
    // solve for the eigen values & vectors
    //
    if (theEigenSOE->solve(numMode, generalized, findSmallest) < 0) {
	opserr << "WARNING TransientDomainDecomposition::eigen() - EigenSOE failed in solve()\n";
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
TransientDomainDecompositionAnalysis::initialize(void)
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
TransientDomainDecompositionAnalysis::domainChanged(void)
{
  int result = 0;
  theAnalysisModel->clearAll();    
  theConstraintHandler->clearAll();

  // now we invoke handle() on the constraint handler which
  // causes the creation of FE_Element and DOF_Group objects
  // and their addition to the AnalysisModel.
  result = theConstraintHandler->handle();
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::handle() - ";
    opserr << "ConstraintHandler::handle() failed";
    return -1;
  }	

  // we now invoke number() on the numberer which causes
  // equation numbers to be assigned to all the DOFs in the
  // AnalysisModel.
  result = theDOF_Numberer->numberDOF();
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::handle() - ";
    opserr << "DOF_Numberer::numberDOF() failed";
    return -2;
  }	    

  result = theConstraintHandler->doneNumberingDOF();
  
  // we invoke setSize() on the LinearSOE which
  // causes that object to determine its size
  Graph &theGraph = theAnalysisModel->getDOFGraph();
  result = theSOE->setSize(theGraph);
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::handle() - ";
    opserr << "LinearSOE::setSize() failed";
    return -3;
  }	    

  if (theEigenSOE != 0) {
    result = theEigenSOE->setSize(theGraph);
    if (result < 0) {
      opserr << "TransientDomainDecompositionAnalysis::handle() - ";
      opserr << "EigenSOE::setSize() failed";
      return -3;
    }	    
  }

  theAnalysisModel->clearDOFGraph();
  
  // finally we invoke domainChanged on the Integrator and Algorithm
  // objects .. informing them that the model has changed

  result = theIntegrator->domainChanged();
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::setAlgorithm() - ";
    opserr << "Integrator::domainChanged() failed";
    return -4;
  }	    
  
  result = theAlgorithm->domainChanged();
  if (result < 0) {
    opserr << "TransientDomainDecompositionAnalysis::setAlgorithm() - ";
    opserr << "Algorithm::domainChanged() failed";
    return -5;
  }	        

  // if get here successful
  return 0;
}    


int  
TransientDomainDecompositionAnalysis::getNumExternalEqn(void)
{
  opserr << "TransientDomainDecompositionAnalysis::getNumExternalEqn() - should never be called\n";
  return 0;
}
int  
TransientDomainDecompositionAnalysis::getNumInternalEqn(void)
{
  opserr << "TransientDomainDecompositionAnalysis::getNumInternalEqn() - should never be called\n";
  return 0;
}

int  
TransientDomainDecompositionAnalysis::analysisStep(double dT) 
{
  this->analyze(dT);
  return 0;
}

int  
TransientDomainDecompositionAnalysis::eigenAnalysis(int numMode, bool generalized, bool findSmallest) 
{
  this->eigen(numMode, generalized, findSmallest);
  return 0;
}


int  
TransientDomainDecompositionAnalysis::computeInternalResponse(void)
{
  opserr << "TransientDomainDecompositionAnalysis::computeInternalResponse() - should never be called\n";
  return 0;
}
int  
TransientDomainDecompositionAnalysis::formTangent(void)
{
  opserr << "TransientDomainDecompositionAnalysis::formTangent() - should never be called\n";
  return 0;
}
int  
TransientDomainDecompositionAnalysis::formResidual(void)
{
  opserr << "TransientDomainDecompositionAnalysis::formResidual() - should never be called\n";
  return 0;
}
int  
TransientDomainDecompositionAnalysis::formTangVectProduct(Vector &force)
{
  opserr << "TransientDomainDecompositionAnalysis::formTangVectProduct() - should never be called\n";
  return 0;
}

const Matrix &
TransientDomainDecompositionAnalysis::getTangent(void)
{
  static Matrix errMatrix;
  opserr << "TransientDomainDecompositionAnalysis::getTangent() - should never be called\n";
  return errMatrix;
}

const Vector &
TransientDomainDecompositionAnalysis::getResidual(void)
{
  static Vector errVector;
  opserr << "TransientDomainDecompositionAnalysis::getResidual() - should never be called\n";
  return errVector;
}

const Vector &
TransientDomainDecompositionAnalysis::getTangVectProduct(void)
{
  static Vector errVector;
  opserr << "TransientDomainDecompositionAnalysis::getTangVectProduct() - should never be called\n";
  return errVector;
}
    
int 
TransientDomainDecompositionAnalysis::sendSelf(int commitTag, Channel &theChannel)
{

  // receive the data identifyng the objects in the aggregation
  int dataTag = this->getDbTag();
  static ID data(8);

  if (theAlgorithm == 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - no objects exist!\n";
    return -1;
  }

  LinearSOESolver *theSolver = theSOE->getSolver();
  
  data(0) = theConstraintHandler->getClassTag();
  data(1) = theDOF_Numberer->getClassTag();
  data(2) = theAnalysisModel->getClassTag();
  data(3) = theAlgorithm->getClassTag();
  data(4) = theSOE->getClassTag();
  data(5) = theSolver->getClassTag();
  data(6) = theIntegrator->getClassTag();

  if (theTest != 0)
    data(7) = theTest->getClassTag();
  else
    data(7) = -1;

  theChannel.sendID(dataTag, commitTag, data);  

  // invoke sendSelf() on all the objects
  if (  theConstraintHandler->sendSelf(commitTag, theChannel) != 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send handler\n";
    return -1;
  }

  if (  theDOF_Numberer->sendSelf(commitTag, theChannel) != 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send numberer\n";
    return -1;
  }

  if (  theAnalysisModel->sendSelf(commitTag, theChannel) != 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send model\n";
    return -1;
  }

  if (  theAlgorithm->sendSelf(commitTag, theChannel) != 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send algorithm\n";
    return -1;
  }

  if (  theSOE->sendSelf(commitTag, theChannel) != 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send SOE\n";
    return -1;
  } else 
    ;

  //    theSOE->setAnalysisModel(*theAnalysisModel);

  if (  theSolver->sendSelf(commitTag, theChannel) != 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send Solver\n";
    return -1;
  }

  if (  theIntegrator->sendSelf(commitTag, theChannel) != 0) {
    opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send integrator\n";
    return -1;
  }

  if (theTest != 0)
    if (theTest->sendSelf(commitTag, theChannel) != 0) {
      opserr << "TransientDomainDecompositionAnalysis::sendSelf() - failed to send integrator\n";
      return -1;
  }

  return 0;
}
int 
TransientDomainDecompositionAnalysis::recvSelf(int commitTag, Channel &theChannel, 
					    FEM_ObjectBroker &theBroker)
{

  myChannel = &theChannel;

  Domain *the_Domain = this->getSubdomainPtr();

  // receive the data identifyng the objects in the aggregation
  static ID data(8);
  int dataTag = this->getDbTag();
  theChannel.recvID(dataTag, commitTag, data);

  // for all objects in the aggregation:
  //  1. make sure objects exist & are of correct type; create new objects if not
  //  2. invoke recvSelf on the object
  if (theConstraintHandler == 0 || theConstraintHandler->getClassTag() != data(0)) {
    if (theConstraintHandler != 0)
      delete theConstraintHandler;
    
    theConstraintHandler = theBroker.getNewConstraintHandler(data(0));
    if (theConstraintHandler == 0) {
      opserr << "TransientDomainDecompositionAnalysis::recvSelf";
      opserr << " - failed to get the ConstraintHandler\n";
      return -1;
    }
  }
  theConstraintHandler->recvSelf(commitTag, theChannel,theBroker);
  
  if (theDOF_Numberer == 0 || theDOF_Numberer->getClassTag() != data(1)) {
    if (theDOF_Numberer != 0)
      delete theDOF_Numberer;
    
    theDOF_Numberer = theBroker.getNewNumberer(data(1));
    if (theDOF_Numberer == 0) {
      opserr << "TransientDomainDecompositionAnalysis::recvSelf";
      opserr << " - failed to get the ConstraintHandler\n";
      return -1;
    }
  }
  theDOF_Numberer->recvSelf(commitTag, theChannel,theBroker);

  if (theAnalysisModel == 0 || theAnalysisModel->getClassTag() != data(2)) {
    if (theAnalysisModel != 0)
      delete theAnalysisModel;
    
    theAnalysisModel = theBroker.getNewAnalysisModel(data(2));
    if (theAnalysisModel == 0) {
      opserr << "TransientDomainDecompositionAnalysis::recvSelf";
      opserr << " - failed to get the Analysis Model\n";
      return -1;
    }
  }
  theAnalysisModel->recvSelf(commitTag, theChannel,theBroker);

  if (theAlgorithm == 0 || theAlgorithm->getClassTag() != data(3)) {
    if (theAlgorithm != 0)
      delete theAlgorithm;
    
    theAlgorithm = theBroker.getNewEquiSolnAlgo(data(3));
    if (theAlgorithm == 0) {
      opserr << "TransientDomainDecompositionAnalysis::recvSelf";
      opserr << " - failed to get the Solution Algorithm\n";
      return -1;
    }
  }
  theAlgorithm->recvSelf(commitTag, theChannel,theBroker);

  if (theSOE == 0 || theSOE->getClassTag() != data(4)) {
    if (theSOE != 0)
      delete theSOE;
    
    theSOE = theBroker.getNewLinearSOE(data(4));
    if (theSOE == 0) {
      opserr << "TransientDomainDecompositionAnalysis::recvSelf";
      opserr << " - failed to get the LinearSOE\n";
      return -1;
    }
  }

  theSOE->recvSelf(commitTag, theChannel, theBroker);
  LinearSOESolver *theSolver = theSOE->getSolver();
  if (theSolver == 0) {
      opserr << "StaticDomainDecompositionAnalysis::recvSelf";
      opserr << " - failed to get the Solver\n";
      return -1;
  }
  theSolver->recvSelf(commitTag, theChannel, theBroker);  
  //  theSOE->setAnalysisModel(*theAnalysisModel);

  if (theIntegrator == 0 || theIntegrator->getClassTag() != data(6)) {
    if (theIntegrator != 0)
      delete theIntegrator;
    
    theIntegrator = theBroker.getNewTransientIntegrator(data(6));
    if (theIntegrator == 0) {
      opserr << "TransientDomainDecompositionAnalysis::recvSelf";
      opserr << " - failed to get the Integrator\n";
      return -1;
    }
  }

  theIntegrator->recvSelf(commitTag, theChannel,theBroker);


  if (theTest == 0 || theTest->getClassTag() != data(7)) {
    if (theTest != 0)
      delete theIntegrator;
    
    if (data(7) != -1) {
      theTest = theBroker.getNewConvergenceTest(data(7));
      if (theTest == 0) {
	opserr << "TransientDomainDecompositionAnalysis::recvSelf";
	opserr << " - failed to get the ConvergenceTest\n";
	return -1;
      }
    }
  }
  if (theTest != 0)
    theTest->recvSelf(commitTag, theChannel,theBroker);

  // set up the links needed by the elements in the aggregation
  theAnalysisModel->setLinks(*the_Domain, *theConstraintHandler);
  theConstraintHandler->setLinks(*the_Domain, *theAnalysisModel, *theIntegrator);
  theDOF_Numberer->setLinks(*theAnalysisModel);
  theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);
  theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
  theAlgorithm->setConvergenceTest(theTest);

  return 0;
}

int 
TransientDomainDecompositionAnalysis::setAlgorithm(EquiSolnAlgo &theNewAlgorithm)
{
  // invoke the destructor on the old one
  if (theAlgorithm != 0)
    delete theAlgorithm;

  // first set the links needed by the Algorithm
  theAlgorithm = &theNewAlgorithm;

  if (theAnalysisModel != 0 && theIntegrator != 0 && theSOE != 0)
    theAlgorithm->setLinks(*theAnalysisModel,*theIntegrator,*theSOE, theTest);

  if (theTest != 0)
    theAlgorithm->setConvergenceTest(theTest);

  // invoke domainChanged() either indirectly or directly
  //  domainStamp = 0;
  if (domainStamp != 0)
    theAlgorithm->domainChanged();  
  
  return 0;
}

int 
TransientDomainDecompositionAnalysis::setIntegrator(IncrementalIntegrator &theNewIntegrator) 
{
  // invoke the destructor on the old one
  if (theIntegrator != 0) {
    delete theIntegrator;
  }
  
  // set the links needed by the other objects in the aggregation
  Domain *the_Domain = this->getDomainPtr();
  
  theIntegrator = (TransientIntegrator *)(&theNewIntegrator);

  if (theIntegrator != 0 && theConstraintHandler != 0 && theAlgorithm != 0 && theAnalysisModel != 0 && theSOE != 0) {
    theIntegrator->setLinks(*theAnalysisModel,*theSOE, theTest);
    theConstraintHandler->setLinks(*the_Domain, *theAnalysisModel, *theIntegrator);
    theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
  }

  // cause domainChanged to be invoked on next analyze
  domainStamp = 0;
  /*
  if (domainStamp != 0)
    theIntegrator->domainChanged();
  */
  return 0;
}


int 
TransientDomainDecompositionAnalysis::setLinearSOE(LinearSOE &theNewSOE)
{
    // invoke the destructor on the old one
    if (theSOE != 0)
	delete theSOE;

    // set the links needed by the other objects in the aggregation
    theSOE = &theNewSOE;
    if (theIntegrator != 0 && theAlgorithm != 0 && theAnalysisModel != 0 && theSOE != 0) {
      theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);
      theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
      theSOE->setLinks(*theAnalysisModel);
    }

    if (theEigenSOE != 0) 
      theEigenSOE->setLinearSOE(*theSOE);

    // cause domainChanged to be invoked on next analyze
    domainStamp = 0;

    return 0;
}

int 
TransientDomainDecompositionAnalysis::setEigenSOE(EigenSOE &theNewSOE)
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
TransientDomainDecompositionAnalysis::setConvergenceTest(ConvergenceTest &theConvergenceTest) 
{
  // invoke the destructor on the old one
  if (theTest != 0) {
    delete theTest;
  }
  theTest = &theConvergenceTest;

  if (theIntegrator != 0)
    theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);

  if (theAlgorithm != 0)
    return theAlgorithm->setConvergenceTest(theTest);

  return 0;
}

