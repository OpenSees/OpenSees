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
                                                                        
//----------------------------------------------------------------------------------------------------------------------------
 // Developed by:
 // Adam M. Knaack (adam.knaack@schaefer-inc.com) 
 // Schaefer-Inc, Cincinnati, Ohio, USA
 // Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
 // Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
 // Yahya C. Kurama (ykurama@nd.edu)
 // Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Created: 2012
 // Last updated: 2019
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Description: This file contains the source code of CreepAnalysis. 
 // CreepAnalysis is an analysis procedure that calculates time-dependent
 // concrete creep and shrinkage strains.
 //---------------------------------------------------------------------------------------------------------------------------- 
 // Detailed descriptions of the model and its implementation can be found in the following:
 // (1) Knaack, A.M., Kurama, Y.C. 2018. Modeling Time-Dependent Deformations: Application for Reinforced Concrete Beams with 
 //     Recycled Concrete Aggregates. ACI Structural J. 115, 175–190. doi:10.14359/51701153
 // (2) Knaack, A.M., 2013. Sustainable concrete structures using recycled concrete aggregate: short-term and long-term behavior
 //     considering material variability. PhD Dissertation, Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, Notre Dame, Indiana, USA, 680 pp.
 // A manual describing the use of the model and sample files can be found at:
 // ***Mendeley Data Link***(will be added later; ntosic)
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Disclaimer: This software is provided “as is”, without any warranties, expressed or implied. In no event shall the developers be liable for any claim, damages, or liability arising from or in connection with this software.
 //----------------------------------------------------------------------------------------------------------------------------

#include <CreepAnalysis.h>
#include <EquiSolnAlgo.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <StaticIntegrator.h>
#include <Domain.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>
#include <Timer.h>
#include <OPS_Globals.h> //Added by AMK

// AddingSensitivity:BEGIN //////////////////////////////////
#ifdef _RELIABILITY
#include <SensitivityAlgorithm.h>
#endif
// AddingSensitivity:END ////////////////////////////////////


// Constructor
//    sets theModel and theSysOFEqn to 0 and the Algorithm to the one supplied

CreepAnalysis::CreepAnalysis(Domain &the_Domain,
			       ConstraintHandler &theHandler,
			       DOF_Numberer &theNumberer,
			       AnalysisModel &theModel,
			       EquiSolnAlgo &theSolnAlgo,		   
			       LinearSOE &theLinSOE,
			       StaticIntegrator &theStaticIntegrator,
			       ConvergenceTest *theConvergenceTest)
:Analysis(the_Domain), theConstraintHandler(&theHandler),
 theDOF_Numberer(&theNumberer), theAnalysisModel(&theModel), 
 theAlgorithm(&theSolnAlgo), theSOE(&theLinSOE), theEigenSOE(0),
 theIntegrator(&theStaticIntegrator), theTest(theConvergenceTest),
 domainStamp(0)
{
    // first we set up the links needed by the elements in the 
    // aggregation
    theAnalysisModel->setLinks(the_Domain, theHandler);
    theConstraintHandler->setLinks(the_Domain, theModel, theStaticIntegrator);
    theDOF_Numberer->setLinks(theModel);

    theIntegrator->setLinks(theModel, theLinSOE, theTest);
    theAlgorithm->setLinks(theModel, theStaticIntegrator, theLinSOE, theTest);

    if (theTest != 0)
      theAlgorithm->setConvergenceTest(theTest);

    // AddingSensitivity:BEGIN ////////////////////////////////////
#ifdef _RELIABILITY
    theSensitivityAlgorithm = 0;
#endif
    // AddingSensitivity:END //////////////////////////////////////
}    


CreepAnalysis::~CreepAnalysis()
{
  // we don't invoke the destructors in case user switching
  // from a static to a direct integration analysis 
  // clearAll() must be invoked if user wishes to invoke destructor
}    

void
CreepAnalysis::clearAll(void)
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
  if (theEigenSOE != 0)
    delete theEigenSOE;

  // now set the pointers to NULL
  theAnalysisModel =0;
  theConstraintHandler =0;
  theDOF_Numberer =0;
  theIntegrator =0;
  theAlgorithm =0;
  theSOE =0;
  theEigenSOE =0;
  theTest = 0;

  // AddingSensitivity:BEGIN ////////////////////////////////////
#ifdef _RELIABILITY
  delete theSensitivityAlgorithm;
  theSensitivityAlgorithm =0;
#endif
  // AddingSensitivity:END //////////////////////////////////////
}    


int 
CreepAnalysis::analyze(int numSteps)
{
    int result = 0;
    Domain *the_Domain = this->getDomainPtr();

	ops_Creep = 1; //Added by AMK
	
    for (int i=0; i<numSteps; i++) {

	result = theAnalysisModel->analysisStep();
	if (result < 0) {
	    opserr << "CreepAnalysis::analyze() - the AnalysisModel failed";
	    opserr << " at iteration: " << i << " with domain at load factor ";
	    opserr << the_Domain->getCurrentTime() << endln;
	    the_Domain->revertToLastCommit();

	    return -2;
	}

	// check for change in Domain since last step. As a change can
	// occur in a commit() in a domaindecomp with load balancing
	// this must now be inside the loop

	int stamp = the_Domain->hasDomainChanged();

	if (stamp != domainStamp) {
	    domainStamp = stamp;

	    result = this->domainChanged();

	    if (result < 0) {
		opserr << "CreepAnalysis::analyze() - domainChanged failed";
		opserr << " at step " << i << " of " << numSteps << endln;
		return -1;
	    }	
	}


	result = theIntegrator->newStep();
	if (result < 0) {
	    opserr << "CreepAnalysis::analyze() - the Integrator failed";
	    opserr << " at iteration: " << i << " with domain at load factor ";
	    opserr << the_Domain->getCurrentTime() << endln;
	    the_Domain->revertToLastCommit();

	    return -2;
	}

		opserr << "\nCreepAnalysis at t = " << the_Domain->getCurrentTime(); //Added by AMK
	       opserr << ": i = "; //Added by AMK
	       
	       
	result = theAlgorithm->solveCurrentStep();
	if (result < 0) {
	    opserr << "CreepAnalysis::analyze() - the Algorithm failed";
	    opserr << " at iteration: " << i << " with domain at load factor ";
	    opserr << the_Domain->getCurrentTime() << endln;
	    the_Domain->revertToLastCommit();	    
	    theIntegrator->revertToLastStep();

	    return -3;
	}    

// AddingSensitivity:BEGIN ////////////////////////////////////
//#ifdef _RELIABILITY
	//if (theSensitivityAlgorithm != 0) {
		//result = theSensitivityAlgorithm->computeSensitivities();
		//if (result < 0) {
			//opserr << "CreepAnalysis::analyze() - the SensitivityAlgorithm failed";
			//opserr << " at iteration: " << i << " with domain at load factor ";
			//opserr << the_Domain->getCurrentTime() << endln;
			//the_Domain->revertToLastCommit();	    
			//theIntegrator->revertToLastStep();
			//return -5;
		//}    
	//}
//#endif
// AddingSensitivity:END //////////////////////////////////////

	result = theIntegrator->commit();
	if (result < 0) {
	    opserr << "CreepAnalysis::analyze() - ";
	    opserr << "the Integrator failed to commit";
	    opserr << " at iteration: " << i << " with domain at load factor ";
	    opserr << the_Domain->getCurrentTime() << endln;
	    the_Domain->revertToLastCommit();	    
	    theIntegrator->revertToLastStep();

	    return -4;
		}
    }
	
	ops_Creep = 0; //Added by AMK
    
    return 0;
}


int 
CreepAnalysis::eigen(int numMode, bool generalized, bool findSmallest)
{

    if (theAnalysisModel == 0 || theEigenSOE == 0) {
      opserr << "WARNING CreepAnalysis::eigen() - no EigenSOE has been set\n";
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
	opserr << "CreepAnalysis::eigen() - domainChanged failed";
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
	opserr << "WARNING CreepAnalysis::eigen() -";
	opserr << " failed in addA for ID " << elePtr->getID();	    
	result = -2;
      }
    }

    //
    // if generalized is true, form M
    //

    if (generalized == true) {
      FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
      while((elePtr = theEles2()) != 0) {     
	elePtr->zeroTangent();
	elePtr->addMtoTang(1.0);
	if (theEigenSOE->addM(elePtr->getTangent(0), elePtr->getID()) < 0) {
	  opserr << "WARNING CreepAnalysis::eigen() -";
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
	  opserr << "WARNING CreepAnalysis::eigen() -";
	  opserr << " failed in addM for ID " << dofPtr->getID();	    
	  result = -3;
	}
      }
    }
    
    // 
    // solve for the eigen values & vectors
    //

    if (theEigenSOE->solve(numMode, generalized) < 0) {
	opserr << "WARNING CreepAnalysis::eigen() - EigenSOE failed in solve()\n";
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
CreepAnalysis::initialize(void)
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
CreepAnalysis::domainChanged(void)
{
    int result = 0;

    Domain *the_Domain = this->getDomainPtr();
    int stamp = the_Domain->hasDomainChanged();
    domainStamp = stamp;

    // Timer theTimer; theTimer.start();
    // opserr << "CreepAnalysis::domainChanged(void)\n";

    theAnalysisModel->clearAll();    
    theConstraintHandler->clearAll();

    // theTimer.pause(); 
    // cout <<  "CreepAnalysis::clearAll() " << theTimer.getReal();
    // cout << theTimer.getCPU() << endln;
    // theTimer.start();    

    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.

    result = theConstraintHandler->handle();
    if (result < 0) {
	opserr << "CreepAnalysis::handle() - ";
	opserr << "ConstraintHandler::handle() failed";
	return -1;
    }

    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.

    result = theDOF_Numberer->numberDOF();
    if (result < 0) {
	opserr << "CreepAnalysis::handle() - ";
	opserr << "DOF_Numberer::numberDOF() failed";
	return -2;
    }	    

    result = theConstraintHandler->doneNumberingDOF();
    if (result < 0) {
	opserr << "CreepAnalysis::handle() - ";
	opserr << "ConstraintHandler::doneNumberingDOF() failed";
	return -2;
    }	    

    // we invoke setSize() on the LinearSOE which
    // causes that object to determine its size

    Graph &theGraph = theAnalysisModel->getDOFGraph();

    result = theSOE->setSize(theGraph);
    if (result < 0) {
	opserr << "CreepAnalysis::handle() - ";
	opserr << "LinearSOE::setSize() failed";
	return -3;
    }	    

    if (theEigenSOE != 0) {
      result = theEigenSOE->setSize(theGraph);
      if (result < 0) {
	opserr << "CreepAnalysis::handle() - ";
	opserr << "EigenSOE::setSize() failed";
	return -3;
      }	    
    }

    theAnalysisModel->clearDOFGraph();

    // finally we invoke domainChanged on the Integrator and Algorithm
    // objects .. informing them that the model has changed

    result = theIntegrator->domainChanged();
    if (result < 0) {
	opserr << "CreepAnalysis::setAlgorithm() - ";
	opserr << "Integrator::domainChanged() failed";
	return -4;
    }	    

    result = theAlgorithm->domainChanged();
    if (result < 0) {
	opserr << "CreepAnalysis::setAlgorithm() - ";
	opserr << "Algorithm::domainChanged() failed";
	return -5;
    }	        

    // if get here successfull
    return 0;
}    

// AddingSensitivity:BEGIN //////////////////////////////
#ifdef _RELIABILITY
int 
CreepAnalysis::setSensitivityAlgorithm(SensitivityAlgorithm *passedSensitivityAlgorithm)
{
    int result = 0;

    // invoke the destructor on the old one
    if (theSensitivityAlgorithm != 0) {
      delete theSensitivityAlgorithm;
    }
    
    theSensitivityAlgorithm = passedSensitivityAlgorithm;
    
    return 0;
}
#endif
// AddingSensitivity:END ///////////////////////////////


int 
CreepAnalysis::setNumberer(DOF_Numberer &theNewNumberer) 
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
CreepAnalysis::setAlgorithm(EquiSolnAlgo &theNewAlgorithm) 
{
    // invoke the destructor on the old one
    if (theAlgorithm != 0)
	delete theAlgorithm;

    // first set the links needed by the Algorithm
    theAlgorithm = &theNewAlgorithm;
    theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
    
    if (theTest != 0)
      theAlgorithm->setConvergenceTest(theTest);
    else   // this else is for backward compatability.
      theTest = theAlgorithm->getConvergenceTest();
    
    // invoke domainChanged() either indirectly or directly
    //    domainStamp = 0;
    if (domainStamp != 0)
      theAlgorithm->domainChanged();

    return 0;
}

int 
CreepAnalysis::setIntegrator(StaticIntegrator &theNewIntegrator)
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
    domainStamp = 0;

    /*
    if (domainStamp != 0)
      theIntegrator->domainChanged();
    */

  return 0;

}

int 
CreepAnalysis::setLinearSOE(LinearSOE &theNewSOE)
{
    // invoke the destructor on the old one
    if (theSOE != 0)
	delete theSOE;

    // set the links needed by the other objects in the aggregation
    theSOE = &theNewSOE;
    theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);
    theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
    theSOE->setLinks(*theAnalysisModel);

    // cause domainChanged to be invoked on next analyze
    /*
    if (domainStamp != 0) {
      Graph &theGraph = theAnalysisModel->getDOFGraph();
      int result = theSOE->setSize(theGraph);
    }
    */
    domainStamp = 0;
    return 0;
}


int 
CreepAnalysis::setEigenSOE(EigenSOE &theNewSOE)
{
    // invoke the destructor on the old one
    if (theEigenSOE != 0)
	delete theEigenSOE;

    // set the links needed by the other objects in the aggregation
    theEigenSOE = &theNewSOE;
    theEigenSOE->setLinks(*theAnalysisModel);

    // cause domainChanged to be invoked on next analyze
	domainStamp = 0;
    /*
	if (domainStamp != 0) {
      Graph &theGraph = theAnalysisModel->getDOFGraph();
      int result = theEigenSOE->setSize(theGraph);
    }
	*/
    return 0;
}


int 
CreepAnalysis::setConvergenceTest(ConvergenceTest &theNewTest)
{
    // invoke the destructor on the old one
    if (theTest != 0)
	delete theTest;

    // set the links needed by the other objects in the aggregation
    theTest = &theNewTest;

    theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);
    return theAlgorithm->setConvergenceTest(theTest);
}

EquiSolnAlgo *
CreepAnalysis::getAlgorithm(void)
{
  return theAlgorithm;
}

StaticIntegrator *
CreepAnalysis::getIntegrator(void)
{
  return theIntegrator;
}


ConvergenceTest *
CreepAnalysis::getConvergenceTest(void)
{
  return theTest;
}










