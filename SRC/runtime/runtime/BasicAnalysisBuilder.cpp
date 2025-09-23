//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: Claudio Perez
//
#include <assert.h>
#include <stdio.h>
#include <unordered_map>

#include "BasicAnalysisBuilder.h"
#include <Domain.h>
#include <G3_Logging.h>
// Abstract classes
#include <EquiSolnAlgo.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <LinearSOE.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <AnalysisModel.h>
#include <TimeSeries.h>
#include <LoadPattern.h>

// For eigen()
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>

// Default concrete analysis classes
#include <Newmark.h>
#include <EigenSOE.h>
#include <SymBandEigenSolver.h>
#include <SymBandEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <ArpackSOE.h>
#include <ProfileSPDLinSOE.h>
#include <NewtonRaphson.h>
#include <RCM.h>
#include <LoadControl.h>
#include <ProfileSPDLinSolver.h>
#include <CTestNormUnbalance.h>
#include <ProfileSPDLinDirectSolver.h>
#include <PlainHandler.h>
#include <TransformationConstraintHandler.h>


static std::unordered_map<int, std::string> SolveFailedMessage {
   {SolutionAlgorithm::BadFormResidual, "Failed to form residual\n"},
   {SolutionAlgorithm::BadFormTangent,  "Failed to form tangent\n"},
   {SolutionAlgorithm::BadLinearSolve,  "Failed to solve system, tangent may be singular\n"},
// {SolutionAlgorithm::TestFailed,      ""},// no output; information will have been printed by the test
   {SolutionAlgorithm::BadTestStart,    "Failed to initialize the convergence test\n"},
};

BasicAnalysisBuilder::BasicAnalysisBuilder(Domain* domain)
:
  theDomain(domain),
  theHandler(nullptr),
  theNumberer(nullptr),
  theAlgorithm(nullptr),
  theSOE(nullptr),
  theEigenSOE(nullptr),
  theStaticIntegrator(nullptr),
  theTransientIntegrator(nullptr),
  theTest(nullptr),
  theVariableTimeStepTransientAnalysis(nullptr),
 CurrentAnalysisFlag(EMPTY_ANALYSIS)
{
  theAnalysisModel = new AnalysisModel();
}

BasicAnalysisBuilder::~BasicAnalysisBuilder()
{
  this->wipe();

  if (theAnalysisModel != nullptr) {
    delete theAnalysisModel;
    theAnalysisModel = nullptr;
  }
}

void
BasicAnalysisBuilder::wipe()
{

  if (theAlgorithm != nullptr) {
      delete theAlgorithm;
      theAlgorithm = nullptr;
  }
  if (theStaticIntegrator != nullptr) {
      delete theStaticIntegrator;
      theStaticIntegrator = nullptr;
  }
  if ((theTransientIntegrator != nullptr) && freeTI) {
      delete theTransientIntegrator;
      theTransientIntegrator = nullptr;
  }
  if ((theSOE != nullptr) && freeSOE) {
      delete theSOE;
      theSOE = nullptr;
  }
  if (theNumberer != nullptr) {
      delete theNumberer;
      theNumberer = nullptr;
  }
  if (theHandler != nullptr) {
      delete theHandler;
      theHandler = nullptr;
  }
  if (theTest != nullptr) {
      delete theTest;
      theTest = nullptr;
  }
  if (theEigenSOE != nullptr) {
      delete theEigenSOE;
      theEigenSOE = nullptr;
  }
  if (theAnalysisModel != nullptr) {
    delete theAnalysisModel;
    theAnalysisModel = new AnalysisModel();
  }
  theVariableTimeStepTransientAnalysis = nullptr;
}

void
BasicAnalysisBuilder::setLinks(CurrentAnalysis flag)
{
  if (theSOE && theAnalysisModel)
    theSOE->setLinks(*theAnalysisModel);

  if (theDomain && theHandler && theAnalysisModel)
    theAnalysisModel->setLinks(*theDomain, *theHandler);

  if (theAnalysisModel && theNumberer)
    theNumberer->setLinks(*theAnalysisModel);

  if (theTest && theAlgorithm)
    theAlgorithm->setConvergenceTest(theTest);


  switch (flag) {
  case EMPTY_ANALYSIS:
    break;

  case TRANSIENT_ANALYSIS:
    if (theDomain && theAnalysisModel && theTransientIntegrator && theHandler)
      theHandler->setLinks(*theDomain, *theAnalysisModel, *theTransientIntegrator);

    if (theAnalysisModel && theTransientIntegrator && theSOE && theTest && theAlgorithm)
      theAlgorithm->setLinks(*theAnalysisModel, *theTransientIntegrator, *theSOE, theTest);

    if (theAnalysisModel && theSOE && theTest && theTransientIntegrator) {
      theTransientIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);
    }
    // if (theTransientIntegrator && domainStamp != 0)
    //   theTransientIntegrator->domainChanged();
      // this->domainChanged();

    // domainStamp  = 0;
    break;

  case STATIC_ANALYSIS:
    // opserr << "setLinks(STATIC)\n";
    if (theDomain && theAnalysisModel && theStaticIntegrator && theHandler)
      theHandler->setLinks(*theDomain, *theAnalysisModel, *theStaticIntegrator);

    if (theAnalysisModel && theSOE && theTest && theStaticIntegrator)
      theStaticIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);

    if (theAnalysisModel && theStaticIntegrator && theSOE && theTest && theAlgorithm)
      theAlgorithm->setLinks(*theAnalysisModel, *theStaticIntegrator, *theSOE, theTest);

    // domainStamp  = 0;
    break;
  }
}

int
BasicAnalysisBuilder::initialize(void)
{
  // check if domain has undergone change
  int stamp = theDomain->hasDomainChanged();
  if (stamp != domainStamp) {
    domainStamp = stamp;
    if (this->domainChanged() < 0) {
      opserr << G3_ERROR_PROMPT << "initialize - domainChanged() failed\n";
      return -1;
    }
  }

  switch (this->CurrentAnalysisFlag) {
    case EMPTY_ANALYSIS:
      break;

    case STATIC_ANALYSIS:
      if (theStaticIntegrator->initialize() < 0) {
          opserr << G3_WARN_PROMPT << "integrator initialize() failed\n";
          return -2;
      } else
        theStaticIntegrator->commit();
      break;

    case TRANSIENT_ANALYSIS:
      if (theTransientIntegrator->initialize() < 0) {
          opserr << "integrator initialize() failed\n";
          return -2;
      } else
        theTransientIntegrator->commit();
      break;
  }
  theDomain->initialize();

  return 0;
}

int
BasicAnalysisBuilder::domainChanged(void)
{
  Domain *domain = this->getDomain();
  int stamp = domain->hasDomainChanged();
  domainStamp = stamp;

  opsdbg << G3_DEBUG_PROMPT << "Domain changed\n";

  theAnalysisModel->clearAll();
  if (theHandler != nullptr) {
    theHandler->clearAll();

    // Invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.
    if (theHandler->handle() < 0) {
      opserr << "BasicAnalysisBuilder::domainChange() - ConstraintHandler::handle() failed\n";
      return -1;
    }
    // Invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.
    if (theNumberer != nullptr && theNumberer->numberDOF() < 0) {
      opserr << "BasicAnalysisBuilder::domainChange() - DOF_Numberer::numberDOF() failed\n";
      return -2;
    }

    if (theHandler->doneNumberingDOF() < 0) {
      opserr << "BasicAnalysisBuilder::domainChange() - ConstraintHandler::doneNumberingDOF() failed\n";
      return -2;
    }
  }

  // Invoke setSize() on the LinearSOE which
  // causes that object to determine its size
  Graph &theGraph = theAnalysisModel->getDOFGraph();

  if (theSOE != nullptr) {
    if (theSOE->setSize(theGraph) < 0) {
      opserr << "BasicAnalysisBuilder::domainChange() - LinearSOE::setSize() failed\n";
      return -3;
    }
  }

  if (theEigenSOE != nullptr) {
    int result = theEigenSOE->setSize(theGraph);
    if (result < 0) {
      return -3;
    }
  }

  theAnalysisModel->clearDOFGraph();

  // finally we invoke domainChanged on the Integrator and Algorithm
  // objects .. informing them that the model has changed
  switch (this->CurrentAnalysisFlag) {

  case STATIC_ANALYSIS:
    if (theStaticIntegrator->domainChanged() < 0) {
      opserr << "BasicAnalysisBuilder::domainChange - Integrator::domainChanged() failed\n";
      return -4;
    }
    break;

  case TRANSIENT_ANALYSIS:

    if (theTransientIntegrator->domainChanged() < 0) {
      opserr << "BasicAnalysisBuilder::domainChange - Integrator::domainChanged() failed\n";
      return -4;
    }
    break;
  default:
    break;
  }

//  if (theAlgorithm && theAlgorithm->domainChanged() < 0) {
//    opserr << "BasicAnalysisBuilder::domainChange - Algorithm::domainChanged failed\n";
//    return -5;
//  }

  return 0;
}

int
BasicAnalysisBuilder::analyze(int num_steps, double size_steps)
{

  switch (this->CurrentAnalysisFlag) {

    case STATIC_ANALYSIS:
      return this->analyzeStatic(num_steps);
      break;

    case TRANSIENT_ANALYSIS: {
      // TODO: Set global timestep variable
      ops_Dt = size_steps;
      return this->analyzeTransient(num_steps, size_steps);
      break;
    }

    default:
      opserr << G3_ERROR_PROMPT << "No Analysis type has been specified \n";
      return -1;
  }
}

int
BasicAnalysisBuilder::analyzeStatic(int numSteps)
{
  int result = 0;

  for (int i=0; i<numSteps; i++) {

      result = theAnalysisModel->analysisStep();

      if (result < 0) {
        opserr << "StaticAnalysis::analyze - the AnalysisModel failed\n";
        opserr << " at step: " << i << " with domain at load factor ";
        opserr << theDomain->getCurrentTime() << "\n";
        theDomain->revertToLastCommit();
        return -2;
      }

      // Check for change in Domain since last step. As a change can
      // occur in a commit() in a domaindecomp with load balancing
      // this must now be inside the loop
      int stamp = theDomain->hasDomainChanged();

      if (stamp != domainStamp) {
        domainStamp = stamp;
        result = this->domainChanged();
        if (result < 0) {
          opserr << "domainChanged failed";
          opserr << " at step " << i << " of " << numSteps << "\n";
          return -1;
        }
      }

      result = theStaticIntegrator->newStep();
      if (result < 0) {
        opserr << "The Integrator failed at step: " << i
               << " with domain at load factor " << theDomain->getCurrentTime() << "\n";
        theDomain->revertToLastCommit();
        theStaticIntegrator->revertToLastStep();
        return -2;
      }

      result = theAlgorithm->solveCurrentStep();
      if (result < 0) {
        // Print error message if we have one
        if (SolveFailedMessage.find(result) != SolveFailedMessage.end()) {
            opserr << OpenSees::PromptAnalysisFailure << SolveFailedMessage[result];
        }
        theDomain->revertToLastCommit();
        theStaticIntegrator->revertToLastStep();
        return -3;
      }

      result = theStaticIntegrator->commit();
      if (result < 0) {
        opserr << "StaticAnalysis::analyze - ";
        opserr << "the Integrator failed to commit";
        opserr << " at step: " << i << " with domain at load factor ";
        opserr << theDomain->getCurrentTime() << "\n";

        theDomain->revertToLastCommit();
        theStaticIntegrator->revertToLastStep();
        return -4;
      }
  }

  return 0;
}

int
BasicAnalysisBuilder::analyzeTransient(int numSteps, double dT)
{
  int result = 0;

  for (int i=0; i<numSteps; i++) {
    result = this->analyzeStep(dT);
    if (result < 0) {
      if (numSubLevels != 0)
        result = this->analyzeSubLevel(1, dT);
      if (result < 0)
        return result;
    }
  }
  return result;
}

int
BasicAnalysisBuilder::analyzeSubLevel(int level, double dT)
{
  int result = 0;
  if (numSubSteps == 0)
    return -1;

  double stepDT = dT/(numSubSteps*1.);

  for (int i=0; i<numSubSteps; i++) {
    result = this->analyzeStep(stepDT);
    if (result < 0) {
      if (level == numSubLevels) {
        return result;
      } else {
        result = this->analyzeSubLevel(level+1, stepDT);
        if (result < 0)
          return result;
      }
    }
  }
  return result;
}

// analyze a transient step
int
BasicAnalysisBuilder::analyzeStep(double dT)
{
  int result = 0;
  if (theAnalysisModel->analysisStep(dT) < 0) {
    opserr << "DirectIntegrationAnalysis::analyze() - the AnalysisModel failed";
    opserr << " at time " << theDomain->getCurrentTime() << "\n";
    theDomain->revertToLastCommit();
    return -2;
  }

  // check if domain has undergone change
  int stamp = theDomain->hasDomainChanged();
  if (stamp != domainStamp) {
    domainStamp = stamp;
    if (this->domainChanged() < 0) {
      opserr << "DirectIntegrationAnalysis::analyze() - domainChanged() failed\n";
      return -1;
    }
  }

  if (theTransientIntegrator->newStep(dT) < 0) {
    opserr << "DirectIntegrationAnalysis::analyze() - the Integrator failed";
    opserr << " at time " << theDomain->getCurrentTime() << "\n";
    theDomain->revertToLastCommit();
    theTransientIntegrator->revertToLastStep();
    return -2;
  }

  result = theAlgorithm->solveCurrentStep();
  if (result < 0) {
    if (SolveFailedMessage.find(result) != SolveFailedMessage.end()) {
        opserr << OpenSees::PromptAnalysisFailure << SolveFailedMessage[result];
    }
    theDomain->revertToLastCommit();
    theTransientIntegrator->revertToLastStep();
    return -3;
  }


  result = theTransientIntegrator->commit();
  if (result < 0) {
    opserr << "DirectIntegrationAnalysis::analyze() - ";
    opserr << "the Integrator failed to commit";
    opserr << " at time " << theDomain->getCurrentTime() << "\n";
    theDomain->revertToLastCommit();
    theTransientIntegrator->revertToLastStep();
    return -4;
  }

  return result;
}



void
BasicAnalysisBuilder::set(ConstraintHandler* obj)
{
  if (theHandler != nullptr)
    delete theHandler;

  theHandler = obj;
}

void
BasicAnalysisBuilder::set(DOF_Numberer* obj)
{
  // free the old numberer
  if (theNumberer != nullptr)
    delete theNumberer;

  // set the links needed by the Algorithm
  theNumberer = obj;
  theNumberer->setLinks(*theAnalysisModel);

  domainStamp = 0;
  return;
}

void
BasicAnalysisBuilder::set(EquiSolnAlgo* obj)
{
  if (theAlgorithm != nullptr)
    delete theAlgorithm;

  theAlgorithm = obj;

  if (theTest != nullptr)
    theAlgorithm->setConvergenceTest(theTest);
  else   // this else is for backward compatibility. // ?
    theTest = theAlgorithm->getConvergenceTest();

  this->setLinks(this->CurrentAnalysisFlag);
}

void
BasicAnalysisBuilder::set(LinearSOE* obj, bool free)
{

  // if free is false then we cant free either
  if ((theSOE != nullptr) && free && freeSOE)
    delete theSOE;

  freeSOE = free;

  theSOE = obj;

  this->setLinks(this->CurrentAnalysisFlag);

  if (theEigenSOE != nullptr)
    theEigenSOE->setLinearSOE(*theSOE);


  domainStamp = 0;
}


LinearSOE*
BasicAnalysisBuilder::getLinearSOE() {
  return theSOE;
}


void
BasicAnalysisBuilder::set(StaticIntegrator& obj)
{
  if (theStaticIntegrator != nullptr)
    delete theStaticIntegrator;

  theStaticIntegrator = &obj;

  this->setLinks(STATIC_ANALYSIS);

  if (domainStamp != 0 && this->CurrentAnalysisFlag != EMPTY_ANALYSIS)
    theStaticIntegrator->domainChanged();

  else
    domainStamp = 0;
}

void
BasicAnalysisBuilder::set(TransientIntegrator& obj, bool free)
{

  if ((theTransientIntegrator != nullptr) && free && freeTI)
    delete theTransientIntegrator;

  freeTI = free;

  theTransientIntegrator = &obj;

  this->setLinks(TRANSIENT_ANALYSIS);

  if (domainStamp != 0  && this->CurrentAnalysisFlag != EMPTY_ANALYSIS)
    theTransientIntegrator->domainChanged();

  else
    domainStamp = 0;
}

void
BasicAnalysisBuilder::set(ConvergenceTest* obj)
{

  if (theTest != nullptr)
    delete theTest;

  theTest = obj;
  this->setLinks(this->CurrentAnalysisFlag);

}

void
BasicAnalysisBuilder::set(EigenSOE &theNewSOE)
{
  // invoke the destructor on the old one if not the same!
  if (theEigenSOE != nullptr) {
    if (theEigenSOE->getClassTag() != theNewSOE.getClassTag()) {
      delete theEigenSOE;
      theEigenSOE = nullptr;
    }
  }

  if (theEigenSOE == nullptr) {
    theEigenSOE = &theNewSOE;
    theEigenSOE->setLinks(*theAnalysisModel);
    theEigenSOE->setLinearSOE(*theSOE);

    domainStamp = 0;
  }

}

void
BasicAnalysisBuilder::fillDefaults(BasicAnalysisBuilder::CurrentAnalysis flag)
{

  switch (flag) {
    case EMPTY_ANALYSIS:
      break;

    case STATIC_ANALYSIS:
      if (theStaticIntegrator == nullptr)
        theStaticIntegrator = new LoadControl(1, 1, 1, 1);
      break;

    case TRANSIENT_ANALYSIS:
      if (theTransientIntegrator == nullptr)
          theTransientIntegrator = new Newmark(0.5,0.25);
      break;
  }

  if (theTest == nullptr)
    theTest = new CTestNormUnbalance(1.0e-6, 25, ConvergenceTest::PrintFailure);

  if (theAlgorithm == nullptr)
    theAlgorithm = new NewtonRaphson(*theTest);


  if (theHandler == nullptr) {
    // Dont show a warning if the user has no constraints
    if (theDomain->getNumMPs() > 0) {
      opserr << G3_WARN_PROMPT << "constraints were used but no ConstraintHandler has been specified; \n";
      opserr << "        PlainHandler default will be used\n";
    }
    theHandler = new PlainHandler();
  }

  if (theNumberer == nullptr)
    theNumberer = new DOF_Numberer(*(new RCM(false)));

  if (theSOE == nullptr)
    // TODO: CHANGE TO MORE GENERAL SOE
      theSOE = new ProfileSPDLinSOE(*(new ProfileSPDLinDirectSolver()));

}


int
BasicAnalysisBuilder::setStaticAnalysis()
{
  domainStamp = 0;
  this->fillDefaults(STATIC_ANALYSIS);
  this->setLinks(STATIC_ANALYSIS);

  this->CurrentAnalysisFlag = STATIC_ANALYSIS;

  return 0;
}

int
BasicAnalysisBuilder::setTransientAnalysis()
{
  domainStamp = 0;
  this->CurrentAnalysisFlag = TRANSIENT_ANALYSIS;
  this->fillDefaults(TRANSIENT_ANALYSIS);
  this->setLinks(TRANSIENT_ANALYSIS);

  return 1;
}

int
BasicAnalysisBuilder::newTransientAnalysis()
{
    assert(theDomain != nullptr);

    this->fillDefaults(TRANSIENT_ANALYSIS);

    return 1;
}


void
BasicAnalysisBuilder::newEigenAnalysis(int typeSolver, double shift)
{
  assert(theAnalysisModel != nullptr);

  if (theHandler == nullptr)
    theHandler = new TransformationConstraintHandler();

  // this->CurrentAnalysisFlag = TRANSIENT_ANALYSIS;
  if (this->CurrentAnalysisFlag == EMPTY_ANALYSIS)
    this->CurrentAnalysisFlag = TRANSIENT_ANALYSIS;
  this->fillDefaults(this->CurrentAnalysisFlag); //TRANSIENT_ANALYSIS);
  this->setLinks(this->CurrentAnalysisFlag); //TRANSIENT_ANALYSIS);

  // create a new eigen system and solver
  if (theEigenSOE != nullptr) {
    if (theEigenSOE->getClassTag() != typeSolver) {
      delete theEigenSOE;
      theEigenSOE = nullptr;
    }
  }

  if (theEigenSOE == nullptr) {
    domainStamp = 0;
    if (typeSolver == EigenSOE_TAGS_SymBandEigenSOE) {
      SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver();
      theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theAnalysisModel);

    } else if (typeSolver == EigenSOE_TAGS_FullGenEigenSOE) {
        FullGenEigenSolver *theEigenSolver = new FullGenEigenSolver();
        theEigenSOE = new FullGenEigenSOE(*theEigenSolver, *theAnalysisModel);

    } else {
        theEigenSOE = new ArpackSOE(shift);
    }

    //
    // set the eigen soe in the system
    //
    theEigenSOE->setLinks(*theAnalysisModel);
    theEigenSOE->setLinearSOE(*theSOE);
  } // theEigenSOE == 0
}

int
BasicAnalysisBuilder::eigen(int numMode, bool generalized, bool findSmallest)
{
  // TODO: merge with newEigenAnalysis

  assert(theAnalysisModel != nullptr);
  assert(     theEigenSOE != nullptr);

  int result = 0;
  Domain *the_Domain = this->getDomain();

  // for parallel processing, want all analysis doing an eigenvalue analysis
  result = theAnalysisModel->eigenAnalysis(numMode, generalized, findSmallest);

  int stamp = the_Domain->hasDomainChanged();

  if (stamp != domainStamp) {
    //domainStamp = stamp; // commented out so domainChanged() gets called with integrator,
                         //  which isnt updated here
//    result = this->domainChanged();

    theAnalysisModel->clearAll();
    theHandler->clearAll();

    // Now invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.
    result = theHandler->handle();

    // Now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.
    result = theNumberer->numberDOF();

    result = theHandler->doneNumberingDOF();

    Graph &theGraph = theAnalysisModel->getDOFGraph();

    result = theSOE->setSize(theGraph);

    result = theEigenSOE->setSize(theGraph);

    theAnalysisModel->clearDOFGraph();

    if (result < 0) {
      opserr << "BasicAnalysisBuilder::eigen() - domainChanged failed\n";
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
  while ((elePtr = theEles()) != nullptr) {
    elePtr->zeroTangent();
    elePtr->addKtToTang(1.0);
    if (theEigenSOE->addA(elePtr->getTangent(0), elePtr->getID()) < 0) {
      opserr << G3_WARN_PROMPT << "eigen -";
      opserr << " failed in addA for ID " << elePtr->getID();
      result = -2;
    }
  }

  //
  // If generalized is true, form M
  //
  if (generalized == true) {
    FE_EleIter &theEles2 = theAnalysisModel->getFEs();
    while ((elePtr = theEles2()) != nullptr) {
      elePtr->zeroTangent();
      elePtr->addMtoTang(1.0);
      if (theEigenSOE->addM(elePtr->getTangent(0), elePtr->getID()) < 0) {
        opserr << "WARNING BasicAnalysisBuilder::eigen -";
        opserr << " failed in addA for ID " << elePtr->getID() << "\n";
        result = -2;
      }
    }

    DOF_Group *dofPtr;
    DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();
    while ((dofPtr = theDofs()) != nullptr) {
      dofPtr->zeroTangent();
      dofPtr->addMtoTang(1.0);
      if (theEigenSOE->addM(dofPtr->getTangent(0), dofPtr->getID()) < 0) {
        opserr << G3_WARN_PROMPT << "theEigenSOE failed in addM for ID " << dofPtr->getID() << "\n";
        result = -3;
      }
    }
  }

  //
  // Solve for the eigen values & vectors
  //
  if (theEigenSOE->solve(numMode, generalized, findSmallest) < 0) {
      opserr << G3_WARN_PROMPT << "EigenSOE failed in solve()\n";
      return -4;
  }

  //
  // Store the eigenvalues and eigenvectors in the model
  //
  theAnalysisModel->setNumEigenvectors(numMode);
  Vector theEigenvalues(numMode);
  for (int i = 1; i <= numMode; i++) {
    theEigenvalues[i-1] = theEigenSOE->getEigenvalue(i);
    theAnalysisModel->setEigenvector(i, theEigenSOE->getEigenvector(i));
  }
  theAnalysisModel->setEigenvalues(theEigenvalues);
  this->numEigen = numMode;

  return 0;
}

Domain*
BasicAnalysisBuilder::getDomain()
{
  return theDomain;
}

EquiSolnAlgo*
BasicAnalysisBuilder::getAlgorithm()
{
  return theAlgorithm;
}

StaticIntegrator*
BasicAnalysisBuilder::getStaticIntegrator()
{
  return theStaticIntegrator;
}

TransientIntegrator*
BasicAnalysisBuilder::getTransientIntegrator() {

  return theTransientIntegrator;
}

ConvergenceTest*
BasicAnalysisBuilder::getConvergenceTest()
{
  return theTest;
}

int
BasicAnalysisBuilder::formUnbalance()
{
    if (theStaticIntegrator != nullptr)
      return theStaticIntegrator->formUnbalance();

    else if (theTransientIntegrator != nullptr)
      return theTransientIntegrator->formUnbalance();

    return -1;
}

