//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
// Revision: A
//
// Purpose: This file contains the class definition for FEM_ObjectBroker.
// FEM_ObjectBroker is is an object broker class for the finite element
// method. All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
// What: "@(#) TclPackageClassBroker.h, revA"

#ifndef TclPackageClassBroker_h
#define TclPackageClassBroker_h

#include <FEM_ObjectBroker.h>

class TclPackageClassBroker : public FEM_ObjectBroker {
public:
  TclPackageClassBroker();
  ~TclPackageClassBroker();

  Actor *getNewActor(int classTag, Channel *theChannel);

  PartitionedModelBuilder *getPtrNewPartitionedModelBuilder(Subdomain &theSub,
                                                            int classTag);

  GraphNumberer *getPtrNewGraphNumberer(int classTag);

  // methods to get new modelling class objects
  Element *getNewElement(int classTag);
  Node *getNewNode(int classTag);
  MP_Constraint *getNewMP(int classTag);
  SP_Constraint *getNewSP(int classTag);
  Pressure_Constraint *getNewPC(int classTag);
  ElementalLoad *getNewElementalLoad(int classTag);

  CrdTransf *getNewCrdTransf(int classTag);

  BeamIntegration *getNewBeamIntegration(int classTag);

  UniaxialMaterial *getNewUniaxialMaterial(int classTag);
  SectionForceDeformation *getNewSection(int classTag);
  NDMaterial *getNewNDMaterial(int classTag);
  Fiber *getNewFiber(int classTag);
  FrictionModel *getNewFrictionModel(int classTag);

  ConvergenceTest *getNewConvergenceTest(int classTag);
  LoadPattern *getNewLoadPattern(int classTag);
  GroundMotion *getNewGroundMotion(int classTag);
  TimeSeries *getNewTimeSeries(int classTag);
  TimeSeriesIntegrator *getNewTimeSeriesIntegrator(int classTag);

  // matrix vector and id objects
  Matrix *getPtrNewMatrix(int classTag, int noRows, int noCols);
  Vector *getPtrNewVector(int classTag, int size);
  ID *getPtrNewID(int classTag, int size);

  // methods for ouput objects
  //    DataOutputHandler *getPtrNewDataOutputHandler(int classTag);
  OPS_Stream *getPtrNewStream(int classTag);
  Recorder *getPtrNewRecorder(int classTag);

  // methods to get new analysis objects
  ConstraintHandler *getNewConstraintHandler(int classTag);
  DOF_Numberer *getNewNumberer(int classTag);
  AnalysisModel *getNewAnalysisModel(int classTag);
  EquiSolnAlgo *getNewEquiSolnAlgo(int classTag);
  Accelerator *getAccelerator(int classTag);
  LineSearch *getLineSearch(int classTag);
  DomainDecompAlgo *getNewDomainDecompAlgo(int classTag);
  StaticIntegrator *getNewStaticIntegrator(int classTag);
  TransientIntegrator *getNewTransientIntegrator(int classTag);
  IncrementalIntegrator *getNewIncrementalIntegrator(int classTag);

  LinearSOE *getNewLinearSOE(int classTagSOE);
  EigenSOE *getNewEigenSOE(int classTagSOE);

  LinearSOE *getPtrNewDDLinearSOE(int classTagSOE, int classTagDDSolver);

  DomainSolver *getNewDomainSolver(void);

  DomainDecompositionAnalysis *getNewDomainDecompAnalysis(int classTag,
                                                          Subdomain &theDomain);

  Subdomain *getSubdomainPtr(int classTag);

  Parameter *getParameter(int classTag);

  int addUniaxialMaterial(int classTag, const char *lib, const char *funcName,
                          UniaxialMaterial *(*)(void));

protected:
private:
  DomainSolver *lastDomainSolver;
};

#endif
