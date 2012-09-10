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
                                                                        
// $Revision: 1.5 $
// $Date: 2009-05-14 22:52:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.h,v $
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for FEM_ObjectBroker.
// FEM_ObjectBroker is is an object broker class for the finite element
// method. All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
// What: "@(#) FEM_ObjectBrokerAllClasses.h, revA"


#ifndef FEM_ObjectBrokerAllClasses_h
#define FEM_ObjectBrokerAllClasses_h

#include<FEM_ObjectBroker.h>

class FEM_ObjectBrokerAllClasses : public FEM_ObjectBroker
{
  public:
    FEM_ObjectBrokerAllClasses();
    ~FEM_ObjectBrokerAllClasses();

    Actor*getNewActor(int classTag, Channel *theChannel);
    
    PartitionedModelBuilder *
      getPtrNewPartitionedModelBuilder(Subdomain &theSub,
				       int classTag);
    
    GraphNumberer *getPtrNewGraphNumberer(int classTag);
    
    // methods to get new modelling class objects
    Element       *getNewElement(int classTag);
    Node          *getNewNode(int classTag);
    MP_Constraint *getNewMP(int classTag);
    SP_Constraint *getNewSP(int classTag);
    Pressure_Constraint *getNewPC(int classTag);
    NodalLoad     *getNewNodalLoad(int classTag);
    ElementalLoad *getNewElementalLoad(int classTag);
    
    CrdTransf *getNewCrdTransf(int classTag);

    BeamIntegration *getNewBeamIntegration(int classTag);

    UniaxialMaterial  *getNewUniaxialMaterial(int classTag);
    SectionForceDeformation  *getNewSection(int classTag);    
    NDMaterial *getNewNDMaterial(int classTag);
    Fiber *getNewFiber(int classTag);
    FrictionModel *getNewFrictionModel(int classTag);

    ConvergenceTest       *getNewConvergenceTest(int classTag);
    LoadPattern           *getNewLoadPattern(int classTag);
    GroundMotion          *getNewGroundMotion(int classTag);
    TimeSeries            *getNewTimeSeries(int classTag);    
    TimeSeriesIntegrator  *getNewTimeSeriesIntegrator(int classTag);    
    
    // matrix vector and id objects
    Matrix	  *getPtrNewMatrix(int classTag, int noRows, int noCols);
    Vector	  *getPtrNewVector(int classTag, int size);
    ID	          *getPtrNewID(int classTag, int size);

    // methods for ouput objects
    //    DataOutputHandler *getPtrNewDataOutputHandler(int classTag);
    OPS_Stream *getPtrNewStream(int classTag);
    Recorder *getPtrNewRecorder(int classTag);
    
    
    // methods to get new analysis objects
    ConstraintHandler     *getNewConstraintHandler(int classTag);
    DOF_Numberer          *getNewNumberer(int classTag);
    AnalysisModel         *getNewAnalysisModel(int classTag);
    EquiSolnAlgo          *getNewEquiSolnAlgo(int classTag);
    Accelerator           *getAccelerator(int classTag);
    LineSearch            *getLineSearch(int classTag);
    DomainDecompAlgo      *getNewDomainDecompAlgo(int classTag);
    StaticIntegrator      *getNewStaticIntegrator(int classTag);
    TransientIntegrator   *getNewTransientIntegrator(int classTag);
    IncrementalIntegrator *getNewIncrementalIntegrator(int classTag);

    LinearSOE *getNewLinearSOE(int classTagSOE);
    EigenSOE *getNewEigenSOE(int classTagSOE);
    
    LinearSOE *getPtrNewDDLinearSOE(int classTagSOE, 
				    int classTagDDSolver);

    DomainSolver *getNewDomainSolver(void);

    DomainDecompositionAnalysis *
      getNewDomainDecompAnalysis(int classTag, Subdomain &theDomain);

    Subdomain  *getSubdomainPtr(int classTag);

    Parameter *getParameter(int classTag);

    int addUniaxialMaterial(int classTag, const char *lib, const char *funcName, UniaxialMaterial *(*)(void));
    
  protected:
    
  private:
    DomainSolver *lastDomainSolver;
};

#endif


