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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-04 22:37:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/objectBroker/FEM_ObjectBroker.h,v $
                                                                        
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for FEM_ObjectBroker.
// FEM_ObjectBroker is is an object broker class for the finite element
// method. All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
// What: "@(#) FEM_ObjectBroker.h, revA"


#ifndef FEM_ObjectBroker_h
#define FEM_ObjectBroker_h

#include<ObjectBroker.h>

class Element;
class Node;
class MP_Constraint;
class SP_Constraint;
class NodalLoad;
class ElementalLoad;
class LoadPattern;
class TimeSeries;
class TimeSeriesIntegrator;

class Matrix;
class Vector;
class ID;
class Subdomain;
class ConstraintHandler;
class DOF_Numberer;   
class AnalysisModel;    
class EquiSolnAlgo;
class DomainDecompAlgo;
class StaticIntegrator;
class TransientIntegrator;
class IncrementalIntegrator;
class LinearSOE;
class DomainSolver;
class DomainDecompositionAnalysis;
class LinearSOESolver;
class PartitionedModelBuilder;

class CrdTransf2d;
class CrdTransf3d;

class BeamIntegration;

class UniaxialMaterial;
class SectionForceDeformation;
class NDMaterial;
class Fiber;

class ConvergenceTest;
class SectionForceDeformation;
class GroundMotion;
class OPS_Stream;
class Recorder;

class Actor;

class FEM_ObjectBroker
{
  public:
    FEM_ObjectBroker();
    virtual ~FEM_ObjectBroker();

    virtual Actor*getNewActor(int classTag, Channel *theChannel);
    
    virtual PartitionedModelBuilder *
      getPtrNewPartitionedModelBuilder(Subdomain &theSub,
				       int classTag);
    
    virtual GraphNumberer *getPtrNewGraphNumberer(int classTag);
    
    // methods to get new modelling class objects
    virtual Element       *getNewElement(int classTag);
    virtual Node          *getNewNode(int classTag);
    virtual MP_Constraint *getNewMP(int classTag);
    virtual SP_Constraint *getNewSP(int classTag);
    virtual NodalLoad     *getNewNodalLoad(int classTag);
    virtual ElementalLoad *getNewElementalLoad(int classTag);
    
    virtual CrdTransf2d *getNewCrdTransf2d(int classTag);
    virtual CrdTransf3d *getNewCrdTransf3d(int classTag);

    virtual BeamIntegration *getNewBeamIntegration(int classTag);

    virtual UniaxialMaterial  *getNewUniaxialMaterial(int classTag);
    virtual SectionForceDeformation  *getNewSection(int classTag);    
    virtual NDMaterial *getNewNDMaterial(int classTag);
    virtual Fiber *getNewFiber(int classTag);

    virtual ConvergenceTest *getNewConvergenceTest(int classTag);
    virtual LoadPattern *getNewLoadPattern(int classTag);
    virtual GroundMotion *getNewGroundMotion(int classTag);
    virtual TimeSeries  *getNewTimeSeries(int classTag);    
    virtual TimeSeriesIntegrator  *getNewTimeSeriesIntegrator(int classTag);    
    
    // matrix vector and id objects
    virtual Matrix	  *getPtrNewMatrix(int classTag, int noRows, int noCols);
    virtual Vector	  *getPtrNewVector(int classTag, int size);
    virtual ID	          *getPtrNewID(int classTag, int size);

    // methods for ouput objects
    //    virtual DataOutputHandler *getPtrNewDataOutputHandler(int classTag);
    virtual OPS_Stream *getPtrNewStream(int classTag);
    virtual Recorder *getPtrNewRecorder(int classTag);
    
    
    // methods to get new analysis objects
    virtual ConstraintHandler   *getNewConstraintHandler(int classTag);
    virtual DOF_Numberer        *getNewNumberer(int classTag);
    virtual AnalysisModel       *getNewAnalysisModel(int classTag);
    virtual EquiSolnAlgo        *getNewEquiSolnAlgo(int classTag);
    virtual DomainDecompAlgo    *getNewDomainDecompAlgo(int classTag);
    virtual StaticIntegrator    *getNewStaticIntegrator(int classTag);
    virtual TransientIntegrator *getNewTransientIntegrator(int classTag);
    virtual IncrementalIntegrator *getNewIncrementalIntegrator(int classTag);

    virtual LinearSOE *getNewLinearSOE(int classTagSOE, int classTagSolver);
    virtual LinearSOESolver *getNewLinearSolver(void);
    
    virtual LinearSOE *getPtrNewDDLinearSOE(int classTagSOE, 
					    int classTagDDSolver);
    virtual DomainSolver *getNewDomainSolver(void);

    virtual DomainDecompositionAnalysis *
      getNewDomainDecompAnalysis(int classTag, Subdomain &theDomain);

    virtual Subdomain  *getSubdomainPtr(int classTag);


    virtual int addUniaxialMaterial(int classTag, const char *lib, const char *funcName, UniaxialMaterial *(*)(void));
    
  protected:
    
  private:
    LinearSOESolver *lastLinearSolver;
    DomainSolver *lastDomainSolver;

    
};

#endif


