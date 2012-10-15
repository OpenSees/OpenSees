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
                                                                        
// $Revision: 1.48 $
// $Date: 2009-05-14 22:52:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/objectBroker/FEM_ObjectBroker.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for FEM_ObjectBroker.
// FEM_ObjectBroker is is an object broker class for the finite element
// method. All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
// What: "@(#) FEM_ObjectBroker.C, revA"


#include <FEM_ObjectBroker.h>


FEM_ObjectBroker::FEM_ObjectBroker()
:lastDomainSolver(0)
{

}


FEM_ObjectBroker::~FEM_ObjectBroker()
{

}


Actor *
FEM_ObjectBroker::getNewActor(int classTag, Channel *theChannel)
{
  return 0;
}


PartitionedModelBuilder          *
FEM_ObjectBroker::getPtrNewPartitionedModelBuilder(Subdomain &theSubdomain,
						   int classTag)
{
  return 0;
}


GraphNumberer *
FEM_ObjectBroker::getPtrNewGraphNumberer(int classTag)
{
  return 0;
}

/*****************************************
 *
 * METHODS TO GET NEW MODELLING CLASSES
 *
 *****************************************/


Element       *
FEM_ObjectBroker::getNewElement(int classTag)
{
  return 0;
}
				
Node          *
FEM_ObjectBroker::getNewNode(int classTag)
{
  return 0;
}


MP_Constraint *
FEM_ObjectBroker::getNewMP(int classTag)
{
  return 0;
}


SP_Constraint *
FEM_ObjectBroker::getNewSP(int classTag)
{
  return 0;
}

Pressure_Constraint *
FEM_ObjectBroker::getNewPC(int classTag)
{
  return 0;
}

NodalLoad     *
FEM_ObjectBroker::getNewNodalLoad(int classTag)
{
  return 0;
}


ElementalLoad *
FEM_ObjectBroker::getNewElementalLoad(int classTag)
{
  return 0;
}

CrdTransf*
FEM_ObjectBroker::getNewCrdTransf(int classTag)
{
  return 0;
}

BeamIntegration *
FEM_ObjectBroker::getNewBeamIntegration(int classTag)
{
  return 0;
}


UniaxialMaterial *
FEM_ObjectBroker::getNewUniaxialMaterial(int classTag)
{
  return 0;
}

SectionForceDeformation *
FEM_ObjectBroker::getNewSection(int classTag)
{
  return 0;
}

NDMaterial*
FEM_ObjectBroker::getNewNDMaterial(int classTag)
{
  return 0;   
}

Fiber*
FEM_ObjectBroker::getNewFiber(int classTag)
{
  return 0;
}

FrictionModel *
FEM_ObjectBroker::getNewFrictionModel(int classTag)
{
  return 0;
}

ConvergenceTest *
FEM_ObjectBroker::getNewConvergenceTest(int classTag)
{
  return 0;
}


LoadPattern *
FEM_ObjectBroker::getNewLoadPattern(int classTag)
{
  return 0;
}


GroundMotion *
FEM_ObjectBroker::getNewGroundMotion(int classTag)
{
  return 0;
}

TimeSeries *
FEM_ObjectBroker::getNewTimeSeries(int classTag)
{
  return 0;
}

TimeSeriesIntegrator *
FEM_ObjectBroker::getNewTimeSeriesIntegrator(int classTag)
{
  return 0;
}


Matrix	  *
FEM_ObjectBroker::getPtrNewMatrix(int classTag, int noRows, int noCols)
{
  return 0;
}


Vector	  *
FEM_ObjectBroker::getPtrNewVector(int classTag, int size)
{
  return 0;
}


ID	          *
FEM_ObjectBroker::getPtrNewID(int classTag, int size)
{
  return 0;
}

/*****************************************
 *
 * METHODS TO GET NEW OUTPUT CLASS OBJECTS
 *
 *****************************************/

OPS_Stream *
FEM_ObjectBroker::getPtrNewStream(int classTag)
{
  return 0;
}

Recorder *
FEM_ObjectBroker::getPtrNewRecorder(int classTag)
{
  return 0;
}



/*****************************************
 *
 * METHODS TO GET NEW ANALYSIS CLASSES
 *
 *****************************************/

ConstraintHandler   *
FEM_ObjectBroker::getNewConstraintHandler(int classTag)
{
  return 0;
}


DOF_Numberer        *
FEM_ObjectBroker::getNewNumberer(int classTag)
{
  return 0;
}


AnalysisModel       *
FEM_ObjectBroker::getNewAnalysisModel(int classTag)
{
  return 0;
}


EquiSolnAlgo        *
FEM_ObjectBroker::getNewEquiSolnAlgo(int classTag)
{
  return 0;
}

Accelerator        *
FEM_ObjectBroker::getAccelerator(int classTag)
{
  return 0;
}


LineSearch        *
FEM_ObjectBroker::getLineSearch(int classTag)
{
  return 0;
}


DomainDecompAlgo    *
FEM_ObjectBroker::getNewDomainDecompAlgo(int classTag)
{
  return 0;
}


StaticIntegrator    *
FEM_ObjectBroker::getNewStaticIntegrator(int classTag)
{
  return 0;
}


TransientIntegrator *
FEM_ObjectBroker::getNewTransientIntegrator(int classTag)
{
  return 0;
}


IncrementalIntegrator *
FEM_ObjectBroker::getNewIncrementalIntegrator(int classTag)
{
  return 0;
}

LinearSOE *
FEM_ObjectBroker::getNewLinearSOE(int classTagSOE)
{
  return 0;		 
}


EigenSOE *
FEM_ObjectBroker::getNewEigenSOE(int classTagSOE)
{
  return 0;		 
}

DomainSolver *
FEM_ObjectBroker::getNewDomainSolver(void)
{
    return 0;
}
    
LinearSOE *
FEM_ObjectBroker::getPtrNewDDLinearSOE(int classTagSOE, 
				       int classTagDDSolver)
{
  return 0;		 
}


DomainDecompositionAnalysis *
FEM_ObjectBroker::getNewDomainDecompAnalysis(int classTag, 
						Subdomain &theSubdomain)
{
  return 0;
}


Subdomain 	  *
FEM_ObjectBroker::getSubdomainPtr(int classTag)
{
    return 0;
}


int 
FEM_ObjectBroker::addUniaxialMaterial(int classTag, 
				      const char *lib, 
				      const char *funcName, 
				      UniaxialMaterial *(*funcPtr)(void))
{
  return 0;
}


Parameter *
FEM_ObjectBroker::getParameter(int classTag)
{
  return 0;
}

