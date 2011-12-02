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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-12-01 00:07:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/StaticDomainDecompositionAnalysis.h,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for 
// StaticDomainDecompositionAnalysis.StaticDomainDecompositionAnalysis is a subclass 
// of DomainDecompositionAnalysis used to perform a static analysis step on a subdomain.

//
// What: "@(#) StaticDomainDecompositionAnalysis.h, revA"

#ifndef StaticDomainDecompositionAnalysis_h
#define StaticDomainDecompositionAnalysis_h

#ifndef _bool_h
#include <bool.h>
#endif

#include <DomainDecompositionAnalysis.h>
#include <Matrix.h>
#include <Vector.h>
#include <MovableObject.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class StaticIntegrator;
class ConvergenceTest;
class LinearSOE;
class LinearSolver;
class EquiSolnAlgo;
class Subdomain;
class ConvergenceTest;

class StaticDomainDecompositionAnalysis: public DomainDecompositionAnalysis
{
  public:
    StaticDomainDecompositionAnalysis(Subdomain &theDomain);
    StaticDomainDecompositionAnalysis(Subdomain &theDomain,
				      ConstraintHandler &theHandler,
				      DOF_Numberer &theNumberer,
				      AnalysisModel &theModel,
				      EquiSolnAlgo &theSolnAlgo,		   
				      LinearSOE &theSOE,
				      StaticIntegrator &theIntegrator,
				      ConvergenceTest *theTest,
				      bool setLinks = true);

    ~StaticDomainDecompositionAnalysis();
    void clearAll(void);	    
    int initialize(void);
    int domainChanged(void);

    // methods for non standard domain deomposition analysis
    int analyze(double dT);
    bool doesIndependentAnalysis(void);    

    // methods for standard domain deomposition analysis
    // that do some form of condensation to the tangent
    int  getNumExternalEqn(void);
    int  getNumInternalEqn(void);
    int  newStep(double dT);
    int  computeInternalResponse(void);
    int  formTangent(void);
    int  formResidual(void);
    int  formTangVectProduct(Vector &force);
    const Matrix &getTangent(void);
    const Vector &getResidual(void);
    const Vector &getTangVectProduct(void);

    // methods to change the analysis aggregates
    int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    int setIntegrator(IncrementalIntegrator &theIntegrator);
    int setLinearSOE(LinearSOE &theSOE);
    int setConvergenceTest(ConvergenceTest &theTest);

    // methods to send/receive
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected: 
    
  private:
    ConstraintHandler 	*theConstraintHandler;    
    DOF_Numberer 	*theDOF_Numberer;
    AnalysisModel 	*theAnalysisModel;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    StaticIntegrator    *theIntegrator;
    ConvergenceTest    *theTest;
    int domainStamp;			   
};

#endif











