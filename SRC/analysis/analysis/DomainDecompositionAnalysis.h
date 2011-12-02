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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DomainDecompositionAnalysis.h,v $
                                                                        
                                                                        
// File: ~/analysis/method/DomainDecompositionAnalysis.h
// 
// Written: fmk 
// Created: Tue Sept 17 16:34:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// DomainDecompositionAnalysis. DomainDecompositionAnalysis is a subclass 
// of AnalysisAnalysis, it is used when performing a domain decomposition
// analysis. It provides methods which can be invoked by a subdomain to
// perform the numerical computations required.
//
// What: "@(#) DomainDecompositionAnalysis.h, revA"

#ifndef DomainDecompositionAnalysis_h
#define DomainDecompositionAnalysis_h

#ifndef _bool_h
#include <bool.h>
#endif

#include <Analysis.h>
#include <Matrix.h>
#include <Vector.h>
#include <MovableObject.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class IncrementalIntegrator;
class LinearSOE;
class DomainSolver;
class DomainDecompAlgo;
class Subdomain;
class Vector;

class DomainDecompositionAnalysis: public Analysis, public MovableObject
{
  public:
    DomainDecompositionAnalysis(Subdomain &theDomain);

    DomainDecompositionAnalysis(int classTag, 
				Subdomain &theDomain);    
    
    DomainDecompositionAnalysis(Subdomain &theDomain,
		   ConstraintHandler &theHandler,
		   DOF_Numberer &theNumberer,
		   AnalysisModel &theModel,
		   DomainDecompAlgo &theSolnAlgo,		   
		   IncrementalIntegrator &theIntegrator,	
		   LinearSOE &theSOE,
		   DomainSolver &theSolver);


    virtual ~DomainDecompositionAnalysis();
    
    virtual int  analyze(void);
    virtual int domainChanged(void);

    virtual int  getNumExternalEqn(void);
    virtual int  getNumInternalEqn(void);
    virtual int  computeInternalResponse(void);
    virtual int  formTangent(void);
    virtual int  formResidual(void);
    virtual int  formTangVectProduct(Vector &force);
    virtual const Matrix &getTangent(void);
    virtual const Vector &getResidual(void);
    virtual const Vector &getTangVectProduct(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
  protected: 
    Subdomain		*getSubdomainPtr(void) const;
    ConstraintHandler 	*getConstraintHandlerPtr(void) const;
    DOF_Numberer 	*getDOF_NumbererPtr(void) const;
    AnalysisModel 	*getAnalysisModelPtr(void) const;
    DomainDecompAlgo    *getDomainDecompAlgoPtr(void) const;        
    IncrementalIntegrator  *getIncrementalIntegratorPtr(void) const;    
    LinearSOE 		*getLinSOEPtr(void) const;
    DomainSolver        *getDomainSolverPtr(void) const;
    

    
  private:
    Subdomain 		     *theSubdomain;
    ConstraintHandler 	     *theHandler;    
    DOF_Numberer 	     *theNumberer;
    AnalysisModel 	     *theModel;
    DomainDecompAlgo 	     *theAlgorithm;
    IncrementalIntegrator    *theIntegrator;        
    LinearSOE 		     *theSOE;
    DomainSolver	     *theSolver;

    Vector 		     *theResidual;
    int numEqn;
    int numExtEqn;

    // the following 2 variables are used to allow formResidual()
    // and formTangVectProduct() to be called before formTangent()
    // - this must be allowed as typical elements will not have to fromTangent
    // before being asked to form Residual(). 
    bool tangFormed;
    int tangFormedCount; // saves the expense of computing formTangent() 
	               // for same state of Subdomain.
    int domainStamp;			   
};

#endif











