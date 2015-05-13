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
                                                                        
// $Revision: 1.12 $
// $Date: 2009-08-25 23:18:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/AnalysisModel.h,v $
                                                                        
                                                                        
#ifndef AnalysisModel_h
#define AnalysisModel_h

// Written: fmk 
// Created: 9/96
// Revision: A
//
// Description: This file contains the class definition for AnalysisModel.
// AnalysisModel is a container class. This class is responsible for holding
// and providing access to the FE_Element and DOF_Group objects that the 
// ConstraintHandler creates. It is also responsible for updating the 
// response quantities at the DOF_Groups and for triggering methods 
// in the associated Domain.
//
// What: "@(#) AnalysisModel.h, revA"

#include <MovableObject.h>

class TaggedObjectStorage;
class Domain;
class FE_EleIter;
class DOF_GrpIter;
class Graph;
class FE_Element;
class DOF_Group;
class Vector;
class FEM_ObjectBroker;
class ConstraintHandler;

class AnalysisModel: public MovableObject
{
  public:
    AnalysisModel();    
    AnalysisModel(int classTag);
    AnalysisModel(TaggedObjectStorage &theDofStorage,
		  TaggedObjectStorage &theFeStorage);
    virtual ~AnalysisModel();    

    // methods to populate/depopulate the AnalysisModel
    virtual bool addFE_Element(FE_Element *theFE_Ele);
    virtual bool addDOF_Group(DOF_Group *theDOF_Grp);
    virtual void clearAll(void);
    virtual void clearDOFGraph(void);
    virtual void clearDOFGroupGraph(void);
    
    // methods to access the FE_Elements and DOF_Groups and their numbers
    virtual int getNumDOF_Groups(void) const;		
    virtual DOF_Group *getDOF_GroupPtr(int tag);	
    virtual FE_EleIter &getFEs();
    virtual DOF_GrpIter &getDOFs();

    // method to access the connectivity for SysOfEqn to size itself
    virtual void setNumEqn(int) ;	
    virtual int getNumEqn(void) const ; 
    virtual Graph &getDOFGraph(void);
    virtual Graph &getDOFGroupGraph(void);
    
    // methods to update the response quantities at the DOF_Groups,
    // which in turn set the new nodal trial response quantities.
    virtual void setResponse(const Vector &disp, 
			     const Vector &vel, 
			     const Vector &accel);
    virtual void setDisp(const Vector &disp);    
    virtual void setVel(const Vector &vel);        
    virtual void setAccel(const Vector &vel);            

    virtual void incrDisp(const Vector &disp);    
    virtual void incrVel(const Vector &vel);        
    virtual void incrAccel(const Vector &vel);            

    // methods added to store the eigenvalues and vectors in the domain
    virtual void setNumEigenvectors(int numEigenvectors);
    virtual void setEigenvector(int mode, const Vector &);
    virtual void setEigenvalues(const Vector &);    
    virtual const Vector &getEigenvalues(void);    
    const Vector *getModalDampingFactors(void);
    bool inclModalDampingMatrix(void);
    
    // methods which trigger operations in the Domain
    virtual void setLinks(Domain &theDomain, ConstraintHandler &theHandler);
	
    virtual void   applyLoadDomain(double newTime);
    virtual int    updateDomain(void);
    virtual int    updateDomain(double newTime, double dT);
    virtual int    analysisStep(double dT =0.0);
    virtual int    eigenAnalysis(int numMode, bool generalized, bool findSmallest);
    virtual int    commitDomain(void);
    virtual int    revertDomainToLastCommit(void);
    virtual double getCurrentDomainTime(void);
    virtual void   setCurrentDomainTime(double newTime);    
    virtual void   setRayleighDampingFactors(double alphaM, double betaK, double betaKi, double betaKc);    
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    Domain *getDomainPtr(void) const;

  protected:

    
  private:
    Domain *myDomain;
    ConstraintHandler *myHandler;

    Graph *myDOFGraph;
    Graph *myGroupGraph;    
    
    int numFE_Ele;             // number of FE_Elements objects added
    int numDOF_Grp;            // number of DOF_Group objects added
    int numEqn;                // numEqn set by the ConstraintHandler typically

    TaggedObjectStorage  *theFEs;
    TaggedObjectStorage  *theDOFs;
    
    FE_EleIter    *theFEiter;     
    DOF_GrpIter   *theDOFiter;    
};

#endif
