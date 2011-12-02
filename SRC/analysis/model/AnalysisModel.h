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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/AnalysisModel.h,v $
                                                                        
                                                                        
#ifndef AnalysisModel_h
#define AnalysisModel_h

// File: ~/analysis/model/AnalysisModel.h
// 
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

#ifndef _bool_h
#include <bool.h>
#endif
    
#include <SimpleDOF_Iter.h>
#include <SimpleFE_Iter.h>

#include <MovableObject.h>

class Domain;
class FE_EleIter;
class DOF_GrpIter;
class Graph;
class FE_Element;
class DOF_Group;
class Vector;
class DOF_Graph;
class DOF_GroupGraph;
class FEM_ObjectBroker;

class AnalysisModel: public MovableObject
{
  public:
    AnalysisModel();    
    AnalysisModel(int classTag);
    virtual ~AnalysisModel();    

    // methods to populate/depopulate the AnalysisModel
    virtual bool addFE_Element(FE_Element *theFE_Ele);
    virtual bool addDOF_Group(DOF_Group *theDOF_Grp);
    virtual void clearAll(void);
    
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
    
    // methods which trigger operations in the Domain
    virtual void setLinks(Domain &theDomain);
	
    virtual void applyLoadDomain(double pseudoTime);
    virtual void updateDomain(void);
    virtual int commitDomain(void);
    virtual int revertDomainToLastCommit(void);
    virtual double getCurrentDomainTime(void);
    virtual void   setCurrentDomainTime(double newTime);    
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    friend class SimpleFE_Iter;
    friend class SimpleDOF_Iter;    
    
  protected:    
    Domain *getDomainPtr(void) const;
    
  private:
    Domain *myDomain;
    DOF_Graph *myDOFGraph;
    DOF_GroupGraph *myGroupGraph;    
    
    int numFE_Ele;             // number of FE_Elements objects added
    int numDOF_Grp;            // number of DOF_Group objects added
    int numEqn;                // numEqn set by the ConstraintHandler typically

    int sizeEle;               // size of array for FE_Elements
    int sizeDOF;               // size of array for DOF_Groups
    
    FE_Element   **theFEs;     // array of pointers to the FE_Elements
    DOF_Group    **theDOFs;    // array of pointers to the DOF_Groups

    SimpleFE_Iter  theFEiter;     
    SimpleDOF_Iter theDOFiter;    
   
};

#endif
