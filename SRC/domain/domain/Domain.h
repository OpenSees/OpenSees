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
                                                                        
// $Revision: 1.16 $
// $Date: 2005-11-30 23:32:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/Domain.h,v $
                                                                        
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for Domain.
// Domain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, SP_Constraints 
// MP_Constraints, and LoadPatterns.
//
// What: "@(#) Domain.h, revA"

#ifndef Domain_h
#define Domain_h

#include <OPS_Stream.h>
#include <Vector.h>

#ifndef _bool_h
#include <bool.h>
#endif

class Element;
class Node;
class SP_Constraint;
class MP_Constraint;
class NodalLoad;
class ElementalLoad;
class LoadPattern;

class ElementIter;
class NodeIter;
class SP_ConstraintIter;
class MP_ConstraintIter;
class LoadPatternIter;

class SingleDomEleIter;
class SingleDomNodIter;
class SingleDomSP_Iter;
class SingleDomMP_Iter;
class  SingleDomAllSP_Iter;

class MeshRegion;
class Recorder;
class Graph;
class NodeGraph;
class ElementGraph;
class Channel;
class FEM_ObjectBroker;

class TaggedObjectStorage;

class Domain
{
  public:
    Domain();
    Domain(int numNodes, int numElements, int numSPs, int numMPs,
	   int numLoadPatterns);
    
    Domain(TaggedObjectStorage &theNodesStorage,
	   TaggedObjectStorage &theElementsStorage,
	   TaggedObjectStorage &theMPsStorage,
	   TaggedObjectStorage &theSPsStorage,
	   TaggedObjectStorage &theLoadPatternsStorage);

    Domain(TaggedObjectStorage &theStorageType);
    
    virtual ~Domain();    

    // methods to populate a domain
    virtual  bool addElement(Element *);
    virtual  bool addNode(Node *);
    virtual  bool addSP_Constraint(SP_Constraint *);
    virtual  bool addMP_Constraint(MP_Constraint *); 
    virtual  bool addLoadPattern(LoadPattern *);            
    
    // methods to add components to a LoadPattern object
    virtual  bool addSP_Constraint(SP_Constraint *, int loadPatternTag); 
    virtual  bool addNodalLoad(NodalLoad *, int loadPatternTag);
    virtual  bool addElementalLoad(ElementalLoad *, int loadPatternTag);
    
    // methods to remove the components 
    virtual void clearAll(void);	
    virtual Element       *removeElement(int tag);
    virtual Node          *removeNode(int tag);    
    virtual SP_Constraint *removeSP_Constraint(int tag);
    virtual MP_Constraint *removeMP_Constraint(int tag);    
    virtual LoadPattern   *removeLoadPattern(int loadTag);

    virtual NodalLoad     *removeNodalLoad(int tag, int loadPattern);
    virtual ElementalLoad *removeElementalLoad(int tag, int loadPattern);
    virtual SP_Constraint *removeSP_Constraint(int tag, int loadPattern);
    
    // methods to access the components of a domain
    virtual  ElementIter       &getElements();
    virtual  NodeIter          &getNodes();
    virtual  SP_ConstraintIter &getSPs();
    virtual  MP_ConstraintIter &getMPs();
    virtual  LoadPatternIter   &getLoadPatterns();
    virtual  SP_ConstraintIter &getDomainAndLoadPatternSPs();
    
    virtual  Element       *getElement(int tag);
    virtual  Node          *getNode(int tag);
    virtual  SP_Constraint *getSP_Constraint(int tag);    
    virtual  MP_Constraint *getMP_Constraint(int tag);    
    virtual  LoadPattern   *getLoadPattern(int tag);        

    // methods to query the state of the domain
    virtual double  getCurrentTime(void) const;
    virtual int     getCommitTag(void) const;    	
    virtual int getNumElements(void) const;
    virtual int getNumNodes(void) const;
    virtual int getNumSPs(void) const;
    virtual int getNumMPs(void) const;
    virtual int getNumLoadPatterns(void) const;            
    virtual const Vector &getPhysicalBounds(void); 


    // methods to get element and node graphs
    virtual  Graph  &getElementGraph(void);
    virtual  Graph  &getNodeGraph(void);
    
    // methods to update the domain
    virtual  void setCommitTag(int newTag);    	
    virtual  void setCurrentTime(double newTime);    
    virtual  void setCommittedTime(double newTime);        
    virtual  void applyLoad(double pseudoTime);
    virtual  void setLoadConstant(void);    
    virtual  int  initialize(void);    
    virtual  int  setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);
    
    virtual  int  commit(void);
    virtual  int  revertToLastCommit(void);
    virtual  int  revertToStart(void);    
    virtual  int  update(void);
    virtual  int  update(double newTime, double dT);
    virtual  int  newStep(double dT);

    
    // methods for eigenvalue analysis
    virtual int setEigenvalues(const Vector &theEigenvalues);
    virtual const Vector &getEigenvalues(void);
    virtual double getTimeEigenvaluesSet(void);
    
    // methods for other objects to determine if model has changed
    virtual void domainChange(void);    
    virtual int hasDomainChanged(void);
    virtual void setDomainChangeStamp(int newStamp);
    
    // methods for output
    virtual int  addRecorder(Recorder &theRecorder);    	
    virtual int  removeRecorders(void);

    virtual int  addRegion(MeshRegion &theRegion);    	
    virtual MeshRegion *getRegion(int region);    	

    virtual void Print(OPS_Stream &s, int flag =0);
    friend OPS_Stream &operator<<(OPS_Stream &s, Domain &M);    

    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);    

    // nodal methods required in domain interface for parallel interprter
    virtual double getNodeDisp(int nodeTag, int dof, int &errorFlag);
    virtual int setMass(const Matrix &mass, int nodeTag);

    virtual int calculateNodalReactions(bool inclInertia);

  protected:    
    virtual int buildEleGraph(Graph *theEleGraph);
    virtual int buildNodeGraph(Graph *theNodeGraph);

    Recorder **theRecorders;
    int numRecorders;    

  private:
    double currentTime;               // current pseudo time
    double committedTime;             // the committed pseudo time
    double dT;                        // difference between committed and current time
    int	   currentGeoTag;             // an integer used to mark if domain has changed
    bool   hasDomainChangedFlag;      // a bool flag used to indicate if GeoTag needs to be ++
    int    theDbTag;                   // the Domains unique database tag == 0
    int    lastGeoSendTag;            // the value of currentGeoTag when sendSelf was last invoked
    int dbEle, dbNod, dbSPs, dbMPs, dbLPs; // database tags for storing info

    bool eleGraphBuiltFlag;
    bool nodeGraphBuiltFlag;
    
    Graph *theNodeGraph;
    Graph *theElementGraph;

    TaggedObjectStorage  *theElements;
    TaggedObjectStorage  *theNodes;
    TaggedObjectStorage  *theSPs;    
    TaggedObjectStorage  *theMPs;    
    TaggedObjectStorage  *theLoadPatterns;        

    SingleDomEleIter      *theEleIter;
    SingleDomNodIter  	  *theNodIter;
    SingleDomSP_Iter      *theSP_Iter;
    SingleDomMP_Iter      *theMP_Iter;
    LoadPatternIter       *theLoadPatternIter;        
    SingleDomAllSP_Iter   *allSP_Iter;
    
    MeshRegion **theRegions;
    int numRegions;    

    int commitTag;
    
    Vector theBounds;
    
    Vector *theEigenvalues;
    double theEigenvalueSetTime;

    int lastChannel;
};

#endif


