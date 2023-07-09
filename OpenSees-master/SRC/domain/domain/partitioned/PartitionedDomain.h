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
                                                                        
// $Revision: 1.14 $
// $Date: 2010-09-16 00:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomain.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: Wed Sep 25 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for PartitionedDomain.
// PartitionedDomain is an abstract class. The class is responsible for holding
// and providing access to the Elements, Nodes, SP_Constraints 
// and MP_Constraints just like a normal domain. In addition the domain provides
// a method to partition the domain into Subdomains.
//
// ModelBuilder. There are no partitions in a PartitionedDomain.
//
// What: "@(#) PartitionedDomain.h, revA"

#ifndef PartitionedDomain_h
#define PartitionedDomain_h

#include <Domain.h>

class DomainPartitioner;
class Subdomain;
class SubdomainIter; 
class ArrayOfTaggedObjects;
class PartitionedDomainSubIter;
class PartitionedDomainEleIter;
class SingleDomEleIter;
class Parameter;
class GraphPartitioner;

class PartitionedDomain: public Domain
{
  public:
    PartitionedDomain();    
    PartitionedDomain(DomainPartitioner &thePartitioner);    

    PartitionedDomain(int numNodes, int numElements, 
		      int numSPs, int numMPs, int numLoadPatterns,
		      int numSubdomains,
		      DomainPartitioner &thePartitioner);
    
    virtual  ~PartitionedDomain();    

    // public methods to populate a domain	
    virtual  bool addElement(Element *elePtr);
    virtual  bool addNode(Node *nodePtr);

    virtual  bool addLoadPattern(LoadPattern *);            
    virtual  bool addSP_Constraint(SP_Constraint *); 
    virtual  int  addSP_Constraint(int axisDirn, double axisValue, 
				   const ID &fixityCodes, double tol=1e-10);
    virtual  bool addSP_Constraint(SP_Constraint *, int loadPatternTag); 
    virtual  bool addMP_Constraint(MP_Constraint *); 

    virtual  bool addNodalLoad(NodalLoad *, int loadPatternTag);
    virtual  bool addElementalLoad(ElementalLoad *, int loadPatternTag);

    // methods to remove the components     
    virtual void clearAll(void);
    virtual Element *removeElement(int tag);
    virtual Node *removeNode(int tag);        
    virtual SP_Constraint *removeSP_Constraint(int tag);
    virtual int removeSP_Constraint(int nodeTag, int dof, int loadPatternTag);
    virtual MP_Constraint *removeMP_Constraint(int tag);
    virtual int removeMP_Constraints(int tag);
    virtual LoadPattern   *removeLoadPattern(int loadTag);
    
    // methods to access the elements
    virtual  ElementIter       &getElements();
    virtual  Element           *getElement(int tag);
    virtual  int 		getNumElements(void) const;

    // public methods to update the domain
    virtual int hasDomainChanged(void);

    virtual  void setCommitTag(int newTag);    	
    virtual  void setCurrentTime(double newTime);    
    virtual  void setCommittedTime(double newTime);        
    virtual  void applyLoad(double pseudoTime);
    virtual  void setLoadConstant(void);    
    virtual  int  setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);

    virtual  int commit(void);    
    virtual  int revertToLastCommit(void);        
    virtual  int revertToStart(void);    

    virtual  int update(void);        
    virtual  int update(double newTime, double dT);

    virtual  int analysisStep(double dT);
    virtual  int eigenAnalysis(int, bool, bool);

    virtual int  record(bool fromAnalysis = true);    
    virtual int  addRecorder(Recorder &theRecorder);    	
    virtual int  removeRecorders(void);
    virtual int  removeRecorder(int tag);
    
    virtual  void Print(OPS_Stream &s, int flag =0);    
    virtual void Print(OPS_Stream &s, ID *nodeTags, ID *eleTags, int flag =0);

    // methods dealing with Parameters
    virtual  bool           addParameter(Parameter *);            
    virtual  Parameter     *removeParameter(int tag);
    virtual  int  updateParameter(int tag, int value);
    virtual  int  updateParameter(int tag, double value);    

    // public member functions in addition to the standard domain
    virtual int setPartitioner(DomainPartitioner *thePartitioner);
    virtual int partition(int numPartitions, bool usingMain = false, int mainPartitionID = 0, int specialElementTag = 0);
    virtual int repartition(int numPartitions, bool usingMain = false, int mainPartitionID = 0, int specialElementTag = 0);
			
    virtual bool addSubdomain(Subdomain *theSubdomain);
    virtual int getNumSubdomains(void);
    virtual Subdomain *getSubdomainPtr(int tag);
    virtual SubdomainIter &getSubdomains(void);
    virtual Node *removeExternalNode(int tag);        
    virtual Graph &getSubdomainGraph(void);

    // nodal methods required in domain interface for parallel interprter
    virtual const Vector *getNodeResponse(int nodeTag, NodeResponseType); 
    virtual const Vector *getElementResponse(int eleTag, const char **argv, int argc); 

    virtual double getNodeDisp(int nodeTag, int dof, int &errorFlag);
    virtual int setMass(const Matrix &mass, int nodeTag);

    virtual int calculateNodalReactions(bool inclInertia);
    
    virtual int activateElements(const ID& elementList);
    virtual int deactivateElements(const ID& elementList);

    // friend classes
    friend class PartitionedDomainEleIter;
    
  protected:    
    int barrierCheck(int result);        
    GraphPartitioner *getGraphPartitioner(void);
    DomainPartitioner *getPartitioner(void) const;
    virtual int buildEleGraph(Graph *theEleGraph);
    virtual TaggedObjectStorage* getElementsStorage();
    
    
    TaggedObjectStorage  *elements;    
    ArrayOfTaggedObjects *theSubdomains;
    DomainPartitioner    *theDomainPartitioner;

    SingleDomEleIter	       *mainEleIter;  // for ele that belong to elements
    PartitionedDomainSubIter   *theSubdomainIter;
    PartitionedDomainEleIter   *theEleIter;
    
    Graph *mySubdomainGraph;    // a graph of subdomain connectivity

    bool has_sent_yet;
};

#endif


