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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:00:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomain.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/partitioned/PartitionedDomain.h
// 
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
class  ArrayOfTaggedObjects;
class  PartitionedDomainSubIter;
class  PartitionedDomainEleIter;
class SingleDomEleIter;

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

    // methods to remove the components     
    virtual void clearAll(void);
    virtual Element *removeElement(int tag);
    virtual Node *removeNode(int tag);        
    virtual SP_Constraint *removeSP_Constraint(int tag);
    virtual MP_Constraint *removeMP_Constraint(int tag);
    
    // methods to access the elements
    virtual  ElementIter       &getElements();
    virtual  Element           *getElement(int tag);
    
    virtual  int 		getNumElements(void) const;

    // public methods to update the domain
    virtual  void setCommitTag(int newTag);    	
    virtual  void setCurrentTime(double newTime);    
    virtual  void setCommittedTime(double newTime);        
    virtual  void applyLoad(double pseudoTime);
    virtual  void setLoadConstant(void);    

    virtual  int commit(void);    
    virtual  int revertToLastCommit(void);        
    virtual  int revertToStart(void);    

    virtual  int update(void);        
    virtual  int update(double newTime, double dT);

    
    virtual  void Print(OPS_Stream &s, int flag =0);    

    // public member functions in addition to the standard domain
    virtual int partition(int numPartitions);
			
    virtual bool addSubdomain(Subdomain *theSubdomain);
    virtual int getNumSubdomains(void);
    virtual Subdomain *getSubdomainPtr(int tag);
    virtual SubdomainIter &getSubdomains(void);
    virtual Node *removeExternalNode(int tag);        
    virtual Graph &getSubdomainGraph(void);
    
    
    // friend classes
    friend class PartitionedDomainEleIter;
    
  protected:    
    DomainPartitioner *getPartitioner(void) const;
    virtual int buildEleGraph(Graph *theEleGraph);
    
  private:
    TaggedObjectStorage  *elements;    
    ArrayOfTaggedObjects *theSubdomains;
    DomainPartitioner *theDomainPartitioner;

    SingleDomEleIter	       *mainEleIter;  // for ele that belong to elements
    PartitionedDomainSubIter   *theSubdomainIter;
    PartitionedDomainEleIter   *theEleIter;
    
    Graph *mySubdomainGraph;    // a graph of subdomain connectivity
};

#endif


