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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-11-18 01:59:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/ShadowSubdomain.h,v $
                                                                        
                                                                        
#ifndef ShadowSubdomain_h
#define ShadowSubdomain_h

// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for ShadowSubdomain.
// ShadowSubdomain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints that have been added to the ShadowSubdomain.
//
// What: "@(#) ShadowSubdomain.h, revA"

#include <Subdomain.h>
#include <Shadow.h>
#include <remote.h>

class ShadowSubdomain: public Shadow, public Subdomain
{
  public:
    ShadowSubdomain(int tag,
		    MachineBroker &theMachineBroker,
		    FEM_ObjectBroker &theObjectBroker);

    ShadowSubdomain(int tag, 
		    Channel &theChannel, 
		    FEM_ObjectBroker &theObjectBroker);

    virtual  ~ShadowSubdomain();    

    // method added for parallel domain generation
    virtual int buildSubdomain(int numSubdomains, 
			       PartitionedModelBuilder &theBuilder);
    virtual int getRemoteData(void);

    // Methods inherited from Domain, Subdomain and Element
    // which must be rewritten

    virtual  bool addElement(Element *);
    virtual  bool addNode(Node *);
    virtual  bool addExternalNode(Node *);
    virtual  bool addSP_Constraint(SP_Constraint *);
    virtual  bool addMP_Constraint(MP_Constraint *);    
    virtual  bool addLoadPattern(LoadPattern *);            
    virtual  bool addNodalLoad(NodalLoad *, int loadPattern);
    virtual  bool addElementalLoad(ElementalLoad *, int loadPattern);
    virtual  bool addSP_Constraint(SP_Constraint *, int loadPattern);    

    virtual Element 	  *removeElement(int tag);
    virtual Node 	  *removeNode(int tag);    
    virtual SP_Constraint *removeSP_Constraint(int tag);
    virtual MP_Constraint *removeMP_Constraint(int tag);
    virtual LoadPattern   *removeLoadPattern(int tag);
    virtual NodalLoad     *removeNodalLoad(int tag, int loadPattern);
    virtual ElementalLoad *removeElementalLoad(int tag, int loadPattern);
    virtual SP_Constraint *removeSP_Constraint(int tag, int loadPattern);
    
    virtual  ElementIter       &getElements();
    virtual  NodeIter          &getNodes();
    virtual  NodeIter	       &getInternalNodeIter(void);
    virtual  NodeIter	       &getExternalNodeIter(void);    
    
    virtual  Element       *getElementPtr(int tag);
    virtual  Node          *getNodePtr(int tag);

    virtual int getNumElements(void) const;
    virtual int getNumNodes(void) const;
    virtual int getNumSPs(void) const;
    virtual int getNumMPs(void) const;
    virtual int getNumLoadPatterns(void) const;
    
    virtual  Graph &getElementGraph(void);
    virtual  Graph &getNodeGraph(void);

    // methods to update the domain
    virtual  void setCommitTag(int newTag);    	
    virtual  void setCurrentTime(double newTime);    
    virtual  void setCommittedTime(double newTime);        
    virtual  void applyLoad(double pseudoTime);
    virtual  void setLoadConstant(void);    
    
    virtual  int update(void);    
    virtual  int commit(void);
    virtual  int revertToLastCommit(void);    
    virtual  int revertToStart(void);    

    virtual void setDomainDecompAnalysis(DomainDecompositionAnalysis &theAnalysis);
    virtual int setAnalysisAlgorithm(EquiSolnAlgo &theAlgorithm);
    virtual int setAnalysisIntegrator(IncrementalIntegrator &theIntegrator);
    virtual int setAnalysisLinearSOE(LinearSOE &theSOE);
    virtual void clearAnalysis(void);
    virtual void domainChange(void);
    
    virtual int 	getNumExternalNodes(void) const;    
    virtual const ID   &getExternalNodes(void);
    virtual int 	getNumDOF(void);

    virtual const Matrix &getTang(void);    
    virtual const Vector &getResistingForce(void);    

    virtual int  computeTang(void);
    virtual int  computeResidual(void);

    const Vector &getLastExternalSysResponse(void);
    virtual int computeNodalResponse(void);    
    virtual int newStep(double deltaT);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);    

    virtual double getCost(void);
    
    virtual  void Print(OPS_Stream &s, int flag =0);
    
  protected:    
    virtual int buildMap(void);
    virtual int buildEleGraph(Graph *theEleGraph);
    virtual int buildNodeGraph(Graph *theNodeGraph);    
    
  private:
    ID msgData;
    ID theElements;
    ID theNodes;
    ID theExternalNodes;    
    ID theLoadCases;
    
    int numDOF;
    int numElements;
    int numNodes;
    int numExternalNodes;
    int numSPs;
    int numMPs;
    int numLoadPatterns;

    bool buildRemote;
    bool gotRemoteData;
    
    FE_Element *theFEele;

    Vector *theVector; // for storing residual info
    Matrix *theMatrix; // for storing tangent info
    
    static char *shadowSubdomainProgram;

    static int count;  // MHS
    static int numShadowSubdomains;
    static ShadowSubdomain **theShadowSubdomains;
};

#endif
