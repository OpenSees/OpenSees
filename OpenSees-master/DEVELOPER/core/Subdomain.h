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
                                                                        
// $Revision: 1.13 $
// $Date: 2009-08-26 20:33:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/Subdomain.h,v $
                                                                        
                                                                        
#ifndef Subdomain_h
#define Subdomain_h

// Written: fmk 
// Created:  11/96
// Revision: A
// Revision: B 03/98 - revised to allow parallel model generation
//
// Description: This file contains the class definition for Subdomain.
// Subdomain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints that have been added to the subdomain.
//
// What: "@(#) Subdomain.h, revA"

#include <Domain.h>
#include <Element.h>

class Node;
class ID;
class TaggedObjectStorage;
class DomainDecompositionAnalysis;
//class PartitionedModelBuilder;
class EquiSolnAlgo;
class IncrementalIntegrator;
class LinearSOE;
class EigenSOE;
class ConvergenceTest;
class FE_Element;

#include <SubdomainNodIter.h>

class Subdomain: public Element, public Domain
{
  public:
    Subdomain(int tag);

    Subdomain(int tag, 
	      TaggedObjectStorage &theInternalNodeStorage,
	      TaggedObjectStorage &theExternalNodeStorage,
	      TaggedObjectStorage &theElementsStorage,
	      TaggedObjectStorage &theLoadPatternsStorage,	      
	      TaggedObjectStorage &theMPsStorage,
	      TaggedObjectStorage &theSPsStorage);
    
    virtual  ~Subdomain();    

    // method added for parallel domain generation
    //    virtual int buildSubdomain(int numSubdomains, 
    //			       PartitionedModelBuilder &theBuilder); 

    // Domain methods which must be rewritten
    virtual void clearAll(void);
    virtual bool addNode(Node *);	
    virtual Node *removeNode(int tag);        
    virtual NodeIter &getNodes(void);    
    virtual Node *getNode(int tag);            
    virtual Node **getNodePtrs(void);            

    virtual bool hasNode(int tag);
    virtual bool hasElement(int tag);

    virtual int getNumNodes(void) const;    
    virtual int commit(void);
    virtual int revertToLastCommit(void);    
    virtual int revertToStart(void);        
    virtual int update(void);
    virtual int update(double newTime, double dT);

//#ifdef _PARALLEL_PROCESSING
    virtual  int barrierCheckIN(void) {return 0;};
    virtual  int barrierCheckOUT(int) {return 0;};
//#endif
   
    virtual  void Print(OPS_Stream &s, int flag =0);
    virtual void Print(OPS_Stream &s, ID *nodeTags, ID *eleTags, int flag =0);
    
    // Domain type methods unique to a Subdomain
    virtual NodeIter &getInternalNodeIter(void);
    virtual NodeIter &getExternalNodeIter(void);
    virtual bool addExternalNode(Node *);

    virtual void wipeAnalysis(void);
    virtual void setDomainDecompAnalysis(DomainDecompositionAnalysis &theAnalysis);
    virtual int setAnalysisAlgorithm(EquiSolnAlgo &theAlgorithm);
    virtual int setAnalysisIntegrator(IncrementalIntegrator &theIntegrator);
    virtual int setAnalysisLinearSOE(LinearSOE &theSOE);
    virtual int setAnalysisEigenSOE(EigenSOE &theSOE);
    virtual int setAnalysisConvergenceTest(ConvergenceTest &theTest);
    virtual int invokeChangeOnAnalysis(void);
    
    // Element methods which must be written
    virtual int getNumExternalNodes(void) const;    
    virtual const ID &getExternalNodes(void);
    virtual int getNumDOF(void);

    virtual int commitState(void);    
    
    virtual const Matrix &getTangentStiff(void);
    virtual const Matrix &getInitialStiff(void);    
    virtual const Matrix &getDamp(void);    
    virtual const Matrix &getMass(void);    

    virtual void  zeroLoad(void);
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor);
    virtual int addInertiaLoadToUnbalance(const Vector &accel);

    virtual const Vector &getResistingForce(void);    
    virtual const Vector &getResistingForceIncInertia(void);        
    virtual bool isSubdomain(void);    
    virtual int setRayleighDampingFactors(double alphaM, 
					  double betaK, 
					  double betaK0, 
					  double betaKc);

    virtual  int  updateParameter(int tag, int value);
    virtual  int  updateParameter(int tag, double value);    

    // Element type methods unique to a subdomain
    virtual int computeTang(void);
    virtual int computeResidual(void);
    virtual const Matrix &getTang(void);    

    void setFE_ElementPtr(FE_Element *theFE_Ele);
    virtual const Vector &getLastExternalSysResponse(void);
    virtual int computeNodalResponse(void);    
    virtual int analysisStep(double deltaT);
    virtual int eigenAnalysis(int numMode, bool generalized, bool findSmallest);
    virtual bool doesIndependentAnalysis(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    virtual double getCost(void);
    virtual int addResistingForceToNodalReaction(bool inclInertia);
    
  protected:    
    virtual int buildMap(void);
    bool mapBuilt;
    ID *map;
    Vector *mappedVect;
    Matrix *mappedMatrix;


    FE_Element *getFE_ElementPtr(void);
    TaggedObjectStorage  *internalNodes;
    TaggedObjectStorage  *externalNodes;    

    DomainDecompositionAnalysis *getDDAnalysis(void);

  private:
    double realCost;
    double cpuCost;
    int pageCost;
    DomainDecompositionAnalysis *theAnalysis;
    ID *extNodes;
    FE_Element *theFEele;
    
    //    TaggedObjectStorage  *realExternalNodes;        

    SingleDomNodIter   *internalNodeIter;
    SingleDomNodIter   *externalNodeIter;    
    SubdomainNodIter   *theNodIter;

    //    PartitionedModelBuilder *thePartitionedModelBuilder;
    static Matrix badResult;
};

#endif


