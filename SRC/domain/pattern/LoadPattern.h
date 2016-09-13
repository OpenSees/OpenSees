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
// $Date: 2008-02-29 20:45:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/LoadPattern.h,v $
                                                                        
                                                                        
#ifndef LoadPattern_h
#define LoadPattern_h

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for LoadPattern.
// LoadPattern is a concrete class. A LoadPattern object is used to 
// to store reference loads and single point constraints and a TimeSeries function
// which is used to determine the load factor given the pseudo-time
// to the model. 
//
// What: "@(#) LoadPattern.h, revA"

#include <DomainComponent.h>
#include <Vector.h>

class NodalLoad;
class TimeSeries;
class ElementalLoad;
class SP_Constraint;
class NodalLoadIter;
class ElementalLoadIter;
class SingleDomSP_Iter;
class SP_ConstraintIter;
class TaggedObjectStorage;
class GroundMotion;

class LoadPattern : public DomainComponent    
{
  public:
    // constructors
    LoadPattern(int tag, double fact = 1.0);    
    LoadPattern();                                         // for FEM_ObjectBroker
    LoadPattern(int tag, int classTag, double fact = 1.0); // for subclasses
    
    // destructor
    virtual ~LoadPattern();

    // method to set the associated TimeSeries and Domain
    virtual void setTimeSeries(TimeSeries *theSeries);
    virtual void setDomain(Domain *theDomain);

    // methods to add loads
    virtual bool addSP_Constraint(SP_Constraint *);
    virtual bool addNodalLoad(NodalLoad *);
    virtual bool addElementalLoad(ElementalLoad *);
    virtual NodalLoadIter     &getNodalLoads(void);
    virtual ElementalLoadIter &getElementalLoads(void);    
    virtual SP_ConstraintIter &getSPs(void);        
    
    // methods to remove loads
    virtual void clearAll(void);
    virtual NodalLoad *removeNodalLoad(int tag);
    virtual ElementalLoad *removeElementalLoad(int tag);
    virtual SP_Constraint *removeSP_Constraint(int tag);

    // methods to apply loads
    virtual void applyLoad(double pseudoTime = 0.0);
    virtual void setLoadConstant(void);
	virtual void unsetLoadConstant(void);
    virtual double getLoadFactor(void);

    // methods for o/p
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    virtual void Print(OPS_Stream &s, int flag =0);        

    // method to obtain a blank copy of the LoadPattern
    virtual LoadPattern *getCopy(void);

    virtual int addMotion(GroundMotion &theMotion, int tag);    
    virtual GroundMotion *getMotion(int tag);        

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual void applyLoadSensitivity(double pseudoTime = 0.0);
    virtual int setParameter(const char **argv, int argc, Parameter &param);
    virtual int  updateParameter(int parameterID, Information &info);
    virtual int  activateParameter(int parameterID);
    virtual const Vector & getExternalForceSensitivity(int gradNumber);

    virtual int saveLoadFactorSensitivity(double dlambdadh, int gradIndex, int numGrads);
    virtual double getLoadFactorSensitivity(int gradIndex);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    int    isConstant;     // to indictae whether setConstant has been called
	
  private:
    double loadFactor;     // current load factor
    double scaleFactor;    // factor to scale load factor from time series

    TimeSeries *theSeries; // pointer to associated TimeSeries

    int	   currentGeoTag;
    int    lastGeoSendTag;
    int    dbSPs, dbNod, dbEle; // database tags for storing info about components
    
    // storage objects for the loads and constraints
    TaggedObjectStorage  *theNodalLoads;
    TaggedObjectStorage  *theElementalLoads;
    TaggedObjectStorage  *theSPs; 	  

    // iterator objects for the objects added to the storage objects
    NodalLoadIter       *theNodIter;
    ElementalLoadIter   *theEleIter;
    SingleDomSP_Iter    *theSpIter;    

    // AddingSensitivity:BEGIN //////////////////////////////////////
    Vector *randomLoads;
    bool RVisRandomProcessDiscretizer;
    Vector *dLambdadh;
    // AddingSensitivity:END ////////////////////////////////////////

    int lastChannel; 
};

#endif







