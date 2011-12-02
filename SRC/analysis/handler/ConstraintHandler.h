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
// $Date: 2005-08-31 17:28:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/ConstraintHandler.h,v $
                                                                        
                                                                        
#ifndef ConstraintHandler_h
#define ConstraintHandler_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for ConstraintHandler.
// ConstraintHandler is an abstract base class, i.e. no objects of it's
// type can be created. ConstraintHandlers enforce the single and multi point 
// constraints that exist in the domain by creating the appropriate FE_Element
// and DOF_Group objects.
//
// What: "@(#) ConstraintHandler.h, revA"

#include <MovableObject.h>

class AnalysisMethod;
class ID;
class Domain;
class AnalysisModel;
class Integrator;
class FEM_ObjectBroker;

class ConstraintHandler : public MovableObject
{
  public:
    ConstraintHandler(int classTag);
    virtual ~ConstraintHandler();

    void setLinks(Domain &theDomain, 
		  AnalysisModel &theModel,
		  Integrator &theIntegrator);

    // pure virtual functions
    virtual int handle(const ID *nodesNumberedLast =0) =0;
    virtual int update(void);
    virtual int applyLoad(void);
    virtual int doneNumberingDOF(void);
    virtual void clearAll(void) =0;    

  protected:
    Domain *getDomainPtr(void) const;
    AnalysisModel *getAnalysisModelPtr(void) const;
    Integrator *getIntegratorPtr(void) const;
    
  private:
    Domain *theDomainPtr;
    AnalysisModel *theAnalysisModelPtr;
    Integrator *theIntegratorPtr;
};

#endif

