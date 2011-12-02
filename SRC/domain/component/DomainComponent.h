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
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/DomainComponent.h,v $
                                                                        
                                                                        
#ifndef DomainComponent_h
#define DomainComponent_h

// File: ~/domain/component/DomainComponent.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for DomainComponent.
// The DomainComponent class is an abstract class, example subclasses include
// Element, Node, SP_Constraint, MP_Constraint, NodalLoad, ElementalLoad. 
// Each of these objects forms part of Domain and has methods to set and obtain
// the associated Domain object.
//
// What: "@(#) DomainComponent.h, revA"

#include <TaggedObject.h>
#include <MovableObject.h>

class Domain;
class Renderer;
class Information;

class DomainComponent: public TaggedObject, public MovableObject
{
  public:
    virtual ~DomainComponent();

    virtual void setDomain(Domain *theDomain);
    virtual Domain *getDomain(void) const;

    // Method for visualisation, default does nothing
    virtual int displaySelf(Renderer &, int displayMode, float fact);

    // methods for sensitivity studies
    virtual int setParameter(char **argv, int argc, Information &eleInformation);
    virtual int updateParameter(int responseID, Information &eleInformation);	
    
  protected:
    DomainComponent(int tag, int classTag);
    
  private:    
    Domain *theDomain; // a pointer to the enclosing Domain object
};

#endif


