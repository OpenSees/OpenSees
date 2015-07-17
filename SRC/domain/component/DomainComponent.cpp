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
                                                                        
// $Revision: 1.5 $
// $Date: 2006-09-05 19:47:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/DomainComponent.cpp,v $
                                                                        
                                                                        
// File: ~/domain/component//DomainComponent.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class implementation for DomainComponent
//
// What: "@(#) DomainComponent.h, revA"


#include <Domain.h>
#include <DomainComponent.h>

DomainComponent::DomainComponent(int tag, int clasTag)
  :TaggedObject(tag), MovableObject(clasTag), theDomain(0)
{
    // does nothing else
}


DomainComponent::~DomainComponent()
{
    // does nothing
}


void
DomainComponent::setDomain(Domain *model)
{
    // sets the pointer 
    theDomain = model;
}


Domain *
DomainComponent::getDomain(void) const
{
    // returns the current pointer
    return theDomain;
}

int 
DomainComponent::displaySelf(Renderer &theViewer, int mode, float fact, const char **displayModes, int numModes)
{
  return 0;
}
