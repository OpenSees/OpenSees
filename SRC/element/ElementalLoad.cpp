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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-08-26 16:52:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ElementalLoad.cpp,v $
                                                                        
                                                                        
// Written: fmk 11/95
//          modified 11/01 for new design

// Purpose: This file contains the methods for class ElementalLoad.

#include <ElementalLoad.h>
#include <Element.h>
#include <Domain.h>

ElementalLoad::ElementalLoad(int tag, int cTag, int theEleTag)
  :Load(tag, cTag), eleTag(theEleTag), theElement(0)
{

}

ElementalLoad::ElementalLoad(int tag, int cTag)
  :Load(tag, cTag), eleTag(0), theElement(0)
{

}


// provided for the FEM_Object broker; the tag and elementTag need
// to be supplied in recvSelf();
ElementalLoad::ElementalLoad(int cTag)
:Load(0, cTag), eleTag(0), theElement(0)
{

}


ElementalLoad::~ElementalLoad()
{

}


void
ElementalLoad::setDomain(Domain *theDomain)
{
  this->DomainComponent::setDomain(theDomain);

  if (theDomain == 0) {
    theElement = 0;
    return;
  }

  theElement = theDomain->getElement(eleTag);
  if (theElement == 0) {
    opserr << "WARNING - ElementalLoad::setDomain - no ele with tag ";
    opserr << eleTag << " exists in the domain\n";
  }
}

void 
ElementalLoad::applyLoad(double loadFactor) 
{
  if (theElement != 0)
    theElement->addLoad(this, loadFactor);
}

void 
ElementalLoad::applyLoad(const Vector &loadFactors) 
{
  if (theElement != 0)
    theElement->addLoad(this, loadFactors);
}

const Vector&
ElementalLoad::getSensitivityData(int gradIndex)
{
  static Vector trash(10);

  return trash;
}

int
ElementalLoad::getElementTag(void) 
{
   return eleTag;
}



