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
// $Date: 2003-02-14 23:01:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ElementalLoad.cpp,v $
                                                                        
                                                                        
// Written: fmk 11/95
//          modified 11/01 for new design

// Purpose: This file contains the methods for class ElementalLoad.

#include "ElementalLoad.h"
#include <Element.h>
#include <Domain.h>

ElementalLoad::ElementalLoad(int tag, int cTag, const ID &theEleTags)
  :Load(tag, cTag), theElementTags(0), theElements(0), numElements(0)
{
  theElementTags = new ID(theEleTags);
  if (theElementTags != 0) 
    numElements = theElementTags->Size();
}

ElementalLoad::ElementalLoad(int tag, int cTag)
  :Load(tag, cTag), theElementTags(0), theElements(0), numElements(0)
{

}

// provided for the FEM_Object broker; the tag and elementTag need
// to be supplied in recvSelf();
ElementalLoad::ElementalLoad(int cTag)
:Load(0, cTag), theElementTags(0), theElements(0), numElements(0)
{

}

ElementalLoad::~ElementalLoad()
{
  if (theElementTags != 0)
    delete theElementTags;

  if (theElements != 0)
    delete [] theElements;
}


void
ElementalLoad::setDomain(Domain *theDomain)
{
  this->DomainComponent::setDomain(theDomain);

  int size = 0;
  if (theElementTags != 0) {
    size = theElementTags->Size();
    if (theElements == 0)
      delete [] theElements;
    theElements = new Element *[size];
    
    for (int i=0; i<size; i++) {
      theElements[i] = theDomain->getElement((*theElementTags)(i));
      if (theElements[i] == 0) {
	opserr << "WARNING - ElementalLoad::setDomain - ele with tag ";
	opserr << (*theElementTags)(i) << " does not exist in the domain\n";
      }
    }
  }
}

void 
ElementalLoad::applyLoad(double loadFactor) 
{
  for (int i=0; i<numElements; i++)
    if (theElements[i] != 0)
      theElements[i]->addLoad(this, loadFactor);
}


const ID &
ElementalLoad::getElementTags(void) 
{
   return *theElementTags;
}


int 
ElementalLoad::removeElement(int tag) 
{
  int loc = theElementTags->getLocation(tag);
  if (loc >= 0) {
    (*theElementTags)(loc) = -1;
    for (int i=loc; i<numElements-1; i++)
      theElements[i] = theElements[i+1];
    numElements--;
  }

  return numElements;
}

int 
ElementalLoad::setElementTags(const ID &theEleTags)
{
  if (theElementTags != 0) {
    delete theElementTags;
    numElements = 0;
  }

  theElementTags = new ID(theEleTags);
  if (theElementTags != 0) 
    numElements = theElementTags->Size();

  return 0;
}
  

