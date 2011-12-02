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
// $Date: 2000-09-15 08:23:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Renderer.cpp,v $
                                                                        
                                                                        
// File: ~/renderer/Renderer.C
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class interface for Renderer.
// Renderer is an abstract base class. An Renderer object is used
// to create an image of the domain.
//
// What: "@(#) Renderer.h, revA"

#include <Renderer.h>
#include <ColorMap.h>
#include <Domain.h>


Renderer::Renderer(ColorMap &_theMap)
  :theMap(&_theMap)
{

}

Renderer::~Renderer()
{

}

/*
int 
Renderer::displayModel(int eleFlag, int nodeFlag, 
		       float fact)
{
  // loop over the elements getting each to display itself
  // using this and displayTag as arguments.
  // first clear the image
  this->startImage();
  int res = 0;

  if (eleFlag >= 0) {
      ElementIter &theElements = theDomain->getElements();
      Element *theEle;
      while ((theEle = theElements()) != 0) {
	  res = theEle->displaySelf(*this, eleFlag, fact);
	  if (res < 0) {
	      cerr << "Renderer::displayModel() - Element: ";
	      cerr << theEle->getTag() << " failed to display itself\n";
	  }
      }
  }
  
  if (nodeFlag >= 0) {
      NodeIter &theNodes = theDomain->getNodes();
      Node *theNode;
      while ((theNode = theNodes()) != 0) {
	  res = theNode->displaySelf(*this, nodeFlag, fact);
	  if (res < 0) {
	      cerr << "Renderer::displayModel() - Node: ";
	      cerr << theNode->getTag() << " failed to display itself\n";
	  }
      }
  }  

  // now mark the image has having been completed
  this->doneImage();
  return res;
}
*/

int
Renderer::drawVector(const Vector &position, const Vector &value)
{
    return 0;
}

int
Renderer::drawGText(const Vector &posGlobal, char *string, int length)
{
    return 0;
}

int
Renderer::drawLText(const Vector &posLocal, char *string, int length)
{
    return 0;
}

void
Renderer::setColorMap(ColorMap &map)
{
    theMap = &map;
}














