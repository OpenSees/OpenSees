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
// $Date: 2003-05-15 23:20:47 $
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


int        Renderer::numRenderers(0);
char     **Renderer::theTitles =0;
Renderer **Renderer::theRenderers =0;

Renderer::Renderer(ColorMap &_theMap)
  :theMap(&_theMap)
{

}


Renderer::Renderer(const char *title, ColorMap &_theMap)
  :theMap(&_theMap)
{
  int loc = -1;

  // look for an empty slot
  for (int i=0; i<numRenderers; i++)
    if (theRenderers[i] == 0) {
      loc = i;
      i = numRenderers;
    }

  // if no space or not already there add
  if (loc == -1) {
    Renderer **theNewRenderers = new Renderer *[numRenderers+1];
    char **theNewTitles = new char *[numRenderers+1];

    for (int i=0; i<numRenderers; i++) {
      theNewRenderers[i] = theRenderers[i];
      theNewTitles[i] = theTitles[i];
    }

    loc = numRenderers;
    numRenderers++;
    
    if (theRenderers != 0) 
      delete [] theRenderers;
    if (theTitles != 0)
      delete [] theTitles;

    theRenderers = theNewRenderers;
    theTitles = theNewTitles;
  }

  // set this in current slot
  theRenderers[loc] = this;
  char *titleCopy = new char [strlen(title)+1];
  strcpy(titleCopy, title);
  theTitles[loc] = titleCopy;
}

Renderer::~Renderer()
{
  for (int i=0; i<numRenderers; i++)
    if (theRenderers[i] == this) {
      theRenderers[i] = 0;
      delete [] (theTitles[i]);
      theTitles[i] = 0;
    }
}

int
Renderer::saveImage(const char *fileName)
{
  opserr << "Renderer::saveImage - no default implementation provided\n";
  return 0;
}


int
Renderer::saveImage(const char *rendererTitle, const char *fileName)
{
  for (int i=0; i<numRenderers; i++)
    if (theRenderers[i] != 0) 
      if (strcmp(rendererTitle, theTitles[i]) == 0)
	return theRenderers[i]->saveImage(fileName);

  opserr << "Renderer::saveImage - no renderer with title: " << rendererTitle << " found\n";
  return 0;
}

int
Renderer::drawVector(const Vector &position, const Vector &value, double factor)
{
    return 0;
}


void
Renderer::setColorMap(ColorMap &map)
{
    theMap = &map;
}














