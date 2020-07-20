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
                                                                        
// $Revision: 1.9 $
// $Date: 2008-02-15 23:37:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/AlgorithmIncrements.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class implementation for AlgorithmIncrements
// AlgorithmIncrements is a class for building a Plane Frame model in an interpreted
// enviroment. The constructor is used to add new commands to the interpreter,
// these commands are also defined in this file.
//
// What: "@(#) AlgorithmIncrements.C, revA"


#include <AlgorithmIncrements.h>
#include "equiSolnAlgo/EquiSolnAlgo.h"
#include <LinearSOE.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <ID.h>
#include <Vector.h>

#include <fstream>
using std::ifstream;

#include <iomanip>
using std::ios;

#ifdef _WGL
#include <OpenGLRenderer.h>
#elif _GLX
#include <OpenGLRenderer.h>
#elif _AGL
#include <OpenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif
#include <PlainMap.h>
#include <Vector.h>



AlgorithmIncrements::AlgorithmIncrements(EquiSolnAlgo *theEquiAlgo,
					 const char *windowTitle, 
					 int xLoc, int yLoc, int width, int height,
					 bool displayOnRecord, const char *theFileName)
  :Recorder(RECORDER_TAGS_AlgorithmIncrements),
   theMap(0), theRenderer(0), numRecord(0), displayRecord(displayOnRecord), fileName(0)
{
  
  theAlgo = theEquiAlgo;

  // create the window in which we plot on the screen
  theMap = new PlainMap();
#ifdef _WGL
  theRenderer = new OpenGLRenderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#elif _GLX
  theRenderer = new OpenGLRenderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#elif _AGL
  theRenderer = new OpenGLRenderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#else
  theRenderer = new X11Renderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#endif
  // set up for standard xy-plot - rendering in the xy plane at z =0
  theRenderer->setVRP(0.0, 0.0, 0.0); 
  theRenderer->setVPN(0.0, 0.0, 1.0);
  theRenderer->setVUP(0.0, 1.0, 0.0);
  theRenderer->setFillMode("wire");             // wire mode
  theRenderer->setPlaneDist(1.0, 100.0);
  theRenderer->setPRP(0.0, 0.0, 10.0);

  if (theFileName != 0) {
    fileName = new char[strlen(theFileName) + 1];
    strcpy(fileName, theFileName);        
  }
}

AlgorithmIncrements::~AlgorithmIncrements()
{
  // may possibly invoke Tcl_DeleteCommand() later
  // for moment just invoke destructor on map and renderer
  // and set pointer to NULL
  delete theMap;
  delete theRenderer;
  
  if (theFile)
    theFile.close();
  if (fileName != 0) 
    delete [] fileName;
    
}
    
int 
AlgorithmIncrements::record(int cTag, double timeStamp)
{
  LinearSOE *theSOE = theAlgo->getLinearSOEptr();
  if (theSOE == 0)
    return 0;


  const Vector &B = theSOE->getB();
  const Vector &X = theSOE->getX();

  if (fileName != 0) {
    if (cTag == 0) {

      if (theFile)
	theFile.close();

      numRecord = 0;

      theFile.open(fileName, ios::out);
      if (!theFile) {
	opserr << "WARNING - AlgorithmIncrements::record()";
	opserr << " - could not open file " << fileName << endln;
      } 
    }
    if (theFile) {
      numRecord ++;
      int i;
      for (i=0; X.Size(); i++) theFile << X(i);
      for (i=0; X.Size(); i++) theFile << B(i);
    }
  }

  if (displayRecord == true)
    return this->plotData(X, B);

  return 0;
}

int 
AlgorithmIncrements::playback(int cTag)
{

  if (fileName != 0) {

    LinearSOE *theSOE = theAlgo->getLinearSOEptr();
    if (theSOE == 0)
      return 0;

    Vector X(theSOE->getNumEqn());
    Vector B(theSOE->getNumEqn());

    if (theFile) {
      theFile.close();
    }

    ifstream aFile;
    aFile.open(fileName, ios::in);
    if (!aFile) {
      opserr << "WARNING - AlgorithmIncrements::playback()";
      opserr << " - could not open file " << fileName << endln;
    }
    for (int i = 0; i<numRecord; i++) {
      int ii;
      for (ii=0; X.Size(); ii++) theFile << X(ii);
      for (ii=0; X.Size(); ii++) theFile << B(ii);

      this->plotData(X,B);
      //      char c = getchar();
    }
  }

  return 0;
}

int
AlgorithmIncrements::restart(void)
{
  return 0;
}

int
AlgorithmIncrements::plotData(const Vector &X, const Vector &B)
{

  int size = X.Size();
  if (size < 2)
    return 0;

  /*
   * first pass plot the unbalance (B in the SOE):
   *             1) set the port window to be upper half
   *             2) determine bounds
   *             3) draw axis and bounds values
   *             4) draw points
   */

  theRenderer->clearImage();
  theRenderer->setPortWindow(-1.0, 1.0, 0.0, 1.0);  // use the whole window

   double xMin, xMax, yMin, yMax;
   xMin = 0.0;
   xMax = B.Size();
   yMin = 0.0;
   yMax = 0.0;
	int i;

   for (i=0; i<size; i++) {
     double value = B(i);
     if (value < yMin) yMin = value;
     if (value > yMax) yMax = value;	  
   } 

   if (-yMin > yMax) 
     yMax = -yMin;
   else
     yMin = -yMax;


   // set the window bounds NOTE small border around the edges
   double xBnd = (xMax-xMin)/10;
   double yBnd = (yMax-yMin)/10;
    
   theRenderer->setViewWindow(xMin-xBnd,xMax+xBnd,yMin-yBnd,yMax+yBnd);
   theRenderer->startImage();

   static Vector pt1(3); 
   static Vector pt2(3);
   pt1(2) = 0.0;  pt2(2) = 0.0;

   // draw the x axis
   pt1(0) = xMin; pt2(0) = xMax;
   pt1(1) = 0.0;  pt2(1) = 0.0;
   theRenderer->drawLine(pt1, pt2, 0.0, 0.0);    

   // draw the y axis
   pt1(0) = 0.0; pt2(0) = 0.0;
   pt1(1) = yMin;  pt2(1) = yMax;
   theRenderer->drawLine(pt1, pt2, 0.0, 0.0);        


   static char theText[20];
   if (yMin != 0.0 && -100 *yMin > yMax) {
     sprintf(theText,"%.2e",yMin);
     theRenderer->drawText(pt1, theText, strlen(theText));
   }
   if (yMax != 0.0) {
     sprintf(theText,"%.2e",yMax);
     theRenderer->drawText(pt2, theText, strlen(theText));
   }

   pt1(0) = 1;
   pt1(1) = B(0);

   for (i=1; i<size; i++) {
     pt2(0) = i+1; 
     pt2(1) = B(i);
     theRenderer->drawLine(pt1, pt2, 1.0, 1.0);
     pt1(0) = pt2(0);
     pt1(1) = pt2(1);
   }

   theRenderer->doneImage();

  /*
   * first pass plot the unbalance (B in the SOE):
   *             1) set the port window to be upper half
   *             2) determine bounds
   *             3) draw axis and bounds values
   *             4) draw points
   */

   //  theRenderer->clearImage();
  theRenderer->setPortWindow(-1.0, 1.0, -1.0, 0.0);  // use the whole window

   xMin = 0.0;
   xMax = X.Size();
   yMin = 0.0;
   yMax = 0.0;

   for (i=0; i<size; i++) {
     double value = X(i);
     if (value < yMin) yMin = value;
     if (value > yMax) yMax = value;	  
   } 

   if (-yMin > yMax) 
     yMax = -yMin;
   else
     yMin = -yMax;


   // set the window bounds NOTE small border around the edges
   xBnd = (xMax-xMin)/10;
   yBnd = (yMax-yMin)/10;
    
   theRenderer->setViewWindow(xMin-xBnd,xMax+xBnd,yMin-yBnd,yMax+yBnd);
   theRenderer->startImage();

   pt1(2) = 0.0;  pt2(2) = 0.0;

   // draw the x axis
   pt1(0) = xMin; pt2(0) = xMax;
   pt1(1) = 0.0;  pt2(1) = 0.0;
   theRenderer->drawLine(pt1, pt2, 0.0, 0.0);    

   // draw the y axis
   pt1(0) = 0.0; pt2(0) = 0.0;
   pt1(1) = yMin;  pt2(1) = yMax;
   theRenderer->drawLine(pt1, pt2, 0.0, 0.0);        

   if (yMin != 0.0 && -100 *yMin > yMax) {
     sprintf(theText,"%.2e",yMin);
     theRenderer->drawText(pt1, theText, strlen(theText));
   }
   if (yMax != 0.0) {
     sprintf(theText,"%.2e",yMax);
     theRenderer->drawText(pt2, theText, strlen(theText));
   }

   pt1(0) = 1;
   pt1(1) = X(0);

   for (i=1; i<size; i++) {
     pt2(0) = i+1; 
     pt2(1) = X(i);
     theRenderer->drawLine(pt1, pt2, 1.0, 1.0);
     pt1(0) = pt2(0);
     pt1(1) = pt2(1);
   }

   theRenderer->doneImage();
   
    return 0;
}




