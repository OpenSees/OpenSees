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
                                                                        
// $Revision: 1.8 $
// $Date: 2003-02-25 23:34:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/FilePlotter.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/99
//
// Description: This file contains the class implementation for FilePlotter
// FilePlotter is a class for building a Plane Frame model in an interpreted
// enviroment. The constructor is used to add new commands to the interpreter,
// these commands are also defined in this file.
//
// What: "@(#) FilePlotter.C, revA"


#include <stdlib.h>
#include <string.h>

#include <fstream>
using std::ifstream;

#include <iomanip>
using std::ios;

#include <ctype.h>
#include <ID.h>

#ifdef _WGL
#include <OpenGLRenderer.h>
#elif _GLX
#include <OpenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif

#include <PlainMap.h>
#include <Vector.h>

#include "FilePlotter.h"

FilePlotter::FilePlotter(const char *_fileName, 
			 const char *windowTitle, 
			 int xLoc, int yLoc, int width, int height, double dT)
  :theMap(0), theRenderer(0), cols(0), deltaT(dT), nextTimeStampToRecord(0.0)
{

  // create the window in which we plot on the screen
  theMap = new PlainMap();
#ifdef _WGL
  theRenderer = new OpenGLRenderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#elif _GLX
  theRenderer = new OpenGLRenderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#else
  theRenderer = new X11Renderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#endif
  // set up for standard xy-plot - rendering in the xy plane at z =0
  theRenderer->setVRP(0.0, 0.0, 0.0); 
  theRenderer->setVPN(0.0, 0.0, 1.0);
  theRenderer->setVUP(0.0, 1.0, 0.0);
  theRenderer->setFillMode("wire");             // wire mode
  theRenderer->setProjectionMode("parallel");  // wire mode
  theRenderer->setPlaneDist(1.0, -1.0);
  theRenderer->setPRP(0.0, 0.0, 10.0);
  theRenderer->setPortWindow(-1.0, 1.0, -1.0, 1.0);  // use the whole window

  // copy the file name
  fileName = new char[strlen(_fileName)+1];
  if (fileName == 0) {
    opserr << "FilePlotter::FilePlotter -- out of memory copying fileName " << endln;
    exit(-1);
  }
  strcpy(fileName, _fileName);    

}

FilePlotter::~FilePlotter()
{
  // may possibly invoke Tcl_DeleteCommand() later
  // for moment just invoke destructor on map and renderer
  // and set pointer to NULL
  delete theMap;
  delete theRenderer;

  if (cols != 0)
    delete cols;

  if (fileName != 0)
    delete [] fileName;
}
    
int 
FilePlotter::record(int cTag, double timeStamp)
{

  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    return this->plotFile();

  } else
    return 0;
}

int 
FilePlotter::playback(int cTag)
{
  //  this->plotFile();
  return 0;
}

void
FilePlotter::restart(void)
{
    
}

int
FilePlotter::setFile(char *newFileName)
{
  strcpy(fileName, newFileName);    
  return 0;
}

int
FilePlotter::plotFile(void)
{
    /*
     * first pass: 1) open file
     *             2) determine number of entries on first line ASSSSSUME same through file
     *             3) determine number of lines and bounds [xmin, xmax, ymin, ymax]
     *             4) close the file
     */
  
     // open file
    ifstream theFile; 
    theFile.open(fileName, ios::in);

    if (theFile.bad()) {
	opserr << "WARNING - FilePlotter::FilePlotter()";
	opserr << " - could not open file " << fileName << endln;
	return -1;
    }    

    double xMin, xMax, yMin, yMax;
    xMin = 0; xMax = 0; yMin =0; yMax =0;
    double xValue, yValue;

    // determine number of elements in each line 
    //  NOTE ASSUMES ALL LINES HAVE THE SAME NUMBER OF ELEMENTS
    char c;
    int numLineEntries = 0;
    int numLines = 0;

    while (theFile.get(c) && (c!= EOF) && c != '\n') {
      if (!isspace(c)) {
	theFile.putback(c);

	theFile >> xValue;
	for (int i=0; i<cols->Size(); i++) {
	  if (i%2 == 0) { // an xValue if numLineEntries == colX
	    if (numLineEntries == (*cols)(i)) {
	      if (xValue < xMin) xMin = xValue;
	      if (xValue > xMax) xMax = xValue;
	    }
	  } else { // a y value if (numLineEntries == colY) {
	    if (numLineEntries == (*cols)(i)) {
	      if (xValue < yMin) yMin = xValue;
	      if (xValue > yMax) yMax = xValue;	  
	    } 
	  }
	}
	

	numLineEntries ++;
      }
      numLines =1;
    }

    // check colX and colY for valid entries
    if (numLines > 0) {
      if (cols == 0) {
	opserr << "FilePLotter::plotFile() - no valid columns have been set\n";

      } else {

	// parse through file checking the bounds
	Vector data(numLineEntries);
	while (theFile >> data(0)) {
	  for (int j=1; j<numLineEntries; j++)
	    theFile >> data(j);

	  for (int i=0; i<cols->Size(); i += 2) {
	    xValue = data((*cols)(i));
	    yValue = data((*cols)(i+1));
	    if (xValue < xMin) xMin = xValue;
	    if (xValue > xMax) xMax = xValue;
	    if (yValue < yMin) yMin = yValue;
	    if (yValue > yMax) yMax = yValue;
	    numLines++;
	  }    
	}

	// set the window bounds NOTE small border around the edges
	double xBnd = (xMax-xMin)/8;
	double yBnd = (yMax-yMin)/8;
    
	theRenderer->setViewWindow(xMin-xBnd,xMax+xBnd,yMin-yBnd,yMax+yBnd);
      }
    }

    // close the file
    theFile.close();

    /*
     * second pass: 1) open file
     *              2) get the renderer ready for drawing virgin image and then draw the x and y axis
     *              3) parse throgh the file and connect the dots
     *              4) close the file
     */

    if (numLines > 1) {

      static Vector pt1(3); 
      static Vector pt2(3);
      static Vector rgb(3);

      // clear the present image and get renderer ready to process data
      theRenderer->clearImage();
      theRenderer->startImage();

      // draw the x axis
      pt1(0) = xMin; pt2(0) = xMax;
      pt1(1) = 0.0;  pt2(1) = 0.0;
      theRenderer->drawLine(pt1, pt2, rgb, rgb);    

      static char theText[20];
      if (xMin != 0.0 && -100*xMin > xMax) {
	sprintf(theText,"%.2e",xMin);
	theRenderer->drawText(pt1, theText, strlen(theText), 'l', 'b');
      }
      if (xMax != 0.0) {
	sprintf(theText,"%.2e",xMax);
	theRenderer->drawText(pt2, theText, strlen(theText), 'r', 'b');
      }

      // draw the y axis
      pt1(0) = 0.0; pt2(0) = 0.0;
      pt1(1) = yMin;  pt2(1) = yMax;
      theRenderer->drawLine(pt1, pt2, rgb, rgb);        

      if (yMin != 0.0 && -100 *yMin > yMax) {
	sprintf(theText,"%.2e",yMin);
	theRenderer->drawText(pt1, theText, strlen(theText), 'c', 't');
      }
      if (yMax != 0.0) {
	sprintf(theText,"%.2e",yMax);
	theRenderer->drawText(pt2, theText, strlen(theText), 'c', 'b');
      }

      Vector data1(numLineEntries);
      Vector data2(numLineEntries);

      // open the file again, read through and connect the dots
      ifstream theFile1; 
      theFile1.open(fileName, ios::in);
      if (theFile1.bad()) {
	opserr << "WARNING - FilePlotter::FilePlotter()";
	opserr << " - could not open file " << fileName << endln;
	return -1;
      }    

      for (int ii=0; ii< numLineEntries; ii++)
	theFile1 >> data1(ii);

      for (int i=1; i<numLines; i++) {
	// read the data
	for (int ii=0; ii< numLineEntries; ii++)
	  theFile1 >> data2(ii);

	// plot the lines
	for (int j=0; j<cols->Size(); j+=2) {
	  pt1(0) = data1((*cols)(j)); 
	  pt1(1) = data1((*cols)(j+1));
	  pt2(0) = data2((*cols)(j)); 
	  pt2(1) = data2((*cols)(j+1));
	  theRenderer->drawLine(pt1, pt2, rgb, rgb);
	}
	
	data1 = data2;

      }
      
      theRenderer->doneImage();

      // close the file
      theFile1.close();

    }
    return 0;
}


int
FilePlotter::setCol(const ID &theCols)
{
  if (theCols.Size()%2 != 0) {
    opserr << "FilePlotter::setCol() - the size of the cols ID " << theCols.Size() << " is not a multiple of 2\n";
    return -1;
  }

  for (int i=0; i<theCols.Size(); i++) {
    if (theCols(i) < 1) {
      opserr << "FilePlotter::FilePlotter() - a value of the cols " << theCols(i) << " is < 1\n";
      return -2;
    }
  }
  // check colX is valid, i.e. >= 1
  // if valid set colX using c indexing

  if (cols != 0) {
    if (cols->Size() != theCols.Size()) {
      delete cols;
      cols = 0;
    } else
      *cols = theCols;
  }

  if (cols == 0)
    cols = new ID(theCols);

  for (int j=0; j<cols->Size(); j++)
    (*cols)(j) -= 1;
    
  return 0;
}







