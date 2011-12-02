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
// $Date: 2000-09-15 08:23:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/FilePlotter.cpp,v $
                                                                        
                                                                        
// File: ~/recorder/FilePlotter
// 
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class implementation for FilePlotter
// FilePlotter is a class for building a Plane Frame model in an interpreted
// enviroment. The constructor is used to add new commands to the interpreter,
// these commands are also defined in this file.
//
// What: "@(#) FilePlotter.C, revA"


#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <ctype.h>
#ifdef _WIN32
#include <OpenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif
#include <PlainMap.h>
#include <Vector.h>

#include "FilePlotter.h"

FilePlotter::FilePlotter(char *_fileName, 
			 char *windowTitle, 
			 int xLoc, int yLoc, int width, int height)
  :theMap(0), theRenderer(0), colX(0), colY(1)
{

  // create the window in which we plot on the screen
  theMap = new PlainMap();
#ifdef _WIN32
  theRenderer = new OpenGLRenderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#else
  theRenderer = new X11Renderer(windowTitle, xLoc, yLoc, width, height, *theMap);
#endif
  // set up for standard xy-plot - rendering in the xy plane at z =0
  theRenderer->setVRP(0.0, 0.0, 0.0); 
  theRenderer->setVPN(0.0, 0.0, 1.0);
  theRenderer->setVUP(0.0, 1.0, 0.0);
  theRenderer->setFillMode(1);             // wire mode
  theRenderer->setPlaneDist(1.0, 100.0);
  theRenderer->setPRP(0.0, 0.0, 10.0);
  theRenderer->setPortWindow(-1.0, 1.0, -1.0, 1.0);  // use the whole window

  // copy the file name
  strcpy(fileName, _fileName);    
}

FilePlotter::~FilePlotter()
{
  // may possibly invoke Tcl_DeleteCommand() later
  // for moment just invoke destructor on map and renderer
  // and set pointer to NULL
  delete theMap;
  delete theRenderer;
}
    
int 
FilePlotter::record(int cTag)
{
  return this->plotFile();
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
    if (!theFile) {
	cerr << "WARNING - FileNodeDispRecorder::FileNodeDispRecorder()";
	cerr << " - could not open file " << fileName << endl;
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

	if (numLineEntries == colX) {
	  theFile >> xValue;
	  if (xValue < xMin) xMin = xValue;
	  if (xValue > xMax) xMax = xValue;
	} else if (numLineEntries == colY) {
	  theFile >> yValue;
	  if (yValue < yMin) yMin = yValue;
	  if (yValue > yMax) yMax = yValue;	  
	} else
	  theFile >> xValue;

	numLineEntries ++;
      }
      numLines =1;
    }

    // check colX and colY for valid entries
    if (numLines > 0) {
      if (colX >= numLineEntries || colY >= numLineEntries) {
	g3ErrorHandler->warning("FilePLotter::plotFile() colX %d or colY %d >= numLineEntries %d\n",
				colX, colY, numLineEntries);

      } else {

	// parse through file checking the bounds
	Vector data(numLineEntries);
	while (theFile >> data(0)) {
	  for (int i=1; i<numLineEntries; i++)
	    theFile >> data(i);
	  xValue = data(colX);
	  yValue = data(colY);
	  if (xValue < xMin) xMin = xValue;
	  if (xValue > xMax) xMax = xValue;
	  if (yValue < yMin) yMin = yValue;
	  if (yValue > yMax) yMax = yValue;
	  numLines++;
	}    

	// set the window bounds NOTE small border around the edges
	double xBnd = (xMax-xMin)/10;
	double yBnd = (yMax-yMin)/10;
    
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

      // open the file
      theFile.open(fileName, ios::in);
      if (!theFile) {
	cerr << "WARNING - FileNodeDispRecorder::FileNodeDispRecorder()";
	cerr << " - could not open file " << fileName << endl;
	return -1;
      }    

      Vector pt1(3); Vector pt2(3);
      pt1(2) = 0.0;  pt2(2) = 0.0;

      // clear the present image and get renderer ready to process data
      theRenderer->clearImage();
      theRenderer->startImage();

      // draw the x axis
      pt1(0) = xMin; pt2(0) = xMax;
      pt1(1) = 0.0;  pt2(1) = 0.0;
      theRenderer->drawLine(pt1, pt2, 0.0, 0.0);    

      // draw the y axis
      pt1(0) = 0.0; pt2(0) = 0.0;
      pt1(1) = yMin;  pt2(1) = yMax;
      theRenderer->drawLine(pt1, pt2, 0.0, 0.0);        

      Vector data(numLineEntries);
      theFile >> data;
      pt1(0) = data(colX); 
      pt1(1) = data(colY);

      for (int i=1; i<numLines; i++) {
	theFile >> data;
	pt2(0) = data(colX); 
	pt2(1) = data(colY);
	theRenderer->drawLine(pt1, pt2, 1.0, 1.0);
	pt1(0) = pt2(0);
	pt1(1) = pt2(1);
      }
      
      theRenderer->doneImage();

      // close the file
      theFile.close();

    }
    return 0;
}


int
FilePlotter::setCol(int _colX, int _colY)
{
  // check colX is valid, i.e. >= 1
  // if valid set colX using c indexing
  if (_colX > 0) 
    colX = _colX-1;
  else {
    g3ErrorHandler->warning("WARNING FilePlotter::setFile() colX %d must be >= 1\n",_colX);
    return -1;
  }

  // check colY is valid, i.e. >= 1
  // if valid set colY using c indexing
  if (_colY > 0) 
    colY = _colY-1;
  else {
    g3ErrorHandler->warning("WARNING FilePlotter::setFile() colY %d must be >= 1\n",_colX);
    return -1;
  }

  return 0;
}







