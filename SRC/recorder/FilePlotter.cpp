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
                                                                        
// $Revision: 1.13 $
// $Date: 2008-02-15 23:37:35 $
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
#elif _AGL
#include <OpenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif

#include <PlainMap.h>
#include <Vector.h>

#include "FilePlotter.h"

FilePlotter::FilePlotter(const char *_fileName1, 
			 const char *windowTitle, 
			 int xLoc, int yLoc, int width, int height, double dT)
  :Recorder(RECORDER_TAGS_FilePlotter), 
   theMap(0), theRenderer(0), cols(0), deltaT(dT), nextTimeStampToRecord(0.0),
   data1a(0), data1b(0), data2a(0), data2b(0)
{

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
  theRenderer->setProjectionMode("parallel");  // wire mode
  theRenderer->setPlaneDist(1.0, -1.0);
  theRenderer->setPRP(0.0, 0.0, 10.0);
  theRenderer->setPortWindow(-1.0, 1.0, -1.0, 1.0);  // use the whole window

  // copy the file name
  fileName1 = new char[strlen(_fileName1)+1];
  if (fileName1 == 0) {
    opserr << "FilePlotter::FilePlotter -- out of memory copying fileName1 " << endln;
  }
  strcpy(fileName1, _fileName1);    

  fileName2 = 0;
}



FilePlotter::FilePlotter(const char *_fileName1, 
			 const char *_fileName2, 
			 const char *windowTitle, 
			 int xLoc, int yLoc, int width, int height, double dT)
  :Recorder(RECORDER_TAGS_FilePlotter), 
   theMap(0), theRenderer(0), cols(0), deltaT(dT), nextTimeStampToRecord(0.0),
   data1a(0), data1b(0), data2a(0), data2b(0)
{

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
  theRenderer->setProjectionMode("parallel");  // wire mode
  theRenderer->setPlaneDist(1.0, -1.0);
  theRenderer->setPRP(0.0, 0.0, 10.0);
  theRenderer->setPortWindow(-1.0, 1.0, -1.0, 1.0);  // use the whole window

  // copy the file name
  fileName1 = new char[strlen(_fileName1)+1];
  fileName2 = new char[strlen(_fileName2)+1];
  if (fileName1 == 0 || fileName2 == 0) {
    opserr << "FilePlotter::FilePlotter -- out of memory copying fileName1 " << endln;
  }
  strcpy(fileName1, _fileName1);    
  strcpy(fileName2, _fileName2);    
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

  if (fileName1 != 0)
    delete [] fileName1;

  if (fileName2 != 0)
    delete [] fileName2;

  if (data1a != 0)
    delete data1a;
  if (data1b != 0)
    delete data1b;
  if (data2a != 0)
    delete data2a;
  if (data2b != 0)
    delete data2b;
}
    
int 
FilePlotter::record(int cTag, double timeStamp)
{

  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    if (fileName1 != 0 && fileName2 == 0)
      return this->plotFile();
    else if (fileName1 != 0 && fileName2 != 0)
      return this->plotFiles();      
    else 
      return 0;

  } else
    return 0;
}

int 
FilePlotter::playback(int cTag)
{
  //  this->plotFile();
  return 0;
}

int
FilePlotter::restart(void)
{
  return 0;
}

int
FilePlotter::setFile(char *newFileName)
{
  strcpy(fileName1, newFileName);    
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
    theFile.open(fileName1, ios::in);

    if (theFile.bad()) {
	opserr << "WARNING - FilePlotter::FilePlotter()";
	opserr << " - could not open file " << fileName1 << endln;
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


    if (data1a == 0 || data1b == 0 || 
	data1a->Size() != numLineEntries || data1b->Size() != numLineEntries) {
      if (data1a != 0) delete data1a;
      if (data1b != 0) delete data1b;
      
      data1a = new Vector(numLineEntries);
      data1b = new Vector(numLineEntries);
    }

    // check colX and colY for valid entries
    if (numLines > 0) {
      if (cols == 0) {
	opserr << "FilePLotter::plotFile() - no valid columns have been set\n";

      } else {

	// parse through file checking the bounds
	Vector data(numLineEntries);
	while (theFile >> (*data1a)(0)) {
	  for (int j=1; j<numLineEntries; j++)
	    theFile >> (*data1a)(j);

	  for (int i=0; i<cols->Size(); i += 2) {
	    xValue = (*data1a)((*cols)(i));
	    yValue = (*data1a)((*cols)(i+1));
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

      // open the file again, read through and connect the dots
      ifstream theFile1; 
      theFile1.open(fileName1, ios::in);
      if (theFile1.bad()) {
	opserr << "WARNING - FilePlotter::FilePlotter()";
	opserr << " - could not open file " << fileName1 << endln;
	return -1;
      }    

      for (int ii=0; ii< numLineEntries; ii++)
	theFile1 >> (*data1a)(ii);

      for (int i=1; i<numLines; i++) {
	// read the data
	for (int ii=0; ii< numLineEntries; ii++)
	  theFile1 >> (*data1b)(ii);

	// plot the lines
	for (int j=0; j<cols->Size(); j+=2) {
	  pt1(0) = (*data1a)((*cols)(j)); 
	  pt1(1) = (*data1a)((*cols)(j+1));
	  pt2(0) = (*data1b)((*cols)(j)); 
	  pt2(1) = (*data1b)((*cols)(j+1));
	  theRenderer->drawLine(pt1, pt2, rgb, rgb);
	}
	
	*data1a = *data1b;

      }
      
      theRenderer->doneImage();

      // close the file
      theFile1.close();

    }
    return 0;
}



int
FilePlotter::plotFiles(void)
{
    /*
     * first pass: 1) open files
     *             2) determine number of entries on first line ASSSSSUME same through file
     *             3) determine number of lines and bounds [xmin, xmax, ymin, ymax]
     *             4) close the file
     */
  
     // open file
    ifstream theXfile; 
    theXfile.open(fileName1, ios::in);

    ifstream theYfile; 
    theYfile.open(fileName2, ios::in);

    if (theXfile.bad() || theYfile.bad()) {
	opserr << "WARNING - FilePlotter::FilePlotter()";
	opserr << " - could not open files " << fileName1 << " " << fileName2 << endln;
	return -1;
    }    

    double xMin, xMax, yMin, yMax;
    xMin = 0; xMax = 0; yMin =0; yMax =0;
    double xValue, yValue;

    // determine number of elements in each line 
    //  NOTE ASSUMES ALL LINES HAVE THE SAME NUMBER OF ELEMENTS
    char c;
    int numLineEntriesX = 0;
    int numLineEntriesY = 0;
    int numLinesX = 0;
    int numLinesY = 0;

    while (theXfile.get(c) && (c!= EOF) && c != '\n') {
      if (!isspace(c)) {
	theXfile.putback(c);

	theXfile >> xValue;
	for (int i=0; i<cols->Size(); i++) {
	  if (i%2 == 0) { // an xValue if numLineEntries == colX
	    if (numLineEntriesX == (*cols)(i)) {
	      if (xValue < xMin) xMin = xValue;
	      if (xValue > xMax) xMax = xValue;
	    }
	  }
	}
	numLineEntriesX ++;
      }
      numLinesX =1;
    }

    while (theYfile.get(c) && (c!= EOF) && c != '\n') {
      if (!isspace(c)) {
	theYfile.putback(c);

	theYfile >> yValue;
	for (int i=0; i<cols->Size(); i++) {
	  if (i%2 != 0) { // a yValue
	    if (numLineEntriesY == (*cols)(i)) {
	      if (yValue < yMin) yMin = yValue;
	      if (yValue > yMax) yMax = yValue;
	    }
	  }
	}
	numLineEntriesY ++;
      }
      numLinesY =1;
    }


    if (data1a == 0 || data1b == 0 || 
	data1a->Size() != numLineEntriesX || data1b->Size() != numLineEntriesX ||
	data2a == 0 || data2b == 0 || 
	data2a->Size() != numLineEntriesY || data2b->Size() != numLineEntriesY) {

      if (data1a != 0) delete data1a;
      if (data1b != 0) delete data1b;
      if (data2a != 0) delete data2a;
      if (data2b != 0) delete data2b;
      
      data1a = new Vector(numLineEntriesX);
      data1b = new Vector(numLineEntriesX);
      data2a = new Vector(numLineEntriesY);
      data2b = new Vector(numLineEntriesY);
    }

    // check colX and colY for valid entries
    if (numLinesX > 0 && numLinesY > 0) {
      if (cols == 0) {
	opserr << "FilePLotter::plotFile() - no valid columns have been set\n";
	return 0;
      } else {

	// parse through file checking the bounds
	while (theXfile >> (*data1a)(0)) {
	  for (int j=1; j<numLineEntriesX; j++)
	    theXfile >> (*data1a)(j);
	  
	  for (int i=0; i<cols->Size(); i += 2) {
	    xValue = (*data1a)((*cols)(i));
	    if (xValue < xMin) xMin = xValue;
	    if (xValue > xMax) xMax = xValue;
	    numLinesX++;
	  }    
	}

	while (theYfile >> (*data2a)(0)) {
	  for (int j=1; j<numLineEntriesY; j++)
	    theYfile >> (*data2a)(j);
	  
	  for (int i=0; i<cols->Size(); i += 2) {
	    yValue = (*data2a)((*cols)(i+1));
	    if (yValue < yMin) yMin = yValue;
	    if (yValue > yMax) yMax = yValue;
	    numLinesY++;
	  }    
	}

	// set the window bounds NOTE small border around the edges
	double xBnd = (xMax-xMin)/8;
	double yBnd = (yMax-yMin)/8;
    
	theRenderer->setViewWindow(xMin-xBnd,xMax+xBnd,yMin-yBnd,yMax+yBnd);
      }
    }

    // close the file
    theXfile.close();
    theYfile.close();


    if (numLinesX != numLinesY) {
	opserr << "FilePLotter::plotFile() - two files have different line counts\n";
	return -1;
    }

    /*
     * second pass: 1) open file
     *              2) get the renderer ready for drawing virgin image and then draw the x and y axis
     *              3) parse throgh the file and connect the dots
     *              4) close the file
     */

    if (numLinesX > 1) {

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

      // open the file again, read through and connect the dots
      ifstream theXfile1; 
      theXfile1.open(fileName1, ios::in);
      ifstream theYfile1; 
      theYfile1.open(fileName1, ios::in);

      if (theXfile1.bad() || theYfile1.bad()) {
	opserr << "WARNING - FilePlotter::FilePlotter()";
	opserr << " - could not open files " << fileName1 << " " << fileName2 << endln;
	return -1;
      }    

      for (int ii=0; ii< numLineEntriesX; ii++)
	theXfile1 >> (*data1a)(ii);

      for (int jj=0; jj< numLineEntriesY; jj++)
	theYfile1 >> (*data2a)(jj);

      for (int i=1; i<numLinesX; i++) {
	// read the data
	for (int ii=0; ii< numLineEntriesX; ii++)
	  theXfile1 >> (*data1b)(ii);
	for (int kk=0; kk< numLineEntriesY; kk++)
	  theYfile1 >> (*data2b)(kk);

	// plot the lines
	for (int j=0; j<cols->Size(); j+=2) {
	  pt1(0) = (*data1a)((*cols)(j)); 
	  pt1(1) = (*data2a)((*cols)(j+1));
	  pt2(0) = (*data1b)((*cols)(j)); 
	  pt2(1) = (*data2b)((*cols)(j+1));
	  theRenderer->drawLine(pt1, pt2, rgb, rgb);
	}
	
	*data1a = *data1b;
	*data2a = *data2b;
      }
      
      theRenderer->doneImage();

      // close the file
      theXfile1.close();
      theYfile1.close();

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







