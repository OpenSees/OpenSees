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
// $Date: 2008-02-15 23:46:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/TclVideoPlayer.cpp,v $
                                                                        
                                                                        
// File: ~/tcl/TclVideoPlayer.C
// 
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class implementation for TclVideoPlayer
// TclVideoPlayer is a class for building a Plane Frame model in an interpreted
// environment. The constructor is used to add new commands to the interpreter,
// these commands are also defined in this file.
//
// What: "@(#) TclVideoPlayer.C, revA"


#include <stdlib.h>
#include <string.h>

#include <Renderer.h>

#ifndef _NOGRAPHICS

#ifdef _WGL
#include <OpenGLRenderer.h>
#elif _GLX
#include <OpenGLRenderer.h>
#elif _AGL
#include <OpenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif

#endif

#include <PlainMap.h>
#include "TclVideoPlayer.h"
#include <Vector.h>

#include <iomanip>
using std::ios;

//
// some static variables used in the functions
//


extern TclVideoPlayer *theTclVideoPlayer;


int
TclVideoPlayer_play(ClientData clientData, Tcl_Interp *interp, int argc, 
		    TCL_Char **argv);

//
// the class constructor, destructor and methods
//

TclVideoPlayer::TclVideoPlayer(const char *title, const char *fileName, const char *imageName,
			       Tcl_Interp *interp, const char *offsetFileName, double fact)
  :theMap(0), theRenderer(0), theOffsetFileName(0), factor(fact)
{
    // open the file
    if (fileName != 0) {

      // make room for copy of filename string and copy it
      theFileName = new char[strlen(fileName)]; 
      if (theFileName == 0) {
	opserr << "FATAL TclVideoPlayer::TclVideoPlayer() - out of memory copy " << fileName << endln;
	exit(-1);
      }

      strcpy(theFileName, fileName);    
      // test we can open the file for when we get to play - 
      theFile.open(fileName, ios::in);
      if (!theFile) {
	opserr << "WARNING TclVideoPlayer::TclVideoPlayer() - could not open file " << fileName << endln;
	exit(-1);
      } else {

	// read in the window properties from the file
	int xLoc, yLoc, width, height;
	char newTitle[50];
	theFile >> newTitle;
	theFile >> xLoc >> yLoc >> width >> height;

	// create the map and renderer
	theMap = new PlainMap();

#ifndef _NOGRAPHICS

#ifdef _WGL
	if (imageName == 0)
	  theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
	else
	  theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap,0,imageName);
#elif _GLX
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
#elif _AGL
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
#else
	theRenderer = new X11Renderer(title, xLoc, yLoc, width, height, *theMap);
#endif

#endif
	if (theMap == 0 || theRenderer == 0) 
	  opserr << "WARNING TclVideoPlayer::TclVideoPlayer() - could not create renderer\n";

	theFile.close();
      }
    }



    // see if a file containing rigid offsets has been supplied,
    // if so - make space, copy the string and try opening the file
    if (offsetFileName != 0) {
      // make copy of offset file name
      theOffsetFileName = new char[strlen(offsetFileName)]; 
      if (theOffsetFileName == 0) {
	opserr << "FATAL TclVideoPlayer::TclVideoPlayer() - out of memory copy " << offsetFileName << endln;
	exit(-1);
      }
      strcpy(theOffsetFileName, offsetFileName);    

      theOffsetFile.open(offsetFileName, ios::in);
      if (!theOffsetFile) {
	opserr << "WARNING TclVideoPlayer::TclVideoPlayer() - could not open file " << offsetFileName << endln;
      }
      theOffsetFile.close();
    }

    // call Tcl_CreateCommand for class specific commands
    Tcl_CreateCommand(interp, "play", TclVideoPlayer_play,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
}


TclVideoPlayer::~TclVideoPlayer() {
  if (theMap != 0)
    delete theMap;

  if (theRenderer != 0)
    delete theRenderer;
}

    
int 
TclVideoPlayer::play(void)
{
#ifdef _NOGRAPHICS
  return 0;
#else
  
  // rigid body offsets

  double dX, dY, dZ;
  dX = 0.0; dY = 0.0; dZ = 0.0;

  // open the file
  theFile.open(theFileName, ios::in);
  if (!theFile) 
    return -1;  // no error printed - print once in constructor


  // open the offset file
  if (theOffsetFileName != 0) {
    theOffsetFile.open(theOffsetFileName, ios::in);
    if (!theOffsetFile) 
      return -1;  // no error printed - print once in constructor
  }

  int xLoc, yLoc, width, height;
  char newTitle[50];
  theFile >> newTitle;
  theFile >> xLoc >> yLoc >> width >> height;

  // till we are done
  float a,b,c,d;
  int mode;
  Vector pt(3);
  Vector rgb(3);
  Vector pt2(3);
  Vector rgb2(3);

  int count =1;
  char drivel[50];
  while (theFile >> drivel) { // reads in StartImage
    // read in VRP data and set in renderer

    theFile >> drivel; 
    theFile >> a >> b >> c;
    theRenderer->setVRP(a,b,c);

    // read in VPN data and set in renderer
    theFile >> drivel; 
    theFile >> a >> b >> c;
    theRenderer->setVPN(a,b,c);

    // read in VUV data and set in renderer
    theFile >> drivel; 
    theFile >> a >> b >> c;
    theRenderer->setVUP(a,b,c);

    // read in COP data and set in renderer
    theFile >> drivel; 
    theFile >> a >> b >> c;
    theRenderer->setPRP(a,b,c);

    // read in projection mode data and set in renderer
    theFile >> drivel; 
    theFile >> drivel;
    theRenderer->setProjectionMode(drivel);

    // read in view data and set in renderer
    theFile >> drivel; 
    theFile >> a >> b >> c >> d;
    theRenderer->setViewWindow(a,b,c,d);

    // read in clipping plane distances and set in renderer
    theFile >> drivel; 
    theFile >> a >> b;
    theRenderer->setPlaneDist(a,b);

    // read in viewport data and set in renderer
    theFile >> drivel; 
    theFile >> a >> b >> c >> d;
    theRenderer->setPortWindow(a,b,c,d);

    theRenderer->startImage();
    theRenderer->clearImage();

    // read in the rigid body offsets
    if (theOffsetFileName != 0) 
      theOffsetFile >> dX >> dY >> dZ;

    int endImage = 0;
    while (endImage == 0 && theFile >> drivel) { 
      if (strcmp(drivel,"DoneImage") == 0) {
	theRenderer->doneImage();
	endImage = 1;
      }
      else if (strcmp(drivel,"Line") == 0) {
	theFile >> pt(0) >> pt(1) >> pt(2) >> rgb(0) >> rgb(1) >> rgb(2);
	theFile >> pt2(0) >> pt2(1) >> pt2(2) >> rgb2(0) >> rgb2(1) >> rgb2(2);

	// add rigid body offsets
	if (theOffsetFileName != 0) {
	  pt(0) += factor * dX; pt(1) += factor * dY; pt(2) += factor * dZ;
	  pt2(0) += factor * dX; pt2(1) += factor * dY; pt2(2) += factor * dZ;
	}

	theRenderer->drawLine(pt,pt2,rgb,rgb2);
      }
    }
  }

  // close the files
  theFile.close();
  if (theOffsetFileName != 0) 
    theOffsetFile.close();	

  return 0;
#endif
}

int
TclVideoPlayer_play(ClientData clientData, Tcl_Interp *interp, int argc, 
		    TCL_Char **argv)
{
#ifdef _NOGRAPHICS
  
#else
  if (theTclVideoPlayer != 0)
    theTclVideoPlayer->play();
#endif
  return TCL_OK;
}

