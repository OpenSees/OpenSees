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
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/TclVideoPlayer.cpp,v $
                                                                        
                                                                        
// File: ~/tcl/TclVideoPlayer.C
// 
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class implementation for TclVideoPlayer
// TclVideoPlayer is a class for building a Plane Frame model in an interpreted
// enviroment. The constructor is used to add new commands to the interpreter,
// these commands are also defined in this file.
//
// What: "@(#) TclVideoPlayer.C, revA"


#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#ifdef _WIN32
#include <OPenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif


#include <PlainMap.h>
#include "TclVideoPlayer.h"
#include <Vector.h>


//
// some static variables used in the functions
//

extern TclVideoPlayer *theTclVideoPlayer;

int
TclVideoPlayer_play(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv);

//
// the class constructor, destructor and methods
//

TclVideoPlayer::TclVideoPlayer(char *title, char *fileName, char *imageName,
			       Tcl_Interp *interp)
  :theMap(0), theRenderer(0)
{
    // open the file
    if (fileName != 0) {

      // make room for copy of filename string and copy it
      theFileName = new char[strlen(fileName)]; 
      if (theFileName == 0) 
	g3ErrorHandler->fatal("FATAL TclVideoPlayer::TclVideoPlayer() - out of memory copy %s\n",
			      fileName);
      strcpy(theFileName, fileName);    
      // test we can open the file for when we get to play - 
      theFile.open(fileName, ios::in);
      if (!theFile) {
	g3ErrorHandler->warning("WARNING TclVideoPlayer::TclVideoPlayer() - could not open file %s\n",
				  fileName);	  
      } else {

	// read in the window properties from the file
	int xLoc, yLoc, width, height;
	char newTitle[50];
	theFile >> newTitle;
	theFile >> xLoc >> yLoc >> width >> height;

	// create the map and renderer
	theMap = new PlainMap();
#ifdef _WIN32
	if (imageName == 0)
	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
	else

	theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap,0,imageName);
#else
	theRenderer = new X11Renderer(title, xLoc, yLoc, width, height, *theMap);

#endif
	if (theMap == 0 || theRenderer == 0) 
	  g3ErrorHandler->warning("WARNING TclVideoPlayer::TclVideoPlayer() - could not create renderer\n");

	theFile.close();
      }
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
  // open the file
  theFile.open(theFileName, ios::in);
  if (!theFile) 
    return -1;  // no error printed - print once in constructor


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
    theFile >> mode;
    theRenderer->setProjectionMode(mode);

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

    
    int endImage = 0;
    while (endImage == 0 && theFile >> drivel) { 
      if (strcmp(drivel,"DoneImage") == 0) {
	theRenderer->doneImage();
	endImage = 1;
      }
      else if (strcmp(drivel,"Line") == 0) {
	theFile >> pt >> rgb;
	theFile >> pt2 >> rgb2;
	theRenderer->drawLine(pt,pt2,rgb,rgb2);
      }
    }
  }
	
  theFile.close();
  return 0;
}

int
TclVideoPlayer_play(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv)
{
  if (theTclVideoPlayer != 0)
	theTclVideoPlayer->play();
  return TCL_OK;
}

