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
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/TclFeViewer.cpp,v $
                                                                        
                                                                        
// File: ~/tcl/TclFeViewer.C
// 
// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the class implementation for TclFeViewer
//
// What: "@(#) TclFeViewer.C, revA"


#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>

#ifdef _WIN32
#include <OpenGLRenderer.h>
#else
#include <X11Renderer.h>
#endif

#include <PlainMap.h>

#include "TclFeViewer.h"

//
// some static variables used in the functions
//

static TclFeViewer *theTclFeViewer = 0;

// 
// the functions that will be invoked by the interpreter while building the model
//

int
TclFeViewer_setVRP(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv);
int
TclFeViewer_setVPN(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv);		   
int
TclFeViewer_setVUP(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv);
int
TclFeViewer_setViewWindow(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv);
int
TclFeViewer_setPlaneDist(ClientData clientData, Tcl_Interp *interp, int argc, 
			 char **argv);		   
int
TclFeViewer_setProjectionMode(ClientData clientData, Tcl_Interp *interp, int argc, 
			      char **argv);		   			 
int
TclFeViewer_setFillMode(ClientData clientData, Tcl_Interp *interp, int argc, 
			char **argv);		   			 			      
int
TclFeViewer_setPRP(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv);		   
int
TclFeViewer_setPortWindow(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv);		   			 
int
TclFeViewer_displayModel(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv);		   			 
			  
int
TclFeViewer_clearImage(ClientData clientData, Tcl_Interp *interp, int argc, 
		  char **argv);		 
//
// the class constructor, destructor and methods
//

TclFeViewer::TclFeViewer(char *title, int xLoc, int yLoc, int width, int height,
			 Domain &_theDomain, int WipeFlag,
			 Tcl_Interp *interp)
  :theMap(0),theRenderer(0), theDomain(&_theDomain), wipeFlag(WipeFlag),
   theEleMode(-1), theNodeMode(-1), theDisplayFact(1)
{

  // set the static pointer used in the class
  theTclFeViewer = this;
  theMap = new PlainMap();

#ifdef _WIN32
  theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap);
#else
  theRenderer = new X11Renderer(title, xLoc, yLoc, width, height, *theMap);
#endif

  
  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "vrp", TclFeViewer_setVRP,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "vpn", TclFeViewer_setVPN,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "vup", TclFeViewer_setVUP,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "viewWindow", TclFeViewer_setViewWindow,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  
  Tcl_CreateCommand(interp, "plane", TclFeViewer_setPlaneDist,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
  
  Tcl_CreateCommand(interp, "projection", TclFeViewer_setProjectionMode,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
  
  Tcl_CreateCommand(interp, "fill", TclFeViewer_setFillMode,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);      
  
  Tcl_CreateCommand(interp, "prp", TclFeViewer_setPRP,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "port", TclFeViewer_setPortWindow,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  

  Tcl_CreateCommand(interp, "display", TclFeViewer_displayModel,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
  
  Tcl_CreateCommand(interp, "clearImage", TclFeViewer_clearImage,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);      
}



TclFeViewer::TclFeViewer(char *title, int xLoc, int yLoc, int width, int height, 
						 char *fileName,
			  Domain &_theDomain,
			 Tcl_Interp *interp)
  :theMap(0),theRenderer(0), theDomain(&_theDomain),
   theEleMode(-1), theNodeMode(-1), theDisplayFact(1), wipeFlag(1), count(0)
{

  // set the static pointer used in the class
  theTclFeViewer = this;
  theMap = new PlainMap();
#ifdef _WIN32
  cerr << "HERE I AM\n";
  theRenderer = new OpenGLRenderer(title, xLoc, yLoc, width, height, *theMap, 0, fileName);
#else
  theRenderer = new X11Renderer(title, xLoc, yLoc, width, height, *theMap, fileName);
#endif
  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "vrp", TclFeViewer_setVRP,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "vpn", TclFeViewer_setVPN,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "vup", TclFeViewer_setVUP,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "viewWindow", TclFeViewer_setViewWindow,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  
  Tcl_CreateCommand(interp, "plane", TclFeViewer_setPlaneDist,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
  
  Tcl_CreateCommand(interp, "projection", TclFeViewer_setProjectionMode,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
  
  Tcl_CreateCommand(interp, "fill", TclFeViewer_setFillMode,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);      
  
  Tcl_CreateCommand(interp, "prp", TclFeViewer_setPRP,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "port", TclFeViewer_setPortWindow,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  

  Tcl_CreateCommand(interp, "display", TclFeViewer_displayModel,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
  
  Tcl_CreateCommand(interp, "clearImage", TclFeViewer_clearImage,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);      
}

TclFeViewer::~TclFeViewer()
{
  // may possibly invoke Tcl_DeleteCommand() later
  // for moment just invoke destructor on map and renderer
  // and set pointer to NULL
  delete theMap;
  delete theRenderer;
  theTclFeViewer = 0;  
}
    
int 
TclFeViewer::record(int cTag)
{
  //  theRenderer->displayModel(thetheEleMode, theNodeMode, theDisplayFact);

  // loop over the elements getting each to display itself
  // using theRenderer and displayTag as arguments.
  // first clear the image
		
  if (wipeFlag == 1)
		theRenderer->clearImage();

  theRenderer->startImage();
  int res = 0;

  if (theEleMode >= 0) {
      ElementIter &theElements = theDomain->getElements();
      Element *theEle;
      while ((theEle = theElements()) != 0) {
	  res = theEle->displaySelf(*theRenderer, theEleMode, theDisplayFact);
	  if (res < 0) {
	      cerr << "Renderer::displayModel() - Element: ";
	      cerr << theEle->getTag() << " failed to display itself\n";
	  }
      }
  }
  
  if (theNodeMode >= 0) {
      NodeIter &theNodes = theDomain->getNodes();
      Node *theNode;
      while ((theNode = theNodes()) != 0) {
	  res = theNode->displaySelf(*theRenderer, theNodeMode, theDisplayFact);
	  if (res < 0) {
	      cerr << "Renderer::displayModel() - Node: ";
	      cerr << theNode->getTag() << " failed to display itself\n";
	  }
      }
  }  

  // now mark the image has having been completed
  theRenderer->doneImage();

  return res;
}

int 
TclFeViewer::playback(int cTag)
{
  return 0;
}

void
TclFeViewer::restart(void)
{
    
}


int
TclFeViewer::setVRP(float xLoc, float yLoc , float zLoc)
{
    // point on view plane    
    return theRenderer->setVRP(xLoc, yLoc, zLoc);
}

int
TclFeViewer::setVPN(float xdirn, float ydirn, float zdirn)
{
    // view plane normal
    return theRenderer->setVPN(xdirn, ydirn, zdirn);
}

int
TclFeViewer::setVUP(float xdirn, float ydirn, float zdirn)
{
    // view-up vector
    return theRenderer->setVUP(xdirn, ydirn, zdirn);
}    

int
TclFeViewer::setViewWindow(float uMin, float uMax, float vMin, float vMax)
{
    // view bounds  umin, umax, vmin, vmax
    return theRenderer->setViewWindow(uMin, uMax, vMin, vMax);
}        

int
TclFeViewer::setPlaneDist(float anear, float afar)
{
    // location of near, view & far clipping planes
    return theRenderer->setPlaneDist(anear,afar);
}            

int
TclFeViewer::setProjectionMode(int mode)
{
    return theRenderer->setProjectionMode(mode);
}                

int
TclFeViewer::setFillMode(int mode)
{
    // 1 = wire, otherwise fill
    return theRenderer->setFillMode(mode);
}                

int
TclFeViewer::setPRP(float uLoc, float vLoc , float nLoc)
{
    // eye location -- global coordinates
    return theRenderer->setPRP(uLoc, vLoc, nLoc);
}
    
    
int
TclFeViewer::setPortWindow(float left, float right, float bottom, float top)
{     
    // view port left, right, bottom, top [-1,1,-1,1]
    return theRenderer->setPortWindow(left,right,bottom,top);    
}

int
TclFeViewer::displayModel(int eleFlag, int nodeFlag, float displayFact)
{
    // methods invoked on the FE_Viewer
    theEleMode = eleFlag;
    theNodeMode = nodeFlag;    
    theDisplayFact = displayFact;    
    return this->record(0);
}

int
TclFeViewer::clearImage(void)
{
    return theRenderer->clearImage();
}
    



int
TclFeViewer_setVRP(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // ensure corrcet num arguments
  if (argc < 4) {
      interp->result = "WARNING args incorrect - vrp xloc yloc zloc ";
      return TCL_ERROR;
  }    

  // get the xLoc, yLoc and zLoc
  double xLoc, yLoc, zLoc;
  if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
      interp->result = "WARNING invalid x_loc - vrp x_loc y_loc z_loc";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &yLoc) != TCL_OK) {
      interp->result = "WARNING invalid y_loc - vrp x_loc y_loc z_loc";
      return TCL_ERROR;
  }  
  if (Tcl_GetDouble(interp, argv[3], &zLoc) != TCL_OK) {
      interp->result = "WARNING invalid z_loc - vrp x_loc y_loc z_loc";
      return TCL_ERROR;
  }    
  
  theTclFeViewer->setVRP(xLoc,yLoc,zLoc);
  return TCL_OK;  
}

int
TclFeViewer_setVPN(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv)
{  
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
      interp->result = "WARNING args incorrect - vpn xdirn ydirn zdirn ";
      return TCL_ERROR;
  }    

  // get the id, x_dirn and y_dirn
  double xDirn, yDirn, zDirn;
  if (Tcl_GetDouble(interp, argv[1], &xDirn) != TCL_OK) {
      interp->result = "WARNING invalid x_dirn - vpn x_dirn y_dirn z_dirn";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &yDirn) != TCL_OK) {
      interp->result = "WARNING invalid y_dirn - vpn x_dirn y_dirn z_dirn";
      return TCL_ERROR;
  }  
  if (Tcl_GetDouble(interp, argv[3], &zDirn) != TCL_OK) {
      interp->result = "WARNING invalid z_dirn - vpn x_dirn y_dirn z_dirn";
      return TCL_ERROR;
  }    
  
  theTclFeViewer->setVPN(xDirn,yDirn,zDirn);
  return 0;
}

int
TclFeViewer_setVUP(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv)
{  
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
      interp->result = "WARNING args incorrect - vup xdirn ydirn zdirn ";
      return TCL_ERROR;
  }    

  // get the id, x_dirn and y_dirn
  double xDirn, yDirn, zDirn;
  if (Tcl_GetDouble(interp, argv[1], &xDirn) != TCL_OK) {
      interp->result = "WARNING invalid x_dirn - vup x_dirn y_dirn z_dirn";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &yDirn) != TCL_OK) {
      interp->result = "WARNING invalid y_dirn - vup x_dirn y_dirn z_dirn";
      return TCL_ERROR;
  }  
  if (Tcl_GetDouble(interp, argv[3], &zDirn) != TCL_OK) {
      interp->result = "WARNING invalid z_dirn - vup x_dirn y_dirn z_dirn";
      return TCL_ERROR;
  }    
  
  theTclFeViewer->setVUP(xDirn,yDirn,zDirn);
  return TCL_OK;  
}

int
TclFeViewer_setViewWindow(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv)
{  
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // check num args
  if (argc < 5) {
      interp->result = "WARNING args incorrect - prp uMin uMax vMin vMax ";
      return TCL_ERROR;
  }    

  double uMin, uMax, vMin, vMax;
  if (Tcl_GetDouble(interp, argv[1], &uMin) != TCL_OK) {
      interp->result = "WARNING invalid uMin - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &uMax) != TCL_OK) {
      interp->result = "WARNING invalid uMax - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }  
  if (Tcl_GetDouble(interp, argv[3], &vMin) != TCL_OK) {
      interp->result = "WARNING invalid vMin - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }    
  if (Tcl_GetDouble(interp, argv[4], &vMax) != TCL_OK) {
      interp->result = "WARNING invalid vMin - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }      
  
  theTclFeViewer->setViewWindow(uMin,uMax,vMin,vMax);
  return TCL_OK;    
}
int
TclFeViewer_setPlaneDist(ClientData clientData, Tcl_Interp *interp, int argc, 
			 char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // check num args
  if (argc < 3) {
      interp->result = "WARNING args incorrect - dist near far ";
      return TCL_ERROR;
  }    
  // get distnces to near view and far planes
  double anear, afar;
  if (Tcl_GetDouble(interp, argv[1], &anear) != TCL_OK) {
      interp->result = "WARNING invalid near - vup near far";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &afar) != TCL_OK) {
      interp->result = "WARNING invalid far - vup near far";
      return TCL_ERROR;
  }  

  theTclFeViewer->setPlaneDist(anear, afar);    
  return TCL_OK;  
}

int
TclFeViewer_setProjectionMode(ClientData clientData, Tcl_Interp *interp, int argc, 
			      char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // ensure corrcet num arguments
  if (argc < 2) {
      interp->result = "WARNING args incorrect - projection modeID ";
      return TCL_ERROR;
  }    

  // get the mode
  int i;
  if (Tcl_GetInt(interp, argv[1], &i) != TCL_OK) {
      interp->result = "WARNING invalid modeID - projection modeID";
      return TCL_ERROR;
  }
  
  theTclFeViewer->setProjectionMode(i);    
  return TCL_OK;  
}

int
TclFeViewer_setFillMode(ClientData clientData, Tcl_Interp *interp, int argc, 
			char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // ensure corrcet num arguments
  if (argc < 2) {
      interp->result = "WARNING args incorrect - fill modeID ";
      return TCL_ERROR;
  }    

  // get the mode
  int i;
  if (Tcl_GetInt(interp, argv[1], &i) != TCL_OK) {
      interp->result = "WARNING invalid modeID - fill modeID";
      return TCL_ERROR;
  }
  
  theTclFeViewer->setFillMode(i);    
  return TCL_OK;  
}

int
TclFeViewer_setPRP(ClientData clientData, Tcl_Interp *interp, int argc, 
		   char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // ensure corrcet num arguments
  if (argc < 4) {
      interp->result = "WARNING args incorrect - cop xloc yloc zloc ";
      return TCL_ERROR;
  }    

  // get the xLoc, yLoc and zLoc
  double xLoc, yLoc, zLoc;
  if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
      interp->result = "WARNING invalid x_loc - cop x_loc y_loc z_loc";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &yLoc) != TCL_OK) {
      interp->result = "WARNING invalid y_loc - cop x_loc y_loc z_loc";
      return TCL_ERROR;
  }  
  if (Tcl_GetDouble(interp, argv[3], &zLoc) != TCL_OK) {
      interp->result = "WARNING invalid z_loc - cop x_loc y_loc z_loc";
      return TCL_ERROR;
  }    
  
  theTclFeViewer->setPRP(xLoc,yLoc,zLoc);
  return TCL_OK;  
}

int
TclFeViewer_setPortWindow(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // check num args  
  if (argc < 5) {
      interp->result = "WARNING args incorrect - prp uMin uMax vMin vMax ";
      return TCL_ERROR;
  }    

  double uMin, uMax, vMin, vMax;
  if (Tcl_GetDouble(interp, argv[1], &uMin) != TCL_OK) {
      interp->result = "WARNING invalid uMin - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &uMax) != TCL_OK) {
      interp->result = "WARNING invalid uMax - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }  
  if (Tcl_GetDouble(interp, argv[3], &vMin) != TCL_OK) {
      interp->result = "WARNING invalid vMin - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }    
  if (Tcl_GetDouble(interp, argv[4], &vMax) != TCL_OK) {
      interp->result = "WARNING invalid vMin - vup uMin uMax vMin vMax";
      return TCL_ERROR;
  }      

  theTclFeViewer->setPortWindow(uMin,uMax,vMin,vMax);    
  return TCL_OK;
}

int
TclFeViewer_displayModel(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  // check number of args  
  if (argc != 3 && argc != 4) {
      interp->result = "WARNING args incorrect - display eleMode <nodeMode> fact";
      return TCL_ERROR;
  }    

  if (argc == 3) {
      int displayMode;
      double displayFact;
      if (Tcl_GetInt(interp, argv[1], &displayMode) != TCL_OK) {
	  interp->result = "WARNING invalid displayMode - display displayMode displayFact";
	  return TCL_ERROR;
	      }
      if (Tcl_GetDouble(interp, argv[2], &displayFact) != TCL_OK) {
	  interp->result = "WARNING invalid displayFact - display displayMode displayFact";
	  return TCL_ERROR;
      }  

      theTclFeViewer->displayModel(displayMode, -1, displayFact);
      return TCL_OK;    
  } else {
      
      int eleFlag, nodeFlag;
      double displayFact;
      if (Tcl_GetInt(interp, argv[1], &eleFlag) != TCL_OK) {
	  interp->result = "WARNING invalid displayMode - display eleFlag nodeFlag displayFact";
	  return TCL_ERROR;
	      }
      if (Tcl_GetInt(interp, argv[2], &nodeFlag) != TCL_OK) {
	  interp->result = "WARNING invalid displayMode - display eleFlag nodeFlahg displayFact";
	  return TCL_ERROR;
	      }      
      if (Tcl_GetDouble(interp, argv[3], &displayFact) != TCL_OK) {
	  interp->result = "WARNING invalid displayMode - display eleFlag nodeFlahg displayFact";
	  return TCL_ERROR;
      }  

      theTclFeViewer->displayModel(eleFlag, nodeFlag, displayFact);
      return TCL_OK;    
  }
}



int
TclFeViewer_clearImage(ClientData clientData, Tcl_Interp *interp, int argc, 
			  char **argv)
{
  // check destructor has not been called
  if (theTclFeViewer == 0)
      return TCL_OK;    
  
  theTclFeViewer->clearImage();
  return TCL_OK;
}




