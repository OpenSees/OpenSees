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
// $Date: 2007-06-26 20:13:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/TclFeViewer.h,v $
                                                                        
                                                                        
// File: ~/modelbuilder/tcl/TclFeViewer.h.h
// 
// Written: fmk 
// Created: 4/99
// Revision: A
//
// Description: This file contains the class definition for TclFeViewer.
// A TclFeViewer adds commands to the interpreter for displaying the model.
//
// What: "@(#) ModelBuilder.h, revA"

#ifndef TclFeViewer_h
#define TclFeViewer_h

#include <Recorder.h>
class Renderer;
class ColorMap;

extern "C" {
#include <tcl.h>
}

class TclFeViewer : public Recorder
{
  public:
    TclFeViewer(const char *title, int xLoc, int yLoc, int width, int height,
		Domain &theDomain, int wipeFlag, Tcl_Interp *interp,
        double deltaT = 0.0);

    TclFeViewer(const char *title, int xLoc, int yLoc, int width, int height,
        const char *fileName, Domain &theDomain, Tcl_Interp *interp,
        double deltaT = 0.0);

    TclFeViewer();
    ~TclFeViewer();

    int buildFE_Model(void);
    
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);
    int restart(void);    

    // methods invoked on the ViewingSystem
    int setVRP(float, float, float); // point on view plane    
    int setVPN(float, float, float); // view plane normal
    int setVUP(float, float, float); // view-up vector
    int setViewWindow(float, float, float, float); // view bounds
                               // umin, umax, vmin, vmax

    int setPlaneDist(float, float); // location of
                               // near, view & far clipping planes

    int setProjectionMode(const char *); // 
    int setFillMode(const char *);    // 1 = wire, otherwise fill
    
    int setPRP(float, float, float); // eye location, global coords

    int setPortWindow(float, float, float, float); // view port
                              // left, right, bottom, top [-1,1,-1,1]


    // methods invoked on the FE_Viewer
    int displayModel(int eleFlag, int nodeFlag, float displayFact);
    int clearImage(void);
    int saveImage(const char *fileName);
    int saveImage(const char *imageName, const char *fileName);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

  protected:

  private:
    ColorMap *theMap;
    Renderer *theRenderer;
    Domain *theDomain;
    int theEleMode;
    int theNodeMode;    
    double theDisplayFact;
    double deltaT;
    double nextTimeStampToRecord;
    int wipeFlag;
    int vrpSet;
    int vpwindowSet;
    int clippingPlaneDistancesSet;
};

#endif
