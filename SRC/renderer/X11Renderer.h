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
// $Date: 2003-02-25 23:34:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/X11Renderer.h,v $
                                                                        
                                                                        
// File: ~/renderer/X11Renderer.h
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for X11Renderer.
// X11Renderer is an abstract base class. An X11Renderer object is used
// to create an image of the domain.
//
// What: "@(#) X11Renderer.h, revA"

#ifndef X11Renderer_h
#define X11Renderer_h

#include <Renderer.h>

class View;
class Projection;
class Clipping;
class Viewport;
class ScanLineConverter;
class Device;
class FACE;

#include <fstream>
using std::ofstream;

class X11Renderer : public Renderer
{
 public:
    X11Renderer(const char *title, int xLoc, int yLoc, int width, int height,
		ColorMap &theMap);
    X11Renderer(const char *title, int xLoc, int yLoc, int width, int height,
		ColorMap &theMap, const char *fileName);

    virtual ~X11Renderer();

    virtual int clearImage(void);    
    virtual int startImage(void);
    virtual int doneImage(void);

    virtual int drawPoint(const Vector &, float V1, int width = 1);
    virtual int drawPoint(const Vector &, const Vector &rgb1, int width = 1);
    
    virtual int drawLine(const Vector &, const Vector &, 
			 float V1, float V2, int width = 1, int style = 1);

    virtual int drawLine(const Vector &end1, const Vector &end2, 
			 const Vector &rgb1, const Vector &rgb2,
			 int width = 1, int style = 1);
    
    virtual int drawPolygon(const Matrix &points, const Vector &values);
    virtual int drawPolygon(const Matrix &points, const Matrix &rgbValues);

    // 
    // the following are for setting up the vieing system
    //

    // the following are in world coordinates & define view coord system
    virtual int setVRP(float x, float y, float z); // point on view plane    
    virtual int setVPN(float x, float y, float z); // view plane normal
    virtual int setVUP(float x, float y, float z); // view-up vector
	
    // the following are in view coordinates	
    virtual int setViewWindow(float, float, float, float); // view bounds
                               // umin, umax, vmin, vmax

    virtual int setPlaneDist(float, float); // location of
                               // near and far clipping planes

    virtual int setProjectionMode(const char *mode); // 
    virtual int setFillMode(const char *);    // 1 = wire, otherwise fill
    
    virtual int setPRP(float u, float v, float n); // eye location if 
	                         // perspective, dirn to +ViewPlane if parallel

    // the following are in normalized coordinates
    virtual int setPortWindow(float, float, float, float); // view port
                              // left, right, bottom, top [-1,1,-1,1]
				  
    virtual int drawText(const Vector &posGlobal, char *string, int length, 
			 char horizontalJustify = 'l', char verticalJustify = 'b');    
    
 protected:

    int displayFace(FACE &);    

    View *theView;
    Projection *theProjection;
    Clipping *theClipping;
    Viewport *theViewport;
    ScanLineConverter *theScan;
    Device *theDevice;

  private:
    int aFile;                              // int flag indicating if data to be sent to file or not
    char *theFileName;                      // file name  
    ofstream theFile; 	                    // output stream    
};

#endif

