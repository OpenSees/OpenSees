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
                                                                        
// $Revision: 1.12 $
// $Date: 2008-04-17 19:33:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGLRenderer.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for OpenGLRenderer.
// OpenGLRenderer is an abstract base class. An OpenGLRenderer object is used
// to create an image of the domain.
//
// What: "@(#) OpenGLRenderer.h, revA"

#ifndef OpenGLRenderer_h
#define OpenGLRenderer_h

#include <Renderer.h>
#include <fstream>
using std::ofstream;

#include <OpenGlDevice.h>

#include <db.H>
#include <Matrix.h>
#include <Vector.h>

class OpenGLRenderer : public Renderer
{
  public:
    OpenGLRenderer(const char *title, int xLoc, int yLoc, int width, int height,
		   ColorMap &theMap);
    OpenGLRenderer(const char *title, int xLoc, int yLoc, int width, int height,
		   ColorMap &theMap, const char *texFileName, const char *bmpFileName);	

    virtual ~OpenGLRenderer();

    virtual int clearImage(void);    
    virtual int saveImage(const char *imageName);    
    virtual int startImage(void);
    virtual int doneImage(void);

    virtual int drawPoint(const Vector &, float V1, int tag = 0, int mode =0, int width = 1);
    virtual int drawPoint(const Vector &, const Vector &rgb1, int tag = 0, int mode=0, int width = 1);    

    virtual int drawLine(const Vector &, const Vector &, 
			 float V1, float V2, int tag = 0, int mode = 0, int width = 1, int style = 1);
    virtual int drawLine(const Vector &end1, const Vector &end2, 
			 const Vector &rgb1, const Vector &rgb2,
			 int tag = 0, int mode = 0, int width = 1, int style = 1);
   
    virtual int drawPolygon(const Matrix &points, const Vector &values, int tag = 0, int mode = 0);
    virtual int drawPolygon(const Matrix &points, const Matrix &rgbValues, int tag = 0, int mode = 0);

    virtual int drawText(const Vector &posGlobal, char *string, int length, 
			 char horizontalJustify = 'l', char verticalJustify = 'b');    

    // the following are in world coordinates & define view coord system
    virtual int setVRP(float x, float y, float z); // point on view plane    
    virtual int setVPN(float x, float y, float z); // view plane normal
    virtual int setVUP(float x, float y, float z); // view-up vector
	
    // the following are in view coordinates	
    virtual int setViewWindow(float, float, float, float); // view bounds
                               // umin, umax, vmin, vmax

    virtual int setPlaneDist(float, float); // location of
                               // near and far clipping planes

    virtual int setProjectionMode(const char *mode); // parallel or perspective
    virtual int setFillMode(const char *mode);    // wire or  fill
    
    virtual int setPRP(float u, float v, float n); // eye location if 
	                         // perspective, dirn to +ViewPlane if parallel

    // the following are in normalized coordinates
    virtual int setPortWindow(float, float, float, float); // view port
                              // left, right, bottom, top [-1,1,-1,1]
				  
 protected:

 private:
    char *windowTitle; // title name of the window
    int height;        // current height of window in pixels
    int width;         // current width of window in pixels
    int xLoc;          // upper xLocation of window
    int yLoc;          // upper yLocation of window

    int count;	               // number of times done image has been invoked
    ofstream theFile; 	       // output stream if saving drawing commands
    char *theOutputFileName;   // file name for output stream

    OpenGlDevice *theDevice;

    // viewing 
    Vector vrp;  // point on the view plane - global coords
    Vector vuv;  // vector defining the view up vector, 
    Vector vpn;  // vector defining the view plane normal
    Vector cop;  // eye location - NOW IN GLOBAL COORDINATES
    Matrix ViewMat;

    // projection
    int projectionMode;        // flag indicating projection mode
    Vector vpWindow;           // view window bounds - local window coordinates (u,v)
    double clippingPlanes[2];  // distance to front and back clipping planes FROM THE PLANE (n)
    Matrix ProjMat;         

    // viewport
    Vector portWindow;  // mapping to window - port window coords [-1,-1] to [1,1]

    int fillMode;        // flag indicating fill mode

    float viewData[16];
    float projData[16];
};

#endif

