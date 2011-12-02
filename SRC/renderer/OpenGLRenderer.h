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
// $Date: 2000-09-15 08:23:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGLRenderer.h,v $
                                                                        
                                                                        
// File: ~/renderer/OpenGLRenderer.h
//
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
#include <G3Globals.h>
#include <fstream.h>

#ifdef _UNIX

#else 
#include <windows.h>
#include <gl\gl.h>
#include <gl\glaux.h>
// include open gl stuff for win32
#endif

#include <db.H>
#include <Matrix.h>
#include <Vector.h>

class OpenGLRenderer : public Renderer
{
 public:
    OpenGLRenderer(char *title, int xLoc, int yLoc, int width, int height,
		   ColorMap &theMap);
    OpenGLRenderer(char *title, int xLoc, int yLoc, int width, int height,
		   ColorMap &theMap, char *texFileName, char *bmpFileName);	

    virtual ~OpenGLRenderer();

    virtual int clearImage(void);    
    virtual int startImage(void);
    virtual int doneImage(void);
    
    virtual int drawLine(const Vector &, const Vector &, 
			 float V1, float V2);
    
    virtual int drawLine(const Vector &end1, const Vector &end2, 
			 const Vector &rgb1, const Vector &rgb2);
   
    virtual int drawTriangle(const Vector &, const Vector &, const Vector &,
			     float V1, float V2, float V3);

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

    virtual int setProjectionMode(int); // 
    virtual int setFillMode(int);    // 1 = wire, otherwise fill
    
    virtual int setPRP(float u, float v, float n); // eye location if 
	                         // perspective, dirn to +ViewPlane if parallel

    // the following are in normalized coordinates
    virtual int setPortWindow(float, float, float, float); // view port
                              // left, right, bottom, top [-1,1,-1,1]
				  
    virtual int drawGText(const Vector &posGlobal, char *string, int length);    
    virtual int drawLText(const Vector &posLocal, char *string, int length);     
    
 protected:
    int saveBmpImage(void);  // to save the current image into a .BMP file

  // view
  VECTOR vrp;
  VECTOR vuv;
  VECTOR vpn;
  MATRIX ViewMat;
  
  // projection
  int projection_mode;
  VECTOR vpwindow;
  VECTOR planedist;
  VECTOR cop;
  MATRIX ProjMat;

  // clipping
  float X, Y, Zfar, Znear;

  float viewData[16];
  float projData[16];

  // viewport
  VECTOR portwindow;

  private:
#ifdef _UNIX
 
#else
  // win32 stuff
  HDC   theHDC;        // device context
  HGLRC theHRC;        // openGL context
  HWND  theWND;        // the window
#endif
  int winOpen;
  int height;       // current height of window in pixels
  int width;        // current width of window in pixels
  int numPoints;	
  int drawingPolygon;
  int xLoc;
  int yLoc;

  char title[50];
  int aFile;                              // int flag indicating if data to be 
                                          // sent to file or not
					      
  char theFileName[MAX_FILENAMELENGTH];   // tex file name  
  char theBmpFileName[MAX_FILENAMELENGTH];// bmp file name        
  int count;	                          // number of times done image has been invoked
  ofstream theFile; 	                  // output stream
  BITMAPINFO	info;
  GLint	        viewport[4];
  GLubyte       *bits;
  long	        currentBitSize;
  HBITMAP		theBitmap;
  
};

#endif

