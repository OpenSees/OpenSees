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
// $Date: 2000-09-15 08:23:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGlDevice.h,v $
                                                                        
                                                                        
#ifndef Device_H
#define Device_H

#ifdef XWINDOWS
#include <X11/Xlib.h>
#include <X11/X.h>
#else
#include <gl.h>
#endif

#define MAX_NUM_POINTS_FOR_POLYGON 256

class Device
{
 public:
  Device();
  virtual ~Device();  

  // Specify a 2D point to the hardware  
  virtual void V2F(float x, float y);

  // Specify a color to the hardware
  virtual void C3F(float r, float g, float b);

  // Gets the width of the current window
  virtual int GetWidth();

  // Gets the height of the current window
  virtual int GetHeight();

  // Call when about to begin/end drawing a polygon. All V2F calls
  // from then on until ENDPOLYGON will be interpreted as 
  // vertices of the polygon
  virtual void BGNPOLYGON();
  virtual void ENDPOLYGON();


  // Same as BGNPOLYGON but for wireframe polygons
  virtual void BGNCLOSEDLINE();
  virtual void ENDCLOSEDLINE();

  // Call when about to begin drawing a set of points. All V2F
  // calls from then on until ENDPOINT will be interpreted as
  // points.
  virtual void BGNPOINT();
  virtual void ENDPOINT();

  // Necessary when operating in XWINDOWS mode since the drawn
  // image is buffered until this call is made.
  virtual void ENDIMAGE();

  // Opens a window of the specified width & height.
  virtual void WINOPEN(int width, int height);

  // Clears the currently opened window
  virtual void CLEAR();


 private:
#ifdef XWINDOWS
  int numPoints;
  int drawingPolygon;
  XPoint polygonPointArray[MAX_NUM_POINTS_FOR_POLYGON+1]; // +1 for wireframe polygons
  Display *theDisplay;
  Window theWindow;
  GC theGC;
  Colormap cmap;

#else
  int window_id;
#endif
  int width, height;		// Width and height of our window
};

#endif



  
