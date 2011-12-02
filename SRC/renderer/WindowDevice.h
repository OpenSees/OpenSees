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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-26 18:56:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/WindowDevice.h,v $
                                                                        
#ifndef WindowDevice_H
#define WindowDevice_H

#include <Device.h>

#ifdef _UNIX
// include the X11 stuff
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/X.h>
#include <X11/Xatom.h>
#else 
// include open gl stuff for win32
#include <windows.h>
#include <gl\gl.h>
#include <gl\glaux.h>
#endif

#define MAX_NUM_POINTS_FOR_POLYGON 64
#define X11_MAX_COLORS 256

class WindowDevice: public  Device
{
 public:
  WindowDevice();
  ~WindowDevice();
  
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

  // draw text
  virtual void drawText(float x, float y, char *text, int length);  
  
  // Necessary when operating in XWINDOWS mode since the drawn
  // image is buffered until this call is made.
  virtual void STARTIMAGE();
  virtual void ENDIMAGE();

  // Opens a window of the specified width & height.
  virtual void WINOPEN(const char *title, int xLoc, int yLoc, int width, int height);

  // Clears the currently opened window
  virtual void CLEAR();

 private:
  void initWindow(void); // procedure called on construction of 1st Window
#ifdef _UNIX
  // X11 stuff

  Window theWindow; // the window associated with the Window
  GC theGC;         // the graphic context associated with the Window
  XSizeHints hints; // conatins the infor about where window is and its size
  static unsigned long foreground, background;  
  XStandardColormap theMap;  
  XEvent theEvent;
  static Display *theDisplay;  // the display all Window objecs display on
  static Colormap cmap;        // the colormap all X11 Window objects share   
  static int theScreen;        // the screen 
  static unsigned long pixels[X11_MAX_COLORS]; // pixels obtained from default
  static XColor colors[X11_MAX_COLORS];   // colors we define for our pixels
  static int colorFlag;        // flag indicating num of colors in colormap
  XPoint polygonPointArray[MAX_NUM_POINTS_FOR_POLYGON+1]; // +1 for wireframe
#else
  // win32 stuff
  HDC   theHDC;        // device context
  HGLRC theHRC;        // openGL context
  HWND  theWND;       // the window
#endif

  int winOpen;
  int height;       // current height of window in pixels
  int width;        // current width of window in pixels
  int numPoints;	
  int drawingPolygon;
  int xLoc;
  int yLoc;

  char title[50];   
  
  // class variables
  static int numWindowDevice;     // the number of WindowDevice objects in app  
};

#endif



  
