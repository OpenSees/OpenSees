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
// $Date: 2000-09-15 08:23:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/X11Device.h,v $
                                                                        
                                                                        
#ifndef X11Device_H
#define X11Device_H

#include <Device.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/X.h>
#include <X11/Xatom.h>

#define MAX_NUM_POINTS_FOR_POLYGON 64
#define X11_MAX_COLORS 256

class X11Device: public  Device
{
 public:
  X11Device();
  ~X11Device();
  
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
  virtual void WINOPEN(int width, int height);

  // Clears the currently opened window
  virtual void CLEAR();

 private:
  void initX11(void); // procedure called once on construction of first X11Device

   Window theWindow; // the window associated with the X11Device
  GC theGC;         // the graphic context associated with the X11Device
  int winOpen;
  unsigned long foreground, background;  
  XSizeHints hints; // conatins the infor about where window is and its size
  int height;       // current height of the window

  XStandardColormap theMap;  
  XEvent theEvent;
  
  int numPoints;	
  int drawingPolygon;
  XPoint polygonPointArray[MAX_NUM_POINTS_FOR_POLYGON+1]; // +1 for wireframe 
  
  // class variables
  static int numX11Device;     // the number of X11Device objects in app  
  static Display *theDisplay;  // the display all X11Devices display on
  static Colormap cmap;        // the colormap all X11Device objects share      
  static int theScreen;        // the screen 
  static unsigned long pixels[X11_MAX_COLORS]; // the pixels we obtain from default
  static XColor colors[X11_MAX_COLORS];   // the colors we define for our pixels
  static int colorFlag;        // flag indicating num of colors in colormap

};

#endif



  

