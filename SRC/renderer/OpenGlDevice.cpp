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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGlDevice.cpp,v $
                                                                        
                                                                        
#include "Device.h"
#include <iostream.h>

Device::Device()
{
#ifdef XWINDOWS
  drawingPolygon = 0;
  numPoints = 0;
#else
  window_id = 0;
#endif
}

Device::~Device()
{
#ifdef XWINDOWS
  drawingPolygon = 0;
  numPoints = 0;
#else
  window_id = 0;
#endif
}

void
Device::WINOPEN(int _width, int _height)
{
  width = _width;
  height = _height;
#ifdef XWINDOWS  
  theDisplay = XOpenDisplay("");
  /*
  theWindow = XCreateSimpleWindow(theDisplay,RootWindow(theDisplay,0),50,50,
				  width,height,4,
				  WhitePixel(theDisplay,0),
				  BlackPixel(theDisplay,0));
				  */
  theWindow = XCreateWindow(theDisplay,RootWindow(theDisplay,0),50,50,
			    width,height,4,
			    CopyFromParent, CopyFromParent, 
			    CopyFromParent, 0, 0);
  
  XGCValues values;
  theGC = XCreateGC(theDisplay, theWindow, 0, &values);
  cmap = XCreateColormap(theDisplay,theWindow,
				  DefaultVisual(theDisplay,0),AllocAll);
  XSetWindowColormap(theDisplay, theWindow, cmap);

  XColor Color;
  Color.flags = DoRed | DoGreen | DoBlue;
  for (int red = 0; red < 8; red++)
    {
      for (int green = 0; green < 8; green++)
	{
	  for (int blue = 0; blue < 4; blue++)
	    {
	      Color.pixel = 32*red + 4*green + blue;
	      Color.red = (65536/7)*red;
	      Color.green = (65536/7)*green;
	      Color.blue = (65536/3)*blue;
	      XStoreColor(theDisplay, cmap, &Color);
	    }
	}
    }

	      
  XMapWindow(theDisplay,theWindow);
  XFlush(theDisplay);
  
#else
  prefsize(width,height);
  foreground();
  window_id = winopen("Test");
  RGBmode();
  shademodel(GOURAUD);
  gconfig();
  CLEAR();
#endif
}

void

Device::CLEAR()
{
cerr << "Device::CLEAR()    \n";
#ifdef XWINDOWS
  XSetBackground(theDisplay, theGC, 0);  
  XSetForeground(theDisplay, theGC, 0);  
  XClearWindow(theDisplay, theWindow);
  XFlush(theDisplay);
#else
  if (window_id)
    {
      winset(window_id);
      RGBcolor(0,0,0);
      clear();
    }
#endif
}

void
Device::C3F(float r, float g, float b)
{
#ifdef XWINDOWS
  int index = (((int)((r * 7.0)+.5))*32 + ((int)((g * 7.0)+.5))*4 +
	       ((int)((b * 3.0)+.5)));  
  XSetForeground(theDisplay, theGC, index);

#else
  float argToGl[3];
  argToGl[0] = r;
  argToGl[1] = g;
  argToGl[2] = b;
  c3f(argToGl);
#endif
}


void
Device::V2F(float x, float y)
{
#ifdef XWINDOWS
  // Flip the Y-Coordinate because X goes from 0->height as we 
  // go top->bottom while GL goes from height->0 as we go top->bottom. 
  y = height-y;	
  if (drawingPolygon)
  {
    if (numPoints == MAX_NUM_POINTS_FOR_POLYGON)
      {
	cerr << "ERROR: Maximum number of points has been exceeded" << endl;
	return;
      }
    polygonPointArray[numPoints].x = (int)x;
    polygonPointArray[numPoints].y = (int)y;
    numPoints++;
  }
  else
  {
    XDrawPoint(theDisplay, theWindow, theGC, (int) x, (int) y);
  }
#else
  float argToGl[2];
  argToGl[0] = x;
  argToGl[1] = y;
  v2f(argToGl);
#endif
}

void
Device::ENDIMAGE()
{
#ifdef XWINDOWS
  // Copy the image from our internal 
  // buffer (theImage) onto the display.
  XFlush(theDisplay);		// Update the XServer
#else
				// GL does not need to do anything 
				// after an image is ready to draw
				// as the image is drawn "on the fly".
#endif
}
  



void 
Device::BGNPOLYGON()
{
#ifdef XWINDOWS
  numPoints = 0;
  drawingPolygon = 1;
#else
  bgnpolygon();
#endif
}

void 
Device::ENDPOLYGON()
{
#ifdef XWINDOWS
  drawingPolygon = 0;
  // Draw the polygon with the GCs color
  cerr << " NUMPOINTS: "<< numPoints << endl;

  XFillPolygon(theDisplay, theWindow, theGC,		
	       polygonPointArray, numPoints, Complex, CoordModeOrigin);

#else
  endpolygon();
#endif
}

void 
Device::BGNCLOSEDLINE()
{
#ifdef XWINDOWS
  numPoints = 0;
  drawingPolygon = 1;
#else
  bgnclosedline();
#endif
}

void 
Device::ENDCLOSEDLINE()
{
#ifdef XWINDOWS
  drawingPolygon = 0;
  // Draw the polygon with the GCs color
  polygonPointArray[numPoints] = polygonPointArray[0]; // Close the loop
  XDrawLines(theDisplay, theWindow, theGC,		
	       polygonPointArray, numPoints+1, CoordModeOrigin);
#else
  endclosedline();
#endif
}

void
Device::BGNPOINT()
{
#ifdef XWINDOWS
#else
  bgnpoint();
#endif
}

void
Device::ENDPOINT()
{
#ifdef XWINDOWS
#else
  endpoint();
#endif
}
  
int
Device::GetWidth()
{
#ifdef XWINDOWS
  return width;
#else
  winset(window_id);
  long width, height;
  getsize(&width,&height);
  return width;
#endif
}

int
Device::GetHeight()
{
#ifdef XWINDOWS
  return height;
#else
  winset(window_id);
  long width, height;
  getsize(&width,&height);
  return height;
#endif
}
  

