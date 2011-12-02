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
// $Date: 2003-02-18 23:38:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/X11Device.cpp,v $
                                                                        
                                                                        
#include "X11Device.h"
#include <OPS_Globals.h>
#include <stdlib.h>


int X11Device::numX11Device(0);
Display *X11Device::theDisplay;  
Colormap X11Device::cmap;        
int X11Device::theScreen;        
unsigned long X11Device::pixels[X11_MAX_COLORS];
XColor X11Device::colors[X11_MAX_COLORS];
int X11Device::colorFlag;  


X11Device::X11Device()
:winOpen(1)
{
    hints.x = 50;
    hints.y = 50;
    hints.width = 0; 
    hints.height = 0;
    drawingPolygon = 0;
    numPoints = 0;

    // call the initX11 method if this is the first object
    if (numX11Device == 0) {
	this->initX11();
	numX11Device++;
    }
}

X11Device::~X11Device()
{
    numX11Device--;
    if (winOpen == 0) { // we must close the old window
	XFreeGC(theDisplay, theGC);
	XDestroyWindow(theDisplay, theWindow); 
    }
    
    if (numX11Device == 0) {
	if (colorFlag == 0) 
	    XFreeColors(theDisplay, cmap, pixels, 256, 0);
	else if (colorFlag == 1)
	    XFreeColors(theDisplay, cmap, pixels, 192, 0);
	else if (colorFlag == 2)
	    XFreeColors(theDisplay, cmap, pixels, 64, 0);
	
	XFreeColormap(theDisplay, cmap);
	XCloseDisplay(theDisplay);
    }
}

void
X11Device::WINOPEN(int _width, int _height)
{
    if (winOpen == 0) { // we must close the old window
	XFreeGC(theDisplay, theGC);
	XDestroyWindow(theDisplay, theWindow); 
    }

    // define the position and size of the window - only hints
    hints.x = 50;
    hints.y = 50;
    hints.width = _width;
    hints.height = _height;
    height = _height;

    // now open a window
    theWindow = XCreateSimpleWindow(theDisplay,RootWindow(theDisplay,0),
				    hints.x, hints.y,
				    hints.width,hints.height,4,
				    foreground, background);
    if (theWindow == 0) {
	opserr << "X11Device::WINOPEN() - could not open a window\n";
	exit(-1);
    }	
    
    XSetStandardProperties(theDisplay, theWindow, "g3", "g3", None, 0, 0, &hints);
    
    // create a graphical context
    theGC = XCreateGC(theDisplay, theWindow, 0, 0);
    XSetBackground(theDisplay, theGC, background);
    XSetForeground(theDisplay, theGC, foreground);

    if (colorFlag == 3) {
	cmap = XCreateColormap(theDisplay,theWindow,
			   DefaultVisual(theDisplay,0),AllocAll);
	if (cmap == 0) {
	    opserr << "X11Device::initX11() - could not get a new color table\n";
	    exit(-1);
	}	    

	// we are going to try to allocate 256 new colors -- need 8 planes for this
	int depth = DefaultDepth(theDisplay, theScreen);
	if (depth < 8) {
	    opserr << "X11Device::initX11() - needed at least 8 planes\n";
	    exit(-1);
	}	    
	int cnt = 0;
	for (int red = 0; red < 8; red++) {
	    for (int green = 0; green < 8; green++) {
		for (int blue = 0; blue < 4; blue++) {
		    colors[cnt].pixel = pixels[32*red + 4*green + blue];
		    colors[cnt].red = (65536/7)*red;
		    colors[cnt].green = (65536/7)*green;
		    colors[cnt].blue = (65536/3)*blue;
		    colors[cnt].flags = DoRed | DoGreen | DoBlue;
		    cnt++;
		}			
	    }
	}
	XStoreColors(theDisplay, cmap, colors, cnt);			    

	XSetWindowColormap(theDisplay, theWindow, cmap);    
    }

    XMapWindow(theDisplay,theWindow);
    XClearWindow(theDisplay, theWindow);      
    XFlush(theDisplay);
    
    winOpen = 0;
}

void

X11Device::CLEAR()
{
  XSetBackground(theDisplay, theGC, background);
  XClearWindow(theDisplay, theWindow);  
  XFlush(theDisplay);
}

void
X11Device::C3F(float r, float g, float b)
{
    int index, val;
    // check range of rgb values
    if (r<0 || r>1.0 || g<0 || g>1.0 || b<0 || b>1.0) {
	opserr << "X11Device::X11Device::C3F() rgb val out of range ";
	opserr << r << " " << g << " " << b << endln;
	return;
    }
    
    if (colorFlag == 0 || colorFlag == 3) {
	val = (((int)((r * 7.0)+.5))*32 + ((int)((g * 7.0)+.5))*4 +
		 ((int)((b * 3.0)+.5)));  	
    } else if (colorFlag == 1) {
	val = (((int)((r * 7.0)+.5))*24 + ((int)((g * 5.0)+.5))*4 +
		 ((int)((b * 3.0)+.5)));  
    } else if (colorFlag == 2) {
	val = ((int)((r * 3.0)+.5))*16 + ((int)((g * 3.0)+.5))*4 +
	    ((int)((b * 3.0)+.5));
    } else
      val = 0; // should never be called

    index = pixels[val];
    XSetForeground(theDisplay, theGC, index);
}


void
X11Device::V2F(float x, float y)
{
  // Flip the Y-Coordinate because X goes from 0->height as we 
  // go top->bottom while GL goes from height->0 as we go top->bottom. 
  y = height-y;	
  if (drawingPolygon)
  {
    if (numPoints == MAX_NUM_POINTS_FOR_POLYGON)
      {
	opserr << "ERROR: Maximum number of points has been exceeded" << endln;
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
}

void
X11Device::ENDIMAGE()
{
  // Copy the image from our internal 
  // buffer (theImage) onto the display.
  XFlush(theDisplay);		// Update the XServer
}


void
X11Device::STARTIMAGE()
{
  // Copy the image from our internal 
  // buffer (theImage) onto the display.
  XFlush(theDisplay);		// Update the XServer
}
  



void 
X11Device::BGNPOLYGON()
{
  numPoints = 0;
  drawingPolygon = 1;
}

void 
X11Device::ENDPOLYGON()
{
  drawingPolygon = 0;
  // Draw the polygon with the GCs color
  opserr << " NUMPOINTS: "<< numPoints << endln;

  XFillPolygon(theDisplay, theWindow, theGC,		
	       polygonPointArray, numPoints, Complex, CoordModeOrigin);

}

void 
X11Device::BGNCLOSEDLINE()
{
  numPoints = 0;
  drawingPolygon = 1;
}

void 
X11Device::ENDCLOSEDLINE()
{
  drawingPolygon = 0;
  // Draw the polygon with the GCs color
  polygonPointArray[numPoints] = polygonPointArray[0]; // Close the loop
  XDrawLines(theDisplay, theWindow, theGC,		
	       polygonPointArray, numPoints+1, CoordModeOrigin);
}

void
X11Device::BGNPOINT()
{

}

void
X11Device::ENDPOINT()
{

}




void
X11Device::drawText(float x, float y, char *text, int length)
{
  y = height-y;	
  XDrawString(theDisplay, theWindow, theGC, (int) x, (int) y, text, length);
}


int
X11Device::GetWidth()
{
  int x,y;
  unsigned int width, h, borderWidth, depth;
  
  XGetGeometry(theDisplay, theWindow, &RootWindow(theDisplay,0),
	       &x, &y, &width, &h, &borderWidth, &depth);

  hints.width = width;
  hints.height = h;
  height = h;
  return width;
}

int
X11Device::GetHeight()
{
  int x,y;
  unsigned int width, h, borderWidth, depth;


  XGetGeometry(theDisplay, theWindow, &RootWindow(theDisplay,0),
	       &x, &y, &width, &h, &borderWidth, &depth);

  hints.width = width;
  hints.height = h;
  height = h;
  return height;    
}
  




void
X11Device::initX11(void) {
    // set the display and screen variables
    theDisplay = XOpenDisplay("");	// init a display connection
    if (theDisplay == 0) {              // and check we got one
	opserr << "X11Device::initX11() - could not connect to display\n";
	exit(-1);
    }

    theScreen = DefaultScreen(theDisplay);
    
    // set the defualt foreground and background colors
    foreground = BlackPixel(theDisplay, theScreen);
    background = WhitePixel(theDisplay, theScreen);    

    // lets try using the default colormap
    cmap = DefaultColormap(theDisplay, theScreen);
    XVisualInfo vinfo;    
    if (XMatchVisualInfo(theDisplay, theScreen, 8, PseudoColor, &vinfo) == false)
	opserr << "X11Device::initX11 - could not get info\n";

    // we now try to allocate some color cells from the colormap
    // we start by tring to obtain 256 colors, then 192, finally 64
    // if we can't get these (and as a last resort) we create a new color map

    if (XAllocColorCells(theDisplay, cmap, false, NULL, 0, pixels, 256) != 0) {
	// we were able to allocate 256 colors from the table for our use
	colorFlag = 0;
	int cnt = 0;
	for (int red =0; red <8; red++) {
	    for (int green = 0; green<8; green++) {
		for (int blue =0; blue<4; blue++) {
		    colors[cnt].pixel = pixels[32*red + 4*green + blue];
		    colors[cnt].red = (65536/7)*red;
		    colors[cnt].green = (65536/7)*green;
		    colors[cnt].blue = (65536/3)*blue;
		    colors[cnt].flags = DoRed | DoGreen | DoBlue;
		    cnt++;
		}			
	    }
	}
	XStoreColors(theDisplay, cmap, colors, cnt);
	
    } else if (XAllocColorCells(theDisplay, cmap, false, NULL, 0, pixels, 192) != 0) {
	// we were able to allocate 192 colors from the table for our use	
	colorFlag = 1;	
	int cnt = 0;
	for (int red =0; red <8; red++) {
	    for (int green = 0; green<6; green++) {
		for (int blue =0; blue<4; blue++) {
		    colors[cnt].pixel = pixels[24*red + 4*green + blue];
		    colors[cnt].red = (65536/7)*red;
		    colors[cnt].green = (65536/5)*green;
		    colors[cnt].blue = (65536/3)*blue;
		    colors[cnt].flags = DoRed | DoGreen | DoBlue;
		    cnt++;
		}			
	    }
	}
	XStoreColors(theDisplay, cmap, colors, cnt);
    } else if (XAllocColorCells(theDisplay, cmap, false, NULL, 0, pixels, 64) != 0) {
	colorFlag = 2;	
	int cnt = 0;
	for (int red =0; red <4; red++) {
	    for (int green = 0; green<4; green++) {
		for (int blue =0; blue<4; blue++) {
		    colors[cnt].pixel = pixels[16*red + 4*green + blue];
		    colors[cnt].red = (65536/3)*red;
		    colors[cnt].green = (65536/3)*green;
		    colors[cnt].blue = (65536/3)*blue;
		    colors[cnt].flags = DoRed | DoGreen | DoBlue;
		    cnt++;
		}			
	    }
	}
	XStoreColors(theDisplay, cmap, colors, cnt);
    } else {
	colorFlag = 3;
	// lets create our own color table - 
	// problem with this is that screen colors change as we enter
	opserr << "X11Device::initX11() - could not any colors from the default\n";
	opserr << "colormap - have to create our own for the app - windows will\n";
	opserr << "windows will change color as move mouse from one window to another\n\n";	
    }
}    
	
    

	
