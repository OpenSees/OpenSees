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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-04-15 18:01:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/WindowDevice.cpp,v $
                                                                        
                                                                        
                                                                        
#include "WindowDevice.h"
#include <OPS_Globals.h>
#include <stdlib.h>
#include <string.h>

int WindowDevice::numWindowDevice(0);

#ifdef _UNIX
Display *WindowDevice::theDisplay;  
Colormap WindowDevice::cmap;        
int WindowDevice::theScreen;        
unsigned long WindowDevice::pixels[X11_MAX_COLORS];
XColor WindowDevice::colors[X11_MAX_COLORS];
int WindowDevice::colorFlag; 
unsigned long WindowDevice::foreground(0);
unsigned long WindowDevice::background(0);
 

#else
/* WindowProc()
 *  Minimum Window Procedure
 */
LONG WINAPI WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{ 
    LONG        lRet = 1;
    PAINTSTRUCT ps;

    switch(uMsg) {
    case WM_CREATE:
        break; 

    case WM_DESTROY:
        break; 

    case WM_PAINT: 
        BeginPaint(hWnd, &ps); 
        EndPaint(hWnd, &ps); 
        break; 

    default: 
        lRet = DefWindowProc (hWnd, uMsg, wParam, lParam); 
        break; 
    }

    return lRet;
}


/* oglPixelFormat()
 *  Sets the pixel format for the context
 */
int oglSetPixelFormat(HDC hDC, BYTE type, DWORD flags)
{
    int pf;
    PIXELFORMATDESCRIPTOR pfd;

    /* fill in the pixel format descriptor */
    pfd.nSize        = sizeof(PIXELFORMATDESCRIPTOR);
    pfd.nVersion     = 1;		    /* version (should be 1) */
    pfd.dwFlags      = flags | /* draw to window (not bitmap) */
                       PFD_SUPPORT_OPENGL;  /* draw using opengl */
    pfd.iPixelType   = type;                /* PFD_TYPE_RGBA or COLORINDEX */
    pfd.cColorBits   = 24;
    pfd.cRedBits = 8;
    pfd.cGreenBits = 8;
    pfd.cBlueBits = 8;
    pfd.cDepthBits = 16;
    /* other criteria here */
    
    /* get the appropriate pixel format */
    pf = ChoosePixelFormat(hDC, &pfd);
    if (pf == 0) {
       MessageBox(NULL,
		  "ChoosePixelFormat() failed:  Cannot find format specified.",
		  "Error", MB_OK); 
       return 0;
    } 
 
    /* set the pixel format */
    if (SetPixelFormat(hDC, pf, &pfd) == FALSE) {
	MessageBox(NULL,
		   "SetPixelFormat() failed:  Cannot set format specified.",
		   "Error", MB_OK);
        return 0;
    } 

    return pf;
}    


/* oglCreateWindow
 *  Create a window suitable for OpenGL rendering
 */
HWND oglCreateWindow(char* title, int x, int y, int width, int height,
					 HGLRC *hRC, HDC *hDC)
{
    WNDCLASS  wc;
    HWND hWnd;
    HINSTANCE hInstance;

    /* get this modules instance */
    hInstance = GetModuleHandle(NULL);

    /* fill in the window class structure */
    wc.style         = CS_HREDRAW | CS_VREDRAW;  // to redraw if moved
    wc.lpfnWndProc   = (WNDPROC)WndProc;         /* event handler */
    wc.cbClsExtra    = 0;                           /* no extra class data */
    wc.cbWndExtra    = 0;                           /* no extra window data */
    wc.hInstance     = hInstance;                   /* instance */
    wc.hIcon         = LoadIcon(NULL, IDI_WINLOGO); /* load a default icon */
    wc.hCursor       = LoadCursor(NULL, IDC_ARROW); /* load a default cursor */
    wc.hbrBackground = NULL;                        /* redraw our own bg */
    wc.lpszMenuName  = NULL;                        /* no menu */
    wc.lpszClassName = title;                       /* use a special class */

    /* register the window class */
    if (!RegisterClass(&wc)) {
      MessageBox(NULL, 
		   "RegisterClass() failed:  Cannot register window class,",
		   "Error", MB_OK);
	return NULL;
    }

    /* create a window */
    hWnd = CreateWindow(title,          /* class */
			title,          /* title (caption) */
			WS_CLIPSIBLINGS | WS_CLIPCHILDREN,  /* style */
			x, y, width, height, /* dimensions */
			NULL,		/* no parent */
			NULL,		/* no menu */
			hInstance,	/* instance */
			NULL);		/* don't pass anything to WM_CREATE */

    /* make sure we got a window */
    if (hWnd == NULL) {
	MessageBox(NULL,
		   "CreateWindow() failed:  Cannot create a window.",
		   "Error", MB_OK);
	return NULL;
    }

    /* show the window (map it) */
    ShowWindow(hWnd, SW_SHOW);

    /* send an initial WM_PAINT message (expose) */
    UpdateWindow(hWnd);

	    /* get the device context */
    *hDC = GetDC(hWnd);
	
    /* set the pixel format */
    if (oglSetPixelFormat(*hDC, PFD_TYPE_RGBA, 
				PFD_DRAW_TO_WINDOW | PFD_DOUBLEBUFFER) == 0)
      exit(1);

    /* create an OpenGL context */
    *hRC = wglCreateContext(*hDC);
	wglMakeCurrent(*hDC, *hRC);
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0, 0, (GLsizei)width, (GLsizei)height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    gluOrtho2D(0.0, (GLdouble)width, 0.0, (GLdouble)height);
    glFlush();

   
    return hWnd;
}


/* oglCreateWindow
 *  Create a window suitable for OpenGL rendering
 */
int oglDestroyWindow(char* title, HWND hWnd, HGLRC hRC, HDC hDC)
{
    HINSTANCE hInstance;

    /* get this modules instance */
    hInstance = GetModuleHandle(NULL);

    /*
	 * now release the device context, destroy the rendering context
	 * and destroy the window.
	 */
	wglMakeCurrent(NULL, NULL);	//make the gl context 'un-'current 
    ReleaseDC(hWnd, hDC);		//release handle to DC 
    wglDeleteContext(hRC);		//delete the rendering context 
    DestroyWindow(hWnd);	    // destroy the window
	
	/* unregister the window class - so can use window name again*/    
	if (!UnregisterClass(title, hInstance)) {
      MessageBox(NULL, 
		   "UnregisterClass() failed:  Cannot unregister window class,",
		   "Error", MB_OK);
	return -1;
    }
	
	return 0;
}


/* oglCreateWindow
 *  Create a window suitable for OpenGL rendering
 */
int oglCreateBitmap(int width, int height, HGLRC *hRC, HDC *hDC, 
						HBITMAP *theBitmap, BITMAPINFO *info, GLubyte **bits)
{

 	*hDC = CreateCompatibleDC(NULL);


//	memset(&info, 0, sizeof(BITMAPINFO));
	
	GLint Width = width;
	GLint Height = height;

	info->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	info->bmiHeader.biPlanes = 1;
	info->bmiHeader.biBitCount = 24;
	info->bmiHeader.biCompression = BI_RGB;
	info->bmiHeader.biXPelsPerMeter = 11808;
	info->bmiHeader.biYPelsPerMeter = 11808;
	info->bmiHeader.biClrUsed = 0;
	info->bmiHeader.biClrImportant = 0;
	info->bmiHeader.biWidth = width;
	info->bmiHeader.biHeight = height;

	/*
	if ((*bits = (GLubyte *)(calloc(width*height, 1))) == 0){
		opserr << "BITS ZERO\n";
		return -1;
	}
	*/
	void *theBits = *bits;
//	void **theBitsPtr = &theBits;

	*theBitmap = CreateDIBSection(*hDC, info, DIB_RGB_COLORS, &theBits, NULL, 0);
	*bits = (GLubyte *)theBits;

	SelectObject(*hDC, *theBitmap);
   /* set the pixel format */
    if (oglSetPixelFormat(*hDC, PFD_TYPE_RGBA, PFD_DRAW_TO_BITMAP) == 0)
      exit(1);

    /* create an OpenGL context */
    *hRC = wglCreateContext(*hDC);
  
    return 0;
}


/* oglCreateWindow
 *  Create a window suitable for OpenGL rendering
 */
int oglDestroyBitmap(HBITMAP *theBitmap, HGLRC hRC, HDC hDC)
{
    /*
	 * now release the device context, destroy the rendering context
	 * and destroy the bitmap.
	 */
	wglMakeCurrent(NULL, NULL);	//make the gl context 'un-'current 
    wglDeleteContext(hRC);		//delete the rendering context 
    DeleteObject(theBitmap);// destroy the window
	

	return 0;
}

#endif


WindowDevice::WindowDevice()
:winOpen(1),height(0),width(0),numPoints(0),drawingPolygon(0)
{
#ifdef _UNIX
    hints.x = 50;
    hints.y = 50;
    hints.width = 0; 
    hints.height = 0;
#else

#endif
    // call the initX11 method if this is the first object
    if (numWindowDevice == 0) {
	this->initWindow();
    }
    numWindowDevice++;
}

WindowDevice::~WindowDevice()
{
   numWindowDevice--;
#ifdef _UNIX
    if (winOpen == 0) { // we must close the old window
	XFreeGC(theDisplay, theGC);
	XDestroyWindow(theDisplay, theWindow); 
    }

    if (numWindowDevice == 0) {
	if (colorFlag == 0 ) 
	  XFreeColors(theDisplay, cmap, pixels, 256, 0);
	else if (colorFlag == 1)
	    XFreeColors(theDisplay, cmap, pixels, 192, 0);
	else if (colorFlag == 2)
	    XFreeColors(theDisplay, cmap, pixels, 64, 0);
	
       	XFreeColormap(theDisplay, cmap);
	XCloseDisplay(theDisplay);
    }

#else
    if (winOpen == 0) { // we must close the window
      oglDestroyWindow(title,theWND, theHRC, theHDC);
    }
#endif
}

void
WindowDevice::WINOPEN(const char *_title, int _xLoc, int _yLoc, int _width, int _height)
{
    // set the WindowDevices title, height, wdth, xLoc and yLoc
    strcpy(title, _title);

    height = _height;
    width = _width;  
    xLoc = _xLoc;
    yLoc = _yLoc;

#ifdef _UNIX
    if (winOpen == 0) { // we must close the old window
	XFreeGC(theDisplay, theGC);
	XDestroyWindow(theDisplay, theWindow); 
    }

    // define the position and size of the window - only hints
    hints.x = _xLoc;
    hints.y = _yLoc;
    hints.width = _width;
    hints.height = _height;
    hints.flags = PPosition | PSize;

    // set the defualt foreground and background colors
    XVisualInfo visual; 
    visual.visual = 0;
    int depth = DefaultDepth(theDisplay, theScreen);

    if (background == 0) {
      if (XMatchVisualInfo(theDisplay, theScreen, depth, PseudoColor, &visual) == 0) {
	  foreground = BlackPixel(theDisplay, theScreen);
	  background = WhitePixel(theDisplay, theScreen);    

      } else {
	foreground = 0;
	background = 255;
      }
    }

    // now open a window
    theWindow = XCreateSimpleWindow(theDisplay,RootWindow(theDisplay,0),
				    hints.x, hints.y,
				    hints.width,hints.height,4,
				    foreground, background);

    if (theWindow == 0) {
	opserr << "WindowDevice::WINOPEN() - could not open a window\n";
	exit(-1);
    }	
    
    XSetStandardProperties(theDisplay, theWindow, title, title, None, 0, 0, &hints);
    
    // create a graphical context
    theGC = XCreateGC(theDisplay, theWindow, 0, 0);

    // if we were unable to get space for our colors
    // we must create and use our own colormap
    if (colorFlag == 3 ) {

      // create the colormap if the 1st window
      if (numWindowDevice == 1) {
	int fail = false;
	//	XMatchVisualInfo(theDisplay, theScreen, depth, PseudoColor, &visual);
	if (XMatchVisualInfo(theDisplay, theScreen, depth, PseudoColor, &visual) == 0) {
	  opserr << "WindowDevice::initX11() - could not get a visual for PseudoColor\n";
	  opserr << "Colors diplayed will be all over the place\n";
	  cmap = DefaultColormap(theDisplay, theScreen);
	  fail = true;
        } else {
	  opserr << "WindowDevice::WINOPEN have created our own colormap, \n";
	  opserr << "windows may change color as move mouse from one window to\n";
	  opserr << "another - depends on your video card to use another colormap\n\n";	

	  cmap = XCreateColormap(theDisplay,theWindow,
				 visual.visual, AllocAll);
	}


	/*
	cmap = XCreateColormap(theDisplay,theWindow,
			   DefaultVisual(theDisplay,0),AllocAll);
	*/

	if (cmap == 0) {
	    opserr << "WindowDevice::initX11() - could not get a new color table\n";
	    exit(-1);
	}	    

	// we are going to try to allocate 256 new colors -- need 8 planes for this
	depth = DefaultDepth(theDisplay, theScreen);
	if (depth < 8) {
	    opserr << "WindowDevice::initX11() - needed at least 8 planes\n";
	    exit(-1);
	}	    
	if (fail == false) {
	  int cnt = 0;
	  for (int red = 0; red < 8; red++) {
	    for (int green = 0; green < 8; green++) {
		for (int blue = 0; blue < 4; blue++) {
		  pixels[32*red + 4*green + blue] = cnt;
		  colors[cnt].pixel = pixels[32*red + 4*green + blue];
		  colors[cnt].red = (65536/7)*red;
		  colors[cnt].green = (65536/7)*green;
		  colors[cnt].blue = (65536/3)*blue;
		  colors[cnt].flags = DoRed | DoGreen | DoBlue;
		  cnt++;
		}			
	    }
	  }
	  background = 0; //pixels[0];
	  foreground = 255; // pixels[255];
	  XStoreColors(theDisplay, cmap, colors, cnt);			    
	}
      }

      // now set the windows to use the colormap
      XSetWindowColormap(theDisplay, theWindow, cmap);    
    
    }

    XSetBackground(theDisplay, theGC, background);
    XSetForeground(theDisplay, theGC, foreground);

    XMapWindow(theDisplay,theWindow);
    XClearWindow(theDisplay, theWindow);      
    XFlush(theDisplay);

#else
    //    auxInitDisplayMode(AUX_SINGLE | AUX_RGBA);
    //    auxInitPosition(100,100,_width,_height);
    //    auxInitWindow("G3");

    if (winOpen == 0)
      oglDestroyWindow(title,theWND, theHRC, theHDC);      

    theWND = oglCreateWindow(title, xLoc, yLoc, width, height, &theHRC, &theHDC);
    if (theWND == NULL)
      exit(1);
    winOpen = 0;

    wglMakeCurrent(theHDC, theHRC);
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0, 0, (GLsizei)width, (GLsizei)height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    gluOrtho2D(0.0, (GLdouble)width, 0.0, (GLdouble)height);
    glFlush();

#endif

    winOpen = 0;
}

void

WindowDevice::CLEAR()
{
#ifdef _UNIX
  XSetBackground(theDisplay, theGC, background);
  XClearWindow(theDisplay, theWindow);  
  XFlush(theDisplay);
#else
   wglMakeCurrent(theHDC, theHRC);
   glClearColor(1.0f,1.0f,1.0f,1.0f);
   glClear(GL_COLOR_BUFFER_BIT);
 
   glFlush();
#endif
}

void
WindowDevice::C3F(float r, float g, float b)
{
    
    // check range of rgb values
    if (r<0 || r>1.0 || g<0 || g>1.0 || b<0 || b>1.0) {
	opserr << "WindowDevice::WindowDevice::C3F() rgb val out of range ";
	opserr << r << " " << g << " " << b << endln;
	return;
    }

#ifdef _UNIX 
    int index, val;
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
      val = 0;
    
    index = pixels[val];
    XSetForeground(theDisplay, theGC, index);
#else
	
	glColor3f(r,g,b);

 
#endif
}


void
WindowDevice::V2F(float x, float y)
{
#ifdef _UNIX
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
#else
		glVertex2f(x,y);
#endif
}

void
WindowDevice::ENDIMAGE()
{
#ifdef _UNIX
  // Copy the image from our internal 
  // buffer (theImage) onto the display.
  XFlush(theDisplay);		// Update the XServer

#else
				// GL does not need to do anything 
				// after an image is ready to draw
				// as the image is drawn "on the fly".
  SwapBuffers(theHDC);
  glFlush();
#endif
}

void
WindowDevice::STARTIMAGE()
{
#ifdef _UNIX
  // Copy the image from our internal 
  // buffer (theImage) onto the display.
  XFlush(theDisplay);		// Update the XServer

#else
   
  // ensure that the gl contex is the one for this object & do a flush
  wglMakeCurrent(theHDC, theHRC);
 // glFlush();
#endif
}
  



void 
WindowDevice::BGNPOLYGON()
{
  numPoints = 0;
  drawingPolygon = 1;
}

void 
WindowDevice::ENDPOLYGON()
{
  drawingPolygon = 0;
  // Draw the polygon with the GCs color
  #ifdef _UNIX
  XFillPolygon(theDisplay, theWindow, theGC,		
	       polygonPointArray, numPoints, Complex, CoordModeOrigin);
  #else

  #endif
}

void 
WindowDevice::BGNCLOSEDLINE()
{
  numPoints = 0;
  drawingPolygon = 1;
#ifdef _UNIX

#else
  glBegin(GL_LINES);
#endif
}

void 
WindowDevice::ENDCLOSEDLINE()
{
  drawingPolygon = 0;
  // Draw the polygon with the GCs color
  
#ifdef _UNIX
  polygonPointArray[numPoints] = polygonPointArray[0]; // Close the loop
  XDrawLines(theDisplay, theWindow, theGC,		
	       polygonPointArray, numPoints+1, CoordModeOrigin);
#else
	glEnd();
#endif
}

void
WindowDevice::BGNPOINT()
{

}

void
WindowDevice::ENDPOINT()
{

}




void
WindowDevice::drawText(float x, float y, char *text, int length)
{
#ifdef _UNIX
  y = height-y;	
  XDrawString(theDisplay, theWindow, theGC, (int) x, (int) y, text, length);
#else

#endif
}


int
WindowDevice::GetWidth()
{
  

#ifdef _UNIX
  int x,y;
  unsigned int borderWidth, depth;
  unsigned int w, h;
  XGetGeometry(theDisplay, theWindow, &RootWindow(theDisplay,0),
	       &x, &y, &w, &h, &borderWidth, &depth);
  width = w;
  height = h;
  hints.width = width;
  hints.height = h;
  height = h;
#else

#endif

  return width;
}

int
WindowDevice::GetHeight()
{
#ifdef _UNIX
  unsigned int borderWidth, depth;
  int x,y;
  unsigned int w, h;

  XGetGeometry(theDisplay, theWindow, &RootWindow(theDisplay,0),
	       &x, &y, &w, &h, &borderWidth, &depth);
  width = w;
  height = h;
  hints.width = width;
  hints.height = height;
#else

#endif

  return height;    
}
  




void
WindowDevice::initWindow(void) {
    // set the display and screen variables
#ifdef _UNIX
    theDisplay = XOpenDisplay("");	// init a display connection
    if (theDisplay == 0) {              // and check we got one
	opserr << "WindowDevice::initX11() - could not connect to display\n";
	exit(-1);
    }

    theScreen = DefaultScreen(theDisplay);
    
    // set the defualt foreground and background colors
    //    foreground = BlackPixel(theDisplay, theScreen);
    // background = WhitePixel(theDisplay, theScreen);    

    // lets try using the default colormap
    cmap = DefaultColormap(theDisplay, theScreen);

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
	foreground = pixels[0];
	background = pixels[255];
	
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
	foreground = pixels[0];
	background = pixels[191];
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
	foreground = pixels[0];
	background = pixels[63];
    } else {
	colorFlag = 3;
	// lets create our own color table - 
	// problem with this is that screen colors change as we enter
	opserr << "WindowDevice::initWindow() - could not add any colors to the\n";
	opserr << "existing colormap - will try to create our own colormap\n";
    }
#else


#endif
}    
	
    

	
