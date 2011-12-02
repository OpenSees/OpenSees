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
                                                                        
// $Revision: 1.14 $
// $Date: 2004-06-07 23:09:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGlDevice.cpp,v $
                                                                        
                                                                        
#include <OpenGlDevice.h>
#include <OPS_Globals.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _GLX
#define _PNG
#include <png.h>
#endif


int OpenGlDevice::numWindows(0);
// GLuint OpenGlDevice::FontBase(0);


#ifdef _WGL
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, (GLsizei)width, (GLsizei)height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, (GLdouble)width, 0.0, (GLdouble)height);
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

#elif _GLX

static int attributeListSgl[] = {
	GLX_RGBA,
	GLX_RED_SIZE,   1, /*get the deepest buffer with 1 red bit*/
	GLX_GREEN_SIZE, 1,
	GLX_BLUE_SIZE,  1,
	None };

static int attributeListDbl[] = {
	GLX_RGBA,
	GLX_DOUBLEBUFFER, /*In case single buffering is not supported*/
	GLX_RED_SIZE,   1,
	GLX_GREEN_SIZE, 1,
	GLX_BLUE_SIZE,  1,
	None };

static Bool WaitForNotify(Display *d, XEvent *e, char *arg) {
	return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
}

static const char *FontName = "fixed";

//XFontStruct *OpenGlDevice::fontInfo(0);

#else

#endif










OpenGlDevice::OpenGlDevice()
  :FontBase(0), winOpen(1), width(0), height(0), windowTitle(0)
{

#ifdef _WGL

#elif _GLX
  fontInfo = 0;
#else

#endif

  // if this is the first window to open call init windows
  if (numWindows == 0)
    this->initWindow();

  numWindows++;
}

OpenGlDevice::~OpenGlDevice()
{
  numWindows--;

#ifdef _WGL
  if (winOpen == 0) { // we must close the window
      oglDestroyWindow(windowTitle,theWND, theHRC, theHDC);
  }

  if (FontBase != 0) {
    glDeleteLists(FontBase, 256);
  }
#elif _GLX

  if (fontInfo != 0) {
    Font id;
    unsigned int first, last;    

    id = fontInfo->fid;
    first = fontInfo->min_char_or_byte2;
    last = fontInfo->max_char_or_byte2;
  
    glDeleteLists(FontBase, (GLuint) last + 1);
    XFreeFont(theDisplay, fontInfo);
  }

  if (winOpen == 0) { // we must close the old window
    XFreeGC(theDisplay, theGC);
    XDestroyWindow(theDisplay, theWindow); 
    XCloseDisplay(theDisplay);
  }

  if (windowTitle != 0)
    delete [] windowTitle;

#else

#endif

}

void
OpenGlDevice::WINOPEN(const char *_title, int _xLoc, int _yLoc, int _width, int _height)
{
  if (windowTitle != 0)
    delete [] windowTitle;
  
  windowTitle = new char[strlen(_title)+1];
  strcpy(windowTitle, _title);

  width = _width;
  height = _height;

#ifdef _WGL
  xLoc = _xLoc;
  yLoc = _yLoc;

  if (winOpen == 0)
    oglDestroyWindow(windowTitle,theWND, theHRC, theHDC);      
 
  theWND = oglCreateWindow(windowTitle, _xLoc, _yLoc, _width, _height, &theHRC, &theHDC);
  if (theWND == NULL)
    exit(1);
  winOpen = 0;
  wglMakeCurrent(theHDC, theHRC);

  FontBase = glGenLists((GLuint) 256);
  if (!FontBase) {
    printf("Error: unable to allocate display lists\n");
    exit(0);
  }

  wglUseFontBitmaps(theHDC, 0, 256, FontBase);

#elif _GLX

  theDisplay = XOpenDisplay("");      // init a display connection
  if (theDisplay == 0) {              // and check we got one
    opserr << "OpenGlDevice::initX11() - could not connect to display\n";
    exit(-1);
  }

  theScreen = DefaultScreen(theDisplay);

  //XVisualInfo *visual;
  XSetWindowAttributes swa;
  swap_flag = GL_FALSE;
  //Colormap cmap; //this overrides the other color map
  //GLXContext cx;
  //XEvent event;

  if (winOpen == 0) { // we must close the old window
    XFreeGC(theDisplay, theGC);
    XDestroyWindow(theDisplay, theWindow);
  }

  /* get an appropriate visual*/
  visual = glXChooseVisual(theDisplay, theScreen, attributeListSgl);
  if (visual == NULL) {
    visual = glXChooseVisual(theDisplay, theScreen, attributeListDbl);
    swap_flag = GL_TRUE;
  }

  if(visual == NULL) {
    opserr << "OpenGlDevice::WINOPEN - unable to get visual\n";
    exit(2);
  }

  /* create a color map */
  //  cmap = XCreateColormap(theDisplay, RootWindow(theDisplay, visual->screen),
  //		 visual->visual, AllocNone);
  
  /* create a window */
  swa.colormap = XCreateColormap(theDisplay, RootWindow(theDisplay, visual->screen),
				 visual->visual, AllocNone);

  swa.border_pixel = 0;
  swa.event_mask = StructureNotifyMask;

  unsigned long mask;
  mask =  CWBackPixel|CWBorderPixel|CWColormap|CWEventMask;
  theWindow = XCreateWindow(theDisplay, RootWindow(theDisplay, visual->screen), 
			    _xLoc, _yLoc, _width, _height,
			    0, visual->depth, InputOutput, visual->visual,
			    mask, &swa);

  if (theWindow == 0) {
    opserr << "OpenGlDevice::WINOPEN() - could not open a window\n";
    exit(-1);
  }
    
  // define the position and size of the window - only hints
  hints.x = _xLoc;
  hints.y = _yLoc;
  hints.width = _width;
  hints.height = _height;
  hints.flags = PPosition | PSize;
  XSetNormalHints(theDisplay, theWindow, &hints);
  XSetStandardProperties(theDisplay, theWindow, windowTitle, windowTitle, None, 0, 0, &hints);

  /* create a GLX context */
  cx = glXCreateContext(theDisplay, visual, NULL, GL_TRUE);
  if (cx == 0) {
    opserr << "OpenGlDevice::WINOPEN() - could not create a glx context\n";
    exit(-1);
  }    

  
  XMapWindow(theDisplay, theWindow);
    
  /* connect the context to the window */
  glXMakeCurrent(theDisplay, theWindow, cx);

  // create a graphical context
  theGC = XCreateGC(theDisplay, theWindow, 0, 0);
  winOpen = 0;

  XVisualInfo visual; 
  visual.visual = 0;
  int depth = DefaultDepth(theDisplay, theScreen);

  XMapWindow(theDisplay,theWindow);
  XClearWindow(theDisplay, theWindow);      
  XFlush(theDisplay);

  Font id;
  unsigned int first, last;

  fontInfo = XLoadQueryFont(theDisplay, FontName);
  if (!fontInfo) {
    printf("Error: font %s not found\n", FontName);
    exit(0);
  }

  id = fontInfo->fid;
  first = fontInfo->min_char_or_byte2;
  last = fontInfo->max_char_or_byte2;

  FontBase = glGenLists((GLuint) last + 1);
  if (!FontBase) {
    printf("Error: unable to allocate display lists\n");
    exit(0);
  }

  glXUseXFont(id, first, last - first + 1, FontBase + first);

#else

#endif


}




void

OpenGlDevice::CLEAR()
{
#ifdef _WGL
    wglMakeCurrent(theHDC, theHRC);

#elif _GLX
    glXMakeCurrent(theDisplay, theWindow, cx);
#else

#endif

}

void
OpenGlDevice::STARTIMAGE()
{
#ifdef _WGL
    wglMakeCurrent(theHDC, theHRC);
#elif _GLX
    glXMakeCurrent(theDisplay, theWindow, cx);
#else

#endif
}

void
OpenGlDevice::ENDIMAGE()
{
#ifdef _WGL
    wglMakeCurrent(theHDC, theHRC);
    SwapBuffers(theHDC);
#elif _GLX
    if (swap_flag == GL_TRUE)
      glXSwapBuffers(theDisplay, theWindow);
    glXMakeCurrent(theDisplay, theWindow, cx);
#else

#endif
}


int
OpenGlDevice::GetWidth()
{
  return width;
}

int
OpenGlDevice::GetHeight()
{
  return height;
}
  


void
OpenGlDevice::initWindow(void) {
  // set the display and screen variables

#ifdef _WGL

#elif _GLX

#else

#endif
}    
	
    
void
OpenGlDevice::drawText(float x, float y, float z, char *text, int length,
		       char horizontalJustify, char verticalJustify)
{

#ifdef _WGL
  glColor3f(0.0,0.0,0.0);
  glRasterPos3f(x,y,z);
  glListBase(FontBase);
  glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
#elif _GLX

  // set color and raster position
  glColor3f(0.0, 0.0 , 0.0);
  glRasterPos3f(x,y,z);

  // now we move the raster position based on text justification
  int moveX = 0;
  int moveY = 0;

  int height = fontInfo->ascent;  
  if (horizontalJustify != 'l') {
    char *s;
    int width = 0;
    for (s=text; *s; s++) {
      width += fontInfo->per_char[*s].width;
    }
    if (horizontalJustify == 'r')
      moveX = -width;
    else
      moveX = -width/2;
  }

  if (verticalJustify != 'b') {
    if (verticalJustify == 't')
      moveY = -height - 1;
    else
      moveY = -height/2;
  } else
    moveY = +2;

  glBitmap(0, 0, 0, 0, (GLfloat)moveX, (GLfloat)moveY, NULL);

  // finally draw the text
  glListBase(FontBase);
  glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
#else

#endif
}


int 
OpenGlDevice::saveImage(const char *fileName, int type)
{
  // make the context current
#ifdef _WGL
  return this->saveImageAsBMP(fileName);
#elif _GLX
  return this->saveImageAsPNG(fileName);
#else

#endif

  return 0;
}

int 
OpenGlDevice::saveImageAsBMP(const char *fileName)
{
  // make the context current
#ifdef _WGL
    wglMakeCurrent(theHDC, theHRC);
#elif _GLX
    glXMakeCurrent(theDisplay, theWindow, cx);
#else

#endif

    // create the file name 'bmpFileName$count.BMP'
    char *newFileName = new char[strlen(fileName)+4];
    if (newFileName == 0) {
      opserr << "OpenGlDevice::saveImageAsBMP() failed to open file: " << fileName << endln;
      return -1;
    }	
    strcpy(newFileName, fileName);
    strcat(newFileName,".BMP");

    // open the file
    FILE *fp;
    if ((fp = fopen(newFileName,"wb")) == NULL) {
      opserr << "OpenGLDevice::saveBmpImage() - could not open file named" <<  newFileName << endln;
	return -1;
    }	

#ifdef _WGL	
    int bitsize = (info.bmiHeader.biWidth *
        	   info.bmiHeader.biBitCount + 7) / 8 *
		  abs(info.bmiHeader.biHeight);
    int infosize;
    infosize = sizeof(BITMAPINFOHEADER);
    switch (info.bmiHeader.biCompression) {
	case BI_BITFIELDS :
            infosize += 12; /* Add 3 RGB doubleword masks */
            if (info.bmiHeader.biClrUsed == 0)
	      break;
	case BI_RGB :
            if (info.bmiHeader.biBitCount > 8 &&
        	info.bmiHeader.biClrUsed == 0)
	      break;
	case BI_RLE8 :
	case BI_RLE4 :
            if (info.bmiHeader.biClrUsed == 0)
              infosize += (1 << info.bmiHeader.biBitCount) * 4;
	    else
              infosize += info.bmiHeader.biClrUsed * 4;
	    break;
    }

    int size = sizeof(BITMAPFILEHEADER) + infosize + bitsize;

    // check the bit map header info has been created, if not create one
/*    if (info == 0) {
	if ((info = (BITMAPINFO *)malloc(sizeof(BITMAPINFOHEADER))) < 0) {
	    opserr << "OpenGLDevice::saveBmpImage() - %s\n",
				    "out of memory creating BITMAPINFO struct");
	    count = -1;
	    return -2;
	}
	info->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	info->bmiHeader.biPlanes = 1;
	info->bmiHeader.biBitCount = 24;
	info->bmiHeader.biCompression = BI_RGB;
	info->bmiHeader.biXPelsPerMeter = 2952;
	info->bmiHeader.biYPelsPerMeter = 2952;
	info->bmiHeader.biClrUsed = 0;
	info->bmiHeader.biClrImportant = 0;
    }
    
    // determine the number of bits needed to save the image
    width = viewport[2]*3;
    width = (width+3) & ~3;
    bitsize = width * viewport[3];

    // check the bits pointeris of correct size, if not
    // delete the old and create a new one
    if (bitsize != currentBitSize) {
	if (currentBitSize != 0)
	    free (bits);
	if ((bits = (GLubyte *)calloc(bitsize, 1)) == 0) {
	    opserr << "OpenGLDevice::saveBmpImage() - %s\n",
				    "out of memory creating BITMAPINFO struct");
	    count = -1;
	    free (info);
	    info = 0;
	    return -3;	    
	}
	currentBitSize = bitsize;
    }

    // set the info for the bit map header
    info->bmiHeader.biWidth = viewport[2];
    info->bmiHeader.biHeight = viewport[3];    
    info->bmiHeader.biSizeImage = currentBitSize;        
*/
    // read the pixels from the frame buffer
    glFinish();
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glPixelStorei(GL_PACK_ROW_LENGTH, 0);
    glPixelStorei(GL_PACK_SKIP_ROWS, 0);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 0);

    if (bits == 0) opserr << "BITS ZERO\n";

    glReadPixels(0, 0, info.bmiHeader.biWidth, info.bmiHeader.biHeight,
		GL_BGR_EXT, GL_UNSIGNED_BYTE, bits);
	

    currentBitSize = info.bmiHeader.biWidth * info.bmiHeader.biHeight;
    // create a header for the BMP file
    BITMAPFILEHEADER header;
    header.bfType      = 'MB'; /* Non-portable... sigh */
    header.bfSize      = size;
    header.bfReserved1 = 0;
    header.bfReserved2 = 0;
    header.bfOffBits   = sizeof(BITMAPFILEHEADER) + infosize;

    if (fwrite(&header, 1, sizeof(BITMAPFILEHEADER), fp) < sizeof(BITMAPFILEHEADER))
	{
    // write the header to the file
//    if (fwrite(&header, 1, sizeof(BITMAPFILEHEADER), fp) < sizeof(BITMAPFILEHEADER)) {
	    opserr << "OpenGLDevice::saveBmpImage() - failed to write BITMAPHEADER" << endln;
	    fclose(fp);
	    return -4;
	}
    if (fwrite(&info, 1, infosize, fp) < infosize)
        {
    // write the bit map information to the file
//    if (fwrite(&info, 1, sizeof(BITMAPINFOHEADER), fp) < sizeof(BITMAPINFOHEADER)) {
	    opserr << "OpenGLDevice::saveBmpImage() - failed to write BITMAPINFOHEADER" << endln;
	    fclose(fp);
	    return -5;	    
	}    
    if (fwrite(bits, 1, bitsize, fp) < bitsize)
        {
    // now we write the bits
    //if (fwrite(bits, 1, currentBitSize, fp) < currentBitSize) {
	opserr << "OpenGLDevice::saveBmpImage() - failed to write BITMAPINFOHEADER" << endln;
	fclose(fp);
	return -6;	
    }        
    // if get here we are done .. close file and return ok
#endif
    
    delete [] newFileName;
    fclose(fp);
    return 0;
}



int 
OpenGlDevice::saveImageAsPNG(const char *fileName)
{

  // make the context current
#ifdef _WGL
    wglMakeCurrent(theHDC, theHRC);
#elif _GLX
    glXMakeCurrent(theDisplay, theWindow, cx);
#else

#endif

  //
  // first we read the image from the buffer
  //

  // create some memory to store the image
  char *image;
  image = new char[3*width*height];
  if (image == 0) {
    opserr << "OpenGlDevice::failed to allocate memory for image\n";
    return(-1);
  }
  
  // read into the buffer
  glReadBuffer(GL_BACK_LEFT);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);


  //
  // now we write the file to a png
  //   .. code for this from Greg Roelofs book: PNG: The Definitive Guide published by O'Reilly.
  // 

  // open the file
  FILE *fp;
  char *newFileName = new char (strlen(fileName+4));
  if (newFileName == 0) {
    opserr << "OpenGlDevice::saveImageAsPNG() failed to open file: " << fileName << endln;
    delete [] image;
    return -1;
  }
    
  strcpy(newFileName, fileName);
  strcat(newFileName, ".png");

  if((fp = fopen(newFileName, "wb"))==NULL) {
    opserr << "OpenGlDevice::saveImageAsPNG() failed to open file: " << fileName << endln;
    delete [] image;
    return -1;
  }

#ifdef _PNG  
  // Allocate write & info structures
  png_structp png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if(!png_ptr) {
    delete [] image;
    fclose(fp);
    opserr << "OpenGlDevice::saveImageAsPNG() - out of memery creating write structure\n";
    return -1;
  }
  
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if(!info_ptr) {
    png_destroy_write_struct(&png_ptr,
			     (png_infopp)NULL);
    delete [] image;
    fclose(fp);
    opserr << "OpenGlDevice::saveImageAsPNG() - out of memery creating info structure\n";

    return -1;
  }

  // setjmp() must be called in every function that calls a PNG-writing
  // libpng function, unless an alternate error handler was installed--
  // but compatible error handlers must either use longjmp() themselves
  // (as in this program) or exit immediately, so here we go: */

  if(setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    opserr << "OpenGlDevice::saveImageAsPNG() - setjmp problem\n";
    return 2;
  }

  // make sure outfile is (re)opened in BINARY mode 
  png_init_io(png_ptr, fp);

  // set the compression levels--in general, always want to leave filtering
  // turned on (except for palette images) and allow all of the filters,
  // which is the default; want 32K zlib window, unless entire image buffer
  // is 16K or smaller (unknown here)--also the default; usually want max
  // compression (NOT the default); and remaining compression flags should
  // be left alone

  png_set_compression_level(png_ptr, Z_BEST_COMPRESSION);
  //
  // this is default for no filtering; Z_FILTERED is default otherwise:
  // png_set_compression_strategy(png_ptr, Z_DEFAULT_STRATEGY);
  //  these are all defaults:
  //   png_set_compression_mem_level(png_ptr, 8);
  //   png_set_compression_window_bits(png_ptr, 15);
  //   png_set_compression_method(png_ptr, 8);


  // Set some options: color_type, interlace_type
  int color_type, interlace_type, numChannels;

  //  color_type = PNG_COLOR_TYPE_GRAY;
  //  color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
  color_type = PNG_COLOR_TYPE_RGB; numChannels = 3;
  // color_type = PNG_COLOR_TYPE_RGB_ALPHA;

  interlace_type =  PNG_INTERLACE_NONE;
  // interlace_type = PNG_INTERLACE_ADAM7;
  

  int bit_depth = 8;
  png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, 
	       color_type,
	       interlace_type,
	       PNG_COMPRESSION_TYPE_BASE, 
	       PNG_FILTER_TYPE_BASE);
  
  // Optional gamma chunk is strongly suggested if you have any guess
  // as to the correct gamma of the image. (we don't have a guess)
  //
  // png_set_gAMA(png_ptr, info_ptr, image_gamma);
  
  // write all chunks up to (but not including) first IDAT 
  png_write_info(png_ptr, info_ptr);
  

  // set up the row pointers for the image so we can use png_write_image

  png_bytep* row_pointers = new png_bytep[height];
  if (row_pointers == 0) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    delete [] image;
    fclose(fp);
    opserr << "OpenGlDevice::failed to allocate memory for row pointers\n";
    return(-1);
  }

  unsigned int row_stride = width*numChannels;
  unsigned char *rowptr = (unsigned char*) image;
  for (int row = height-1; row >=0 ; row--) {
    row_pointers[row] = rowptr;
    rowptr += row_stride;
  }

  // now we just write the whole image; libpng takes care of interlacing for us
  png_write_image(png_ptr, row_pointers);
  
  // since that's it, we also close out the end of the PNG file now--if we
  // had any text or time info to write after the IDATs, second argument
  // would be info_ptr, but we optimize slightly by sending NULL pointer: */

  png_write_end(png_ptr, info_ptr);
  
  //
  // clean up after the write
  //    free any memory allocated & close the file
  //
  png_destroy_write_struct(&png_ptr, &info_ptr);
  
  delete [] row_pointers;

#endif

  delete [] image;
  delete [] newFileName;
  fclose(fp);
  
  return 0;
}



