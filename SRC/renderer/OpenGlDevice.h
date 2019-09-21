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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-04-15 18:01:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGlDevice.h,v $
                                                                        
                                                                        
#ifndef OpenGlDevice_H
#define OpenGlDevice_H

#if _DLL

class OpenGlDevice // dummy
{
public:
	OpenGlDevice() {};
	virtual ~OpenGlDevice() {};

	// Gets the width of the current window
	virtual int GetWidth() { return 0; };

	// Gets the height of the current window
	virtual int GetHeight() { return 0; };

	// Necessary when operating since the drawn
	// image is buffered until this call is made.
	virtual void STARTIMAGE() { return ; };
	virtual void ENDIMAGE() { return ; };

	// Opens a window of the specified width & height.
	virtual void WINOPEN(const char* title, int xLoc, int yLoc, int width, int height) { return ; };;

	// Clears the currently opened window
	virtual void CLEAR() { return ; };;

	virtual void drawText(float x, float y, float z, char* text, int length,
		char horizontalJustify, char verticalJustify) {
		return ;
	};;

	// to save the current image to a file of a specific type
	int saveImage(const char* fileName, int type) { return 0; };;

private:
	void initWindow(void); // procedure called on construction of 1st Window

	// save image methods for specific file formats
	int saveImageAsBMP(const char* fileName) { return 0; };;
	int saveImageAsPNG(const char* fileName) { return 0; };;
};

#else 


#include <Device.h>

#ifdef _GLX

#include <GL/gl.h>
#include <GL/glx.h>

#include <X11/Xlib.h>
#include <X11/X.h>
#include <X11/Xutil.h>
#include <X11/X.h>
#include <X11/Xatom.h>

#define X11_MAX_COLORS 256

#elif _WGL

#include <windows.h>
#include <gl\gl.h>
//#include <gl\glaux.h>

#elif _AGL

class AGL_Device;

#endif

class OpenGlDevice
{
 public:
  OpenGlDevice();
  virtual ~OpenGlDevice();  

  // Gets the width of the current window
  virtual int GetWidth();

  // Gets the height of the current window
  virtual int GetHeight();

  // Necessary when operating since the drawn
  // image is buffered until this call is made.
  virtual void STARTIMAGE();
  virtual void ENDIMAGE();

  // Opens a window of the specified width & height.
  virtual void WINOPEN(const char *title, int xLoc, int yLoc, int width, int height);

  // Clears the currently opened window
  virtual void CLEAR();

  virtual void drawText(float x, float y, float z, char *text, int length, 
			char horizontalJustify, char verticalJustify); 

  // to save the current image to a file of a specific type
  int saveImage(const char *fileName, int type);  

 private:
  void initWindow(void); // procedure called on construction of 1st Window
  
  // save image methods for specific file formats
  int saveImageAsBMP(const char *fileName);  
  int saveImageAsPNG(const char *fileName);  

#ifdef _GLX

  // glx utility toolkit
  Display *theDisplay;  // the display all Window objecs display on
  Window theWindow;
  int theScreen;        // the screen 
  Colormap cmap;        // the colormap all X11 Window objects share   
  GC theGC;
  GLXContext cx;    
  XSizeHints hints; // conatins the infor about where window is and its size
                    //  static unsigned long foreground, background;
  XEvent theEvent;
  XVisualInfo *visual;
  int swap_flag;

  //  static XFontStruct *fontInfo;
  XFontStruct *fontInfo;
  //  static GLuint FontBase;
  GLuint FontBase;

#elif _WGL

  // win32 stuff using wgl toolkit
  HDC         theHDC;        // device context
  HGLRC       theHRC;        // openGL context
  HWND        theWND;        // the window
  BITMAPINFO  info;
  HBITMAP     theBitmap;
  GLint	      viewport[4];
  GLubyte     *bits;
  long	      currentBitSize;
  int xLoc;
  int yLoc;
  //  static GLuint FontBase;
  GLuint FontBase;

#elif _AGL

  AGL_Device *theDevice; 

#else

#endif



  static int numWindows;
  int winOpen;
  int width, height;		// Width and height of our window
  char *windowTitle;
};
#endif
#endif



  
