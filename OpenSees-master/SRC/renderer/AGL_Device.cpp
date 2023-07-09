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
                                                                        
// $Revision: 1.2 $
// $Date: 2010-02-04 01:21:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/AGL_Device.cpp,v $
                                                                        
#include <AGL_Device.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int AGL_Device::numWindows(0);
// GLuint AGL_Device::FontBase(0);

  //  AGL_WINDOW,

static GLint ATTRIBUTES[] =
{
  AGL_RGBA,
  AGL_BUFFER_SIZE, 24,
  AGL_DEPTH_SIZE, 24,
  AGL_DOUBLEBUFFER,
  AGL_NONE
};


AGL_Device::AGL_Device()
  :FontBase(0), winOpen(1), width(0), height(0), windowTitle(0)
{
  context = NULL;
  
  // if this is the first window to open call init windows
  if (numWindows == 0)
    this->initWindow();
  
  numWindows++;
}

AGL_Device::~AGL_Device()
{
  numWindows--;
  
  aglSetCurrentContext(NULL);
  aglSetDrawable(context, NULL);
  aglDestroyContext(context);
  //  ReleaseWindow(window);
}

void
AGL_Device::WINOPEN(const char *_title, int _xLoc, int _yLoc, int _width, int _height)
{
  if (windowTitle != 0)
    delete [] windowTitle;
  
  windowTitle = new char[strlen(_title)+1];
  strcpy(windowTitle, _title);
  
  width = _width;
  height = _height;
  
  //Rect WINDOW_BOUNDS = {_xLoc, _yLoc, _width,_height};
  Rect WINDOW_BOUNDS;
  SetRect(&WINDOW_BOUNDS, _xLoc, _yLoc, _xLoc+_width,_yLoc+_height);
  
  OSStatus carbonError =
    CreateNewWindow(kDocumentWindowClass,
		    kWindowStandardDocumentAttributes |
		    kWindowStandardHandlerAttribute,
		    &WINDOW_BOUNDS,
		    &window);
  if (carbonError != noErr || window == NULL) {
    fprintf(stderr, "Can't create window (%u)\n", (unsigned)carbonError);
  }
  
  screen = GetGWorldDevice(GetWindowPort(window));
  if (screen == 0) {
    fprintf(stderr, "Can't get GDevice for window\n");
    //    ReleaseWindow(window);
  }
  
  //    pixelFormat = aglChoosePixelFormat(&screen, 1, ATTRIBUTES);
  pixelFormat = aglChoosePixelFormat(NULL, 1, ATTRIBUTES);
  if (pixelFormat == NULL) {
    fprintf(stderr, "Can't choose pixel format\n");
    //    ReleaseWindow(window);
  }
  
  context = aglCreateContext(pixelFormat, NULL);
  if (context == NULL) {
    fprintf(stderr, "Can't create context\n");
    aglDestroyPixelFormat(pixelFormat);
  }
  
  // MAC OSX 10.5
  //    aglSetWindowRef(context, window);
  
  // aglDestroyPixelFormat(pixelFormat);
  
  if (!aglSetDrawable(context, GetWindowPort(window))) {
      fprintf(stderr, "Can't set context's drawable\n");
      aglDestroyContext(context);
      //      ReleaseWindow(window);
    }
  
  if (!aglSetCurrentContext(context)) {
      fprintf(stderr, "Can't make context current\n");
      aglDestroyContext(context);
      //      ReleaseWindow(window);
    }
  
    //    aglSetWindowRef(context, window);
    ShowWindow(window);

    aglSetCurrentContext(NULL);
}


void
AGL_Device::CLEAR()
{
  if (!aglSetCurrentContext(context)) {
    fprintf(stderr, "Can't make context current\n");
    aglDestroyContext(context);
    //    ReleaseWindow(window);
  }
}

void
AGL_Device::STARTIMAGE()
{
  /*    if (!aglSetDrawable(context, GetWindowPort(window)))
    {
        fprintf(stderr, "Can't set context's drawable\n");
        aglDestroyContext(context);
        ReleaseWindow(window);

    }
  */

    aglSetCurrentContext(context);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void
AGL_Device::ENDIMAGE()
{
  aglSwapBuffers(context);
  aglSetCurrentContext(NULL);
}


int
AGL_Device::GetWidth()
{
  return width;
}

int
AGL_Device::GetHeight()
{
  return height;
}
  


void
AGL_Device::initWindow(void) {

}    
	
    
void
AGL_Device::drawText(float x, float y, float z, char *text, int length,
		       char horizontalJustify, char verticalJustify)
{

}


int 
AGL_Device::saveImage(const char *fileName, int type)
{
  // make the context current
  return 0;
}

int 
AGL_Device::saveImageAsBMP(const char *fileName)
{
  return -1;
}



int 
AGL_Device::saveImageAsPNG(const char *fileName)
{
  return -1;;
}



