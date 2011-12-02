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
                                                                        
// $Revision: 1.1 $
// $Date: 2008-02-15 23:40:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/AGL_Device.h,v $
                                                                        
                                                                        
#ifndef AGL_Device_H
#define AGL_Device_H

#include <AGL/agl.h>


class AGL_Device
{
 public:
  AGL_Device();
  virtual ~AGL_Device();  

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

  void initWindow(void); // procedure called on construction of 1st Window

 private:

  
  // save image methods for specific file formats
  int saveImageAsBMP(const char *fileName);  
  int saveImageAsPNG(const char *fileName);  

  WindowRef window;
  GDHandle screen;

  AGLPixelFormat pixelFormat;
  AGLContext context;
  
  //  static GLuint FontBase;
  GLuint FontBase;

  static int numWindows;
  int winOpen;
  int width, height;		// Width and height of our window
  char *windowTitle;
};

#endif



  
