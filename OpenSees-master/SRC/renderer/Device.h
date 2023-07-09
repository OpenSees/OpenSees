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
// $Date: 2003-02-26 18:56:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Device.h,v $
                                                                        
                                                                        
#ifndef Device_H
#define Device_H

class Device
{
 public:
  Device();
  virtual ~Device();  

  // Specify a 2D point to the hardware  
  virtual void V2F(float x, float y) =0;

  // Specify a color to the hardware
  virtual void C3F(float r, float g, float b) =0;

  // Gets the width of the current window
  virtual int GetWidth() =0;

  // Gets the height of the current window
  virtual int GetHeight() =0;

  // Call when about to begin/end drawing a polygon. All V2F calls
  // from then on until ENDPOLYGON will be interpreted as 
  // vertices of the polygon
  virtual void BGNPOLYGON() =0;
  virtual void ENDPOLYGON() =0;

  // draw text
  virtual void drawText(float x, float y, char *text, int length) =0;  

  // Same as BGNPOLYGON but for wireframe polygons
  virtual void BGNCLOSEDLINE() =0;
  virtual void ENDCLOSEDLINE() =0;

  // Call when about to begin drawing a set of points. All V2F
  // calls from then on until ENDPOINT will be interpreted as
  // points.
  virtual void BGNPOINT() =0;
  virtual void ENDPOINT() =0;

  // Necessary when operating in XWINDOWS mode since the drawn
  // image is buffered until this call is made.
  virtual void STARTIMAGE() =0;
  virtual void ENDIMAGE() =0;

  // to open a window
  virtual void WINOPEN(const char *title, int xLoc, int yLoc, int width, int height) =0;

  // Clears the currently opened window
  virtual void CLEAR() =0;


 private:

};

#endif



  
