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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Scan.h,v $
                                                                        
                                                                        
#ifndef Scan_H
#define Scan_H

#include "db.H"
class Device;
class Projection;

// This module transforms a list of polygons in device coordinates
// into pixels that represent these polygons.

#define EVEN_ODD 1
#define WINDING  2

class ScanLineConverter
{
 public:
  ScanLineConverter();
  ~ScanLineConverter();  
  
  int update(void);
  void scanLine(FACE &FACE); 
  void scanPolygon(FACE &FACE);   
  void setFillMode(int mode);
  void setScanMode(int mode);  

  void setDevice(Device &theDevice);
  void setProjection(Projection &theProjection);  
  
 private:
  int width;
  int length;
  int fillMode;      // wire == 1, otherwise fill 
  int scanMode;      // WINDING or EVEN_ODD
  float *zBuffer;
  int sizeBuffer;
  
  void DrawWire(FACE &face);	     // Draws a single face in wireframe 
  void DrawFill(FACE &face);         // Draws a face using fill
      
  Device *theDevice;
  Projection *theProjection;  
};

#endif


