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
// $Date: 2001-07-26 00:56:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/PlainMap.h,v $
                                                                        
                                                                        
// File: ~/graphics/PlainMap.h
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for PlainMap.
// PlainMap is an abstract base class. An PlainMap object is used
// to determine the r,g,b values given an input value.
//
// What: "@(#) PlainMap.h, revA"

#ifndef PlainMap_h
#define PlainMap_h

#include <ColorMap.h>

class PlainMap: public ColorMap
{
  public:
    PlainMap();
    float getRed(float value);
    float getGreen(float value);
    float getBlue(float value);
    int   getRGB(float value, float &red, float &green, float &blue);
    int   startImage();

  protected:
    
  private:
    float max, min;
    float maxLast, minLast;
    float *data;
    int sizeData;
};


#endif

