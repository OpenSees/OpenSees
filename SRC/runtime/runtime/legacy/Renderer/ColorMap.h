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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/ColorMap.h,v $
                                                                        
                                                                        
// File: ~/graphics/ColorMap.h
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for ColorMap.
// ColorMap is an abstract base class. An ColorMap object is used
// to determine the r,g,b values given an input value.
//
// What: "@(#) ColorMap.h, revA"

#ifndef ColorMap_h
#define ColorMap_h

class ColorMap
{
  public:
    ColorMap() {};
    virtual ~ColorMap() {};
    virtual float getRed(float value) =0;
    virtual float getGreen(float value) =0;
    virtual float getBlue(float value) =0;
    virtual int   getRGB(float value, float &red, float &green, float &blue) =0;
    virtual int   startImage() =0;
  protected:
    
  private:
};


#endif

