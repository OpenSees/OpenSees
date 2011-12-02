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
// $Date: 2000-09-15 08:23:25 $
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
    
  protected:
    
  private:
};


#endif

