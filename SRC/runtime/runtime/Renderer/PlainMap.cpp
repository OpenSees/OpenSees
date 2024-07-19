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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/PlainMap.cpp,v $
                                                                        
// Written: fmk 
// Created: 10/98
//
// Description: This file contains the class definition for PlainMap.
// PlainMap is an abstract base class. An PlainMap object is used
// to determine the r,g,b values given an input value.
//
// What: "@(#) PlainMap.h, revA"

#include "PlainMap.h"
#include <OPS_Stream.h>
#include <stdlib.h>
#include <math.h>

static float jet[64*3] = {
  0,         0,    0.5625,
  0,         0,    0.6250,
  0,         0,    0.6875,
  0,         0,    0.7500,
  0,         0,    0.8125,
  0,         0,    0.8750,
  0,         0,    0.9375,
  0,         0,    1.0000,
  0,    0.0625,    1.0000,
  0,    0.1250,    1.0000,
  0,    0.1875,    1.0000,
  0,    0.2500,    1.0000,
  0,    0.3125,    1.0000,
  0,    0.3750,    1.0000,
  0,    0.4375,    1.0000,
  0,    0.5000,    1.0000,
  0,    0.5625,    1.0000,
  0,    0.6250,    1.0000,
  0,    0.6875,    1.0000,
  0,    0.7500,    1.0000,
  0,    0.8125,    1.0000,
  0,    0.8750,    1.0000,
  0,    0.9375,    1.0000,
  0,    1.0000,    1.0000,
  0.0625,    1.0000,    1.0000,
  0.1250,    1.0000,    0.9375,
  0.1875,    1.0000,    0.8750,
  0.2500,    1.0000,    0.8125,
  0.3125,    1.0000,    0.7500,
  0.3750,    1.0000,    0.6875,
  0.4375,    1.0000,    0.6250,
  0.5000,    1.0000,    0.5625,
  0.5625,    1.0000,    0.5000,
  0.6250,    1.0000,    0.4375,
  0.6875,    1.0000,    0.3750,
  0.7500,    1.0000,    0.3125,
  0.8125,    1.0000,    0.2500,
  0.8750,    1.0000,    0.1875,
  0.9375,    1.0000,    0.1250,
  1.0000,    1.0000,    0.0625,
  1.0000,    1.0000,         0,
  1.0000,    0.9375,         0,
  1.0000,    0.8750,         0,
  1.0000,    0.8125,         0,
  1.0000,    0.7500,         0,
  1.0000,    0.6875,         0,
  1.0000,    0.6250,         0,
  1.0000,    0.5625,         0,
  1.0000,    0.5000,         0,
  1.0000,    0.4375,         0,
  1.0000,    0.3750,         0,
  1.0000,    0.3125,         0,
  1.0000,    0.2500,         0,
  1.0000,    0.1875,         0,
  1.0000,    0.1250,         0,
  1.0000,    0.0625,         0,
  1.0000,         0,         0,
  0.9375,         0,         0,
  0.8750,         0,         0,
  0.8125,         0,         0,
  0.7500,         0,         0,
  0.6875,         0,         0,
  0.6250,         0,         0,
  0.5625,         0,         0};



PlainMap::PlainMap()
  :max(0.0), min(0.0), maxLast(0.0), minLast(0.0)
{
  data = jet;
  sizeData = 64;
}


float 
PlainMap::getRed(float value){
  if (value > max)
    max = value;
  else if (value < min)
    min = value;

  if (maxLast == minLast) {
    int index = sizeData/2;
    return data[index*3-3];

  } else if (value > maxLast) 
    return data[sizeData*3-3];
  else if (value < minLast)
    return data[0];
  else {
    int index = (floor)((value-minLast)*sizeData/((maxLast-minLast)));
    return data[index*3-3];
  }
    
  return 0.0;
}
      

float 
PlainMap::getGreen(float value) {
  if (value > max)
    max = value;
  else if (value < min)
    min = value;

  if (maxLast == minLast) {
    int index = sizeData/2;
    return data[index*3-2];

  } else if (value > maxLast) 
    return data[sizeData*3-2];
  else if (value < minLast)
    return data[1];
  else {
    int index = (floor)((value-minLast)*sizeData/((maxLast-minLast)));
    return data[index*3-2];
  }
    
  return 0.0;
}

float 
PlainMap::getBlue(float value) {
  if (value > max)
    max = value;
  else if (value < min)
    min = value;

  if (maxLast == minLast) {
    int index = sizeData/2;
    return data[index*3-1];

  } else if (value > maxLast) 
    return data[sizeData*3-1];
  else if (value < minLast)
    return data[2];
  else {
    int index = (floor)((value-minLast)*sizeData/((maxLast-minLast)));
    return data[index*3-1];
  }
    
  return 0.0;
}

int
PlainMap::getRGB(float value, float &red, float &blue, float &green) {
  red = this->getRed(value);
  green = this->getGreen(value);
  blue = this->getBlue(value);
  return 0;
}
    


int
PlainMap::startImage()
{
  maxLast = max;
  minLast = min;
  max = 0;
  min = 0;

  if (maxLast > -minLast)
    maxLast = -minLast;
  if (minLast < -maxLast)
    minLast = -maxLast;
  return 0;
}
