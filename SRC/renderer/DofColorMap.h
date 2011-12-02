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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/DofColorMap.h,v $
                                                                        
                                                                        
// File: ~/graphics/DofColorMap.h
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for DofColorMap.
// DofColorMap is an abstract base class. An DofColorMap object is used
// to determine the r,g,b values given an input value.
//
// What: "@(#) DofColorMap.h, revA"

#ifndef DofColorMap_h
#define DofColorMap_h

#include <ColorMap.h>

class DofColorMap: public ColorMap
{
  public:
    DofColorMap(int numEqn, int numPartitions);
    float getRed(float value);
    float getGreen(float value);
    float getBlue(float value);
    
  protected:
    
  private:
    int numEqn, numPartitions;
    float nPerP;
    float r[8], g[8], b[8];
};


#endif

