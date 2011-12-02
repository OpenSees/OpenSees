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
// $Date: 2000-09-15 08:23:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/PlainMap.cpp,v $
                                                                        
                                                                        
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

#include <PlainMap.h>
#include <iostream.h>
#include <stdlib.h>


PlainMap::PlainMap()
{
}


float 
PlainMap::getRed(float value){
    if (value > 0) 
	return 1.0;
    else
	return 0.0;
}
      

float 
PlainMap::getGreen(float value) {
    return 0.0;
}

float 
PlainMap::getBlue(float value) {
    if (value < 0) 
	return 1.0;
    else
	return 0.0;
}
    

