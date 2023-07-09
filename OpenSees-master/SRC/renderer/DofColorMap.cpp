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
// $Date: 2003-02-14 23:01:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/DofColorMap.cpp,v $
                                                                        
                                                                        
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

#include <DofColorMap.h>
#include <OPS_Globals.h>
#include <stdlib.h>


DofColorMap::DofColorMap(int _numEqn, int _numP)
  :numEqn(_numEqn),numPartitions(_numP)
{
  nPerP=numEqn/numPartitions;

  r[0] = 0; g[0] = 1; b[0] = 1;
  r[1] = 1; g[1] = 1; b[1] = 0;
  r[2] = 1; g[2] = 0; b[2] = 1;
  r[3] = 0; g[3] = 1; b[3] = 0;
  r[4] = 0; g[4] = 0; b[4] = 1;
  r[5] = 0; g[5] = 0; b[5] = 0;
  r[6] = 0.5; g[6] = 0.5; b[6] = 0.5;
  r[7] = 1; g[7] = 0; b[7] = 0;

  if (_numP > 8) {
    opserr << "DofColorMap- can't handle numP > 8\n";
    exit(-1);
  }

}


float 
DofColorMap::getRed(float value){
  int j = value/nPerP;
  return r[j];

}
      

float 
DofColorMap::getGreen(float value) {
  int j = value/nPerP;
  return g[j];
}

float 
DofColorMap::getBlue(float value) {
  int j = value/nPerP;
  return b[j];
}
    

