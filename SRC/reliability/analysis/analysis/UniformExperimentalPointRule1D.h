/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//
 
#ifndef UniformExperimentalPointRule1D_H
#define UniformExperimentalPointRule1D_H

#include "ExperimentalPointRule1D.h"

class UniformExperimentalPointRule1D : public ExperimentalPointRule1D  
{
public:
  int getPointClosestToOrigin();
  UniformExperimentalPointRule1D();
  UniformExperimentalPointRule1D(ExperimentalPointRule1D * pExperimentalPointRule1D);
  virtual ~UniformExperimentalPointRule1D();
  
  Vector * getPointCoordinates();
  double getPointCoordinate(int i);
  char * getType();
  
 private:
  
  char type[40];
  Vector * tmp;
  
};

#endif 
