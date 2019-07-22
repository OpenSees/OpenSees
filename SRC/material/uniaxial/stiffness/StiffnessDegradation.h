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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: August 2000
//
// Description: This file contains the interface for StiffnessDegradation,
// which models hysteretic stiffness degradation.

#ifndef StiffnessDegradation_h
#define StiffnessDegradation_h

#include <MaterialState.h>

class UniaxialMaterial;

class StiffnessDegradation : public MaterialState
{
 public:
  StiffnessDegradation(int tag, int classTag);    
  virtual ~StiffnessDegradation();
  
  virtual const char *getMeasure(void) = 0;
  virtual int setTrialMeasure(double measure) = 0;
  virtual double getValue(void) = 0;
  
  virtual void setNegative(bool flag = false) {return;}
  
  virtual int commitState(void) = 0;
  virtual int revertToLastCommit(void) = 0;
  virtual int revertToStart(void) = 0;
  
  virtual StiffnessDegradation *getCopy(void) = 0;
  virtual StiffnessDegradation *getCopy(UniaxialMaterial *u);
  
 protected:
  
 private:
};

extern bool OPS_addStiffnessDegradation(StiffnessDegradation *newComponent);
extern StiffnessDegradation *OPS_getStiffnessDegradation(int tag);
extern void OPS_clearAllStiffnessDegradation(void);

#endif
