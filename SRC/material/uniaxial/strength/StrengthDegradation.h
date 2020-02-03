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
// Description: This file contains the interface for StrengthDegradation,
// which models hysteretic strength degradation.

#ifndef StrengthDegradation_h
#define StrengthDegradation_h

#include <MaterialState.h>

class UniaxialMaterial;

class StrengthDegradation : public MaterialState
{
 public:
  StrengthDegradation(int tag, int classTag);    
  virtual ~StrengthDegradation();
  
  virtual const char *getMeasure(void) = 0;
  virtual int setTrialMeasure(double measure) = 0;
  virtual double getValue(void) = 0;
  
  virtual void setNegative(bool flag = false) {return;}
  
  virtual int commitState(void) = 0;
  virtual int revertToLastCommit(void) = 0;
  virtual int revertToStart(void) = 0;
  
  virtual StrengthDegradation *getCopy(void) = 0;
  virtual StrengthDegradation *getCopy(UniaxialMaterial *u);
  
 protected:
  
 private:

};

extern bool OPS_addStrengthDegradation(StrengthDegradation *newComponent);
extern StrengthDegradation *OPS_getStrengthDegradation(int tag);
extern void OPS_clearAllStrengthDegradation(void);

#endif
