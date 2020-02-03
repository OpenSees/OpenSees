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
// Created: July 2000
//
// Description: This file contains the implementation of 
// ACIStrengthDegradation, which models hysteretic shear strength
// degradation as a function of beam section curvature ductility.

#ifndef DuctitlityDegradation_h
#define DuctitlityDegradation_h

#include <StrengthDegradation.h>

class ACIStrengthDegradation : public StrengthDegradation
{
 public:
  ACIStrengthDegradation(int tag, double ky, double d1, double V2, double d2);
  ACIStrengthDegradation();
  virtual ~ACIStrengthDegradation();
  
  const char *getMeasure(void);
  int setTrialMeasure(double measure);
  double getValue(void);
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  StrengthDegradation *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag = 0);
  
 protected:
  
 private:
  double V2;
  double d1, d2;
  double oneOverKy;
  double Tductility;
  double Cductility;
};

#endif
