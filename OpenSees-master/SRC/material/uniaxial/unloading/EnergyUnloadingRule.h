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
// Description: This file contains the implementation of 
// EnergyStrengthDegradation, which models hysteretic strength
// degradation with the deterioration parameter developed by
// Rahnama and Krawinkler (1993).

#ifndef EnergyUnloadingRule_h
#define EnergyUnloadingRule_h

#include <UnloadingRule.h>

class EnergyUnloadingRule : public UnloadingRule
{
 public:
  EnergyUnloadingRule(int tag, double Et, double c);
  EnergyUnloadingRule();
  virtual ~EnergyUnloadingRule();
  
  const char *getMeasure(void);
  int setTrialMeasure(double measure);
  double getValue(void);
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  UnloadingRule *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag = 0);
  
 protected:
  
 private:
  double Et;
  double c;
  
  double energyExcursion;
  
  double TenergySum;
  double CenergySum;
  
  double Tfactor;
  double Cfactor;
};

#endif
