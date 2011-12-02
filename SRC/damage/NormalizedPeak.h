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
                                                                        
#ifndef NormalizedPeak_h
#define NormalizedPeak_h         

// $Revision: 1.2 $
// $Date: 2008-04-14 22:38:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/NormalizedPeak.h,v $

// Written: Arash Altoontash, Gregory Deierlein
// Created: 08/02
// Revision: AA
//
// Description: This file contains the class definition for NormalizedPeak
// Damage model. This is a maximum based, non hysteretic damage model. 
// It gets the maximum positive and maximum negative values as initial
// parameters and calculates the damage index based on the maximum and minimum
// values occured.


#include <ErrorHandler.h>
#include <DamageModel.h>

class DamageResponse;
class Element;


class NormalizedPeak : public DamageModel
{
  public:
  NormalizedPeak(int tag, double maxVal, double minVal , const char *argv);
  NormalizedPeak();  
  ~NormalizedPeak();
  
  int setTrial (const Vector &trialVector);
  int setTrial () { return -1; }
  
  double getDamage(void);
  double getPosDamage (void);
  double getNegDamage (void);
  
  int commitState(void);
  int revertToLastCommit (void);
  int revertToStart (void);
  
  DamageModel *getCopy (void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  
 private:
  char damagename[80];
  DamageType damagetype;
  
  // Model parameters
  double MaxValue, MinValue;
  
  // Trial step
  double TrialScalar;
  double TrialDmg;
  Vector TrialVector;
  
  // Commited state
  double CommitScalar;
  double CommitDmg;
  Vector CommitVector;
  
  // Last commit
  double LCommitScalar;
  double LCommitDmg;
  Vector LCommitVector;
  
};

#endif
