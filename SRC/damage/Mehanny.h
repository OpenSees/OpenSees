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
                                                                        
#ifndef Mehanny_h
#define Mehanny_h         

// $Revision: 1.3 $
// $Date: 2008-04-14 22:38:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/Mehanny.h,v $


// Written: Arash Altoontash, Gregory Deierlein
// Created: 10/02
// Revision: AA
//
// Description: This file contains the class definition for 
// Cimbined damage model. It is a subclass od DamageModel
//

#include <ErrorHandler.h>
#include <DamageModel.h>

class DamageResponse;


class Mehanny : public DamageModel
{
public:
	
  Mehanny(int tag, double alpha, double beta, double gamma ,
	  double ultimatePosValue ,  double ultimateNegValue, double abstol, double reltol, double posmodifier, double negmodifier);
  Mehanny();  
  ~Mehanny();
  
  int setTrial (const Vector &trialVector );
  int setTrial () { return -1; } 
  
  int setInputResponse ( Element *elem , const char **argv , int argc, int ndof );
  
  double getDamage(void);
  double getPosDamage (void);
  double getNegDamage (void);
  
  int commitState(void);
  int revertToLastCommit (void);
  int revertToStart (void);
  
  DamageModel *getCopy (void);
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &info);
  int getResponse(int responseID, Information &info);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);    
  
  // method for this damage model to update itself according to its new parameters
  void update(void) {return;}
  
 protected:
  
 private:
  int processData (double PDefo);
  // Model parameters
  
  double Alpha , Beta , Gamma , UltimatePosValue , UltimateNegValue;
  double PosModifier, NegModifier;
  double AbsTol, RelTol;
  
  // Trial step
  double TrialPlasticDefo;
  double TrialDefoIncr;
  double TrialTempPDefo;
  double TrialPosCycle;
  double TrialNegCycle;	
  double TrialSumPosFHC;
  double TrialPosPHC;
  double TrialSumNegFHC;
  double TrialNegPHC;
  double TrialDamage;
  
  // Committed state
  double CommPlasticDefo;
  double CommDefoIncr;
  double CommTempPDefo;
  double CommPosCycle;
  double CommNegCycle;	
  double CommSumPosFHC;
  double CommPosPHC;
  double CommSumNegFHC;
  double CommNegPHC;
  double CommDamage;
  
  // Last commit
  double LCommPlasticDefo;
  double LCommDefoIncr;
  double LCommTempPDefo;
  double LCommPosCycle;
  double LCommNegCycle;	
  double LCommSumPosFHC;
  double LCommPosPHC;
  double LCommSumNegFHC;
  double LCommNegPHC;
  double LCommDamage;
  
  //  FILE *OutputFile;		// For debugging
};


#endif

