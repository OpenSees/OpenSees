
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
  

                                                                        
// $Revision: 1.3 $
// $Date: 2008-04-14 22:38:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/Kratzig.h,v $
                                                                        
#ifndef Kratzig_h
#define Kratzig_h         
                                                               
// Written: AA,GGD
// Created: 08/02
// Revision: AA
//
// Description: This file contains the class definition for 
// Damage. Damage is a base class and 
// thus no objects of it's type can be instantiated. It has pure virtual 
// functions which must be implemented in it's derived classes. 


#include <ErrorHandler.h>
#include <DamageModel.h>

class DamageResponse;


class Kratzig : public DamageModel
{
  public:
    Kratzig(int tag, double ultimatePosVal ,  double ultimateNegVal);

	Kratzig();  
    
	~Kratzig();

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
	
	// Model parameters
	double UltimatePosValue , UltimateNegValue;

	double TrialInfo[10];
	double CommitInfo[10];
	double LastCommitInfo[10];
	
};


#endif
