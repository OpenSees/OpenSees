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
  

                                                                        
// $Revision: 1.1 $
// $Date: 2004-09-01 03:54:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/ParkAng.h,v $
                                                                        
#ifndef ParkAng_h
#define ParkAng_h         
                                                               
// Written: AA,GGD
// Created: 10/02
// Revision: AA
//
// Description: This file contains the class definition for 
// Cimbined damage model. It is a subclass od DamageModel
//

#include <ErrorHandler.h>
#include <DamageModel.h>

class DamageResponse;


class ParkAng : public DamageModel
{
  public:
    ParkAng(int tag, double deltaU , double beta , double sigmaY );
	ParkAng();  
    ~ParkAng();

    int setTrial (double scalar, double scalarRate = 0.0 );
	int setTrial (Vector trialVector );
	int setTrial (); 

	int setInputResponse ( Element *elem , const char **argv , int argc, int ndof );

    double getDamage(void);
	double getPosDamage (void);
	double getNegDamage (void);
    
    int commitState(void);
    int revertToLastCommit (void);
    int revertToStart (void);

    DamageModel *getCopy (void);

    int setVariable(const char *argv);
    int getVariable(int variableID, double &info);

    int setParameter(char **argv, int argc, Information &eleInformation);
    int updateParameter(int responseID, Information &eleInformation);	

    Response *setResponse(char **argv, int argc, Information &info);
    int getResponse(int responseID, Information &info);

	int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag =0);    

    // method for this damage model to update itself according to its new parameters
    void update(void) {return;}

  protected:
    
  private:
	
	double DeltaU , Beta , SigmaY;

	double TrialInfo[6];
	double CommitInfo[6];
	double LastCommitInfo[6];

};


#endif
