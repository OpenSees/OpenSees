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
// $Source: /usr/local/cvs/OpenSees/SRC/damage/HystereticEnergy.h,v $
                                                                        
#ifndef HystereticEnergy_h
#define HystereticEnergy_h         
                                                               
// Written: Arash Altoontash, Gregory Deierlein
// Created: 06/03
// Revision: AA
//
// Description: This file contains the class definition for 
// Cimbined damage model. It is a subclass of DamageModel
//

#include <ErrorHandler.h>
#include <DamageModel.h>

class DamageResponse;


class HystereticEnergy : public DamageModel
{
  public:
    HystereticEnergy(int tag, double Etot , double Cpow);
    HystereticEnergy();  
    ~HystereticEnergy();
    
    int setTrial (const Vector &trialVector );
    int setTrial (); 
    
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
    
    double Etotal,Cpower;
    
    // Trial step
    double TrialInfo[8];
    double CommitInfo[8];
    double LastCommitInfo[8];
};


#endif
