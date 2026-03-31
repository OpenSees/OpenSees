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
                                                                        
// $Date: 2018/11/25 
                                                   
#ifndef Ratchet_h
#define Ratchet_h

// Written: Yi Xiao from Tongji University/University of Washington 
// Modified from GNGMaterial by Jook
//
// Description: This file contains the class definition for 
// Ratchet.
//
// What: "@(#) Ratchet.h, revA"

#include <UniaxialMaterial.h>

class Ratchet : public UniaxialMaterial
{
  public:
    Ratchet(int tag, double E, double fTravel, double fTravelInitial, int RatType);    
    Ratchet();    

    ~Ratchet();

    const char *getClassType(void) const {return "Ratchet";}

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) ;

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
    double commitStrain; //
    double commitStress;
    double commitTangent;
    double commitEngageStrain;
    double trialStrain;	// 
    double trialStress; 
    double trialTangent;
    
    double E;            //d
    double freeTravel;    //d
    double fTravelInitial;    //d
    double engageStrain; //
    double currentStrain;//
   
	  int RatType;
    int nratchet; //ratchet count
	  int commitNratchet;
   
};


#endif



