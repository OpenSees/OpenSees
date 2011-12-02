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
// $Date: 2005-06-14 18:49:58 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/MyElasticPP.h,v $
                                                                        
#ifndef MyElasticPP_h
#define MyElasticPP_h

// Written: fmk 
//
// Description: This file contains the class definition for 
// MyElasticPP. MyElasticPP provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) MyElasticPP.h, revA"

#include <UniaxialMaterial.h>

#define MAT_TAG_MyElasticPP 5002

class MyElasticPP : public UniaxialMaterial
{
  public:
    MyElasticPP(int tag, double E, double eyp);    
    MyElasticPP();    

    ~MyElasticPP();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return E;};

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
    double fyp, fyn;	// positive and negative yield stress
    double ezero;	// initial strain
    double E;		// elastic modulus
    double trialStrain;	// trial strain
    double ep;		// plastic strain at last commit

    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
};


#endif



