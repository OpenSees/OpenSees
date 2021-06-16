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

// $Revision:
// $Date: 2021/6/16 $
// $Source: C:\Users\Edbert Lumbantobing\Documents\GitHub\OpenSees\SRC\GyroMassMaterial.h $

// Written: Edbert Rainer Lumbantobing (Tokyo Institute of Technology)
// Created: June 2021
// Revision:
//
// Description: This file contains the class declaration for GyroMassMaterial that can be used to model an inerter.
// Inerter is a mechanical element which output force is proportional to the relative acceleration between its terminals (between two nodes). 
// F = b (a1 - a2)
// b = inertance
// a1 = acceleration at node 1
// a2 = acceleration at node 2

// This code has been partly published by Hessabi (2017)
// Reference: Hessabi R. (2017). Application of real-time hybrid simulation method in experimental identification of gyromass dampers (Doctoral thesis). University of Toronto, Toronto.
       

#ifndef GyroMassMaterial_h
#define GyroMassMaterial_h

#include <UniaxialMaterial.h>

class GyroMassMaterial : public UniaxialMaterial
{
  public:
    GyroMassMaterial(int tag, double b);    
    GyroMassMaterial();    

    ~GyroMassMaterial();

    const char* getClassType(void) const { return "GyroMassMaterial"; };

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    
    double getTangent(void);
    double getInitialTangent(void);
    double getDampTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    // function that models the Newmark approximation of acceleration
    int newmark(double gama, double beta, double dt, double trialStrain, double commitStrain, double& trialStrainRate, double commitStrainRate, double& trialAccel, double commitAccel);

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
    double trialStrain;         // displacement at next step
    double trialStress;         // inertial force at next step
    double commitStrain;        // displacement at the current step
    double commitStress;        // inertial force at current step
    double trialStrainRate;     // velocity at next step
    double trialAccel;          // acceleration at next step
    double b;                   // equivalent mass provided for the GMD
    double currentTime;         // time at the next step
    double commitStrainRate;    // velocity at current step
    double commitAccel;         // acceleration at current step
    double commitTime;          // time at the current step
};


#endif



