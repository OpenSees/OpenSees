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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HardeningMaterial.h,v $
                                                                        
                                                                        
#ifndef HardeningMaterial_h
#define HardeningMaterial_h

// File: ~/material/HardeningMaterial.h
//
// Written: MHS
// Created: May 2000
// Revision: A
//
// Description: This file contains the class definition for 
// HardeningMaterial.  HardeningMaterial provides the abstraction
// for a one-dimensional rate-independent plasticity model
// with combined isotropic and kinematic hardening.
//
//
// What: "@(#) HardeningMaterial.h, revA"

#include <UniaxialMaterial.h>

class HardeningMaterial : public UniaxialMaterial
{
  public:
    HardeningMaterial(int tag, double E, double sigmaY, double K, double H);
	HardeningMaterial();
    ~HardeningMaterial();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getSecant (void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(ostream &s, int flag =0);
    
  protected:
    
  private:
    // Material parameters
    double E;	// Elastic modulus
    double sigmaY;	// Yield stress
    double K;	// Isotropic hardening parameter
    double H;	// Kinematic hardening parameter
	
    // Trial state variables
    double Tstrain;	// Trial strain
    double Tstress;	// Trial stress
    double Ttangent;	// Trial tangent
    double TplasticStrain;	// Trial plastic strain
    double TbackStress;	// Trial back stress
    double Thardening;	// Trial internal hardening variable
	
    // Committed history variables
    double Cstrain;	// Committed strain
    double CplasticStrain;	// Committed plastic strain
    double CbackStress;	// Committed back stress;
    double Chardening;	// Committed internal hardening variable
};


#endif

