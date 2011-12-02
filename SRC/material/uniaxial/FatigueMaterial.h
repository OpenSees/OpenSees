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
// $Date: 2003-08-14 20:23:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FatigueMaterial.h,v $
                                                      
// Written: Patxi
// Created: Aug 2003
//
// Description: This file contains the class definition for 
// FatigueMaterial.  FatigueMaterial wraps a UniaxialMaterial
// and imposes fatigue limits.

#ifndef FatigueMaterial_h
#define FatigueMaterial_h

#include <UniaxialMaterial.h>

class FatigueMaterial : public UniaxialMaterial
{
  public:
    FatigueMaterial(int tag, UniaxialMaterial &material, 
		    double Dmax = 1.0,
		    double Nf  = 0.5E6,
		    double E0  = 14.5/29000.0,
		    double FE  = -1.0/2.31);
    FatigueMaterial();
    ~FatigueMaterial();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStrainRate(void);
    double getStress(void);
    double getTangent(void);
    double getDampTangent(void);
    double getInitialTangent(void) {return theMaterial->getInitialTangent();}

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
    UniaxialMaterial *theMaterial;

    double D;   //; % Damage index
    double X;   //; % Range in consideration
    double Y;   //; % Previous Adjacent Range
    double A;   //; % Strain at first  cycle peak/valley
    double B;   //; % Strain at second cycle peak/valley
    int R1F;    // % Flag for first  cycle count
    int R2F;    // % Flag for second cycle count
    double CS;  // % Current Slope
    double PS;  // % Previous slope
    double EP;  // % Previous Strain
    int FF;     // % Failure Flag
    int SF;     // % Start Flag - for initializing the very first strain
    int PF;     // % Peak Flag --> Did we reach a peak/valley at current strain?
    
    double Dmax; // = 1.0; % Maximum Damage index, normally 1.0
    double Nf;   //  = 0.5E6; % Number of cycles to failure at strain E0
    double E0;   //  = 14.5/29000;
    double FE;   //  = -1/2.31; % Log log slope of fracture curve
    double b;    //   =  log10(E0) - log10(Nf)*FE; %Theoretical Intercept
    
    bool Tfailed;
    bool Cfailed;
};


#endif

