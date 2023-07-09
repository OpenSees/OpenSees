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

// $Revision: 1.2 $
// $Date: 2008-04-14 21:27:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SelfCenteringMaterial.h,v $

#ifndef SelfCenteringMaterial_h
#define SelfCenteringMaterial_h

// Written: JAE
// Created: Oct 2007
//
// Description: This file contains the class definition for 
// Self-Centering Material.  SelfCenteringMaterial provides 
// the abstraction for a one-dimensional rate-independent 
// flag-shaped hysteresis with the possibility of permanent plastic
// deformation at large strains.

#include <UniaxialMaterial.h>
//#include <Matrix.h>

class SelfCenteringMaterial : public UniaxialMaterial
{
  public:
    SelfCenteringMaterial(int tag, double k1, double k2,
		      double ActF, double beta, double SlipDef, 
			  double BearDef, double rBear);
    SelfCenteringMaterial();
    ~SelfCenteringMaterial();

    const char *getClassType(void) const {return "SelfCenteringMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return k1;};

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
    // Material parameters (from input)
    double k1;		// Precompression Stiffness
    double k2;		// Prestress Stiffness
    double ActF;	// Activation Stress/Force
    double beta;	// Flag-Shape Parameter
    double rBear;	// Bearing Stiffness Factor
    double SlipDef;	// Slip Strain/Deformation
    double BearDef;	// Bearing Strain/Deformation
    
    // Extra Calculated Material Parameters
    double SlipF;	// External Fuse Slip Stress/Force
    double ActDef;	// Actvation Strain/Deformation
    double BearF;	// Bearing Stress/Force
    
    double diffStrain;		// Difference of strain from last step
    double noSlipStrain;	// Tstrain minus the current Slip Strain
    
    // Committed history variables
    double CactivStrainPos;	// Committed activation strain (Pos Quad)
    double CactivStrainNeg;	// Committed activation strain (Neg Quad)
    double CslipStrain;		// Committed slip strain
    double CupperStrainPos;	// Committed upper activation strain (Pos Quad)
    double ClowerStrainPos;	// Committed lower activation strain (Pos Quad)
    double CupperStressPos;	// Committed upper activation stress (Pos Quad)
    double ClowerStressPos;	// Committed lower activation stress (Pos Quad)
    double CupperStrainNeg;	// Committed upper activation strain (Neg Quad)
    double ClowerStrainNeg;	// Committed lower activation strain (Neg Quad)
    double CupperStressNeg;	// Committed upper activation stress (Neg Quad)
    double ClowerStressNeg;	// Committed lower activation stress (Neg Quad)
    
    // Trial history variables
    double TactivStrainPos;	// Trial activation strain (Pos Quad)
    double TactivStrainNeg;	// Trial activation strain (Neg Quad)
    double TslipStrain;		// Trial slip strain
    double TupperStrainPos;	// Trial upper activation strain (Pos Quad)
    double TlowerStrainPos;	// Trial lower activation strain (Pos Quad)
    double TupperStressPos;	// Trial upper activation stress (Pos Quad)
    double TlowerStressPos;	// Trial lower activation stress (Pos Quad)
    double TupperStrainNeg;	// Trial upper activation strain (Neg Quad)
    double TlowerStrainNeg;	// Trial lower activation strain (Neg Quad)
    double TupperStressNeg;	// Trial upper activation stress (Neg Quad)
    double TlowerStressNeg;	// Trial lower activation stress (Neg Quad)
    
    // Trial state variables
    double Tstrain;		// Trial strain
    double Tstress;		// Trial stress
    double Ttangent;	// Trial tangent
    
    // Committed State Variables
    double Cstrain;		// Committed Strain
    double Cstress;		// Committed Strain
    double Ctangent;		// Committed Strain
};


#endif

