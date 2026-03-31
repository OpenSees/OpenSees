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

// $Revision: 1.0 $
// $Date: 2020-09-01 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASD_SMA_3K.cpp,v $

// Written: Luca Aceto
// Created: Aug
// Revision: A
//
// Description: This file contains the class implementation for ASD_SMA_3K. 

#ifndef ASD_SMA_3K_h
#define ASD_SMA_3K_h

/*
ASD_SMA_3K is written by Eng. Luca Aceto (luca.aceto@unich.it), University of Chieti-Pescara, InGeo department in collaboration with ASDEA Software Technology: https://asdeasoft.net


This material is a modified version of Self Centering Material written by JAE at Oct 2007
With ASD_SMA_3K it is possible to replicate the behavior of SMA material with a different unloading stiffness (k3)

        k1 = Load stiffness
        k2 = Post-Activation Stiffness (0<$k2<$k1)
        k3 = Unload stiffness
        sigAct = Forward Activation Stress/Force
        beta= Ratio of Forward to Reverse Activation Stress/Force


ASD_SMA_3K matTag? k1? k2? k3? sigF? beta?

*/

#include <UniaxialMaterial.h>

class ASD_SMA_3K : public UniaxialMaterial
{
  public:
    ASD_SMA_3K(int tag, double k1, double k2, double k3,
		      double ActF, double beta);
    ASD_SMA_3K();
    ~ASD_SMA_3K();

    const char *getClassType(void) const {return "ASD_SMA_3K";}

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return k1;}

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
    double k3;		// Unloading Stiffness
    double ActF;	// Activation Stress/Force
    double beta;	// Flag-Shape Parameter
    
    // Extra Calculated Material Parameters
    double ActDef;	// Actvation Strain/Deformation
    
    // Committed history variables
    double CactivStrainPos;	// Committed activation strain (Pos Quad)
    double CactivStrainNeg;	// Committed activation strain (Neg Quad)
    double CupperStrainPos;	// Committed upper activation strain (Pos Quad)
    double ClowerStrainPos;	// Committed lower activation strain (Pos Quad)
    double CupperStressPos;	// Committed upper activation stress (Pos Quad)
    double ClowerStressPos;	// Committed lower activation stress (Pos Quad)
    double CupperStrainNeg;	// Committed upper activation strain (Neg Quad)
    double ClowerStrainNeg;	// Committed lower activation strain (Neg Quad)
    double CupperStressNeg;	// Committed upper activation stress (Neg Quad)
    double ClowerStressNeg;	// Committed lower activation stress (Neg Quad)

    double CLastStrain;
    
    // Trial history variables
    double TactivStrainPos;	// Trial activation strain (Pos Quad)
    double TactivStrainNeg;	// Trial activation strain (Neg Quad)
    double TupperStrainPos;	// Trial upper activation strain (Pos Quad)
    double TlowerStrainPos;	// Trial lower activation strain (Pos Quad)
    double TupperStressPos;	// Trial upper activation stress (Pos Quad)
    double TlowerStressPos;	// Trial lower activation stress (Pos Quad)
    double TupperStrainNeg;	// Trial upper activation strain (Neg Quad)
    double TlowerStrainNeg;	// Trial lower activation strain (Neg Quad)
    double TupperStressNeg;	// Trial upper activation stress (Neg Quad)
    double TlowerStressNeg;	// Trial lower activation stress (Neg Quad)

    double TLastStrain;
    
    // Trial state variables
    double Tstrain;		// Trial strain
    double Tstress;		// Trial stress
    double Ttangent;	// Trial tangent
    
    // Committed State Variables
    double Cstrain;		// Committed Strain
    double Cstress;		// Committed Strain
    double Ctangent;		// Committed Strain

    int CNo_k2_Pos;
    int CNo_k2_Neg;
    int CNo_Y_Pos;
    int CNo_Y_Neg;

    int No_k2_Pos;
    int No_k2_Neg;
    int No_Y_Pos;
    int No_Y_Neg;  

};


#endif

