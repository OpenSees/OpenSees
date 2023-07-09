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
// $Date: 2008-08-26 16:34:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SecantConcrete.h,v $

// Written: MHS
// Created: Dec 2001
//
// Description: This file contains the class definition for 
// SecantConcrete.  SecantConcrete provides the abstraction
// for a one-dimensional concrete model with a Kent-Park backbone,
// no tension strength, and secant unloading/reloading.

#ifndef SecantConcrete_h
#define SecantConcrete_h

#include <UniaxialMaterial.h>

class Matrix;

class SecantConcrete : public UniaxialMaterial
{
 public:
  SecantConcrete(int tag, double fc, double epsc, double epsu);
  SecantConcrete();
  ~SecantConcrete();
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain(void);          
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return 2.0*fc/epsc;}
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);        
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
  int setParameter(const char **argv, int argc, Parameter &info);
  int updateParameter(int parameterID, Information &info);
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int activateParameter(int parameterID);
  double getStressSensitivity(int gradIndex, bool conditional);
  double getInitialTangentSensitivity(int gradIndex);
  int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////
  
 protected:
  
 private:
  // Material parameters
  double fc;	// Elastic modulus
  double epsc;	// Yield stress
  double epsu;	// Isotropic hardening parameter
  
  // Committed history variables
  double CminStrain;	// Committed plastic strain
  
  // Trial history variables
  double TminStrain;	// Trial plastic strain
  
  // Trial state variables
  double Tstrain;		// Trial strain
  double Tstress;		// Trial stress
  double Ttangent;	// Trial tangent
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  Matrix *SHVs;
  // AddingSensitivity:END ///////////////////////////////////////////
  
  double getStressGradient(int gradIndex);
  int setStrainGradient(int gradIndex, double depsilondh);
  
  void backbone(double strain, double &stress, double &tangent);
  double backboneCondDeriv(double strain, int gradIndex);
  double backboneUncondDeriv(double strain, int gradIndex, double depsilondh);
};

#endif
