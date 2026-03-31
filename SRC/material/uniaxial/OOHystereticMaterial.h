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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: July 2000
//
// Description: This file contains the class definition for 
// OOHystereticMaterial.  OOHystereticMaterial provides the implementation
// of a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  This
// is a modified implementation of Hyster2.f90 by Filippou.

#ifndef OOHystereticMaterial_h
#define OOHystereticMaterial_h

#include <UniaxialMaterial.h>

class HystereticBackbone;
class StiffnessDegradation;
class UnloadingRule;
class StrengthDegradation;

class SectionForceDeformation;

class OOHystereticMaterial : public UniaxialMaterial
{
 public:
  OOHystereticMaterial(int tag,
		       HystereticBackbone &bb,
		       UnloadingRule &unl,			
		       StiffnessDegradation &stiff,
		       StrengthDegradation &str,
		       double pinchX = 0.0, double pinchY = 1.0);
  OOHystereticMaterial(int tag,
		       HystereticBackbone &posBB,
		       HystereticBackbone &negBB,
		       UnloadingRule &posUnl,
		       UnloadingRule &negUnl,
		       StiffnessDegradation &posStiff,
		       StiffnessDegradation &negStiff,
		       StrengthDegradation &posStr,
		       StrengthDegradation &negStr,
		       double pinchX = 0.0, double pinchY = 1.0);
  OOHystereticMaterial();
  ~OOHystereticMaterial();

  const char *getClassType(void) const {return "OOHystereticMaterial";}

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  UniaxialMaterial *getCopy(void);
  UniaxialMaterial *getCopy(SectionForceDeformation *s);
  
  int setVariable(const char *argv, Information &info);
  int getVariable(int varID, Information &info);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
 protected:
  
 private:
  OOHystereticMaterial(int tag,
		       HystereticBackbone &pos,
		       HystereticBackbone &neg,
		       UnloadingRule &posUnl,
		       UnloadingRule &negUnl,
		       StiffnessDegradation &posStiff,
		       StiffnessDegradation &negStiff,
		       StrengthDegradation &posStr,
		       StrengthDegradation &negStr,
		       double pinchX, double pinchY,
		       SectionForceDeformation *s);
  
  SectionForceDeformation *theSection;
  
  // Pinching parameters
  double pinchX;		// Deformation pinching
  double pinchY;		// Force pinching
  
  // Pointers to backbones
  HystereticBackbone *posEnvelope;
  HystereticBackbone *negEnvelope;
  
  // Initial tangents
  double E1p, E1n;
  
  // Yield strains
  double rot1p, rot1n;
  
  // Unloading from positive backbone
  UnloadingRule *posUnlRule;
  int posUnlRuleID;
  
  // Unloading from negative backbone
  UnloadingRule *negUnlRule;
  int negUnlRuleID;
  
  // Stiffness degradation of positive backbone
  StiffnessDegradation *posStfDegr;
  int posStfDegrID;

  // Stiffness degradation of negative backbone
  StiffnessDegradation *negStfDegr;
  int negStfDegrID;
  
  // Strength degradation of positive backbone
  StrengthDegradation *posStrDegr;
  int posStrDegrID;
  
  // Strength degradation of negative backbone
  StrengthDegradation *negStrDegr;
  int negStrDegrID;
  
  // Trial history variables
  double TrotMax;
  double TrotMin;
  double TtargMin;
  double TtargMax;
  double TrotPu;
  double TrotNu;
  double TenergyD;
  int TloadIndicator;
  
  // Trial state variables
  double Ttangent;
  double Tstress;
  double Tstrain;
  
  // Converged history variables
  double CrotMax;
  double CrotMin;
  double CtargMin;
  double CtargMax;
  double CrotPu;
  double CrotNu;
  double CenergyD;
  int CloadIndicator;
  
  // Converged state variables
  double Cstress;
  double Cstrain;
  
  void positiveIncrement(double dStrain);
  void negativeIncrement(double dStrain);
  
  bool firstIter;
};

#endif
