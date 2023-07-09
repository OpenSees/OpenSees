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
// Description: This file contains the implementation of 
// OOHystereticMaterial.  OOHystereticMaterial is
// a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  This
// is a modified implementation of Hyster2.f90 by Filippou.

#include <OOHystereticMaterial.h>
#include <HystereticBackbone.h>
#include <StiffnessDegradation.h>
#include <UnloadingRule.h>
#include <StrengthDegradation.h>
#include <SectionForceDeformation.h>
#include <G3Globals.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <Information.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

void *
OPS_OOHystereticMaterial(void)
{
  UniaxialMaterial *theMaterial = 0;
  
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial OOHysteretic tag? bTag+? unlRulTag+? stfDegTag+? strDegTag+? "
	   << "<bTag-? unlRulTag-? stfDegTag-? strDegTag-?> <pinchX? pinchY?>" << endln;
    return 0;
  }
  
  int tag;
  int bTagPos, bTagNeg;
  int unlTagPos, unlTagNeg;
  int stfTagPos, stfTagNeg;
  int strTagPos, strTagNeg;
  double pinchX = 0.0;
  double pinchY = 1.0;

  int argc = OPS_GetNumRemainingInputArgs();
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &bTagPos) != 0) {
    opserr << "WARNING invalid bTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &unlTagPos) != 0) {
    opserr << "WARNING invalid unlRulTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &stfTagPos) != 0) {
    opserr << "WARNING invalid stfDegTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &strTagPos) != 0) {
    opserr << "WARNING invalid strDegTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }

  if (argc == 7) {
    if (OPS_GetDoubleInput(&numData, &pinchX) != 0) {
      opserr << "WARNING invalid pinchX\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetDoubleInput(&numData, &pinchY) != 0) {
      opserr << "WARNING invalid pinchY\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
  }

  if (argc > 8) {
    if (OPS_GetIntInput(&numData, &bTagNeg) != 0) {
      opserr << "WARNING invalid bTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numData, &unlTagNeg) != 0) {
      opserr << "WARNING invalid unlRulTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numData, &stfTagNeg) != 0) {
      opserr << "WARNING invalid stfDegTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numData, &strTagNeg) != 0) {
      opserr << "WARNING invalid strDegTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
  }
  if (argc == 11) {
    if (OPS_GetDoubleInput(&numData, &pinchX) != 0) {
      opserr << "WARNING invalid pinchX\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetDoubleInput(&numData, &pinchY) != 0) {
      opserr << "WARNING invalid pinchY\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
  }
  
  HystereticBackbone *posBB = OPS_getHystereticBackbone(bTagPos);
  
  if (posBB == 0) {
    opserr << "WARNING backbone does not exist\n";
    opserr << "backbone: " << bTagPos; 
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }
  UnloadingRule *posUnl = OPS_getUnloadingRule(unlTagPos);
  
  if (posUnl == 0) {
    opserr << "WARNING unloadingRule does not exist\n";
    opserr << "unloadingRule: " << unlTagPos; 
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }
  
  StiffnessDegradation *posStf = OPS_getStiffnessDegradation(stfTagPos);
  
  if (posStf == 0) {
    opserr << "WARNING stiffnessDegradation does not exist\n";
    opserr << "stiffnessDegradation: " << stfTagPos; 
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }

  StrengthDegradation *posStr = OPS_getStrengthDegradation(strTagPos);

  if (posStr == 0) {
    opserr << "WARNING strengthDegradation does not exist\n";
    opserr << "strengthDegradation: " << strTagPos; 
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }
  if (argc > 8) {
    HystereticBackbone *negBB = OPS_getHystereticBackbone(bTagNeg);
    
    if (negBB == 0) {
      opserr << "WARNING backbone does not exist\n";
      opserr << "backbone: " << bTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }
    
    UnloadingRule *negUnl = OPS_getUnloadingRule(unlTagNeg);
    
    if (negUnl == 0) {
      opserr << "WARNING unloadingRule does not exist\n";
      opserr << "unloadingRule: " << unlTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }
    
    StiffnessDegradation *negStf = OPS_getStiffnessDegradation(stfTagNeg);
    
    if (negStf == 0) {
      opserr << "WARNING stiffnessDegradation does not exist\n";
      opserr << "stiffnessDegradation: " << stfTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }
    
    StrengthDegradation *negStr = OPS_getStrengthDegradation(strTagNeg);
    
    if (negStr == 0) {
      opserr << "WARNING strengthDegradation does not exist\n";
      opserr << "strengthDegradation: " << strTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }
    
    theMaterial = 
      new OOHystereticMaterial(tag, *posBB, *negBB, *posUnl, *negUnl, 
			       *posStf, *negStf, *posStr, *negStr,
			       pinchX, pinchY);
  }
  else {
    theMaterial = 
      new OOHystereticMaterial(tag, *posBB, *posUnl, *posStf, *posStr,
			       pinchX, pinchY);
  }

  //opserr << "OOHysteretic " << bTagPos << ' ' << unlTagPos << ' ' << stfTagPos << ' ' << strTagPos << ' ' << pinchX << ' ' << pinchY << endln;
  //opserr << "\t" << bTagNeg << ' ' << unlTagNeg << ' ' << stfTagNeg << ' ' << strTagNeg << ' ' << pinchX << ' ' << pinchY << endln;
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type OOHystereticMaterial\n";
    return 0;
  }

  return theMaterial;
  
}
OOHystereticMaterial::OOHystereticMaterial(int tag,
					   HystereticBackbone &bb,
					   UnloadingRule &unl,
					   StiffnessDegradation &stf,
					   StrengthDegradation &str,
					   double px, double py):
  UniaxialMaterial(tag, MAT_TAG_OOHysteretic),
  theSection(0), pinchX(px), pinchY(py),
  posEnvelope(0), negEnvelope(0), E1p(0.0), E1n(0.0),
  posUnlRule(0), negUnlRule(0),
  posStfDegr(0), negStfDegr(0),
  posStrDegr(0), negStrDegr(0),
  firstIter(true)
{
  posEnvelope = bb.getCopy();
  
  if (posEnvelope == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of positive backbone" << endln;
  
  negEnvelope = bb.getCopy();
  
  if (negEnvelope == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of negative backbone" << endln;
  
  E1p = posEnvelope->getTangent(0.0);
  E1n = negEnvelope->getTangent(0.0);
  
  rot1p = posEnvelope->getYieldStrain();
  rot1n = -negEnvelope->getYieldStrain();
  
  
  posUnlRule = unl.getCopy(this);
  
  if (posUnlRule == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of unloading rule" << endln;
 
  Information tmpInfo; // This is a fix for change in interface -- MHS
 
  posUnlRuleID = this->setVariable(posUnlRule->getMeasure(), tmpInfo);
  
  negUnlRule = unl.getCopy(this);
  
  if (negUnlRule == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of unloading rule" << endln;
  
  negUnlRule->setNegative(true);
  negUnlRuleID = this->setVariable(negUnlRule->getMeasure(), tmpInfo);
  
  
  posStfDegr = stf.getCopy(this);
  
  if (posStfDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of stiffness degradation" << endln;
  
  posStfDegrID = this->setVariable(posStfDegr->getMeasure(), tmpInfo);
  
  negStfDegr = stf.getCopy(this);
  
  if (negStfDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of stiffness degradation" << endln;
  
  negStfDegr->setNegative(true);
  negStfDegrID = this->setVariable(negStfDegr->getMeasure(), tmpInfo);
  
  
  posStrDegr = str.getCopy(this);
  
  if (posStrDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of strength degradation" << endln;
  
  posStrDegrID = this->setVariable(posStrDegr->getMeasure(), tmpInfo);
  
  negStrDegr = str.getCopy(this);
  
  if (negStrDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of strength degradation" << endln;
  
  negStrDegr->setNegative(true);
  negStrDegrID = this->setVariable(negStrDegr->getMeasure(), tmpInfo);
  
  
  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
}

OOHystereticMaterial::OOHystereticMaterial():
  UniaxialMaterial(0, MAT_TAG_OOHysteretic),
  theSection(0), pinchX(0.0), pinchY(0.0),
  posEnvelope(0), negEnvelope(0), E1p(0.0), E1n(0.0), rot1p(0.0), rot1n(0.0),
  posUnlRule(0), negUnlRule(0),
  posStfDegr(0), negStfDegr(0),
  posStrDegr(0), negStrDegr(0)
{
  
}

OOHystereticMaterial::OOHystereticMaterial(int tag,
					   HystereticBackbone &posBB,
					   HystereticBackbone &negBB,
					   UnloadingRule &posUnl,
					   UnloadingRule &negUnl,
					   StiffnessDegradation &posStiff,
					   StiffnessDegradation &negStiff,
					   StrengthDegradation &posStr,
					   StrengthDegradation &negStr,
					   double px, double py):
  UniaxialMaterial(tag, MAT_TAG_OOHysteretic),
  theSection(0), pinchX(px), pinchY(py),
  posEnvelope(0), negEnvelope(0), E1p(0.0), E1n(0.0),
  posUnlRule(0), negUnlRule(0),
  posStfDegr(0), negStfDegr(0),
  posStrDegr(0), negStrDegr(0),
  firstIter(true)
{
  posEnvelope = posBB.getCopy();
  
  if (posEnvelope == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of positive backbone" << endln;
  
  negEnvelope = negBB.getCopy();
  
  if (negEnvelope == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of negative backbone" << endln;
  
  E1p = posEnvelope->getTangent(0.0);
  E1n = negEnvelope->getTangent(0.0);
  
  rot1p = posEnvelope->getYieldStrain();
  rot1n = -negEnvelope->getYieldStrain();
  
  
  posUnlRule = posUnl.getCopy(this);
  
  if (posUnlRule == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of unloading rule" << endln;
  
  Information tmpInfo; // This is a fix for change in interface -- MHS

  posUnlRuleID = this->setVariable(posUnlRule->getMeasure(), tmpInfo);
  
  negUnlRule = negUnl.getCopy(this);
  
  if (negUnlRule == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of unloading rule" << endln;
  
  negUnlRule->setNegative(true);
  negUnlRuleID = this->setVariable(negUnlRule->getMeasure(), tmpInfo);
  
  
  posStfDegr = posStiff.getCopy(this);
  
  if (posStfDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of stiffness degradation" << endln;
  
  posStfDegrID = this->setVariable(posStfDegr->getMeasure(), tmpInfo);
  
  negStfDegr = negStiff.getCopy(this);
  
  if (negStfDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of stiffness degradation" << endln;
  
  negStfDegr->setNegative(true);
  negStfDegrID = this->setVariable(negStfDegr->getMeasure(), tmpInfo);
  
  
  posStrDegr = posStr.getCopy(this);
  
  if (posStrDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of strength degradation" << endln;
  
  posStrDegrID = this->setVariable(posStrDegr->getMeasure(), tmpInfo);
  
  negStrDegr = negStr.getCopy(this);
  
  if (negStrDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of strength degradation" << endln;
  
  negStrDegr->setNegative(true);
  negStrDegrID = this->setVariable(negStrDegr->getMeasure(), tmpInfo);
  
  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
}

OOHystereticMaterial::OOHystereticMaterial(int tag,
					   HystereticBackbone &posBB,
					   HystereticBackbone &negBB,
					   UnloadingRule &posUnl,
					   UnloadingRule &negUnl,
					   StiffnessDegradation &posStiff,
					   StiffnessDegradation &negStiff,
					   StrengthDegradation &posStr,
					   StrengthDegradation &negStr,
					   double px, double py,
					   SectionForceDeformation *s):
  UniaxialMaterial(tag, MAT_TAG_OOHysteretic),
  theSection(s), pinchX(px), pinchY(py),
  posEnvelope(0), negEnvelope(0), E1p(0.0), E1n(0.0),
  posUnlRule(0), negUnlRule(0),
  posStfDegr(0), negStfDegr(0),
  posStrDegr(0), negStrDegr(0),
  firstIter(true)
{
  posEnvelope = posBB.getCopy();
  
  if (posEnvelope == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of positive backbone" << endln;
  
  negEnvelope = negBB.getCopy();
  
  if (negEnvelope == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of negative backbone" << endln;
  
  E1p = posEnvelope->getTangent(0.0);
  E1n = negEnvelope->getTangent(0.0);
  
  rot1p = posEnvelope->getYieldStrain();
  rot1n = -negEnvelope->getYieldStrain();
  
  
  posUnlRule = posUnl.getCopy(this);
  
  if (posUnlRule == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of unloading rule" << endln;
  
  Information tmpInfo; // This is a fix for change in interface -- MHS

  posUnlRuleID = this->setVariable(posUnlRule->getMeasure(), tmpInfo);
  
  negUnlRule = negUnl.getCopy(this);
  
  if (negUnlRule == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of unloading rule" << endln;
  
  negUnlRule->setNegative(true);
  negUnlRuleID = this->setVariable(negUnlRule->getMeasure(), tmpInfo);
    
  posStfDegr = posStiff.getCopy(this);
  
  if (posStfDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of stiffness degradation" << endln;
  
  posStfDegrID = this->setVariable(posStfDegr->getMeasure(), tmpInfo);
  
  negStfDegr = negStiff.getCopy(this);
  
  if (negStfDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of stiffness degradation" << endln;
  
  negStfDegr->setNegative(true);
  negStfDegrID = this->setVariable(negStfDegr->getMeasure(), tmpInfo);
  
  
  posStrDegr = posStr.getCopy(this);
  
  if (posStrDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of strength degradation" << endln;
  
  posStrDegrID = this->setVariable(posStrDegr->getMeasure(), tmpInfo);
  
  negStrDegr = negStr.getCopy(this);
  
  if (negStrDegr == 0)
    opserr << "OOHystereticMaterial::OOHystereticMaterial -- failed to get copy of strength degradation" << endln;
  
  negStrDegr->setNegative(true);
  negStrDegrID = this->setVariable(negStrDegr->getMeasure(), tmpInfo);
  
  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
}

OOHystereticMaterial::~OOHystereticMaterial()
{
  if (posEnvelope != 0)
    delete posEnvelope;
  
  if (negEnvelope != 0)
    delete negEnvelope;
  
  if (posUnlRule != 0)
    delete posUnlRule;
  
  if (negUnlRule != 0)
    delete negUnlRule;
  
  if (posStfDegr != 0)
    delete posStfDegr;
  
  if (negStfDegr != 0)
    delete negStfDegr;
  
  if (posStrDegr != 0)
    delete posStrDegr;
  
  if (negStrDegr != 0)
    delete negStrDegr;
}

int
OOHystereticMaterial::setTrialStrain(double strain, double strainRate)
{
  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TtargMax = CtargMax;
  TtargMin = CtargMin;
  TenergyD = CenergyD;
  TrotNu = CrotNu;
  TrotPu = CrotPu;
  TloadIndicator = CloadIndicator;
  
  Tstrain = strain;
  double dStrain = Tstrain - Cstrain;
  
  if (TloadIndicator == 0)
    TloadIndicator = (dStrain < 0.0) ? 2 : 1;
  
  //if (firstIter) {
  //	Ttangent = E1p;
  //	Tstress = Cstress + Ttangent*dStrain;
  //	firstIter = false;
  //	return 0;
  //}
  
  double tmp;
  
  if (Tstrain > CtargMax) {
    TrotMax = Tstrain;
    TtargMax = Tstrain;
    tmp = posStrDegr->getValue();
    Ttangent = tmp*posEnvelope->getTangent(Tstrain);
    Tstress = tmp*posEnvelope->getStress(Tstrain);
    TenergyD = CenergyD + 0.5*(Cstress+Tstress)*dStrain;
  }
  else if (Tstrain < CtargMin) {
    TrotMin = Tstrain;
    TtargMin = Tstrain;
    tmp = negStrDegr->getValue();
    Ttangent = tmp*negEnvelope->getTangent(-Tstrain);
    Tstress = -tmp*negEnvelope->getStress(-Tstrain);
    TenergyD = CenergyD + 0.5*(Cstress+Tstress)*dStrain;
  }
  else {
    if (dStrain < -DBL_EPSILON)
      negativeIncrement(dStrain);
    else if (dStrain > DBL_EPSILON)
      positiveIncrement(dStrain);
  }
  
  //(dStrain < 0.0) ? negativeIncrement(dStrain) : positiveIncrement(dStrain);
  //cerr << dStrain << ' ' << posStfDegr->getValue() << endl;
  
  return 0;
}

double
OOHystereticMaterial::getStrain(void)
{
  return Tstrain;
}

double
OOHystereticMaterial::getStress(void)
{
  return Tstress;
}

double
OOHystereticMaterial::getTangent(void)
{
  return Ttangent;
}

double
OOHystereticMaterial::getInitialTangent(void)
{
  return E1p;
}

void
OOHystereticMaterial::positiveIncrement(double dStrain)
{
  // Get current degradation values
  double vun = negUnlRule->getValue();
  double vkp = posStfDegr->getValue();
  double vsp = posStrDegr->getValue();
  
  if (TloadIndicator == 2) {
    TloadIndicator = 1;
    
    if (Cstress <= 0.0) {
      
      Information info;
      
      // Unloading rule
      this->getVariable(negUnlRuleID,info);
      negUnlRule->setTrialMeasure(info.theDouble);
      
      // Stiffness degradation
      this->getVariable(posStfDegrID,info);
      posStfDegr->setTrialMeasure(info.theDouble);
      
      // Strength degradation
      this->getVariable(posStrDegrID,info);
      posStrDegr->setTrialMeasure(info.theDouble);
      
      // Get updated degradation values
      vun = negUnlRule->getValue();
      vkp = posStfDegr->getValue();
      vsp = posStrDegr->getValue();
      
      TrotNu = Cstrain - Cstress/(E1n*vun);
      TtargMax = TtargMax*vkp;
      //TtargMax = CtargMax*vkp;
    }
  }
  
  double vup = posUnlRule->getValue();
  
  //TrotMax = (TrotMax > rot1p) ? TrotMax : rot1p;
  if (TrotMax < rot1p) {
    TrotMax = rot1p;
    TtargMax = rot1p;
  }
  
  double maxmom = vsp*posEnvelope->getStress(TtargMax);
  double E = negEnvelope->getTangent(-CrotMin);
  double rotlim = (E < 0) ? CrotMin + negEnvelope->getStress(-CrotMin)/E : NEG_INF_STRAIN;
  if (rotlim > NEG_INF_STRAIN && negEnvelope->getStress(-rotlim) > 0.0)
    rotlim = NEG_INF_STRAIN;
  double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;
  
  double rotmp1 = rotrel + pinchY*(TtargMax-rotrel);
  double rotmp2 = TtargMax - (1.0-pinchY)*maxmom/(vup*E1p);
  double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
  
  double tmpmo1;
  double tmpmo2;
  
  if (Tstrain < TrotNu) {
    Ttangent = E1n*vun;
    Tstress = Cstress + Ttangent*dStrain;
    if (Tstress >= 0.0) {
      Tstress = 0.0;
      Ttangent = 0.0;
    }
  }
  
  else if (Tstrain >= TrotNu && Tstrain < rotch) {
    if (Tstrain <= rotrel) {
      Tstress = 0.0;
      Ttangent = 0.0;
    }
    else {
      Ttangent = maxmom*pinchY/(rotch-rotrel);
      tmpmo1 = Cstress + E1n*vun*dStrain;
      tmpmo2 = (Tstrain-rotrel)*Ttangent;
      if (tmpmo1 < tmpmo2) {
	Tstress = tmpmo1;
	Ttangent = E1n*vun;
      }
      else {
	Tstress = tmpmo2;
      }
    }
  }
  
  else {
    Ttangent = (1.0-pinchY)*maxmom/(TtargMax-rotch);
    tmpmo1 = Cstress + E1p*vup*dStrain;
    tmpmo2 = pinchY*maxmom + (Tstrain-rotch)*Ttangent;
    if (tmpmo1 < tmpmo2) {
      Tstress = tmpmo1;
      Ttangent = E1p*vup;
    }
    else {
      Tstress = tmpmo2;
    }
  }
  
  if (TloadIndicator != CloadIndicator)
    TenergyD = 0.0;
  else
    TenergyD = CenergyD + 0.5*(Cstress+Tstress)*dStrain;
}

void
OOHystereticMaterial::negativeIncrement(double dStrain)
{
  // Get current degradation values
  double vup = posUnlRule->getValue();
  double vkn = negStfDegr->getValue();
  double vsn = negStrDegr->getValue();
  
  //cerr << TtargMin << endl;
  
  if (TloadIndicator == 1) {
    TloadIndicator = 2;
    
    if (Cstress >= 0.0) {
      
      Information info;
      
      // Unloading rule
      this->getVariable(posUnlRuleID,info);
      posUnlRule->setTrialMeasure(info.theDouble);
      
      // Stiffness degradation
      this->getVariable(negStfDegrID,info);
      negStfDegr->setTrialMeasure(info.theDouble);
      
      // Strength degradation
      this->getVariable(negStrDegrID,info);
      negStrDegr->setTrialMeasure(info.theDouble);
      
      // Get updated degradation values
      vup = posUnlRule->getValue();
      vkn = negStfDegr->getValue();
      vsn = negStrDegr->getValue();
      
      TrotPu = Cstrain - Cstress/(E1p*vup);      
      TtargMin = TtargMin*vkn;
      //TtargMin = CtargMin*vkn;
      //cerr << "A: " << TtargMin << endl;
    }
  }
  
  double vun = negUnlRule->getValue();
  
  //TrotMin = (TrotMin < rot1n) ? TrotMin : rot1n;
  if (TrotMin > rot1n) {
    TrotMin = rot1n;
    TtargMin = rot1n;
  }
  
  //cerr << "B: " << TtargMin << endl;
  
  double minmom = -vsn*negEnvelope->getStress(-TtargMin);
  double E = posEnvelope->getTangent(CrotMax);
  double rotlim = (E < 0) ? CrotMax - posEnvelope->getStress(CrotMax)/E : POS_INF_STRAIN;
  if (rotlim < POS_INF_STRAIN && posEnvelope->getStress(rotlim) > 0.0)
    rotlim = POS_INF_STRAIN;
  
  double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;
  
  double rotmp1 = rotrel + pinchY*(TtargMin-rotrel);
  double rotmp2 = TtargMin - (1.0-pinchY)*minmom/(vun*E1n);
  double rotch = rotmp1 + (rotmp2-rotmp1)*pinchX;
  
  double tmpmo1;
  double tmpmo2;
  
  if (Tstrain > TrotPu) {
    Ttangent = E1p*vup;
    Tstress = Cstress + Ttangent*dStrain;
    if (Tstress <= 0.0) {
      Tstress = 0.0;
      Ttangent = 0.0;
    }
  }
  
  else if (Tstrain <= TrotPu && Tstrain > rotch) {
    if (Tstrain >= rotrel) {
      Tstress = 0.0;
      Ttangent = 0.0;
    }
    else {
      Ttangent = minmom*pinchY/(rotch-rotrel);
      tmpmo1 = Cstress + E1p*vup*dStrain;
      tmpmo2 = (Tstrain-rotrel)*Ttangent;
      if (tmpmo1 > tmpmo2) {
	Tstress = tmpmo1;
	Ttangent = E1p*vup;
      }
      else {
	Tstress = tmpmo2;
      }
    }
  }
  
  else {
    Ttangent = (1.0-pinchY)*minmom/(TtargMin-rotch);
    tmpmo1 = Cstress + E1n*vun*dStrain;
    tmpmo2 = pinchY*minmom + (Tstrain-rotch)*Ttangent;
    if (tmpmo1 > tmpmo2) {
      Tstress = tmpmo1;
      Ttangent = E1n*vun;
    }
    else {
      Tstress = tmpmo2;
    }
  }
  
  if (TloadIndicator != CloadIndicator)
    TenergyD = 0.0;
  else
    TenergyD = CenergyD + 0.5*(Cstress+Tstress)*dStrain;
}

int
OOHystereticMaterial::commitState(void)
{
  CrotMax = TrotMax;
  CrotMin = TrotMin;
  CtargMax = TtargMax;
  CtargMin = TtargMin;
  CrotPu = TrotPu;
  CrotNu = TrotNu;
  CenergyD = TenergyD;
  CloadIndicator = TloadIndicator;
  
  Cstress = Tstress;
  Cstrain = Tstrain;
  
  //firstIter = true;
  
  int err = 0;
  
  err += posUnlRule->commitState();
  err += negUnlRule->commitState();
  err += posStfDegr->commitState();
  err += negStfDegr->commitState();
  err += posStrDegr->commitState();
  err += negStrDegr->commitState();

  return err;
}

int
OOHystereticMaterial::revertToLastCommit(void)
{
  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TtargMax = CtargMax;
  TtargMin = CtargMin;
  TrotPu = CrotPu;
  TrotNu = CrotNu;
  TenergyD = CenergyD;
  TloadIndicator = CloadIndicator;

  Tstress = Cstress;
  Tstrain = Cstrain;
  
  firstIter = true;
  
  int err = 0;
  
  err += posUnlRule->revertToLastCommit();
  err += negUnlRule->revertToLastCommit();
  err += posStfDegr->revertToLastCommit();
  err += negStfDegr->revertToLastCommit();
  err += posStrDegr->revertToLastCommit();
  err += negStrDegr->revertToLastCommit();
  
  return err;
}

int
OOHystereticMaterial::revertToStart(void)
{
  CrotMax = 0.0;
  CrotMin = 0.0;
  CtargMax = 0.0;
  CtargMin = 0.0;
  CrotPu = 0.0;
  CrotNu = 0.0;
  CenergyD = 0.0;
  CloadIndicator = 0;
  
  Cstress = 0.0;
  Cstrain = 0.0;
  
  Ttangent = E1p;
  
  firstIter = true;
  
  int err = 0;
  
  err += posUnlRule->revertToStart();
  err += negUnlRule->revertToStart();
  err += posStfDegr->revertToStart();
  err += negStfDegr->revertToStart();
  err += posStrDegr->revertToStart();
  err += negStrDegr->revertToStart();
  
  return err;
}

UniaxialMaterial*
OOHystereticMaterial::getCopy(void)
{
  OOHystereticMaterial *theCopy = 
    new OOHystereticMaterial (this->getTag(),
			      *posEnvelope, *negEnvelope, 
			      *posUnlRule, *negUnlRule,
			      *posStfDegr, *negStfDegr,
			      *posStrDegr, *negStrDegr,
			      pinchX, pinchY);
  
  theCopy->CrotMax = CrotMax;
  theCopy->CrotMin = CrotMin;
  theCopy->CtargMax = CtargMax;
  theCopy->CtargMin = CtargMin;
  theCopy->CrotPu = CrotPu;
  theCopy->CrotNu = CrotNu;
  theCopy->CenergyD = CenergyD;
  theCopy->CloadIndicator = CloadIndicator;
  theCopy->Cstress = Cstress;
  theCopy->Cstrain = Cstrain;
  
  theCopy->Ttangent = Ttangent;
  
  return theCopy;
}

UniaxialMaterial*
OOHystereticMaterial::getCopy(SectionForceDeformation *s)
{
  OOHystereticMaterial *theCopy = 
    new OOHystereticMaterial (this->getTag(),
			      *posEnvelope, *negEnvelope, 
			      *posUnlRule, *negUnlRule,
			      *posStfDegr, *negStfDegr,
			      *posStrDegr, *negStrDegr,
			      pinchX, pinchY, s);
  
  theCopy->CrotMax = CrotMax;
  theCopy->CrotMin = CrotMin;
  theCopy->CtargMax = CtargMax;
  theCopy->CtargMin = CtargMin;
  theCopy->CrotPu = CrotPu;
  theCopy->CrotNu = CrotNu;
  theCopy->CenergyD = CenergyD;
  theCopy->CloadIndicator = CloadIndicator;
  theCopy->Cstress = Cstress;
  theCopy->Cstrain = Cstrain;
  
  theCopy->Ttangent = Ttangent;
  
  return theCopy;
}

int
OOHystereticMaterial::setVariable(const char *argv, Information &info)
{
  if (strcmp(argv,"posDuctility") == 0)
    return 1;
  else if (strcmp(argv,"negDuctility") == 0)
    return 2;
  else if (strcmp(argv,"energyExcursion") == 0)
    return 3;
  else if (strcmp(argv,"yieldEnergy") == 0)
    return 4;
  else if (theSection) {
    int id = theSection->setVariable(argv, info);
    if (id >= 0 && id < 100)
      return id + 100;
    else
      return -1;
  }
  else
    return -1;
}

int
OOHystereticMaterial::getVariable(int varID, Information &info)
{
  switch (varID) {
  case 1:
    info = (Cstrain/rot1p);
    return 0;
  case 2:
    info = (Cstrain/rot1n);
    return 0;
  case 3:
    info = TenergyD;
    return 0;
  case 4:
    info = posEnvelope->getEnergy(rot1p) + negEnvelope->getEnergy(-rot1n);
    return 0;
  default:
    if (varID >= 100 && theSection != 0) {
      //return theSection->getVariable(varID-100, info);
      //return theSection->getVariable(argv, info);
      opserr << "OOHysteretic -- Not calling theSection->getVariable";
      return 0;
    }
    else
      return -1;
  }
}

int
OOHystereticMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
OOHystereticMaterial::recvSelf(int commitTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  return -1;
}
    
void
OOHystereticMaterial::Print(OPS_Stream &s, int flag)
{
  s << "OOHystereticMaterial, tag: " << this->getTag() << endln;
  s << "pinchX: " << pinchX << endln;
  s << "pinchY: " << pinchY << endln;
  s << "positive backbone, tag: " << posEnvelope->getTag() << endln;
  s << "negative backbone, tag: " << negEnvelope->getTag() << endln;
  s << "pos unloading rule, tag: " << posUnlRule->getTag() << endln;
  s << "neg unloading rule, tag: " << negUnlRule->getTag() << endln;
  s << "pos stiffness degradation, tag: " << posStfDegr->getTag() << endln;
  s << "neg stiffness degradation, tag: " << negStfDegr->getTag() << endln;
  s << "pos strength degradation, tag: " << posStrDegr->getTag() << endln;
  s << "neg strength degradation, tag: " << negStrDegr->getTag() << endln;
}
