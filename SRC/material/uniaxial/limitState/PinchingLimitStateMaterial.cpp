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
// $Date: 2012/05/01 01:00:00 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/PinchingLimitStateMaterial.cpp,v $
                                                                        
// Written: MRL 
//
// Description: This file contains the class implementation for PinchingLimitStateMaterial.
//
// What: "@(#) PinchingLimitStateMaterial.C, revN/C"

#include <stdlib.h>

#include <elementAPI.h>
#include "PinchingLimitStateMaterial.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <ArrayOfTaggedObjects.h>
#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <DummyStream.h>
#include <Element.h>
#include <ElementResponse.h>
#include <elementAPI.h>
#include <G3Globals.h>
#include <Parameter.h>

static int numPinchingLimitStateMaterial = 0;

#ifndef fmax
#define fmax(a,b) (((a)>(b)) ? (a) : (b))
#endif

#ifndef fmin
#define fmin(a,b) (((a)<(b)) ? (a) : (b))
#endif
void *
OPS_PinchingLimitState(void)
{
  if (numPinchingLimitStateMaterial == 0) {
    numPinchingLimitStateMaterial++;
    opserr << "PinchingLimitStateMaterial unaxial material - Written by MRL UT Austin Copyright 2012 - Use at Your Peril\n";
  }

  UniaxialMaterial *theMaterial = 0;
  LimitCurve *theCurve = 0;	
  int argc = OPS_GetNumRemainingInputArgs();

  if (!(argc == 32 || argc == 21)) {
    opserr << "WARNING PinchingLimitStateMaterial -- insufficient arguments\n";
    opserr << "For direct input of limit state material want:\n\n";
    opserr << "uniaxialMaterial PinchingLimitStateMaterial matTag?\n";
    opserr << "nodeT? nodeB? driftAxis? Kelas? crvTyp? crvTag?\n";
    opserr << "YpinchUPN? YpinchRPN? XpinchRPN?\n";
    opserr << "YpinchUNP? YpinchRNP? XpinchRNP?\n";
    opserr << "dmgStrsLimE? dmgDispMax?\n?";
    opserr << "dmgE1? dmgE2? dmgE3? dmgE4? dmgELim?\n";
    opserr << "dmgR1? dmgR2? dmgR3? dmgR4? dmgRLim? dmgRCyc?\n";
    opserr << "dmgS1? dmgS2? dmgS3? dmgS4? dmgSLim? dmgSCyc?\n" << endln;
    
    opserr << "OR for calibrated limit state material want:\n\n";
    opserr << "uniaxialMaterial PinchingLimitStateMaterial matTag?\n";
    opserr << "nodeT? nodeB? driftAxis? Kelas? crvTyp? crvTag? eleTag?\n";
    opserr << "b? d? h? a? st? As? Acc? ld? db? rhot? f'c?\n";
    opserr << "fy? fyt?\n" << endln;
    
    return 0;
  }
  
  int iTagData[1]; 		
  int iNodeData[3]; 		
  double dKelasData[1];	
  int iCrvData[2]; 		
  double dpinchPN[3];	
  double dpinchNP[3];	
  double dDmgProp[2];	
  double dDmgEdata[5];
  double dDmgRdata[6];
  double dDmgSdata[6];
  int iEleTag[1];	
  double dPropData[13];
  int numData;
  
  numData = 1;
  if (OPS_GetIntInput(&numData, iTagData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid uniaxialMaterial matTag?\n" << endln;
    return 0;
  }
  numData = 3;
  if (OPS_GetIntInput(&numData, iNodeData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid nodeT? nodeB? driftAxis?\n" << endln;
    return 0;
  }
  Domain *theDomain = 0;
  theDomain = OPS_GetDomain();
  if (theDomain == 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- Pointer to Domain was not returned\n" << endln;
    return 0;
  }
  Node *theNodeT = 0;
  theNodeT = theDomain->getNode(iNodeData[0]);
  if (theNodeT == 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- nodeT with tag " << iNodeData[0] << " does not exist for uniaxialMaterial tag " << iTagData[0] << endln << endln;
    return 0;
  }	
  Node *theNodeB = 0;
  theNodeB = theDomain->getNode(iNodeData[1]);
  if (theNodeB == 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- nodeB with tag " << iNodeData[1] << " does not exist for uniaxialMaterial tag " << iTagData[0] << endln << endln;
    return 0;
  }	
  if (iNodeData[2] < 1 || iNodeData[2] > 3) {
    opserr << "WARNING PinchingLimitStateMaterial -- driftAxis is invalid\n";
    opserr << "driftAxis = 1 -- Drift along the x-axis\n";
    opserr << "driftAxis = 2 -- Drift along the y-axis\n";
    opserr << "driftAxis = 3 -- Drift along the z-axis\n";
    return 0;
  }
  numData = 1;
  if (OPS_GetDoubleInput(&numData, dKelasData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid Kelas?\n";
    return 0;	
  }
  if((dKelasData[0] < -4 || dKelasData[0] == 0) && argc == 23) {
    opserr << "WARNING PinchingLimitStateMaterial -- Kelas? is invalid\n";
    opserr << "Kelas = -4 -- Shear stiffness calculated assuming double curvature and shear springs top and bottom\n";
    opserr << "Kelas = -3 -- Shear stiffness calculated assuming double curvature and a shear spring at the bottom\n";
    opserr << "Kelas = -2 -- Shear stiffness calculated assuming single curvature and shear springs top and bottom\n";
    opserr << "Kelas = -1 -- Shear stiffness calculated assuming single curvature and a shear spring at the bottom\n";
    opserr << "Kelas > 0 -- Shear stiffness is the input value\n";
    return 0;
  } 
  else if (dKelasData[0] <= 0 && argc == 34/*39*/) {
    opserr << "WARNING PinchingLimitStateMaterial -- Kelas? is invalid\n";
    opserr << "Kelas must be greater than zero\n";
    return 0;
  }
  numData = 2;
  if (OPS_GetIntInput(&numData, iCrvData) != 0) {
    opserr << "WARNING PinchingLimitStateMaterial -- invalid crvTyp? crvTag?\n" << endln;
    return 0;
  }
  int crvTyp = iCrvData[0];
  int crvTag = iCrvData[1];
  if (crvTyp == 2) {
    theCurve = OPS_getLimitCurve(crvTag);
    if (theCurve == 0 && crvTyp != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- limit curve with tag " << crvTag << " not found for material tag " << iTagData[0] << endln << endln;
      return 0;
    }
  }
  if (crvTyp < 0 || crvTyp > 2) {
    opserr << "WARNING PinchingLimitStateMaterial --  crvTyp? is invalid\n";
    opserr << "crvType = 0 -- no limit curve\n";
    opserr << "crvType = 1 -- axial limit curve\n";
    opserr << "crvType = 2 -- shear limit curve\n" << endln;
    return 0;
  }
  if (crvTyp == 1) {
    opserr << "WARNING PinchingLimitStateMaterial -- Axial curve has not been implemented\n" << endln;
    return 0;
  }
  if (argc == 32) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dpinchPN) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid YpinchUPN? YpinchRPN? XpinchRPN?\n" << endln;
      return 0;	
    }
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dpinchNP) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid YpinchUNP? YpinchRNP? XpinchRNP?\n" << endln;
      return 0;	
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, dDmgProp) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgStrsLimE? dmgDispMax?\n" << endln;
      return 0;	
    }
    if(dDmgProp[0] < 0.0001)
      dDmgProp[0] = 0.0001;
    numData = 5;
    if (OPS_GetDoubleInput(&numData, dDmgEdata) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgE1? dmgE2? dmgE3? dmgE4? dmgELim?\n" << endln;
      return 0;	
    }
    numData = 6;
    if (OPS_GetDoubleInput(&numData, dDmgRdata) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgR1? dmgR2? dmgR3? dmgR4? dmgRLim? dmgRCyc?\n" << endln;
      return 0;	
    }
    numData = 6;
    if (OPS_GetDoubleInput(&numData, dDmgSdata) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid dmgS1? dmgS2? dmgS3? dmgS4? dmgSLim? dmgSCyc?\n" << endln;
      return 0;	
    }
    theMaterial = new PinchingLimitStateMaterial(iTagData[0], 
						 iNodeData[0], iNodeData[1], iNodeData[2], dKelasData[0], iCrvData[0], iCrvData[1],
						 dpinchPN[0],  dpinchPN[1],  dpinchPN[2],
						 dpinchNP[0],  dpinchNP[1],  dpinchNP[2],
						 dDmgProp[0],  dDmgProp[1],
						 dDmgEdata[0], dDmgEdata[1], dDmgEdata[2], dDmgEdata[3], dDmgEdata[4],
						 0.0         , 0.0           , 0.0         , 0.0         , 0.0         , 
						 dDmgRdata[0], dDmgRdata[1], dDmgRdata[2], dDmgRdata[3], dDmgRdata[4], dDmgRdata[5],
						 dDmgSdata[0], dDmgSdata[1], dDmgSdata[2], dDmgSdata[3], dDmgSdata[4], dDmgSdata[5],
						 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
						 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
						 theDomain, theNodeT, theNodeB, *theCurve, 0);													
    
    if (theMaterial == 0) {
      opserr << "WARNING could not create uniaxialMaterial with PinchinLimitState\n";
      return 0;
    }
    return theMaterial;
    
  } else {
    numData = 1;
    if (OPS_GetIntInput(&numData, iEleTag) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid eleTag?\n" << endln;
      return 0;
    }
    int eleTag = iEleTag[0];
    Element *theElement = 0;
    theElement = theDomain->getElement(eleTag);
    if (theElement == 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- Element with tag " << eleTag << " does not exist for uniaxialMaterial tag " << iTagData[0] << endln << endln;
      return 0;
    }
    numData = 13;
    if (OPS_GetDoubleInput(&numData, dPropData) != 0) {
      opserr << "WARNING PinchingLimitStateMaterial -- invalid b? d? h? a? st? As? Acc? ld? db? rhot? f'c? fy? fyt?\n" << endln;
      return 0;
    }
    theMaterial = new PinchingLimitStateMaterial(iTagData[0], 
						 iNodeData[0], iNodeData[1], iNodeData[2], dKelasData[0], iCrvData[0], iCrvData[1],
						 0.0,  0.0,  0.0,
						 0.0,  0.0,  0.0,
						 0.0,  0.0,
						 0.0, 0.0, 0.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
						 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
						 iEleTag[0],	fabs(dPropData[0]), fabs(dPropData[1]), fabs(dPropData[2]), fabs(dPropData[3]), fabs(dPropData[4]), fabs(dPropData[5]), fabs(dPropData[6]),
						 fabs(dPropData[7]), fabs(dPropData[8]), fabs(dPropData[9]), fabs(dPropData[10]), fabs(dPropData[11]), fabs(dPropData[12]),
						 theDomain, theNodeT, theNodeB, *theCurve, theElement);	
    
    if (theMaterial == 0) {
      opserr << "WARNING could not create uniaxialMaterial PinchingLimitState\n ";
      return 0;
    }
    return theMaterial;
  }  
}

PinchingLimitStateMaterial::PinchingLimitStateMaterial(int matTag,
						       int nodeT, int nodeB, int drftAx, double Kelas, int crvTyp, int crvTag,
						       double YpinchUPN, double YpinchRPN,   double XpinchRPN,
						       double YpinchUNP, double YpinchRNP,   double XpinchRNP,
						       double dmgStrsLimE, double dmgDispMax,
						       double dmgE1, double dmgE2, double dmgE3, double dmgE4, double dmgELim,
						       double dmgU1, double dmgU2, double dmgU3, double dmgU4, double dmgULim,
						       double dmgR1, double dmgR2, double dmgR3, double dmgR4,	double dmgRLim, double dmgRCyc,
						       double dmgS1, double dmgS2, double dmgS3, double dmgS4, double dmgSLim, double dmgSCyc, 
						       int eTag, double B, double D, double H, double A, double S, double AreaS, double AreaCC,
						       double Ldev, double DiaBar, double rhoTrans, double Fc,	double Fy, double Fytrans, 
						       Domain *theDom, Node *theNdT, Node *theNdB, LimitCurve &theCrv, Element *theEle)	
  :UniaxialMaterial(matTag,MAT_TAG_PinchingLimitStateMaterial),
   nodeTop(nodeT), nodeBot(nodeB), driftAxis(drftAx), E1(Kelas), curveType(crvTyp), curveTag(crvTag),
   YpinchUnloadPN(YpinchUPN), YpinchReloadPN(YpinchRPN), XpinchReloadPN(XpinchRPN),
   YpinchUnloadNP(YpinchUNP), YpinchReloadNP(YpinchRNP), XpinchReloadNP(XpinchRNP),
   dmgStressLimE(dmgStrsLimE), dmgDeflMax(dmgDispMax),
   dmgElastic1(dmgE1),  dmgElastic2(dmgE2),  dmgElastic3(dmgE3),  dmgElastic4(dmgE4),  dmgElasticLim(dmgELim),
   dmgUnload1(dmgU1),   dmgUnload2(dmgU2),   dmgUnload3(dmgU3),   dmgUnload4(dmgU4),   dmgUnloadLim(dmgULim),
   dmgReload1(dmgR1),   dmgReload2(dmgR2),   dmgReload3(dmgR3),   dmgReload4(dmgR4),   dmgReloadLim(dmgRLim), dmgReloadCyclic(dmgRCyc),
   dmgStrength1(dmgS1), dmgStrength2(dmgS2), dmgStrength3(dmgS3), dmgStrength4(dmgS4), dmgStrengthLim(dmgSLim), dmgStrengthCyclic(dmgSCyc),
   eleTag(eTag), b(B), d(D), h(H), a(A), st(S), As(AreaS), Acc(AreaCC),
   ld(Ldev), db(DiaBar), rhot(rhoTrans), fc(Fc), fy(Fy), fyt(Fytrans), 
   theDomain(theDom), theNodeT(theNdT), theNodeB(theNdB), theCurve(0), theElement(theEle)																		
{
  if (curveType != 0)
    theCurve = theCrv.getCopy();
  if (theCurve == 0 && curveType != 0) {
    opserr << "FATAL ERROR PinchingLimitStateMaterial -- out of memory, could not get a copy of the Limit Curve\n" << endln;
    exit(-1);
  }
  defineE1();
  revertToStart();
  revertToLastCommit();
}

PinchingLimitStateMaterial::PinchingLimitStateMaterial()
  :UniaxialMaterial(0,MAT_TAG_PinchingLimitStateMaterial), 	
   nodeTop(0), nodeBot(0), driftAxis(0), E1(0.0), curveType(0), curveTag(0),
   YpinchUnloadPN(0.0), YpinchReloadPN(0.0), XpinchReloadPN(0.0),
   YpinchUnloadNP(0.0), YpinchReloadNP(0.0), XpinchReloadNP(0.0),
   dmgStressLimE(0.0), dmgDeflMax(0.0),
   dmgElastic1(0.0),  dmgElastic2(0.0),  dmgElastic3(0.0),  dmgElastic4(0.0),  dmgElasticLim(0.0),
   dmgUnload1(0.0),   dmgUnload2(0.0),   dmgUnload3(0.0),   dmgUnload4(0.0),   dmgUnloadLim(0.0),
   dmgReload1(0.0),   dmgReload2(0.0),   dmgReload3(0.0),   dmgReload4(0.0),   dmgReloadLim(0.0), dmgReloadCyclic(0.0),
   dmgStrength1(0.0), dmgStrength2(0.0), dmgStrength3(0.0), dmgStrength4(0.0), dmgStrengthLim(0.0), dmgStrengthCyclic(0.0),
   eleTag(0), b(0.0), d(0.0), h(0.0), a(0.0), st(0.0), As(0.0), Acc(0.0),
   ld(0.0), db(0.0), rhot(0.0), fc(0.0), fy(0.0), fyt(0.0), 
   theDomain(0), theNodeT(0), theNodeB(0), theCurve(0), theElement(0)
{
  defineE1();
  revertToStart();
  revertToLastCommit();
}

PinchingLimitStateMaterial::~PinchingLimitStateMaterial()
{
  if (theCurve != 0)
    delete theCurve; 
}

int 
PinchingLimitStateMaterial::setTrialStrain(double strain, double strainRate)
{
  Tstrain = strain;
  TstrainFlex = getFlexDisp();
  TstrainGlobal = TstrainFlex + Tstrain;
  TstrainRate = strainRate;
  Tdu = Tstrain - Cstrain;
  if(Tdu == 0.0 || fabs(Tdu) > 1.0)
    return 0;
  
  if(Tstrain > TstrainMax)
    TstrainMax = Tstrain;
  else if(Tstrain < TstrainMin)
    TstrainMin = Tstrain;
  
  if(CstateFlag == 0)	
    {
      updateDamageE();
      Tstress = TdmgElasticE*Tstrain;
      Ttangent = TdmgElasticE;
    }
  else 
    {
      TstateFlag = getStateFlag();
      switch(TstateFlag)
	{
	case 1	:	Tstress = TdmgElasticE*Tstrain;
	  Ttangent = TdmgElasticE;
	  break;
	case -1	:	Tstress = TdmgElasticE*Tstrain;
	  Ttangent = TdmgElasticE;
	  break;
	case 2	:	Ttangent = Kdeg;
	  Tstress = Kdeg*fabs(Tstrain)+ TbKdegDmg;
	  break;
	case -2	:	Ttangent = Kdeg;
	  Tstress = -(Kdeg*fabs(Tstrain)+ TbKdegDmg);
	  break;
	case 3	:	Ttangent = 0.0001;
	  Tstress = Fres;	
	  break;
	case -3	:	Ttangent = 0.0001; 
	  Tstress = -Fres;	
	  break;
	case 4	:	definePinchingPN();
	  Ttangent = TdmgElasticE;
	  Tstress = TdmgElasticE*Tstrain + TbUnloadPN;
	  break;
	case -4	:	definePinchingNP();
	  Ttangent = TdmgElasticE;
	  Tstress = TdmgElasticE*Tstrain + TbUnloadNP;
	  break;
	case 5	:	if(CstateFlag == 6)
	    definePinchingPN();
	  Ttangent = TdmgElasticE;
	  Tstress = TdmgElasticE*Tstrain + TbUnloadPN;
	  break;
	case -5	:	if(CstateFlag == -6)
	    definePinchingNP();
	  Ttangent = TdmgElasticE;
	  Tstress = TdmgElasticE*Tstrain + TbUnloadNP;
	  break;
	case 6	:	if(CstateFlag != 6)	
	    {
	      updateDamageR();
	      TbReloadAfterUnloadPN = Cstress - TdmgReloadE*Cstrain;
	    }
	  Ttangent = TdmgReloadE;
	  Tstress = TdmgReloadE*Tstrain + TbReloadAfterUnloadPN;
	  checkEnvelope();
	  break;
	case -6	:	if(CstateFlag != -6)
	    {
	      updateDamageR();
	      TbReloadAfterUnloadNP = Cstress - TdmgReloadE*Cstrain;
	    }
	  Ttangent = TdmgReloadE;
	  Tstress = TdmgReloadE*Tstrain + TbReloadAfterUnloadNP;
	  checkEnvelope();
	  break;
	case 7	:	Ttangent = TpinchSlopePN;
	  Tstress = TpinchSlopePN*Tstrain + TpinchInterceptPN;
	  break;
	case -7	:	Ttangent = TpinchSlopeNP;
	  Tstress = TpinchSlopeNP*Tstrain + TpinchInterceptNP;
	  break;
	case 8	:	TpinchSlopeNP = (TpinchStressUnloadPN-Cstress)/(TpinchStrainUnloadPN-Cstrain);
	  TpinchInterceptNP = Cstress-TpinchSlopeNP*Cstrain;								
	  TstateFlag = -7;																
	  break;
	case -8	:	TpinchSlopePN = (TpinchStressUnloadNP-Cstress)/(TpinchStrainUnloadNP-Cstrain);	
	  TpinchInterceptPN = Cstress-TpinchSlopePN*Cstrain;							
	  TstateFlag = 7;									
	  break;
	case 9	:	updateDamageR();
	  TreloadInterceptPN = Cstress - TdmgReloadE*Cstrain;
	  Ttangent = TdmgReloadE;
	  Tstress = TdmgReloadE*Tstrain + TreloadInterceptPN;
	  TstateFlag = 10;
	  break;
	case -9	:	updateDamageR();
	  TreloadInterceptNP = Cstress - TdmgReloadE*Cstrain;
	  Ttangent = TdmgReloadE;
	  Tstress = TdmgReloadE*Tstrain + TreloadInterceptNP;
	  TstateFlag = -10;
	  break;
	case 10	:	Ttangent = TdmgReloadE;
	  Tstress = TdmgReloadE*Tstrain + TreloadInterceptPN;
	  checkEnvelope();
	  break;
	case -10 :	Ttangent = TdmgReloadE;
	  Tstress = TdmgReloadE*Tstrain + TreloadInterceptNP;
	  checkEnvelope();
	  break;
	case 11 :	definePinchingNP();
	  Ttangent = TdmgElasticE;
	  Tstress = TdmgElasticE*Tstrain + TbUnloadNP;
	  TstateFlag = -5;
	  break;
	case -11 :	definePinchingPN();
	  Ttangent = TdmgElasticE;
	  Tstress = TdmgElasticE*Tstrain + TbUnloadPN;
	  TstateFlag = 5;
	  break;
	default	:	Tstress = TdmgElasticE*Tstrain;
	  Ttangent = TdmgElasticE;
	}
      updateEnergy();
    }
  return 0;
}

double 
PinchingLimitStateMaterial::getStrain(void)
{
  return Tstrain;
}

double 
PinchingLimitStateMaterial::getStress(void)
{
  return Tstress;
}

double 
PinchingLimitStateMaterial::getTangent(void)
{
  return Ttangent;
}

int 
PinchingLimitStateMaterial::commitState(void)
{
	if((TstateFlag == 2 || TstateFlag == -2) && countGlobalEnv < 2)
	{
		countGlobalEnv++;
		slopeGlobalEnv = (fabs(Tstress) - fabs(Cstress))/(fabs(TstrainGlobal)-fabs(CstrainGlobal));
		interceptGlobalEnv = fabs(Cstress) - slopeGlobalEnv*fabs(CstrainGlobal);
		strainGlobalFresKdeg = (Fres - interceptGlobalEnv)/slopeGlobalEnv;
	}
	if(TstateFlag == 4 || TstateFlag == -4)
		countGlobalEnv = 2;
	if(TstateFlag == 3 || TstateFlag == -3) 
		resFlag = 1;
	if(TstateFlag == 4 || TstateFlag == -4 || (TstateFlag == 5 && CstateFlag == -10) || (TstateFlag == -5 && CstateFlag == 10))
	{
		if(TstateFlag == 4 || TstateFlag == -4)
		{
			strainUn = Cstrain;
			stressUn = Cstress;
			strainFlexRev = CstrainFlex;
		}
		else
		{
			strainUn = -strainUn;
			stressUn = -stressUn;
			strainFlexRev = -strainFlexRev;
		}
		strainUnDmg=strainUn;
		stressUnDmg=stressUn;
		strainFlexRevDmg=strainFlexRev;
	}
	if((TstateFlag == 6 && CstateFlag != 6) || (TstateFlag == -6 && CstateFlag != -6) || (TstateFlag == 10 && CstateFlag != 10) || (TstateFlag == -10 && CstateFlag != -10))
	{
		strainReload = Cstrain;
		stressReload = Cstress;
		Ereload = TdmgReloadE;
		if((stressUnDmg < 0.0 && Cstress < 0.0) || (stressUnDmg > 0.0 && Cstress > 0.0))
			slopeFlexPred = (stressUnDmg-Cstress)/(strainFlexRevDmg-CstrainFlex);
		else
			slopeFlexPred = (-stressUnDmg-Cstress)/(-strainFlexRevDmg-CstrainFlex);
		if(resFlag == 0 && ((stressUn > 0.0 && Cstress < 0.0) || (stressUn < 0.0 && Cstress > 0.0)))
		{	
			interceptGlobalEnv = interceptGlobalEnv-fabs(stressUn)*dmgStrengthCyclic;
			updateDamageS();
			strainGlobalFresKdeg = (Fres - interceptGlobalEnv)/slopeGlobalEnv;
		}
		interceptFlexPred = Cstress - slopeFlexPred*CstrainFlex;
	}
	if(TstateFlag == 6 || TstateFlag == -6 || TstateFlag == 10 || TstateFlag == -10)
	{
		double shiftFlexRe = getFlexShift();
		if(shiftFlexRe < 0.0001)
			shiftFlexRe = 0.0001;
		double slopeFlexRe = Tstress/fabs(Tstress)*(Tstress-stressReload)/(shiftFlexRe); 
		TdmgReloadE = 1/(1/slopeFlexRe+1/Ereload)*(1-dmgReloadCyclic);
		if(TstateFlag == 6)
			TbReloadAfterUnloadPN = Tstress - TdmgReloadE*Tstrain;
		else if(TstateFlag == -6)
			TbReloadAfterUnloadNP = Tstress - TdmgReloadE*Tstrain;
		else if(TstateFlag == 10)
			TreloadInterceptPN = Tstress - TdmgReloadE*Tstrain;
		else if(TstateFlag == -10)
			TreloadInterceptNP = Tstress - TdmgReloadE*Tstrain;

	}
	if(((CstateFlag == 6 && TstateFlag == 2) || (CstateFlag == -6 && TstateFlag == -2) || (CstateFlag == 10 && TstateFlag == -2) || (CstateFlag == -10 && TstateFlag == 2)) && resFlag == 0)
	{	
		TbKdegDmg = fabs(Tstress) - Kdeg*fabs(Tstrain);
		TstrainFresKdegDmg = (Fres - TbKdegDmg)/Kdeg;

	}
	if(curveType != 0 && TstateFlag == 0) 
	{
		TstateFlag = theCurve->checkElementState(Tstress);
		if (TstateFlag == 1)	
		{
			TstateFlag = (Tstress < 0.0) ? -1 : 1;	
			if (b != 0.0 && d != 0.0 && h != 0.0)
				defineTargetVars();
			defineBackbone();
		}
	}

	// Commit Trial Variables to memory
	Cstress = Tstress;
	Cstrain = Tstrain;
	Ctangent = Ttangent;	
	CstrainRate = TstrainRate;
	CstrainMax = TstrainMax;
	CstrainMin = TstrainMin;
	Cdu = Tdu;
	CpinchStressUnloadPN = TpinchStressUnloadPN;
	CpinchStrainUnloadPN = TpinchStrainUnloadPN;
	CpinchStressReloadPN = TpinchStressReloadPN;
	CpinchStrainReloadPN = TpinchStrainReloadPN;
	CpinchStressUnloadNP = TpinchStressUnloadNP;
	CpinchStrainUnloadNP = TpinchStrainUnloadNP;
	CpinchStressReloadNP = TpinchStressReloadNP;
	CpinchStrainReloadNP = TpinchStrainReloadNP;
	CbUnloadPN = TbUnloadPN;
	CbUnloadNP = TbUnloadNP;
	CstateFlag = TstateFlag;
	CdmgElasticE = TdmgElasticE;
	CdmgReloadE = TdmgReloadE;
	CbKdegDmg = TbKdegDmg;
	CpinchSlopePN = TpinchSlopePN;
	CpinchInterceptPN = TpinchInterceptPN;
	CreloadInterceptPN = TreloadInterceptPN;
	CpinchSlopeNP = TpinchSlopeNP;
	CpinchInterceptNP = TpinchInterceptNP;
	CreloadInterceptNP = TreloadInterceptNP;
	Cenergy = Tenergy;
	CstrainFresKdegDmg = TstrainFresKdegDmg;
	CstrainShearFailureDmg = TstrainShearFailureDmg;
	CbReloadAfterUnloadPN = TbReloadAfterUnloadPN;
	CbReloadAfterUnloadNP = TbReloadAfterUnloadNP;
	CstrainFlex = TstrainFlex;
	CstrainGlobal = TstrainGlobal;
	CenergyE = TenergyE;

    return 0;
}	

int 
PinchingLimitStateMaterial::revertToLastCommit(void)
{
	Tstress = Cstress;
	Tstrain = Cstrain;
	Ttangent = Ctangent;	
	TstrainRate = CstrainRate;
	TstrainMax = CstrainMax;
	TstrainMin = CstrainMin;
	Tdu = Cdu;
	TpinchStressUnloadPN = CpinchStressUnloadPN;
	TpinchStrainUnloadPN = CpinchStrainUnloadPN;
	TpinchStressReloadPN = CpinchStressReloadPN;
	TpinchStrainReloadPN = CpinchStrainReloadPN;
	TpinchStressUnloadNP = CpinchStressUnloadNP;
	TpinchStrainUnloadNP = CpinchStrainUnloadNP;
	TpinchStressReloadNP = CpinchStressReloadNP;
	TpinchStrainReloadNP = CpinchStrainReloadNP;
	TbUnloadPN = CbUnloadPN;
	TbUnloadNP = CbUnloadNP;
	TstateFlag = CstateFlag;
	TdmgElasticE = CdmgElasticE;
	TdmgReloadE = CdmgReloadE;
	TbKdegDmg = CbKdegDmg;
	TpinchSlopePN = CpinchSlopePN;
	TpinchInterceptPN = CpinchInterceptPN;
	TreloadInterceptPN = CreloadInterceptPN;
	TpinchSlopeNP = CpinchSlopeNP;
	TpinchInterceptNP = CpinchInterceptNP;
	TreloadInterceptNP = CreloadInterceptNP;
	Tenergy = Cenergy;
	TstrainFresKdegDmg = CstrainFresKdegDmg;
	TstrainShearFailureDmg = CstrainShearFailureDmg;
	TbReloadAfterUnloadPN = CbReloadAfterUnloadPN;
	TbReloadAfterUnloadNP = CbReloadAfterUnloadNP;
	TstrainFlex = CstrainFlex;
	TstrainGlobal = CstrainGlobal;
	TenergyE = CenergyE;

    return 0;
}

int 
PinchingLimitStateMaterial::revertToStart(void)
{
	// Zero Commit Variables
    Cstress = 0.0;
	Cstrain = 0.0;
	Ctangent = 0.0;	
	CstrainRate = 0.0;
	CstrainMax = 0.0;
	CstrainMin = 0.0;
	Cdu = 0.0;
	CpinchStressUnloadPN = 0.0;
	CpinchStrainUnloadPN = 0.0;
	CpinchStressReloadPN = 0.0;
	CpinchStrainReloadPN = 0.0;
	CpinchStressUnloadNP = 0.0;
	CpinchStrainUnloadNP = 0.0;
	CpinchStressReloadNP = 0.0;
	CpinchStrainReloadNP = 0.0;
	CbUnloadPN = 0.0;
	CbUnloadNP = 0.0;
	CstateFlag = 0;
	CdmgElasticE = 0.0;
	CdmgReloadE = 0.0;
	CbKdegDmg = 0.0;
	CpinchSlopePN = 0.0;
	CpinchInterceptPN = 0.0;
	CreloadInterceptPN = 0.0;
	CpinchSlopeNP = 0.0;
	CpinchInterceptNP = 0.0;
	CreloadInterceptNP = 0.0;
	Cenergy = 0.0;
	CstrainFresKdegDmg = 0.0;
	CstrainShearFailureDmg = 0.0;
	CbReloadAfterUnloadPN = 0.0;
	CbReloadAfterUnloadNP = 0.0;
	CstrainFlex = 0.0;
	CstrainGlobal = 0.0;
	CenergyE = 0.0;

	// Zero Static Variables for Limit Curve
	Kdeg = 0.0;
	bKdeg = 0.0;
	Fres = 0.0;
	strainShearFailure = 0.0;
	stressShearFailure = 0.0;
	InelastMonoEnergy = 0.0;
	strainFlexRev = 0.0;
	strainUn = 0.0;
	stressUn = 0.0;
	strainReload = 0.0;
	stressReload = 0.0;
	Ereload = 0.0;
	slopeFlexPred = 0.0;
	interceptFlexPred = 0.0;
	resFlag = 0;
	strainFlexRevFlag = 0;
	stressUnDmg = 0.0;
	strainUnDmg = 0.0;
	slopeGlobalEnv = 0.0;
	interceptGlobalEnv = 0.0;
	strainGlobalFresKdeg = 0.0;
	strainFlexRevDmg = 0.0;
	countGlobalEnv = 0;

	//Initialize initial elastic slope used in setTrial strain
	CdmgElasticE = E1;	
	CdmgReloadE = E1;

    return 0;
}


UniaxialMaterial* PinchingLimitStateMaterial::getCopy(void)
{
	PinchingLimitStateMaterial *theCopy = new PinchingLimitStateMaterial (this->getTag(),
	nodeTop, nodeBot, driftAxis, E1, curveType, curveTag,
	YpinchUnloadPN, YpinchReloadPN, XpinchReloadPN,
	YpinchUnloadNP, YpinchReloadNP, XpinchReloadNP,
	dmgStressLimE, dmgDeflMax,
	dmgElastic1,  dmgElastic2,  dmgElastic3,  dmgElastic4,  dmgElasticLim,
	dmgUnload1,   dmgUnload2,   dmgUnload3,   dmgUnload4,   dmgUnloadLim, 
	dmgReload1,   dmgReload2,   dmgReload3,   dmgReload4,   dmgReloadLim, dmgReloadCyclic,
	dmgStrength1, dmgStrength2, dmgStrength3, dmgStrength4, dmgStrengthLim, dmgStrengthCyclic,
	eleTag, b, d, h, a, st, As, Acc, 
	ld, db, rhot, fc, fy, fyt,
	theDomain, theNodeT, theNodeB, *theCurve, theElement);

	theCopy->Cstress = Cstress;
	theCopy->Cstrain = Cstrain;
	theCopy->Ctangent = Ctangent;
	theCopy->CstrainRate = CstrainRate;
	theCopy->CstrainMax = CstrainMax;
	theCopy->CstrainMin = CstrainMin;
	theCopy->Cdu = Cdu;
	theCopy->CpinchStressUnloadPN = CpinchStressUnloadPN;
	theCopy->CpinchStrainUnloadPN = CpinchStrainUnloadPN;
	theCopy->CpinchStressReloadPN = CpinchStressReloadPN;
	theCopy->CpinchStrainReloadPN = CpinchStrainReloadPN;
	theCopy->CpinchStressUnloadNP = CpinchStressUnloadNP;
	theCopy->CpinchStrainUnloadNP = CpinchStrainUnloadNP;
	theCopy->CpinchStressReloadNP = CpinchStressReloadNP;
	theCopy->CpinchStrainReloadNP = CpinchStrainReloadNP;
	theCopy->CbUnloadPN = CbUnloadPN;
	theCopy->CbUnloadNP = CbUnloadNP;
	theCopy->CstateFlag = CstateFlag;
	theCopy->CdmgElasticE = CdmgElasticE;
	theCopy->CdmgReloadE = CdmgReloadE;
	theCopy->CbKdegDmg = CbKdegDmg;
	theCopy->CpinchSlopePN = CpinchSlopePN;
	theCopy->CpinchSlopeNP = CpinchSlopeNP;
	theCopy->CpinchInterceptPN = CpinchInterceptPN;
	theCopy->CpinchInterceptNP = CpinchInterceptNP;
	theCopy->CreloadInterceptPN = CreloadInterceptPN;
	theCopy->CreloadInterceptNP = CreloadInterceptNP;
	theCopy->Cenergy = Cenergy;
	theCopy->CstrainFresKdegDmg = CstrainFresKdegDmg;
	theCopy->CstrainShearFailureDmg = CstrainShearFailureDmg;
	theCopy->CbReloadAfterUnloadPN = CbReloadAfterUnloadPN;
	theCopy->CbReloadAfterUnloadNP = CbReloadAfterUnloadNP;
	theCopy->CstrainFlex = CstrainFlex;
	theCopy->CstrainGlobal = CstrainGlobal;
	theCopy->CenergyE = CenergyE;
	theCopy->Kdeg = Kdeg;
	theCopy->bKdeg = bKdeg;
	theCopy->Fres = Fres;
	theCopy->strainShearFailure = strainShearFailure;
	theCopy->stressShearFailure = stressShearFailure;
	theCopy->InelastMonoEnergy = InelastMonoEnergy;
	theCopy->strainFlexRev = strainFlexRev;
	theCopy->strainUn = strainUn;
	theCopy->stressUn = stressUn;
	theCopy->strainReload = strainReload;
	theCopy->stressReload = stressReload;
	theCopy->Ereload = Ereload;
	theCopy->slopeFlexPred = slopeFlexPred;
	theCopy->interceptFlexPred = interceptFlexPred;
	theCopy->resFlag = resFlag;
	theCopy->strainFlexRevFlag = strainFlexRevFlag;
	theCopy->stressUnDmg = stressUnDmg;
	theCopy->strainUnDmg = strainUnDmg;
	theCopy->slopeGlobalEnv = slopeGlobalEnv;
	theCopy->interceptGlobalEnv = interceptGlobalEnv;
	theCopy->strainGlobalFresKdeg = strainGlobalFresKdeg;
	theCopy->strainFlexRevDmg = strainFlexRevDmg;
	theCopy->countGlobalEnv = countGlobalEnv;

	return theCopy;
}


int 
PinchingLimitStateMaterial::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	static Vector data(6);
	data(0) = this->getTag();
	data(1) = E1;
	data(2) = dmgElastic1;
	data(3) = dmgElastic2;
	data(4) = dmgElastic3;
	data(5) = dmgElastic4;
	data(6) = dmgElasticLim;
	data(7) = dmgUnload1;
	data(8) = dmgUnload2;
	data(9) = dmgUnload3;
	data(10) = dmgUnload4;
	data(11) = dmgUnloadLim;
	data(12) = dmgReload1;
	data(13) = dmgReload2;
	data(14) = dmgReload3;
	data(15) = dmgReload4;
	data(16) = dmgReloadLim;
	data(17) = dmgStrength1;
	data(18) = dmgStrength2;
	data(19) = dmgStrength3;
	data(20) = dmgStrength4;
	data(21) = dmgStrengthLim;
	data(22) = dmgStressLimE;
	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) 
		opserr << "PinchingLimitStateMaterial::sendSelf() - failed to send data\n";

	return res;
}

int 
PinchingLimitStateMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(6);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);
	if (res < 0) 
		opserr << "PinchingLimitStateMaterial::recvSelf() - failed to recv data\n";
	else {
		this->setTag(data(0));
		E1 = data(1);
		dmgElastic1 = data(2);
		dmgElastic2 = data(3);
		dmgElastic3 = data(4);
		dmgElastic4 = data(5);
		dmgElasticLim = data(6);
		dmgUnload1 = data(7);
		dmgUnload2 = data(8);
		dmgUnload3 = data(9);
		dmgUnload4 = data(10);
	}

  return res;
}

void 
PinchingLimitStateMaterial::Print(OPS_Stream &s, int flag)
{
	s << "PinchingLimitStateMaterial tag: " << this->getTag() << endln;
	s << "nodeT: " << nodeTop << endln;
	s << "nodeB: " << nodeBot << endln;
	s << "driftAxis: " << driftAxis << endln;
	s << "Kelas: " << E1 << endln;
	s << "crvTyp: " << curveType << endln;
	s << "crvTag: " << curveTag << endln;
	s << "eleTag: " << eleTag << endln;
	s << "YpinchUPN: " << YpinchUnloadPN << endln;
	s << "YpinchRPN: " << YpinchReloadPN << endln;
	s << "XpinchRPN: " << XpinchReloadPN << endln;
	s << "YpinchUNP: " << YpinchUnloadNP << endln;
	s << "YpinchRNP: " << YpinchReloadNP << endln;
	s << "XpinchRNP: " << XpinchReloadNP << endln;
	s << "dmgStrsLimE: " << dmgStressLimE << endln;
	s << "dmgDispMax: " << dmgDeflMax << endln;
	s << "dmgE1: " << dmgElastic1 << endln;
	s << "dmgE2: " << dmgElastic2 << endln;
	s << "dmgE3: " << dmgElastic3 << endln;
	s << "dmgE4: " << dmgElastic4 << endln;
	s << "dmgELim: " << dmgElasticLim << endln;
	s << "dmgU1: " << dmgUnload1 << endln;
	s << "dmgU2: " << dmgUnload2 << endln;
	s << "dmgU3: " << dmgUnload3 << endln;
	s << "dmgU4: " << dmgUnload4 << endln;
	s << "dmgULim: " << dmgUnloadLim << endln;
	s << "dmgR1: " << dmgReload1 << endln;
	s << "dmgR2: " << dmgReload2 << endln;
	s << "dmgR3: " << dmgReload3 << endln;
	s << "dmgR4: " << dmgReload4 << endln;
	s << "dmgRLim: " << dmgReloadLim << endln;
	s << "dmgRCyc: " << dmgReloadCyclic << endln;
	s << "dmgS1: " << dmgStrength1 << endln;
	s << "dmgS2: " << dmgStrength2 << endln;
	s << "dmgS3: " << dmgStrength3 << endln;
	s << "dmgS4: " << dmgStrength4 << endln;
	s << "dmgSLim: " << dmgStrengthLim << endln;
	s << "dmgSCyc: " << dmgStrengthCyclic << endln;
	s << "b: " << b << endln;
	s << "d: " << d << endln;
	s << "h: " << h << endln;
	s << "a: " << a << endln;
	s << "st: " << st << endln;
	s << "As: " << As << endln;
	s << "Acc: " << Acc << endln;
	s << "ld: " << ld << endln;
	s << "db: " << db << endln;
	s << "rhot: " << rhot << endln;
	s << "f'c: " << fc << endln;
	s << "fy: " << fy << endln;
	s << "fyt: " << fyt << endln;
}

void 
PinchingLimitStateMaterial::updateDamageE(void)	
{
	double strainLimE = dmgStressLimE/TdmgElasticE;
	double DispRatio = fmax(fabs(TstrainMax/strainLimE), fabs(TstrainMin/strainLimE));
	double EnergyMonotonic = fabs(0.5*strainLimE*dmgStressLimE);
	if((Tdu > 0.0 && Tstrain > 0.0) || (Tdu < 0.0 && Tstrain < 0.0))
		TenergyE = CenergyE + fabs(Tdu)*(fabs(Tstress) + fabs(Cstress))/2;
	if(EnergyMonotonic < 0.0001)
		EnergyMonotonic = 0.0001; 
	double EnergyRatio = TenergyE/EnergyMonotonic;
	double Delastic = dmgElastic1*pow(DispRatio,dmgElastic3)+dmgElastic2*pow(EnergyRatio,dmgElastic4);
	if(Delastic > 1 || Delastic > dmgElasticLim)
		Delastic = fmin(1.0, dmgElasticLim);
	TdmgElasticE = fmin(TdmgElasticE,E1*(1-Delastic));
}

void
PinchingLimitStateMaterial::defineBackbone(void)
{
	Kdeg = theCurve->getDegSlope();
	Fres = theCurve->getResForce();
	if (fabs(Tstress) < Fres)
	{
		opserr << "WARNING PinchingLimitStateMaterial::defineBackbone() - Fres must be less than shear load at failure\n";
		Fres = 0.2*Tstress;
		opserr << "Setting Fres to 0.2*Vmax = " << Fres << endln;
	}

	bKdeg = fabs(Tstress - Kdeg*Tstrain);
	TbKdegDmg = bKdeg;
	TstrainFresKdegDmg = (Fres - bKdeg)/Kdeg;
	strainShearFailure = fabs(Tstrain);
	TstrainShearFailureDmg = fabs(Tstrain);
	stressShearFailure = fabs(Tstress);

	InelastMonoEnergy = 0.5*(TstrainFresKdegDmg - Tstrain)*(Tstress - Fres) + (dmgDeflMax - Tstrain)*Fres;
}

void
PinchingLimitStateMaterial::updateEnergy(void)
{
	if(TstateFlag == 2 || TstateFlag == 3 || TstateFlag == -2 || TstateFlag == -3 || TstateFlag == -10 || TstateFlag == -10)
		Tenergy = Cenergy + Tdu*(Tstress + Cstress)/2;
}

void 
PinchingLimitStateMaterial::updateDamageS(void)
{
	double DispRatio = fmax((fabs(TstrainMax) - strainShearFailure)/dmgDeflMax,(fabs(TstrainMin) - strainShearFailure)/dmgDeflMax);
	double EnergyRatio = Tenergy/InelastMonoEnergy;
	double Dstrength = dmgStrength1*pow(DispRatio,dmgStrength3)+dmgStrength2*pow(EnergyRatio,dmgStrength4);
	if(Dstrength > 1 || Dstrength > dmgStrengthLim)
		Dstrength = fmin(1.0, dmgStrengthLim);	
	interceptGlobalEnv = interceptGlobalEnv*(1-Dstrength);
}

void 
PinchingLimitStateMaterial::updateDamageR(void)
{
	double DispRatio = fmax((fabs(TstrainMax) - strainShearFailure)/dmgDeflMax,(fabs(TstrainMin) - strainShearFailure)/dmgDeflMax);
	double EnergyRatio = Tenergy/InelastMonoEnergy;
	double Dreload = dmgReload1*pow(DispRatio,dmgReload3)+dmgReload2*pow(EnergyRatio,dmgReload4);
	if(Dreload > 1 || Dreload > dmgReloadLim)
		Dreload = fmin(1.0, dmgReloadLim);
	if(TstateFlag == 6 || TstateFlag == -6) 
		TdmgReloadE = fmin(TdmgReloadE,fmin(TdmgReloadE*(1-Dreload),fabs((stressUn-Cstress)/(strainUn-Cstrain))));
	if(TstateFlag == 9 || TstateFlag == -9)
		TdmgReloadE = fmin(TdmgReloadE,fmin(TdmgReloadE*(1-Dreload),fabs((-stressUnDmg-Cstress)/(-strainUnDmg-Cstrain))));
	else
		TdmgReloadE = fmin(TdmgReloadE*(1-Dreload),TdmgReloadE);
}

int
PinchingLimitStateMaterial::getStateFlag(void)
{
	if((CstateFlag == 1 || CstateFlag == 2) && Tdu > 0.0 && Tstrain < TstrainFresKdegDmg)					
		return 2;																					
	else if((CstateFlag == 2 || CstateFlag == 3) && Tdu > 0.0 && (Tstrain >= TstrainFresKdegDmg || (resFlag == 1 && Tstress >= Fres)))	
		return 3;																				
	else if((CstateFlag == 1 || CstateFlag == 2 || CstateFlag == 3) && Tdu < 0.0)	
		return 4;
	else if((CstateFlag == 4 || CstateFlag == 5 || CstateFlag == 6) && Tdu < 0.0 && Tstrain >= TpinchStrainUnloadPN)
		return 5;
	else if((CstateFlag == 4 || CstateFlag == 5 || CstateFlag == 6) && Tdu > 0.0)	
		return 6;
	else if((CstateFlag == 5 || CstateFlag == 7) && Tdu < 0.0 && Tstrain < TpinchStrainUnloadPN && Tstrain >= TpinchStrainReloadPN)
		return 7;
	else if(CstateFlag == 7 && Tdu > 0.0 && Tstrain < TpinchStrainUnloadPN && Tstrain >= TpinchStrainReloadPN)	
		return 8;
	else if(CstateFlag == 7 && Tdu < 0.0 && Tstrain < TpinchStrainReloadPN)	
		return 9;
	else if(CstateFlag == 10 && Tdu < 0.0 && Tstrain < TpinchStrainReloadPN)
		return 10;
	else if(CstateFlag == 10 && Tdu > 0.0 && Cstrain < TpinchStrainReloadPN)	
		return 11;

	else if((CstateFlag == -1 || CstateFlag == -2) && Tdu < 0.0 && Tstrain > -TstrainFresKdegDmg)				
		return -2;																					
	else if((CstateFlag == -2 || CstateFlag == -3) && Tdu < 0.0 && (Tstrain <= -TstrainFresKdegDmg || (resFlag == 1 && Tstress <= Fres)))	
		return -3;																				
	else if((CstateFlag == -1 || CstateFlag == -2 || CstateFlag == -3) && Tdu > 0.0)	
		return -4;
	else if((CstateFlag == -4 || CstateFlag == -5 || CstateFlag == -6) && Tdu > 0.0 && Tstrain <= TpinchStrainUnloadNP)	
		return -5;
	else if((CstateFlag == -4 || CstateFlag == -5 || CstateFlag == -6) && Tdu < 0.0)	
		return -6;
	else if((CstateFlag == -5 || CstateFlag == -7) && Tdu > 0.0 && Tstrain > TpinchStrainUnloadNP && Tstrain <= TpinchStrainReloadNP)	
		return -7;
	else if(CstateFlag == -7 && Tdu < 0.0 && Tstrain > TpinchStrainUnloadNP && Tstrain <= TpinchStrainReloadNP)	
		return -8;
	else if(CstateFlag == -7 && Tdu > 0.0 && Tstrain > TpinchStrainReloadNP)	
		return -9;
	else if(CstateFlag == -10 && Tdu > 0.0 && Tstrain > TpinchStrainReloadNP)	
		return -10;
	else if(CstateFlag == -10 && Tdu < 0.0 && Cstrain > TpinchStrainReloadNP)	
		return -11;
	else
		return 999;
}

void
PinchingLimitStateMaterial::definePinchingPN(void)
{
	TpinchStressUnloadPN = Cstress*YpinchUnloadPN;
	TpinchStrainUnloadPN = (TpinchStressUnloadPN-Cstress)/TdmgElasticE + Cstrain;
	TbUnloadPN = Cstress - TdmgElasticE*Cstrain;
	TpinchStressReloadPN = -Cstress*YpinchReloadPN;
	TpinchStrainReloadPN = -Cstrain*XpinchReloadPN;
	if(TpinchStressUnloadPN < TpinchStressReloadPN)
		TpinchStressReloadPN = TpinchStressUnloadPN;
	if(TpinchStrainUnloadPN < TpinchStrainReloadPN)
		TpinchStrainReloadPN = TpinchStrainUnloadPN-fabs(Cdu);
	TpinchSlopePN = (TpinchStressReloadPN - TpinchStressUnloadPN)/(TpinchStrainReloadPN - TpinchStrainUnloadPN);
	TpinchInterceptPN = TpinchStressUnloadPN - TpinchSlopePN * TpinchStrainUnloadPN;
}

void
PinchingLimitStateMaterial::definePinchingNP(void)
{
	TpinchStressUnloadNP = Cstress*YpinchUnloadNP;
	TpinchStrainUnloadNP = (TpinchStressUnloadNP-Cstress)/TdmgElasticE + Cstrain;
	TbUnloadNP = Cstress - TdmgElasticE*Cstrain;
	TpinchStressReloadNP = -Cstress*YpinchReloadNP;
	TpinchStrainReloadNP = -Cstrain*XpinchReloadNP;
	if(TpinchStressUnloadNP > TpinchStressReloadNP)
		TpinchStressReloadNP = TpinchStressUnloadNP;
	if(TpinchStrainUnloadNP > TpinchStrainReloadNP)
		TpinchStrainReloadNP = TpinchStrainUnloadNP+fabs(Cdu);
	TpinchSlopeNP = (TpinchStressReloadNP - TpinchStressUnloadNP)/(TpinchStrainReloadNP - TpinchStrainUnloadNP);
	TpinchInterceptNP = TpinchStressUnloadNP - TpinchSlopeNP * TpinchStrainUnloadNP;
}

double 
PinchingLimitStateMaterial::getFlexDisp(void)
{
	const Vector &dispT = theNodeT->getTrialDisp();
	const Vector &dispB = theNodeB->getTrialDisp();
	double delta = dispT(driftAxis-1) - dispB(driftAxis-1);
	return delta;
}

double
PinchingLimitStateMaterial::getFlexShift(void)
{
	double strainPred = (Tstress-interceptFlexPred)/slopeFlexPred;
	double shift = Tstress/fabs(Tstress)*(strainPred - TstrainFlex);
	return shift;
}

void
PinchingLimitStateMaterial::checkEnvelope(void)
{
	double shiftFlexRe = getFlexShift();
	double KdegLim = slopeGlobalEnv*fabs(TstrainGlobal)+interceptGlobalEnv;
	if(Tstress >= 0.0 && Tstrain >= 0.0)
	{
		if(Tstress >= KdegLim && TstrainGlobal < strainGlobalFresKdeg && resFlag == 0)
		{
			TstateFlag = 2;
			Ttangent = Kdeg;
			Tstress = KdegLim;
		}
		else if(Tstress >= Fres && TstrainGlobal >= strainGlobalFresKdeg)
		{
			TstateFlag = 3;
			Ttangent = 0.0001; 
			Tstress = Fres;
		}
	}
	else if(Tstress < 0.0 && Tstrain < 0.0)
	{
		if(Tstress <= -KdegLim && TstrainGlobal > -strainGlobalFresKdeg && resFlag == 0)
		{
			TstateFlag = -2;
			Ttangent = Kdeg;
			Tstress = -KdegLim;
		}
		else if(Tstress <= -Fres && TstrainGlobal <= -strainGlobalFresKdeg)
		{
			TstateFlag = -3;
			Ttangent = 0.0001; 
			Tstress = -Fres;
		}
	}
}

void
PinchingLimitStateMaterial::defineTargetVars(void)
{
	double Ag = b*h;
	double V = fabs(Tstress);
	double v = V/(b*d);
	double P = getAxialForce();
	double dmgSCyc = 0.037133 + 0.251204*(fy*As/(fc*Ag)) - 0.354989*(Acc/Ag) + 0.056569*(a/d);
	dmgStrengthCyclic = fmax(dmgSCyc,0);
	double YpinchU = -0.169113 + 0.088820*(v*1000/sqrt(fc*1000)) - 44.375649*(rhot) + 0.189494*(st/d);
	YpinchUnloadPN = fmax(YpinchU,0.0);	
	YpinchUnloadNP = fmax(YpinchU,0.0);
	double XpinchR = -0.589984 + 0.685461*(P/(fc*Ag)) + 0.008966*(ld/db) - 0.209699*(fy*As/(fc*Ag));
	XpinchReloadPN = XpinchR;
	XpinchReloadNP = XpinchR;
	double YpinchR = 0.262867 + 0.761220*(P/(fc*Ag)) - 1.066009*(fy*As/(fc*Ag)) + 0.005967*(ld/db);
	YpinchReloadPN = fmax(YpinchR,0.0);
	YpinchReloadNP = fmax(YpinchR,0.0);
}

double
PinchingLimitStateMaterial::getAxialForce(void)
{
	int trash;	
	const char *forceType2[1] = {"localForce"}; 
	Response *theForces = 0;		
	DummyStream dummy;	
	theForces =  theElement->setResponse(forceType2, 1, dummy);
	trash = theForces->getResponse();
	Information &theInfo = theForces->getInformation();
	Vector *forceVec = (theInfo.theVector);	
	if (forceVec == 0) {
		opserr << "FATAL ERROR RotationShearCurve -- unable to assign force vector\n" << endln;
		exit(-1);
	}
	double P = fabs((*forceVec)(0));
	return P;
}

void
PinchingLimitStateMaterial::defineE1(void)
{
	double Ag = b*h; //in^2
	double Ec = 57*sqrt(fc*1000); //ksi
	double Gv = Ec/(2*(1+0.20)); //ksi
	double L = 2*a;
	if (E1 == -4) 
		E1 = 2*Ag*Gv*5/(6*L);
	else if (E1 == -3) 
		E1 = Ag*Gv*5/(6*L);
	else if (E1 == -2) 
		E1 = 2*Ag*Gv*5/(6*a);
	else if (E1 == -1) 
		E1 = Ag*Gv*5/(6*a);
}
