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

// $Revision: 1.6 $
// $Date: 2006-08-28 17:40:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ReinforcingSteel.cpp,v $

/* ****************************************************************** **
** THIS FILE WAS DEVELOPED AT UC DAVIS                                **
**                                                                    **
** Programmed by: Jon Mohle (jfmohle@ucdavis.edu)                     **
** Supervisor: Sashi Kunnath (skkunnath@ucdavis.edu)                  **
**                                                                    **
********************************************************************* */
// Written: Jon Mohle
// Created: October 2003
// Updated: January 2006
//
// Description: This file contains the class definition for 
// ReinforcingSteel

#include "ReinforcingSteel.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#ifdef HelpDebugMat
  int ReinforcingSteel::classCount = 0;
#endif

#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_ReinforcingSteel)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 7) {
	opserr<<"WARNING insufficient arguments\n";
	opserr<<"uniaxialMaterial ReinforcingSteel ";
	opserr<<"tag? fy? fu? Es? Esh? esh? eult? ";
	opserr<<"<-GABuck?> <-DMBuck?> <-CMFatigue?> <-MPCurveParams?> <-IsoHard?>\n";
	return 0;
    }

    int tag;
    numdata = 1;
    if(OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[6];
    numdata = 6;
    if(OPS_GetDoubleInput(&numdata,data) < 0) {
	opserr << "WARNING invalid double data\n";
	return 0;
    }

    int buckModel = 0;
    double gabuckdata[4] = {0.0, 1.0, 1.0, 0.5};
    double dmbuckdata[2] = {0.0, 1.0};
    double fatiguedata[3] = {0.0, -4.46, 0.0};
    double mpdata[3] = {1.0/3.0,18.0,4.0};
    double isodata[2] = {0.0, 0.01};
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();

	if(strcmp(type,"-GABuck") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 4) {
		opserr << "WARNING insufficient optional arguments for -GABuck\n";
		opserr << "Want: <-GABuck lsr? beta? r? gama?>\n";
		return 0;
	    }
	    buckModel = 1;
	    numdata = 4;
	    if(OPS_GetDoubleInput(&numdata,gabuckdata) < 0) {
		opserr << "WARNING invalid double data\n";
		return 0;
	    }
	    
	} else if(strcmp(type, "-DMBuck") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 2) {
		opserr << "WARNING insufficient optional arguments for -DMBuck\n";
		opserr << "Want: <-DMBuck lsr? alpha?>\n";
		return 0;
	    }
	    buckModel = 2;
	    numdata = 2;
	    
	    if(OPS_GetDoubleInput(&numdata,dmbuckdata) < 0) {
		opserr << "WARNING invalid double data\n";
		return 0;
	    }
	    if(dmbuckdata[1]<0.75||dmbuckdata[1]>1.0) {
		opserr << "WARNING alpha usually is between 0.75 and 1.0\n";
		return 0;
	    }
	    
	} else if(strcmp(type,"-CMFatigue") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 3) {
		opserr << "WARNING insufficient optional arguments for -CMFatigue\n";
		opserr << "Want: <-CMFatigue Cf? alpha? Cd?>\n";
		return 0;
	    }
	    numdata = 3;
	    if(OPS_GetDoubleInput(&numdata,fatiguedata) < 0) {
		opserr << "WARNING invalid double data\n";
		return 0;
	    }

	} else if(strcmp(type, "-MPCurveParams") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if (numdata < 3) {
		opserr << "WARNING insufficient optional arguments for -MPCurveParams\n";
		opserr << "Want: <-CMFatigue R1? R2? R3?>\n";
		return 0;
	    }
	    numdata = 3;
	    if (OPS_GetDoubleInput(&numdata,mpdata)) {
		opserr << "WARNING invalid double data\n";
		return 0;
	    }
	} else if (strcmp(type,"-IsoHard") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if (numdata < 2) {
		opserr << "WARNING insufficient optional arguments for -IsoHard\n";
		opserr << "Want: <-IsoHard a1 limit>\n";
		return 0;
	    }
	    numdata = 2;
	    if (OPS_GetDoubleInput(&numdata,isodata)) {
		opserr << "WARNING invalid double data\n";
		return 0;
	    }
	} else {
	    opserr << "WARNING did not recognize optional flag\n";
	    opserr << "Possible Optional Flags: ";
	    opserr << "<-GABuck?> <-DMBuck?> <-CMFatigue?> <-MPCurveParams?> <-IsoHard?>\n";
	    return 0;
	}
    }

    // Parsing was successful, allocate the material
    double slen = 0.0, beta = 1.0;
    if(buckModel == 1) {
	slen = gabuckdata[0];
	beta = gabuckdata[1];
    } else if(buckModel == 2) {
	slen = dmbuckdata[0];
	beta = dmbuckdata[1];
    }
    UniaxialMaterial* theMaterial = 0;
    theMaterial = new ReinforcingSteel(tag,data[0],data[1],data[2],data[3],data[4],data[5],
				       buckModel,slen,beta,gabuckdata[2],gabuckdata[3],
				       fatiguedata[0],fatiguedata[1],fatiguedata[2],
				       mpdata[0],mpdata[1],mpdata[2],
				       isodata[0],isodata[1]);

    if (theMaterial == 0) {
	opserr << "WARNING could not create uniaxialMaterial of type ReinforcingSteel\n";
	return 0;
    }

    return theMaterial;
    
}

ReinforcingSteel::ReinforcingSteel(int tag, double fy_, double fsu_, double Es_, double Esh_, double esh_, double esu_, 
                                   int buckModel, double slenderness, double alpha, double r, double gama, 
                                   double Fatigue1, double Fatigue2, double Degrade, double rc1, double rc2, double rc3, 
                                   double A1, double HardLim)
  :UniaxialMaterial(tag,MAT_TAG_ReinforcingSteel), fy(fy_), fsu(fsu_), Es(Es_), Esh(Esh_), esh(esh_), esu(esu_), a1(A1), hardLim(HardLim),
  BuckleModel(buckModel),LDratio(slenderness),beta(alpha),fsu_fraction(gama),Fat1(Fatigue1),RC1(rc1),RC2(rc2),RC3(rc3)
{ 
  if((r>=0.0) & (r<=1.0)) 
    reduction=r;
  else
    if(r<=0)
      reduction = 0.0;
    else
      reduction = 1.0;

	if((Fatigue1==0) || (Fatigue2==0)) {
		Fat1=9.9e30;
		Fat2=1.0;
		Deg1=0.0;
  } else {
    Fat2=1.0/Fatigue2;
    if(Degrade != 0.0) 
      Deg1=pow(Fat1/Degrade,Fat2);
    else
      Deg1=0.0;
  }
	
  //initial yield point in natural stress-strain
	eyp=log(1.0+fy/Es);
	fyp=fy*(1.0+fy/Es);
	Esp=fyp/eyp;

  // ultimate strain and slope in natural stress-strain
	esup=log(1.0+esu);
  Esup=fsu*(1.0+esu);
  
  //updateHardeningLoaction(1.0);  done in revert to start
  
  ZeroTol=1.0E-14;
  this->revertToStart();
}

void
ReinforcingSteel::updateHardeningLoaction(double PlasticStrain)
{
  double ep;
  double pBranchStrain_t = Temax - Backbone_f(Temax)/Esp;
  double pBranchStrain_c = Temin + Backbone_f(Temin)/Esp;
  if (pBranchStrain_t > -pBranchStrain_c)
    ep = PlasticStrain - pBranchStrain_t;   
  else
    ep = PlasticStrain + pBranchStrain_c;
  THardFact = 1.0 - a1*ep;
  if (THardFact<hardLim) THardFact = hardLim;
  if (THardFact>1.0) THardFact = 1.0;
  updateHardeningLoactionParams();
}

void
ReinforcingSteel::updateHardeningLoactionParams()
{
  double ey = exp(eyp)-1.0;
  double fy = fyp/(1.0+ey);
  double eshLoc = THardFact*(esh-ey)+ey;
  
  // strain hardened point in natural stress-strain
  eshp=log(1.0+eshLoc);
  fshp=fy*(1.0+eshLoc);
  
  // ultimate stress in natural stress-strain
  fsup=Esup-(esup-eshp)*Esup; 
  
  // strain hardedned slope, yield plateau slope, and intersect
  Eshp=Esh*pow(1.0+eshLoc,2.0)+fshp - Esup;
  Eypp=(fshp-fyp)/(eshp-eyp);
  fint = fyp-Eypp*eyp;
  
  p = Eshp*(esup-eshp)/(fsup-fshp);
  // Set backbone transition variables
  double fTemp = Backbone_fNat(eshp+0.0002);
  Eshpb = Eshp*pow((fsup-fTemp)/(fsup-fshp),1.0-1.0/p);
  eshpa = eshp + 0.0002 - 2.0*(fTemp-fshp)/Eshpb;
}

ReinforcingSteel::ReinforcingSteel(int tag)
  :UniaxialMaterial(tag,MAT_TAG_ReinforcingSteel)
{
#ifdef HelpDebugMat
  thisClassNumber = ++classCount;
  thisClassCommit = 0;
#endif
  ZeroTol=1.0E-14;
}

ReinforcingSteel::~ReinforcingSteel(){

}

/***************** material state determination methods ***********/
int 
ReinforcingSteel::setTrialStrain(double strain, double strainRate) {

  // need to reset values to last committed
  this->revertToLastCommit();

  int res = 0;
  #ifdef HelpDebugMat
  thisClassStep++;
  if (thisClassCommit == 4000 && thisClassStep == 1)
    if (scalefactor()<1.0)
      opserr << scalefactor() << "\n";
#endif
  // Reset Trial History Variables to Last Converged State
  revertToLastCommit();
  
#ifdef _WIN32  
  if(_fpclass(strain)< 8 || _fpclass(strain)==512) {
    opserr << "bad trial strain\n";
    return -1;
  }
#endif
  
  if(strain< -0.95) {
    opserr << "Large trial compressive strain\n";
    return -1;
  } else 
    TStrain = log(1.0 + strain);
  
  if (TStrain == CStrain) return 0;
  
  if (TBranchNum==0){
    if (TStrain>0.0) TBranchNum = 1;
    if (TStrain<0.0) TBranchNum = 2;
  }
  
#ifdef _WIN32
  if(_fpclass(Tfch)< 8 || _fpclass(Tfch)==512 || _fpclass(Tfch)< 8 || _fpclass(Tfch)==512) {
    opserr << "bad stress or tangent\n";
    return -1;
  }
#endif
  if(thisClassNumber==51 && thisClassCommit==781) {
    thisClassCommit = thisClassCommit;
  }
  res = BranchDriver(res);
  
#ifdef _WIN32
  if(_fpclass(TStress)< 8 || _fpclass(TStress)==512 || _fpclass(TTangent)< 8 || _fpclass(TTangent)==512) {
    opserr << "bad stress or tangent\n";
    return -1;
  }
#endif
  
  if (res==0)
    return 0;
  else
    return -1;
}

double 
ReinforcingSteel::getStrain(void) {
  return exp(TStrain)-1.0;
}

double 
ReinforcingSteel::getStress(void) {
  if (theBarFailed) return 0.0;
  double tempstr=TStress;
  switch(BuckleModel) {
  case  1:  tempstr = Buckled_stress_Gomes(TStrain,TStress);
    break;
  case  2:  tempstr = Buckled_stress_Dhakal(TStrain,TStress);
    break;
  }
  double tempOut = tempstr*scalefactor()/exp(TStrain);
  
#ifdef _WIN32
  if(_fpclass(tempOut)< 8 || _fpclass(tempOut)==512)
    opserr << "bad Stress in ReinforcingSteel::getStress\n";
#endif
  
  return tempOut;
}

double 
ReinforcingSteel::getTangent(void) {
  double taTan = TTangent;
  switch(BuckleModel) {
  case  1:  taTan = Buckled_mod_Gomes(TStrain,TStress,TTangent);
    break;
  case  2:  taTan = Buckled_mod_Dhakal(TStrain,TStress,TTangent);
    break;
  }
  double scfact = scalefactor();
  double tempOut = (taTan-TStress)*scfact/pow(exp(TStrain),2.0);
  
#ifdef _WIN32
  if(_fpclass(tempOut)< 8 || _fpclass(tempOut)==512)
    opserr << "bad tangent in ReinforcingSteel::getTangentat\n";
#endif
  
  return tempOut;
}

double 
ReinforcingSteel::getInitialTangent(void) {
  return Esp;
}

/***************** path dependent behavior methods ***********/
int 
ReinforcingSteel::commitState(void) {
#ifdef HelpDebugMat
  thisClassCommit++;
  thisClassStep = 0;
#endif
  
  if(TBranchNum <= 1)
    TBranchMem=0;
  else
    TBranchMem = (TBranchNum+1)/2;
  
  for(int i=0; i<=LastRule_RS/2; i++)
    C_ePlastic[i]=T_ePlastic[i];
  
  CFatDamage       = TFatDamage;
  
  // commit trial history variables
  CBranchNum    = TBranchNum;
  Ceo_p         = Teo_p;
  Ceo_n         = Teo_n;
  Cemax         = Temax;
  Cemin         = Temin;
  CeAbsMax      = TeAbsMax;
  CeAbsMin      = TeAbsMin;
  CeCumPlastic  = TeCumPlastic;
  CHardFact     = THardFact;

  if(TBranchNum > 2) {
    CR[TBranchMem]    = TR;
    Cfch[TBranchMem]  = Tfch;
    CQ[TBranchMem]    = TQ;
    CEsec[TBranchMem] = TEsec;
    Cea[TBranchMem]   = Tea;
    Cfa[TBranchMem]   = Tfa;
    CEa[TBranchMem]   = TEa;
    Ceb[TBranchMem]   = Teb;
    Cfb[TBranchMem]   = Tfb;
    CEb[TBranchMem]   = TEb;
  }
  //by SAJalali
  Energy += 0.5*(TStress + CStress)*(TStrain - CStrain);

  // commit trial state variables
  CStrain    = TStrain;  
  CStress    = TStress;
  CTangent   = TTangent;
  return 0;
}

int 
ReinforcingSteel::revertToLastCommit(void) {
  for(int i=0; i<=LastRule_RS/2; i++)
    T_ePlastic[i]=C_ePlastic[i];
  
  TFatDamage = CFatDamage;
  
  // Reset trial history variables to last committed state
  TBranchNum    = CBranchNum;
  Teo_p         = Ceo_p;
  Teo_n         = Ceo_n;
  Temax         = Cemax;
  Temin         = Cemin;
  TeAbsMax      = CeAbsMax;
  TeAbsMin      = CeAbsMin;
  TeCumPlastic  = CeCumPlastic;
  THardFact     = CHardFact;
  updateHardeningLoactionParams();
  
  if(TBranchNum > 2) SetPastCurve(TBranchNum);
  
  // Reset trial state variables to last committed state
  TStress    = CStress;
  TTangent   = CTangent;

  return 0;
}

int 
ReinforcingSteel::revertToStart(void)
{
	Energy = 0;	// by SAJalali
	theBarFailed = 0;
  
  THardFact = 1.0;
  CHardFact = 1.0;
  updateHardeningLoactionParams();

  CFatDamage = TFatDamage;
  for(int i=0; i<=LastRule_RS/2; i++) {
    C_ePlastic[i]  = 0.0;
    T_ePlastic[i]  = 0.0;
    CR[i]         = 0.0;
    Cfch[i]       = 0.0;
    CQ[i]         = 0.0;
    CEsec[i]      = 0.0;
    Cea[i]        = 0.0;
    Cfa[i]        = 0.0;
    CEa[i]        = 0.0;
    Ceb[i]        = 0.0;
    Cfb[i]        = 0.0;
    CEb[i]        = 0.0;
  }
  TR            = 0.0;
  Tfch          = 0.0;
  TQ            = 0.0;
  TEsec         = 0.0;
  Tea           = 0.0;
  Tfa           = 0.0;
  TEa           = 0.0;
  Teb           = 0.0;
  Tfb           = 0.0;
  TEb           = 0.0;

  // reset trial history variables
  CBranchNum    = 0;
  TBranchNum    = 0;
  Ceo_p         = 0.0;
  Teo_p         = 0.0;
  Ceo_n         = 0.0;
  Teo_n         = 0.0;
  Cemax         = 0.0;
  Temax         = 0.0;
  Cemin         = 0.0;
  Temin         = 0.0;
  CeAbsMax      = 0.0;
  TeAbsMax      = 0.0;
  CeAbsMin      = 0.0;
  TeAbsMin      = 0.0;
  TeCumPlastic  = 0.0;
  CeCumPlastic  = 0.0;

  // reset trial state variables
  CStrain       = 0.0; 
  TStrain       = 0.0;  
  CStrain       = 0.0;
  CStress       = 0.0; 
  TStress       = 0.0;
  CStress       = 0.0;  
  CTangent      = Esp; 
  TTangent      = Esp;

  CFatDamage    = 0.0;
  TFatDamage    = 0.0;

  return 0;
}

UniaxialMaterial *
ReinforcingSteel::getCopy(void)
{
  ReinforcingSteel *theCopy =
  new ReinforcingSteel(this->getTag());
  
  theCopy->reduction    = reduction;
  theCopy->fsu_fraction = fsu_fraction;
  theCopy->beta       = beta;
  theCopy->theBarFailed = theBarFailed;

  // natural stress-strain variables
	theCopy->p            = p;
  theCopy->Esp          = Esp;   // Natural Elastic Modulus
  theCopy->eshp         = eshp;	 // Natural Hardening Strain
  theCopy->fshp         = fshp;  // Natural Hardening Stress
  theCopy->Eshp         = Eshp;  // Natural Hardening Modulus
  theCopy->esup         = esup;  // Natural Strain at Peak Stress
  theCopy->fsup         = fsup;  // Natural Peak Stress
  theCopy->Esup         = Esup;  // Natural Peak Stress Moduus
  theCopy->Eypp         = Eypp;  // Natural Yield Plateau Modulus
  theCopy->fint         = fint;  // Natural yield plateau intersect
  theCopy->eyp          = eyp;   // Natural yield strain
  theCopy->fyp          = fyp;   // Natural yield stress

  theCopy->esh          = esh;   // Engineering hardening strain (user input)
  theCopy->Esh          = Esh;   // Engineering hardening stress (user input)

  theCopy->eshpa        = eshpa; // Curve smoothing Parameters (at SH transition)
  theCopy->Eshpb        = Eshpb; // These are used to eliminate a sudden discontinuity in stiffness

  theCopy->a1           = a1;
  theCopy->hardLim      = hardLim;
  theCopy->THardFact    = THardFact;
  theCopy->CHardFact    = CHardFact;

  theCopy->TFatDamage   = TFatDamage;
  theCopy->CFatDamage   = CFatDamage;
  theCopy->LDratio      = LDratio;
	theCopy->Fat1         = Fat1;
  theCopy->Fat2         = Fat2;
	theCopy->Deg1         = Deg1;
  theCopy->BuckleModel  = BuckleModel;
  theCopy->BackStress   = BackStress;

  theCopy->TBranchMem   = TBranchMem;
  theCopy->TBranchNum   = TBranchNum;
  theCopy->Teo_p        = Teo_p;
  theCopy->Teo_n        = Teo_n;
  theCopy->Temax        = Temax;
  theCopy->Temin        = Temin;
  theCopy->TeAbsMax     = TeAbsMax;
  theCopy->TeAbsMin     = TeAbsMin;
  theCopy->TeCumPlastic = TeCumPlastic;

  theCopy->CBranchNum   = CBranchNum;
  theCopy->Ceo_p        = Ceo_p;
  theCopy->Ceo_n        = Ceo_n;
  theCopy->Cemax        = Cemax;
  theCopy->Cemin        = Cemin;
  theCopy->CeAbsMax     = CeAbsMax;
  theCopy->CeAbsMin     = CeAbsMin;
  theCopy->CeCumPlastic = CeCumPlastic;

  for(int i=0; i<=LastRule_RS/2; i++) {
	  theCopy->C_ePlastic[i] = C_ePlastic[i];
    theCopy->T_ePlastic[i] = T_ePlastic[i];
	  theCopy->CR[i]        = CR[i];
	  theCopy->Cfch[i]      = Cfch[i];
	  theCopy->CQ[i]        = CQ[i];
	  theCopy->CEsec[i]     = CEsec[i];
	  theCopy->Cea[i]       = Cea[i];
	  theCopy->Cfa[i]       = Cfa[i];
	  theCopy->CEa[i]       = CEa[i];
	  theCopy->Ceb[i]       = Ceb[i];
	  theCopy->Cfb[i]       = Cfb[i];
	  theCopy->CEb[i]       = CEb[i];
  }
  theCopy->TR           = TR;
  theCopy->Tfch         = Tfch;
  theCopy->TQ           = TQ;
  theCopy->TEsec        = TEsec;
  theCopy->Tea          = Tea;
  theCopy->Tfa          = Tfa;
  theCopy->TEa          = TEa;
  theCopy->Teb          = Teb;
  theCopy->Tfb          = Tfb;
  theCopy->TEb          = TEb;

  theCopy->re           = re;
  theCopy->rE1          = rE1;
  theCopy->rE2          = rE2;

  theCopy->CStrain      = CStrain;  
  theCopy->CStress      = CStress;
  theCopy->CTangent     = CTangent;
  theCopy->TStrain      = TStrain;  
  theCopy->TStress      = TStress;
  theCopy->TTangent     = TTangent;

  theCopy->RC1          = RC1;
  theCopy->RC2          = RC2;
  theCopy->RC3          = RC3;

  return theCopy;
}

int 
ReinforcingSteel::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  int index =0;
  static Vector data(75+12*(LastRule_RS/2+1));
  
  data(index++) = this->getTag();
  data(index++) = reduction;
  data(index++) = fsu_fraction;
  data(index++) = beta;
  data(index++) = theBarFailed;
  data(index++) = p;
  data(index++) = Esp;
  data(index++) = eshp;
  data(index++) = fshp;
  data(index++) = Eshp;
  data(index++) = esup;
  data(index++) = fsup;
  data(index++) = Esup;
  data(index++) = Eypp;
  data(index++) = fint;
  data(index++) = eyp;
  data(index++) = fyp;
  data(index++) = esh;
  data(index++) = CeCumPlastic;
  data(index++) = TeCumPlastic;
  data(index++) = a1;
  data(index++) = hardLim;
  data(index++) = THardFact;
  data(index++) = CHardFact;
  data(index++) = Esh;
  data(index++) = eshpa;
  data(index++) = Eshpb;
  data(index++) = TFatDamage;
  data(index++) = CFatDamage;
  data(index++) = LDratio;
  data(index++) = Fat1;
  data(index++) = Fat2;
  data(index++) = Deg1;
  data(index++) = BuckleModel;
  data(index++) = TBranchMem;
  data(index++) = TBranchNum;
  data(index++) = Teo_p;
  data(index++) = Teo_n;
  data(index++) = Temax;
  data(index++) = Temin;
  data(index++) = TeAbsMax;
  data(index++) = TeAbsMin;
  data(index++) = CBranchNum;
  data(index++) = Ceo_p;
  data(index++) = Ceo_n;
  data(index++) = Cemax;
  data(index++) = Cemin;
  data(index++) = CeAbsMax;
  data(index++) = CeAbsMin;
  data(index++) = TR;
  data(index++) = Tfch;
  data(index++) = TQ;
  data(index++) = TEsec;
  data(index++) = Tea;
  data(index++) = Tfa;
  data(index++) = TEa;
  data(index++) = Teb;
  data(index++) = Tfb;
  data(index++) = TEb;

  data(index++) = re;
  data(index++) = rE1;
  data(index++) = rE2;
  data(index++) = CStrain;  
  data(index++) = CStress;
  data(index++) = CTangent;
  data(index++) = TStrain;  
  data(index++) = TStress;
  data(index++) = TTangent;
  data(index++) = BackStress;

  data(index++) = RC1;
  data(index++) = RC2;
  data(index++) = RC3;

  for(int i=0; i<=LastRule_RS/2; i++) {
    data(index++) = C_ePlastic[i];
    data(index++) = T_ePlastic[i];
    data(index++) = CR[i];
    data(index++) = Cfch[i];
    data(index++) = CQ[i];
    data(index++) = CEsec[i];
    data(index++) = Cea[i];
    data(index++) = Cfa[i];
    data(index++) = CEa[i];
    data(index++) = Ceb[i];
    data(index++) = Cfb[i];
    data(index++) = CEb[i];
  }
#ifdef _NDEBUG
  if (--index != data.Size())
    opserr << "ReinforcingSteel::sendSelf() wrong vector size\n";
#endif
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ReinforcingSteel::sendSelf() - failed to send data\n";
  
  return res;
}

int 
ReinforcingSteel::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int index =0;
  static Vector data(75+12*(LastRule_RS/2+1));
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ReinforcingSteel::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(index++));
    reduction = data(index++);
    fsu_fraction = data(index++);
    beta = data(index++);
    theBarFailed = static_cast<int>(data(index++));
	  p     = data(index++);
    Esp   = data(index++);
    eshp  = data(index++);
    fshp  = data(index++);
    Eshp  = data(index++);
    esup  = data(index++);
    fsup  = data(index++);
    Esup  = data(index++);
    Eypp  = data(index++);
    fint  = data(index++);
    eyp   = data(index++);
    fyp   = data(index++);
    esh   = data(index++);
    CeCumPlastic = data(index++);
    TeCumPlastic = data(index++);
    a1 = data(index++);
    hardLim = data(index++);
    THardFact = data(index++);
    CHardFact = data(index++);
    Esh   = data(index++);
    eshpa = data(index++);
    Eshpb = data(index++);
    TFatDamage  = data(index++);
    CFatDamage  = data(index++);
    LDratio     = data(index++);
	  Fat1        = data(index++);
    Fat2        = data(index++);
	  Deg1        = data(index++);
    BuckleModel = static_cast<int>(data(index++));
    TBranchMem  = static_cast<int>(data(index++));
    TBranchNum  = static_cast<int>(data(index++));
    Teo_p       = data(index++);
    Teo_n       = data(index++);
    Temax       = data(index++);
    Temin       = data(index++);
    TeAbsMax    = data(index++);
    TeAbsMin    = data(index++);
    CBranchNum  = static_cast<int>(data(index++));
    Ceo_p       = data(index++);
    Ceo_n       = data(index++);
    Cemax       = data(index++);
    Cemin       = data(index++);
    CeAbsMax    = data(index++);
    CeAbsMin    = data(index++);
    TR          = data(index++);
    Tfch        = data(index++);
    TQ          = data(index++);
    TEsec       = data(index++);
    Tea         = data(index++);
    Tfa         = data(index++);
    TEa         = data(index++);
    Teb         = data(index++);
    Tfb         = data(index++);
    TEb         = data(index++);
    re          = data(index++);
    rE1         = data(index++);
    rE2         = data(index++);
    CStrain     = data(index++);  
    CStress     = data(index++);
    CTangent    = data(index++);
    TStrain     = data(index++);  
    TStress     = data(index++);
    TTangent    = data(index++);
    BackStress  = data(index++);

    RC1         = data(index++);
    RC2         = data(index++);
    RC3         = data(index++);
    for(int i=0; i<=LastRule_RS/2; i++) {
	    C_ePlastic[i] = data(index++);
      T_ePlastic[i] = data(index++);
	    CR[i]         = data(index++);
	    Cfch[i]       = data(index++);
	    CQ[i]         = data(index++);
	    CEsec[i]      = data(index++);
	    Cea[i]        = data(index++);
	    Cfa[i]        = data(index++);
	    CEa[i]        = data(index++);
	    Ceb[i]        = data(index++);
	    Cfb[i]        = data(index++);
	    CEb[i]        = data(index++);
    }
    #ifdef _NDEBUG
      if (--index != data.Size())
        opserr << "ReinforcingSteel::sendSelf() wrong vector size\n";
    #endif
  }
  return res;
}

void 
ReinforcingSteel::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "ReinforcingSteel, tag: " << this->getTag() << endln;
        s << "  N2p: " << CFatDamage << endln;
        //s << "  sigmaY: " << sigmaY << endln;
        //s << "  Hiso: " << Hiso << endln;
        //s << "  Hkin: " << Hkin << endln;
        //s << "  eta: " << eta << endln;
    }
    
    if (flag == 3) {
        s << CStrain << "  " << CStress << "  " << CTangent << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ReinforcingSteel\", ";
        s << "\"E\": " << Es << ", ";
        s << "\"Eh\": " << Esh << ", ";
        s << "\"fy\": " << fy << ", ";
        s << "\"fu\": " << fsu << ", ";
        s << "\"epsh\": " << esh << ", ";
        s << "\"epsu\": " << esu << "}";
    }
}

int
ReinforcingSteel::Sign(double x) {
  if(x<0.0)
	return -1;
  else
	return 1;
}

/*****************************************************************************************/
/*************************     Base Menegotto-Pinto Equation       ***********************/
/*****************************************************************************************/
double
ReinforcingSteel::MPfunc(double a)
{                
  if(a>=1.0)
    opserr << "a is one in ReinforcingSteel::MPfunc()\n";
  //double temp1 = pow(a,TR+1);
  //double temp2 = pow(a,TR);
  //double temp3 = TEa*a*(1-temp2)/(1-a);
  //double temp4 = TEsec*(1-temp1)/(1-a);
  //return TEb-temp4+temp3;
  return TEb-TEsec*(1-pow(a,TR+1))/(1-a)+TEa*a*(1-pow(a,TR))/(1-a);
}

int
ReinforcingSteel::SetMP()
{
    double Rmin;
    double a = 0.01;
    double ao, ao_last;
    double da;
    int numIter;
    int maxIter = 50;
    bool converged;

    if (TEb - TEsec == 0.0) {
        TQ = 1.0;
        Tfch = Tfb;
    }
    else {
        if (TEsec != TEa) {
            Rmin = (TEb - TEsec) / (TEsec - TEa);
            if (Rmin < 0.0) {
                opserr << "R is negative in ReinforcingSteel::SetMP()\n";
                Rmin = 0.0;
            }
            if (TR <= Rmin)
                TR = Rmin + 0.01;
            
            numIter = 0;
            converged = false;
            while (converged == false && numIter < maxIter) {
                numIter++;
                if (a > DBL_EPSILON) {
                    if (MPfunc(a)*MPfunc(1.0 - a) > 0.0)
                        a = a / 2.0;
                    else
                        converged = true;
                }
                else
                    converged = true;
            }
            if (numIter >= maxIter) {
                opserr << "WARNING: ReinforcingSteel::SetMP() - did not converge finding a\n";
                return -1;
            }

            ao = Rmin / TR;
            if (ao >= 1.0)
                ao = 0.999999;

            numIter = 0;
            converged = false;
            while (converged == false && numIter < maxIter) {
                numIter++;
                if (a > DBL_EPSILON) {
                    if (MPfunc(ao)*MPfunc(1.0 - a) < 0.0)
                        ao = sqrt(ao);
                    else
                        converged = true;
                }
                else
                    converged = true;

                if (ao > 0.999999)
                    converged = true;
            }
            if (numIter >= maxIter) {
                opserr << "WARNING: ReinforcingSteel::SetMP() - did not converge finding ao\n";
                return -2;
            }
            if (ao >= 1.0)
                ao = 0.999999;

            numIter = 0;
            converged = false;
            while (converged == false && numIter < maxIter) {
                numIter++;
                ao_last = ao;

                da = 0.49*(1 - ao);
                if (da > ao / 10.0)
                    da = ao / 10.0;
                if (ao + da >= 1.0)
                    da = (1.0 - ao) / 10.0;

                double tempdenom = MPfunc(ao + da) - MPfunc(ao - da);
                if (tempdenom != 0.0) {
                    ao = ao - 2 * MPfunc(ao)*da / tempdenom;
                    if (ao > 0.99999999999)
                        ao = 0.99999999999;
                    if (ao < 0.0) {
                        ao = 0.0;
                        converged = true;
                    }
                }
                if (fabs(ao_last - ao) < 1.0E-4)
                    converged = true;
            }
            if (numIter >= maxIter) {
                opserr << "WARNING: ReinforcingSteel::SetMP() - did not converge finding da and ao\n";
                da = da / 100.0;
                ao = ao_last;
                ao = ao - 2 * MPfunc(ao)*da / (MPfunc(ao + da) - MPfunc(ao - da));
                return -3;
            }
            if (ao > 0.99999999)
                ao = 0.99999999;
        }
        else
            ao = 0.99999999;

        TQ = (TEsec / TEa - ao) / (1 - ao);
        double temp1 = pow(ao, TR);
        double temp2 = pow(1.0 - temp1, 1.0 / TR);
        double b = temp2 / ao;
        Tfch = Tfa + TEa / b*(Teb - Tea);
    }

    if (fabs(Teb - Tea) < 1.0E-7)
        TQ = 1.0;

    return 0;
}

double
ReinforcingSteel::MP_f(double e) {
  return Tfa+TEa*(e-Tea)*(TQ-(TQ-1.0)/pow(pow(fabs(TEa*(e-Tea)/(Tfch-Tfa)),TR)+1.0,1/TR));
}

double
ReinforcingSteel::MP_E(double e) {
  if(TR>100.0 || e==Tea) {
    return TEa;
  } else {
    double Esec=(MP_f(e)-Tfa)/(e-Tea);
    return Esec-(Esec-TQ*TEa)/(pow(fabs(TEa*(e-Tea)/(Tfch-Tfa)),-TR)+1.0); 
  }
}

void 
ReinforcingSteel::SetTRp(void)
{
  TR=pow(fyp/Esp,RC1)*RC2*(1.0-RC3*(Teb-Tea));
}

void
ReinforcingSteel::SetTRn(void)
{
  TR=pow(fyp/Esp,RC1)*RC2*(1.0-RC3*(Tea-Teb));
}

void 
ReinforcingSteel::SetTRp1(void)
{
  TR=pow(fyp/Esp,RC1)*RC2*(1.0-RC3*(Teb-Tea));
}

void
ReinforcingSteel::SetTRn1(void)
{
  TR=pow(fyp/Esp,RC1)*RC2*(1.0-RC3*(Tea-Teb));
}

void
ReinforcingSteel::SetPastCurve(int branch)
{
	if(branch == 1)
		TBranchMem=0;
	else
		TBranchMem = (branch+1)/2;
  
	Tea   = Cea[TBranchMem];
  Teb   = Ceb[TBranchMem];
  Tfa   = Cfa[TBranchMem];
  Tfb   = Cfb[TBranchMem];
  TEa   = CEa[TBranchMem];
  TEb   = CEb[TBranchMem];
  TR    = CR [TBranchMem];
  Tfch  = Cfch[TBranchMem];
  TQ    = CQ [TBranchMem];
  TEsec = CEsec[TBranchMem];
}
/*****************************************************************************************/
/***********************        Base Stress-Strain Relations         *********************/
/*****************************************************************************************/
double
ReinforcingSteel::Backbone_f(double e) {
		if(e<0.0)
			return -Backbone_fNat(fabs(e));
		else
			return  Backbone_fNat(fabs(e));
}

double
ReinforcingSteel::Backbone_fNat(double essp) {
	if(essp>eshpa) {
		if(essp>esup)
			return fsup + (essp-eshp)*Esup;
		else {
      if (essp < eshp+0.0002)
				return (Eshpb-Eypp)*pow(essp-eshpa,2.0)/(2*(eshp+0.0002-eshpa))+ essp*Eypp + fint;
      else 
				return fshp + (essp-eshp)*Esup + (fsup-fshp)*(1.0-pow((esup-essp)/(esup-eshp),p));
		}
  } else
    return essp * ((Esp - Eypp) / pow(1 + pow((Esp - Eypp) * essp / fint,10.0),0.1) + Eypp);
}

double
ReinforcingSteel::Backbone_E(double e)
{
	double essp  = fabs(e);
	
	if(essp<=eshpa)
		return (Esp - Eypp) / pow(1.0 + pow((Esp - Eypp) * essp / fint,10.0),1.1) + Eypp;
	else {
		if(essp>esup)
			return Esup;
		else {
      if (essp < eshp+0.0002)
				return (Eshpb-Eypp)*(essp-eshpa)/(eshp+0.0002-eshpa) + Eypp;
      else {
        //double temp1 = (fsup-fshp-(fsup-fshp)*(1.0-pow((esup-essp)/(esup-eshp),p)));
				return Eshp*pow((fsup-fshp-(fsup-fshp)*(1.0-pow((esup-essp)/(esup-eshp),p)))/(fsup-fshp),1.0-1.0/p)+Esup;
      }
		}
	} 
}

/*****************************************************************************************/
/***********************       Buckled Stress-Strain Relations       *********************/
/*****************************************************************************************/

double 
ReinforcingSteel::Buckled_stress_Dhakal(double ess, double fss)
{
  if (LDratio <= 0.0) return fss;

  double aveStress;
  double e_cross = Temax - fsup/Esp;
  double e=ess-e_cross;

  if(e < -eyp) {
	    
		double eStar=55.0-2.3*sqrt(fyp/Esp*2000)*LDratio;
		if (eStar < 7.0) eStar=7.0;
		eStar= -eStar*eyp;
	  
		double fStarL=Backbone_f(eStar);
		double fStar=fStarL*beta*(1.1 - 0.016*sqrt(fyp/Esp*2000)*LDratio);
		if (fStar > -0.2*fyp) fStar= -0.2*fyp;
    if (TBranchNum%4 > 1) {
		  if ((e< -eyp) && (e>=eStar)) {
			  aveStress = fss*(1.0-(1.0-fStar/fStarL)*(e+eyp)/(eStar+eyp));
		  } else if (e<eStar) {
			  aveStress = fss*(fStar-0.02*Esp*(e-eStar))/fStarL;
			  if (aveStress>-0.2*fyp) aveStress=-0.2*fyp;
		  }
		  return aveStress;
    } else {
      if (TBranchNum == 4 || TBranchNum == 5)
        BackStress = MP_f(e_cross-eyp);

      if ((e< -eyp) && (e>=eStar)) {
			  aveStress = Tfa*(1.0-(1.0-fStar/fStarL)*(e+eyp)/(eStar+eyp));
		  } else if (e<eStar) {
			  aveStress = Tfa*(fStar-0.02*Esp*(e-eStar))/fStarL;
			  if (aveStress>-0.2*fyp) aveStress=-0.2*fyp;
		  }
		  return BackStress - (BackStress-fss)*(BackStress-aveStress)/(BackStress-Tfa);
      //return Cfa[0]-aveStress/(BackStress-Cfa[0])*(fss-Cfa[0]);
    }
  } else {
	  return fss;
  }	
}

double
ReinforcingSteel::Buckled_stress_Gomes(double ess, double fss)
{
	if (LDratio <= 0.0) return fss;
	
	double e_cross = Temax - fsup/Esp;
	if (ess>=e_cross) return fss;

	double beta=1.0;
	double gama = 0.1;
	double Dft = 0.25;
		
	double fs_buck = beta*sqrt(32.0/(e_cross-ess))/(9.42477796076938*LDratio);
	double stress_diff=fabs(fs_buck-1.0);
	if (stress_diff <= Dft)	beta = 1-gama*(Dft-stress_diff)/Dft;

	//double factor = ((1.0>fs_buck)?fs_buck:1.0)*beta + reduction)/(1.0+reduction);
	double factor = ((1.0>fs_buck)?fs_buck:1.0)*beta*(1-reduction)+reduction;

	double t_s_out = fsup*fsu_fraction-(factor+fsu_fraction)*(fsup*fsu_fraction-fss)/(1.0+fsu_fraction);
	return t_s_out;
}

double
ReinforcingSteel::Buckled_mod_Gomes(double ess, double fss, double Ess)
{
	double Etmp = Ess + (Buckled_stress_Gomes(ess+0.00005, fss)-Buckled_stress_Gomes(ess-0.00005, fss))/0.0001;
	return Etmp;
}
double
ReinforcingSteel::Buckled_mod_Dhakal(double ess, double fss, double Ess)
{
	double Etmp = Ess + (Buckled_stress_Dhakal(ess+0.00005, fss)-Buckled_stress_Dhakal(ess-0.00005, fss))/0.0001;
	return Etmp;
}
/*****************************************************************************************/
/*************************        Branch Rule Definitions          ***********************/
/*****************************************************************************************/

int
ReinforcingSteel::BranchDriver(int res)
{
  switch(TBranchNum) {
  case  1:	res += Rule1(res);
			break;
  case  2:	res += Rule2(res);
			break;
  case  3:	res += Rule3(res);
			break;
  case  4:	res += Rule4(res);
			break;
  case  5:	res += Rule5(res);
			break;
  case  6:	res += Rule6(res);
			break;
  case  7:	res += Rule7(res);
			break;
  case  8:	res += Rule8(res);
			break;
  case -1:  TStress = 0.0;
			TTangent = Esp/1000000.0;
			break;
  case  0:  TStress = 0.0;
			TTangent = Esp;
			break;
  default:	switch(TBranchNum%4) {
			case 0: res += Rule12(res);
			        break;
			case 1: res += Rule9(res);
			        break;
			case 2: res += Rule10(res);
			        break;
			case 3: res += Rule11(res);
					break;
			}
			break;
  }
  return res;
}

// Rule 1: Tension Envelope Branch
int
ReinforcingSteel::Rule1(int res)
{
  double strain=TStrain-Teo_p;
  // check for load reversal
  if (TStrain-CStrain<0.0) {
	  if (strain - eshp > -ZeroTol) {
	    double emin;
	    // reversal from strain hardening range
	    Tea=CStrain;
	    Temax=Tea-Teo_p;
	    if (CStrain > TeAbsMax) TeAbsMax = CStrain;

	    if(Temin>-eshp)
		    emin=-eshp-1.0E-14;
	    else
		    emin=Temin;

	    double ea   = Teo_p + eshp - fshp/Esp;
	    double eb   = Teo_p + Temax - CStress/Esp;
	    double krev = exp(-Temax/(5000*eyp*eyp));
	    double eon = ea*krev+eb*(1.0-krev);
      if (eon > Teo_n) {
        emin-=(eon-Teo_n);
        Teo_n=eon;
      }
	    Teb=Teo_n+emin;

	    // set stress dependent curve parameters
	    Tfa=CStress;
      Cfa[0]=CStress;
	    TEa=ReturnSlope(Tea-Teo_n-Temin);

  	  updateHardeningLoaction(TeCumPlastic+Tea-emin-(Tfa-Backbone_f(emin))/Esp);
	    Tfb= Backbone_f(emin);
	    TEb= Backbone_E(emin);
		  TEsec = (Tfb-Tfa)/(Teb-Tea);

		  if(TEsec < TEb) {
			  Teo_n = (Tfb-Tfa)/TEb+Tea - emin;
			  Teb=Teo_n+emin;
			  TEsec = (Tfb-Tfa)/(Teb-Tea);
			  opserr << "Adjusted Compressive Curve anchor in ReinforcingSteel::Rule1()\n";
		  }

	    SetTRn();
	    res += SetMP();
  	  
	    T_ePlastic[2]=0.0;
	    TBranchNum=3;
	    Rule3(res);	  
	  } else if (strain - eyp > -ZeroTol) {
	    double emin;
	    // Reversal from Yield Plateau
	    Tea=CStrain;
	    Temax=Tea-Teo_p;
	    if (CStrain > TeAbsMax) TeAbsMax = CStrain;

	    Tfa=CStress;
      Cfa[0]=CStress;
	    TEa=ReturnSlope(Tea-Teo_n-Temin);

	    double pr=(Temax-eyp)/(eshp-eyp);
	    emin=pr*(eyp-eshp)-eyp;
	    Teo_n=Tea-Tfa/Esp;
	    Teb=Teo_n+emin;
  	  
	    // set stress dependent curve parameters
      updateHardeningLoaction(TeCumPlastic+Tea-emin-(Tfa-Backbone_f(emin))/Esp);
	    Tfb=Backbone_f(emin);
	    TEb=(1.0/(1.0/Esp+pr*(1.0/Eshp - 1.0/Esp)));
  	  
	    SetTRn();
	    TEsec = (Tfb-Tfa)/(Teb-Tea);
	    if (TEsec<TEb) TEb=TEsec*0.999;
	    if (TEsec>TEa) TEa=TEsec*1.001;
	    res += SetMP();
  	  
	    T_ePlastic[2]=0.0;
	    TBranchNum=3;
	    Rule3(res);
	  } else if (strain > -ZeroTol) {
	    //if(Temax < strain) Temax=strain;
	    TStress  = Backbone_f(strain);
	    TTangent = Backbone_E(strain);
	  } else {
	    TBranchNum=2;
	    Rule2(res);
	  }
  } else {
    TStress  = Backbone_f(strain);
	  TTangent = Backbone_E(strain);
	  //if(Temin<0.0) {
	    TFatDamage-=damage(T_ePlastic[0]);
      TeCumPlastic -= T_ePlastic[0];
	    T_ePlastic[0]=getPlasticStrain(TStrain-TeAbsMin,TStress-Cfa[1]);
	    TFatDamage+=damage(T_ePlastic[0]);
      TeCumPlastic += T_ePlastic[0];
	  //}
  }
  return res;
}

// Rule 2: Compression Envelope Branch
int
ReinforcingSteel::Rule2(int res)
{
  double strain = TStrain-Teo_n;
 
  // check for load reversal
  if (TStrain-CStrain>0.0) {
	if (strain+eshp< ZeroTol) {
	  double emax;
	  // reversal from strain hardening range
	  Tea=CStrain;
	  Temin=Tea-Teo_n;
	  if (CStrain < TeAbsMin) TeAbsMin = CStrain;
	  
	  if(Temax<eshp) 
		emax=eshp+1.0E-14;
	  else
		emax=Temax;

	  double ea   = Teo_n - eshp + fshp/Esp;
	  double eb   = Teo_n + Temin - CStress/Esp;
	  double krev = exp(Temin/(5000*eyp*eyp));
	  double eop=ea*krev+eb*(1.0-krev);
    if (eop<Teo_p) {
      emax+=(Teo_p-eop);
      Teo_p=eop;
    }
      Teb=Teo_p+emax;
	  
	  // set stress dependent curve parameters 
	  Tfa=CStress;
    Cfa[1]=CStress;
	  TEa=ReturnSlope(Temax + Teo_p -Tea);

    updateHardeningLoaction(TeCumPlastic+emax-Tea-(Backbone_f(emax)-Tfa)/Esp);
	  Tfb= Backbone_f(emax);
	  TEb= Backbone_E(emax);
	  
	  SetTRp();
	  TEsec = (Tfb-Tfa)/(Teb-Tea);
	  res += SetMP();
		
	  T_ePlastic[2]=0.0;
	  TBranchNum=4;
	  Rule4(res);

	} else if (strain+eyp < ZeroTol) {
	  double emax;
	  // reversal form yield plateau
	  Tea=CStrain;
	  Temin=Tea-Teo_n;
	  if (CStrain < TeAbsMin) TeAbsMin = CStrain;

	  Tfa=CStress;
    Cfa[1]=CStress;
	  TEa=ReturnSlope(Temax + Teo_p -Tea);

	  double pr=(Temin+eyp)/(eyp-eshp);
	  emax=eyp+pr*(eshp-eyp);
	  Teo_p=Tea-Tfa/Esp;
	  Teb=Teo_p+emax;
	  
	  // stress dependent curve parameters
    updateHardeningLoaction(TeCumPlastic+emax-Tea-(Backbone_f(emax)-Tfa)/Esp);
	  Tfb= Backbone_f(emax);
	  TEb=1.0/(1.0/Esp+pr*(1.0/Eshp - 1.0/Esp));
	  
	  SetTRp();
	  TEsec = (Tfb-Tfa)/(Teb-Tea);
	  if (TEsec<TEb) TEb=TEsec*0.999;
	  if (TEsec>TEa) TEa=TEsec*1.001;
	  res += SetMP();
	  
	  T_ePlastic[2]=0.0;
	  TBranchNum=4;
	  Rule4(res);

	} else if (strain<ZeroTol) {
	  //if(Temin>strain) Temin=strain;
	  TStress  = Backbone_f(strain);
	  TTangent = Backbone_E(strain);
	} else {
	  TBranchNum=1;
	  Rule1(res);
	}
  } else {
    TStress  = Backbone_f(strain);
    TTangent = Backbone_E(strain);
	  //if(Temax>0.0) {
	    TFatDamage-=damage(T_ePlastic[1]);
      TeCumPlastic -= T_ePlastic[1];
	    T_ePlastic[1]=getPlasticStrain(TeAbsMax-TStrain,Cfa[0]-TStress);
	    TFatDamage+=damage(T_ePlastic[1]);
      TeCumPlastic += T_ePlastic[1];
	  //}
  }
  return res;
}

// Rule 3: Unloading Reversal Branch
int
ReinforcingSteel::Rule3(int res)
{
  if (TStrain-CStrain > 0.0) {	// reversal from branch 
	if(Temin > CStrain-Teo_n) Temin=CStrain-Teo_n;

	Tea=CStrain;
	double dere = Cea[2]-Tea-fyp/(1.2*Esp);
    if (dere<0.0)
	  dere=0.0;
	else if (dere>fyp/3/Esp)
	  dere=fyp/3/Esp;
	Teb=Teo_p+Temax+dere;

	Tfa=CStress;
	TEa=ReturnSlope(Cea[2]-CStrain);

  updateHardeningLoaction(TeCumPlastic+Teb-Tea-(Backbone_f(Teb-Teo_p)-Tfa)/Esp);
	Tfb= Backbone_f(Teb-Teo_p);
	TEb= Backbone_E(Teb-Teo_p);

	SetTRp();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
	#ifdef _WIN32
  if(_fpclass(TStress)< 8 || _fpclass(TStress)==512 || _fpclass(TTangent)< 8 || _fpclass(TTangent)==512) {
    opserr << "bad stress or tangent\n";
    return -1;
  }
#endif
	T_ePlastic[3]=0.0;
	TBranchNum=5;
	Rule5(res);
  } else {
	  if (TStrain - Teb <= ZeroTol) {
	    T_ePlastic[1]=T_ePlastic[2];
	    TBranchNum=2;
	    Rule2(res);
	  } else {
      TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage -=damage(T_ePlastic[2]);
      TeCumPlastic -= T_ePlastic[2];
	    T_ePlastic[2]=getPlasticStrain(TeAbsMax-TStrain,Tfa-TStress);
	    TFatDamage +=damage(T_ePlastic[2]);
      TeCumPlastic += T_ePlastic[2];
	  }
  }
  return res;
}

// Rule 4: Loading Reversal Branch
int
ReinforcingSteel::Rule4(int res)
{
  if (TStrain-CStrain < 0.0) {
	if(Temax<CStrain-Teo_p) Temax=CStrain-Teo_p;

	Tea=CStrain;
	double dere = Cea[2]-Tea+fyp/(1.2*Esp);
    if (dere>0.0)
	  dere=0.0;
	else if (dere<-fyp/3/Esp)
	  dere=-fyp/3/Esp;
	Teb=Teo_n+Temin+dere;

	Tfa=CStress;
	TEa=ReturnSlope(CStrain-Cea[2]); 

  updateHardeningLoaction(TeCumPlastic+Tea-Teb-(Tfa-Backbone_f(Teb-Teo_n))/Esp);
	Tfb= Backbone_f(Teb-Teo_n);
	TEb= Backbone_E(Teb-Teo_n);
	
	SetTRn();
	TEsec = (Tfb-Tfa)/(Teb-Tea); 
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
	
	T_ePlastic[3]=0.0;
	TBranchNum=6;
	Rule6(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
	    T_ePlastic[0]=T_ePlastic[2];
	    TBranchNum=1;
	    Rule1(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=damage(T_ePlastic[2]);
      TeCumPlastic -= T_ePlastic[2];
	    T_ePlastic[2]=getPlasticStrain(TStrain-TeAbsMin,TStress-Tfa);
	    TFatDamage+=damage(T_ePlastic[2]);
      TeCumPlastic += T_ePlastic[2];
	  }
  }
  return res;
}

// Rule 5: Loading Returning Branch
int
ReinforcingSteel::Rule5(int res)
{
  if (TStrain-CStrain < 0.0) {
	rE1=0.0;
	rE2=0.0;
	
	Tea   = Ceb[3]*(CStrain-Cea[3])/(Ceb[3]-Cea[3]) + Cea[2]*(Ceb[3]-CStrain)/(Ceb[3]-Cea[3]);
	Teb   = Ceb[2];
	
  updateHardeningLoaction(TeCumPlastic+CStrain-Tea+(Backbone_f(Tea-Teo_p)-CStress)/Esp);
	Tfa   = Backbone_f(Tea-Teo_p);
	TEa   = CEa[2];
	
  updateHardeningLoaction(TeCumPlastic+CStrain-Teb-(CStress-Backbone_f(Teb-Teo_n))/Esp);
	Tfb   = Backbone_f(Teb-Teo_n);
	TEb   = Backbone_E(Teb-Teo_n);

	SetTRn();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	res += SetMP();

	double fb=MP_f(Cea[3]);
	double Eb=MP_E(Cea[3]);
	
	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(CStrain-Cea[3]);
	Teb=Cea[3];
	Tfb=fb;
	TEb=Eb;

	SetTRn();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
  
	T_ePlastic[4]=0.0;
	TBranchNum=7;
	Rule7(res);	
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
	    TFatDamage-=damage(T_ePlastic[3]);
      TeCumPlastic -= T_ePlastic[3];
      double TempPStrain = getPlasticStrain(Teb-Tea,Tfb-Tfa);
	    TFatDamage+=damage(TempPStrain);
      TeCumPlastic += TempPStrain;
	    TBranchNum=1;
	    Rule1(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=damage(T_ePlastic[3]);
      TeCumPlastic -= T_ePlastic[3];
	    T_ePlastic[3]=getPlasticStrain(TStrain-Tea,TStress-Tfa);
	    TFatDamage+=damage(T_ePlastic[3]);
      TeCumPlastic += T_ePlastic[3];
	  }
  }
  return res;
}

// Rule 6: Unloading Returning Branch
int
ReinforcingSteel::Rule6(int res)
{
  if (TStrain-CStrain > 0.0) {	
	rE1=0.0;
	rE2=0.0;

	Tea   = Ceb[3]*(CStrain-Cea[3])/(Ceb[3]-Cea[3]) + Cea[2]*(Ceb[3]-CStrain)/(Ceb[3]-Cea[3]);
	Teb   = Ceb[2];
	
  updateHardeningLoaction(TeCumPlastic+Tea-CStrain+(CStress-Backbone_f(Tea-Teo_n))/Esp);
	Tfa   = Backbone_f(Tea-Teo_n);
	TEa   = CEa[2];

  updateHardeningLoaction(TeCumPlastic+Teb-CStrain-(Backbone_f(Teb-Teo_p)-CStress)/Esp);
	Tfb   = Backbone_f(Teb-Teo_p);
	TEb   = Backbone_E(Teb-Teo_p);

	SetTRp();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	res += SetMP();

	double fb=MP_f(Cea[3]);
	double Eb=MP_E(Cea[3]);

	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(Cea[3]-CStrain);
	Teb=Cea[3];
	Tfb=fb;
	TEb=Eb;  

	SetTRp();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
	
	T_ePlastic[4]=0.0;
	TBranchNum=8;
	Rule8(res);
  } else {
	  if (TStrain - Teb <= ZeroTol) {
	    TFatDamage-=damage(T_ePlastic[3]);
      TeCumPlastic -= T_ePlastic[3];
      double TempPStrain = getPlasticStrain(Tea-Teb,Tfa-Tfb);
	    TFatDamage+=damage(TempPStrain);
      TeCumPlastic += TempPStrain;
	    TBranchNum=2;
	    Rule2(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=damage(T_ePlastic[3]);
      TeCumPlastic -= T_ePlastic[3];
	    T_ePlastic[3]=getPlasticStrain(Tea-TStrain,Tfa-TStress);
	    TFatDamage+=damage(T_ePlastic[3]);
      TeCumPlastic += T_ePlastic[3];
	  }
  }
  return res;
}

// Rule 7: Unloading First Transition Branch
int
ReinforcingSteel::Rule7(int res)
{
  if (TStrain-CStrain > 0.0) {	
	SetPastCurve(TBranchNum-2);
	
	double fb=MP_f(Cea[4]);
	double Eb=MP_E(Cea[4]);
	
	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(Cea[4]-CStrain);
	Teb=Cea[4];
	Tfb=fb;
	TEb=Eb;

	SetTRp1();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();

	re=Tea;

	T_ePlastic[5]=0.0;
	TBranchNum=9;
	Rule9(res);	
  } else {
	  if (TStrain - Teb <= ZeroTol) {
	    TFatDamage-=damage(T_ePlastic[4]);
      TeCumPlastic -= T_ePlastic[4];
      double TempPStrain = getPlasticStrain(Tea-Teb,Tfa-Tfb);
	    TFatDamage+=damage(TempPStrain);
      TeCumPlastic += TempPStrain;
      double tempTeb = Teb;

	    Tea   = Ceb[3]*(Tea-Cea[3])/(Ceb[3]-Cea[3]) + Cea[2]*(Ceb[3]-Tea)/(Ceb[3]-Cea[3]);
      Teb   = Ceb[2];

      updateHardeningLoaction(TeCumPlastic+tempTeb-Tea+(Backbone_f(Tea-Teo_p)-Tfb)/Esp);
	    Tfa   = Backbone_f(Tea-Teo_p);
	    TEa   = CEa[2];

      updateHardeningLoaction(TeCumPlastic+tempTeb-Teb-(Tfb-Backbone_f(Teb-Teo_n))/Esp);
	    Tfb   = Backbone_f(Teb-Teo_n);
	    TEb   = Backbone_E(Teb-Teo_n);

	    SetTRn();
	    TEsec = (Tfb-Tfa)/(Teb-Tea);
	    res += SetMP();

	    TBranchNum=3;
	    Rule3(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=damage(T_ePlastic[4]);
      TeCumPlastic -= T_ePlastic[4];
	    T_ePlastic[4]=getPlasticStrain(Tea-TStrain,Tfa-TStress);
	    TFatDamage+=damage(T_ePlastic[4]);
      TeCumPlastic += T_ePlastic[4];
	  }
  }
  return res;
}

// Rule 8: Loading First Transition Branch
int
ReinforcingSteel::Rule8(int res)
{
  if (TStrain-CStrain < 0.0) {
	SetPastCurve(TBranchNum-2);
	
	double fb=MP_f(Cea[4]);
	double Eb=MP_E(Cea[4]);

	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(CStrain-Cea[4]);
	Teb=Cea[4];
	Tfb=fb;
	TEb=Eb;  

	SetTRn1();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
	
	re=Tea;

	T_ePlastic[5]=0.0;
	TBranchNum=10;
	Rule10(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
	    TFatDamage-=damage(T_ePlastic[4]);
      TeCumPlastic -= T_ePlastic[4];
      double TempPStrain = getPlasticStrain(Teb-Tea,Tfb-Tfa);
	    TFatDamage+=damage(TempPStrain);
      TeCumPlastic += TempPStrain;
  	  double tempTeb = Teb;

	    Tea   = Ceb[3]*(Tea-Cea[3])/(Ceb[3]-Cea[3]) + Cea[2]*(Ceb[3]-Tea)/(Ceb[3]-Cea[3]);
	    Teb   = Ceb[2];

      updateHardeningLoaction(TeCumPlastic+Tea-tempTeb+(Tfb-Backbone_f(Tea-Teo_n))/Esp);
	    Tfa   = Backbone_f(Tea-Teo_n);
	    TEa   = CEa[2];

      updateHardeningLoaction(TeCumPlastic+Teb-tempTeb-(Backbone_f(Teb-Teo_p)-Tfb)/Esp);
	    Tfb   = Backbone_f(Teb-Teo_p);
	    TEb   = Backbone_E(Teb-Teo_p);

	    SetTRp();
	    TEsec = (Tfb-Tfa)/(Teb-Tea);
	    res += SetMP();
  	  
	    TBranchNum=4;
	    Rule4(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=damage(T_ePlastic[4]);
      TeCumPlastic -= T_ePlastic[4];
	    T_ePlastic[4]=getPlasticStrain(TStrain-Tea,TStress-Tfa);
	    TFatDamage+=damage(T_ePlastic[4]);
      TeCumPlastic += T_ePlastic[4];
	  }
  }
  return res;
}

// Rule 9: Loading Second Transition Branch
int
ReinforcingSteel::Rule9(int res)
{
  if (TStrain-CStrain < 0.0) {
	double eb = Tea;
	if (TBranchNum+4<=LastRule_RS) re=Tea;
	SetPastCurve(TBranchNum-2);
	
	// set new curve 
	double fb=MP_f(re);
	double Eb=MP_E(re);
	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(CStrain-eb);;
	Teb=re;
	Tfb=fb;
	TEb=Eb;
	SetTRn1();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();

	TBranchNum+=2;
	TBranchMem = (TBranchNum+1)/2;
	T_ePlastic[TBranchMem]=0.0;
	Rule11(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage -=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic -= T_ePlastic[TBranchMem];
      double TempPStrain = getPlasticStrain(Teb-Tea,Tfb-Tfa);
	    TFatDamage +=damage(TempPStrain);
      TeCumPlastic += TempPStrain;
		  TBranchNum-=4;
	    SetPastCurve(TBranchNum);
	    if (TBranchNum==5)
		  Rule5(res);
	    else
		  Rule9(res);
	  } else {
      TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage -=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic -= T_ePlastic[TBranchMem];
	    T_ePlastic[TBranchMem]=getPlasticStrain(TStrain-Tea,TStress-Tfa);
	    TFatDamage +=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic += T_ePlastic[TBranchMem];
	  }
  }
  return res;
}

// Rule 10: Unloading Second Transition Branch
int
ReinforcingSteel::Rule10(int res)
{
  if (TStrain-CStrain > 0.0) {	
	double eb = Tea;
	if (TBranchNum+4<=LastRule_RS)
	  re=Tea;
	
	SetPastCurve(TBranchNum-2);

	// set new curve
	double fb=MP_f(re);
	double Eb=MP_E(re);
	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(eb-CStrain); 
	Teb=re;
	Tfb=fb;
	TEb=Eb;  
	SetTRp1();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();

	TBranchNum+=2;
	TBranchMem = (TBranchNum+1)/2;
	T_ePlastic[TBranchMem]=0.0;
	Rule12(res);
  } else {
	  if (TStrain - Teb <= ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic -= T_ePlastic[TBranchMem];
      double TempPStrain = getPlasticStrain(Tea-Teb,Tfa-Tfb);
	    TFatDamage+=damage(TempPStrain);
      TeCumPlastic += TempPStrain;

		  TBranchNum-=4;
	    SetPastCurve(TBranchNum);
	    if (TBranchNum==6)
		  Rule6(res);
	    else
		  Rule10(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage  -=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic -= T_ePlastic[TBranchMem];
	    T_ePlastic[TBranchMem]=getPlasticStrain(Tea-TStrain,Tfa-TStress);
	    TFatDamage  +=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic += T_ePlastic[TBranchMem];
	  }
  }
  return res;
}

// Rule 11: Unloading Third Transition Branch
int
ReinforcingSteel::Rule11(int res)
{
  if (TStrain-CStrain > 0.0) {	
	// reset past curve
	double eb=Tea;
	if(TBranchNum+2>LastRule_RS) {
		TBranchMem = (TBranchNum+1)/2;
	  eb=Cea[TBranchMem-2];
	  SetPastCurve(TBranchNum-6);
	} else {
	  SetPastCurve(TBranchNum-2);
	}
	double fb=MP_f(eb);
	double Eb=MP_E(eb);
	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(eb-CStrain);
	Teb=eb;
	Tfb=fb;
	TEb=Eb;
	SetTRp1();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
	
	if(TBranchNum+2>LastRule_RS)
	  TBranchNum-=2;
	else
	  TBranchNum+=2;

	TBranchMem = (TBranchNum+1)/2;
	T_ePlastic[TBranchMem]=0.0;
	Rule9(res);	
  } else {
	  if (TStrain - Teb <= ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=damage(T_ePlastic[TBranchMem-2]);
      TeCumPlastic -= T_ePlastic[TBranchMem-2];
      double TempPStrain = getPlasticStrain(Tea-Teb,Tfa-Tfb);
	    TFatDamage+=damage(TempPStrain);
      TeCumPlastic += TempPStrain;
		  TBranchNum-=4;
	    SetPastCurve(TBranchNum);
	    if (TBranchNum==7)
		  Rule7(res);
	    else
		  Rule11(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic -= T_ePlastic[TBranchMem];
	    T_ePlastic[TBranchMem]=getPlasticStrain(Tea-TStrain,Tfa-TStress);
	    TFatDamage+=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic += T_ePlastic[TBranchMem];
	  }
  }
  return res;
}

// Rule 12: Loading Third Transition Branch
int
ReinforcingSteel::Rule12(int res)
{
  if (TStrain-CStrain < 0.0) {
	// reset past curve
	double eb=Tea;
	if(TBranchNum+2>LastRule_RS) {
		TBranchMem = (TBranchNum+1)/2;
	  eb=Cea[TBranchMem-2];
	  SetPastCurve(TBranchNum-6);
	} else {
	  SetPastCurve(TBranchNum-2);
	}

	double fb=MP_f(eb);
	double Eb=MP_E(eb);
	Tea=CStrain;
	Tfa=CStress;
	TEa=ReturnSlope(CStrain-eb);
	Teb=eb;
	Tfb=fb;
	TEb=Eb; 
	SetTRn1();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();

	if(TBranchNum+2>LastRule_RS)
	  TBranchNum-=2;
	else
	  TBranchNum+=2;

	TBranchMem = (TBranchNum+1)/2;
	T_ePlastic[TBranchMem]=0.0;
	Rule10(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=damage(T_ePlastic[TBranchMem-2]);
      TeCumPlastic -= T_ePlastic[TBranchMem-2];
      double TempPStrain = getPlasticStrain(Teb-Tea,Tfb-Tfa);
	    TFatDamage+=damage(TempPStrain);
      TeCumPlastic += TempPStrain;
		  TBranchNum-=4;
	    SetPastCurve(TBranchNum);
	    if (TBranchNum==8)
		  Rule8(res);
	    else
		  Rule12(res);
	  } else {
      TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic -= T_ePlastic[TBranchMem];
	    T_ePlastic[TBranchMem]=getPlasticStrain(TStrain-Tea,TStress-Tfa);
	    TFatDamage+=damage(T_ePlastic[TBranchMem]);
      TeCumPlastic += T_ePlastic[TBranchMem];
	  }
  }
  return res;
}



/*****************************************************************************************/
/********        Strength and stiffness degradation including buckling         ***********/
/*****************************************************************************************/
double
ReinforcingSteel::getPlasticStrain(double ehalf, double stressAmp){
  double ehalfPlastic = fabs(ehalf)-fabs(stressAmp/Esp);
  if(ehalfPlastic>0.0)
    return ehalfPlastic;
  else
    return 0.0;
}
double
ReinforcingSteel::damage(double ehalfPlastic){
    return pow(ehalfPlastic/Fat1,Fat2);
}

double
ReinforcingSteel::scalefactor()
{
  if(theBarFailed) return 0.0;

  double sf=1.0-Deg1*TFatDamage;
  if(TFatDamage>1.0) sf-= (TFatDamage-1.0)/0.04;

  if(sf<0.0) {
    theBarFailed=1;
    TBranchNum = -1;
    opserr << "-------------------------Bar failed---------------------------\n";
		return 0.0;
  }	else {
		return sf;
  }
}

double
ReinforcingSteel::ReturnSlope(double dea)
{
  if (TeAbsMax > -TeAbsMin) //Dodd and Cooke
    return Esp*(0.82+1.0/(5.55+1000.0*TeAbsMax));
  else
    return Esp*(0.82+1.0/(5.55-1000.0*TeAbsMin));
  //return Esp;

}
