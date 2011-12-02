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
// $Date: 2005-08-23 17:20:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ReinforcingSteel.cpp,v $

// Written: Jon Mohle
// Created: October 2003
//
// Description: This file contains the class definition for 
// ReinforcingSteel

#include "ReinforcingSteel.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

ReinforcingSteel::ReinforcingSteel(int tag, double fy, double fsu, double Es, double Esh, double esh, double esu, 
					 double slenderness, double alpha, double r, double gama, double Fatigue1, double Fatigue2, double Degrade)
  :UniaxialMaterial(tag,MAT_TAG_ReinforcingSteel),
  LDratio(slenderness),alpha_(alpha),fsu_fraction(gama),Fat1(Fatigue1)
{ 
  if((r>=0.0) & (r<=1.0)) 
    reduction=r;
  else
    if(r<=0)
      reduction = 0.0;
    else
      reduction = 1.0;

	if((Fatigue1==0) || (Fatigue2==0)) {
		Fat1=0.0;
		Fat2=1.0;
		Deg1=0.0;
  } else {
    Fat2=1.0/Fatigue2;
    if(Degrade != 0.0) 
      Deg1=pow(Fat1/Degrade,Fat2);
    else
      Deg1=0.0;
  }

	double estep = 0.01;
	double esg = esu+(esu-esh)*0.1;
	int sgn = 1;
	p = Esh*(esu-esh)/(fsu-fy);
	
	for(int i=0;i<20;i++) {
        fyp = fy + (fsu - fy) * (1 - pow(fabs((esu - esg) / (esu - esh)),p));
        eyp = fyp-Esh*(1.0+esg)*pow((fsu-fyp)/(fsu-fy),(p-1.0)/p);
				if(sgn != Sign(eyp)){
          sgn = -1*sgn;
          estep = -estep/4.0;
				}
        esg += estep;
	}

	eyp=log(1.0+fy/Es);
	fyp=fy*(1.0+fy/Es);
	Esp=fyp/eyp;
	eshp=log(1.0+esh);
	fshp=fy*(1.0+esh);
	Eshp=Esh*pow(1.0+esh,2.0)+fshp;
	Eypp=(fshp-fyp)/(eshp-eyp);
	esup=log(1.0+esg)*1.8;
	fint = fyp-Eypp*eyp;
	fsup=(1.0+esg)*(Es*esg/pow(pow(Es*esg/fy,10.0)+1.0,0.1)+
	     static_cast<double>((Sign(esg-esh)+1)/2)*(fsu-fy)*(1.0-pow(fabs((esu-esg)/(esu-esh)),p)));
	p = Eshp*(esup-eshp)/(fsup-fshp);

	// Set backbone transition variables
	double fTemp = Backbone_fNat(eshp+0.0002);
	Eshpb = Eshp*pow((fsup-fTemp)/(fsup-fshp),1.0-1.0/p);
  eshpa = eshp + 0.0002 - 2.0*(fTemp-fshp)/Eshpb;

  ZeroTol=1.0E-14;
  this->revertToStart();
}

ReinforcingSteel::ReinforcingSteel(int tag)
  :UniaxialMaterial(tag,MAT_TAG_ReinforcingSteel)
{
  ZeroTol=1.0E-14;
}

ReinforcingSteel::~ReinforcingSteel(){

}

/***************** material state determination methods ***********/
int 
ReinforcingSteel::setTrialStrain(double strain, double strainRate) {
  int res = 0;
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
  double tempstr = Buckled_stress_Gomes(TStrain,TStress);
  double tempOut = tempstr*scalefactor()/exp(TStrain);

#ifdef _WIN32
  if(_fpclass(tempOut)< 8 || _fpclass(tempOut)==512)
    opserr << "bad Stress in ReinforcingSteel::getStress\n";
 #endif

  return tempOut;
}

double 
ReinforcingSteel::getTangent(void) {
  double taTan = Buckled_mod_Gomes(TStrain,TStress,TTangent);
  double scfact = scalefactor();
  double tempOut = (taTan+TStress)*scfact/pow(exp(TStrain),2.0);
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

/***************** path dependent bahavior methods ***********/
int 
ReinforcingSteel::commitState(void) {
  if(TBranchNum <= 1)
	TBranchMem=0;
  else
	TBranchMem = (TBranchNum+1)/2;

  for(int i=0; i<=LastRule_RS/2; i++)
    C_FDamage[i]=T_FDamage[i];

  CFatDamage       = TFatDamage;
  
  // commit trial history variables
  CBranchNum = TBranchNum;
  Ceo_p      = Teo_p;
  Ceo_n      = Teo_n;
  Cemax      = Temax;
  Cemin      = Temin;
  CeAbsMax   = TeAbsMax;
  CeAbsMin   = TeAbsMin;

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

  // commit trial state variables
  CStrain    = TStrain;  
  CStress    = TStress;
  CTangent   = TTangent;
  return 0;
}

int 
ReinforcingSteel::revertToLastCommit(void) {
  for(int i=0; i<=LastRule_RS/2; i++)
    T_FDamage[i]=C_FDamage[i];

  TFatDamage = CFatDamage;

  // Reset trial history variables to last committed state
  TBranchNum = CBranchNum;
  Teo_p      = Ceo_p;
  Teo_n      = Ceo_n;
  Temax      = Cemax;
  Temin      = Cemin;
  TeAbsMax   = CeAbsMax;
  TeAbsMin   = CeAbsMin;

  if(TBranchNum > 2) SetPastCurve(TBranchNum);

  // Reset trial state variables to last committed state
  TStress    = CStress;
  TTangent   = CTangent;

  return 0;
}

int 
ReinforcingSteel::revertToStart(void)
{
  theBarFailed = 0;
  
  CFatDamage = TFatDamage;
  for(int i=0; i<=LastRule_RS/2; i++) {
    C_FDamage[i]  = 0.0;
    T_FDamage[i]  = 0.0;
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
  theCopy->alpha_       = alpha_;
  theCopy->theBarFailed = theBarFailed;

  // natural stress-strain variables
	theCopy->p            = p;
  theCopy->Esp          = Esp;   // Elastic Modulus
  theCopy->eshp         = eshp;	// Strain Hardening Strain
  theCopy->fshp         = fshp;  // Strain Hardening Stress
  theCopy->Eshp         = Eshp;  // Strain Hardening Modulus
  theCopy->esup         = esup;  // Strain at Peak Stress
  theCopy->fsup         = fsup;  // Peak Stress
  theCopy->Eypp         = Eypp;  // Yield Plateu Modulus
  theCopy->fint         = fint;  // Stress at yield plateu intersect
  theCopy->eyp          = eyp;
  theCopy->fyp          = fyp;

  theCopy->eshpa        = eshpa;	// Curve smoothing Parameters (at SH transition)
  theCopy->Eshpb        = Eshpb;	// These are used to eliminate a sudden discontinuity in stiffness

  //double Nbf;               // Cyclic Backbone factor used correct backbone proporsional to return strain
  theCopy->TFatDamage   = TFatDamage;
  theCopy->CFatDamage   = CFatDamage;
  theCopy->LDratio      = LDratio;
	theCopy->Fat1         = Fat1;
  theCopy->Fat2         = Fat2;
	theCopy->Deg1         = Deg1;
  theCopy->Deg2         = Deg2;

  theCopy->TBranchMem   =TBranchMem;
  theCopy->TBranchNum   =TBranchNum;
  theCopy->Teo_p        =Teo_p;
  theCopy->Teo_n        =Teo_n;
  theCopy->Temax        =Temax;
  theCopy->Temin        =Temin;
  theCopy->TeAbsMax     =TeAbsMax;
  theCopy->TeAbsMin     =TeAbsMin;

  theCopy->CBranchNum   = CBranchNum;
  theCopy->Ceo_p        = Ceo_p;
  theCopy->Ceo_n        = Ceo_n;
  theCopy->Cemax        = Cemax;
  theCopy->Cemin        = Cemin;
  theCopy->CeAbsMax     = CeAbsMax;
  theCopy->CeAbsMin     = CeAbsMin;

  for(int i=0; i<=LastRule_RS/2; i++) {
	  theCopy->C_FDamage[i] = C_FDamage[i];
    theCopy->T_FDamage[i] = T_FDamage[i];
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
    
  return theCopy;
}

int 
ReinforcingSteel::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(59+12*LastRule_RS/2);
  
  data(0) = this->getTag();
  data(1) = reduction;
  data(2) = fsu_fraction;
  data(3) = alpha_;
  data(4) = theBarFailed;
	data(5) = p;
  data(6) = Esp;
  data(7) = eshp;
  data(8) = fshp;
  data(9) = Eshp;
  data(10) = esup;
  data(11) = fsup;
  data(12) = Eypp;
  data(13) = fint;
  data(14) = eyp;
  data(15) = fyp;
  data(16) = eshpa;
  data(17) = Eshpb;
  data(18) = TFatDamage;
  data(19) = CFatDamage;
  data(20) = LDratio;
	data(21) = Fat1;
  data(22) = Fat2;
	data(23) = Deg1;
  data(24) = Deg2;
  data(25) = TBranchMem;
  data(26) = TBranchNum;
  data(27) = Teo_p;
  data(28) = Teo_n;
  data(29) = Temax;
  data(30) = Temin;
  data(31) = TeAbsMax;
  data(32) = TeAbsMin;
  data(33) = CBranchNum;
  data(34) = Ceo_p;
  data(35) = Ceo_n;
  data(36) = Cemax;
  data(37) = Cemin;
  data(38) = CeAbsMax;
  data(39) = CeAbsMin;
  data(40) = TR;
  data(41) = Tfch;
  data(42) = TQ;
  data(43) = TEsec;
  data(44) = Tea;
  data(45) = Tfa;
  data(46) = TEa;
  data(47) = Teb;
  data(48) = Tfb;
  data(49) = TEb;

  data(50) = re;
  data(51) = rE1;
  data(52) = rE2;

  data(53) = CStrain;  
  data(54) = CStress;
  data(55) = CTangent;
  data(56) = TStrain;  
  data(57) = TStress;
  data(58) = TTangent;

  for(int i=0; i<=LastRule_RS/2; i++) {
	  data(59+i*12) = C_FDamage[i];
    data(60+i*12) = T_FDamage[i];
	  data(61+i*12) = CR[i];
	  data(62+i*12) = Cfch[i];
	  data(63+i*12) = CQ[i];
	  data(64+i*12) = CEsec[i];
	  data(65+i*12) = Cea[i];
	  data(66+i*12) = Cfa[i];
	  data(67+i*12) = CEa[i];
	  data(68+i*12) = Ceb[i];
	  data(69+i*12) = Cfb[i];
	  data(70+i*12) = CEb[i];
  }

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
  
  static Vector data(59+12*LastRule_RS/2);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ReinforcingSteel::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    reduction = data(1);
    fsu_fraction = data(2);
    alpha_ = data(3);
    theBarFailed = static_cast<int>(data(4));
	  p = data(5);
    Esp = data(6);
    eshp = data(7);
    fshp = data(8);
    Eshp = data(9);
    esup = data(10);
    fsup = data(11);
    Eypp = data(12);
    fint = data(13);
    eyp = data(14);
    fyp = data(15);
    eshpa = data(16);
    Eshpb = data(17);
    TFatDamage = data(18);
    CFatDamage = data(19);
    LDratio = data(20);
	  Fat1 = data(21);
    Fat2 = data(22);
	  Deg1 = data(23);
    Deg2 = data(24);
    TBranchMem = static_cast<int>(data(25));
    TBranchNum = static_cast<int>(data(26));
    Teo_p = data(27);
    Teo_n = data(28);
    Temax = data(29);
    Temin = data(30);
    TeAbsMax = data(31);
    TeAbsMin = data(32);
    CBranchNum = static_cast<int>(data(33));
    Ceo_p = data(34);
    Ceo_n = data(35);
    Cemax = data(36);
    Cemin = data(37);
    CeAbsMax = data(38);
    CeAbsMin = data(39);
    TR = data(40);
    Tfch = data(41);
    TQ = data(42);
    TEsec = data(43);
    Tea = data(44);
    Tfa = data(45);
    TEa = data(46);
    Teb = data(47);
    Tfb = data(48);
    TEb = data(49);
    re = data(50);
    rE1 = data(51);
    rE2 = data(52);
    CStrain = data(53);  
    CStress = data(54);
    CTangent = data(55);
    TStrain = data(56);  
    TStress = data(57);
    TTangent = data(58);
    for(int i=0; i<=LastRule_RS/2; i++) {
	    C_FDamage[i] = data(59+i*12);
      T_FDamage[i] = data(60+i*12);
	    CR[i] = data(61+i*12);
	    Cfch[i] = data(62+i*12);
	    CQ[i] = data(63+i*12);
	    CEsec[i] = data(64+i*12);
	    Cea[i] = data(65+i*12);
	    Cfa[i] = data(66+i*12);
	    CEa[i] = data(67+i*12);
	    Ceb[i] = data(68+i*12);
	    Cfb[i] = data(69+i*12);
	    CEb[i] = data(70+i*12);
    }
  }
  return res;
}

void 
ReinforcingSteel::Print(OPS_Stream &s, int flag)
{
  if(flag == 3) {
	s << CStrain << "  " << CStress << "  " << CTangent << endln;
  } else {
    s << "ReinforcingSteel, tag: " << this->getTag() << endln;
    s << "  N2p: " << CFatDamage << endln;
    //s << "  sigmaY: " << sigmaY << endln;
    //s << "  Hiso: " << Hiso << endln;
    //s << "  Hkin: " << Hkin << endln;
    //s << "  eta: " << eta << endln;
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
  return TEb-TEsec*(1-pow(a,TR+1))/(1-a)+TEa*a*(1-pow(a,TR))/(1-a);
}

int
ReinforcingSteel::SetMP()
{
  double Rmin = (TEb-TEsec)/(TEsec-TEa);
  double a=0.01;
  double ao;
  double da;
  bool notConverge(true);
 
  //TEsec = (Tfb-Tfa)/(Teb-Tea); 
	if (Rmin < 0.0) {
		opserr << "R is negative in ReinforcingSteel::SetMP()\n";
		Rmin = 0.0;
	}
	if (Rmin == 0.0) {
	TQ=1.0;
	Tfch=Tfb;
  } else {
	if (TR <= Rmin) TR=Rmin + 0.01;
	while(notConverge) {
	  if (MPfunc(a)*MPfunc(1-a)>0) 
		a=a/2.0;
	  else
		notConverge=false;
	}
	
	ao= Rmin/TR;
	notConverge=true;
	while(notConverge) {
	  if (MPfunc(ao)*MPfunc(1-a)<0) 
		ao=sqrt(ao);
	  else
		notConverge=false;
    if(ao > 0.999999) notConverge=false;
	}

	notConverge=true;
	if (ao >= 1.0) ao=0.999999;
	while(notConverge) {
    double ao_last=ao;

	  da=0.5*(1-ao);
	  if (da>ao/10.0) da=ao/10.0;
	  ao=ao-2*MPfunc(ao)*da/(MPfunc(ao+da)-MPfunc(ao-da));
#ifdef _WIN32
    if(_fpclass(ao)< 8 || _fpclass(ao)==512) {
      opserr << "Stuck in infinite loop, return error, ReinforcingSteel::SetMP()\n";
      da=da/100.0;
      ao=ao_last;
      ao=ao-2*MPfunc(ao)*da/(MPfunc(ao+da)-MPfunc(ao-da));
      return -1;
    }
#endif
    if(fabs(ao_last-ao)<0.0001) notConverge=false;
	}
	TQ=(TEsec/TEa-ao)/(1-ao);
	if (ao>0.99999999) ao=0.99999999;
	
	double b=pow(1.0-pow(ao,TR),1.0/TR)/ao;
	Tfch=Tfa+TEa/b*(Teb-Tea);
  }
  if(fabs(Teb-Tea)<1.0e-7)
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
  TR=pow(fyp/Esp,1.0/3.0)*20.0*(1.0-10.0*(Teb-Tea)); // dfault
  //TR=pow(-fy_/Es_,1.0/3.0)*20.0*(1.0-10.0*(TeAbsMax-Tea));
}

void
ReinforcingSteel::SetTRn(void)
{
  TR=pow(fyp/Esp,1.0/3.0)*16.0*(1.0-5.0*(Tea-Teb));
  //TR=pow(fy/Es,1.0/3.0)*16.0*(1.0-5.0*(Tea-TeAbsMin));
  
}

void 
ReinforcingSteel::SetTRp1(void)
{
  TR=pow(fyp/Esp,1.0/3.0)*15.0*(1.0-10.0*(Teb-Tea)); // dfault
  //TR=pow(-fy_/Es_,1.0/3.0)*20.0*(1.0-10.0*(TeAbsMax-Tea));
}

void
ReinforcingSteel::SetTRn1(void)
{
  TR=pow(fyp/Esp,1.0/3.0)*12.0*(1.0-5.0*(Tea-Teb));
  //TR=pow(fy/Es,1.0/3.0)*16.0*(1.0-5.0*(Tea-TeAbsMin));
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
			return fsup;
		else {
			if (essp < eshp+0.0002) 
				return (Eshpb-Eypp)*pow(essp-eshpa,2.0)/(2*(eshp+0.0002-eshpa))+ essp*Eypp + fint;
			else
				return fshp+(fsup-fshp)*(1.0-pow((esup-essp)/(esup-eshp),p));
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
			return 0.0;
		else {
      if (essp < eshp+0.0002)
				return (Eshpb-Eypp)*(essp-eshpa)/(eshp+0.0002-eshpa) + Eypp;
      else {
        double temp1 = (fsup-Backbone_fNat(essp));
        double temp2 = (fsup-fshp);
        double temp3 = 1.0-1.0/p;
				return Eshp*pow((fsup-Backbone_fNat(essp))/(fsup-fshp),1.0-1.0/p);
      }
		}
	} 
}

/*****************************************************************************************/
/***********************       Buckled Stress-Strain Relations       *********************/
/*****************************************************************************************/

double 
ReinforcingSteel::Buckled_fact(double ess)
{
  if (LDratio <= 0.0) return 1.0;

  double aveStress;
  double e=ess-Teo_n;

  if(e < -eyp) {
		double eStar=55.0-2.3*sqrt(fyp/Esp*2000)*LDratio;
		if (eStar < 7.0) eStar=7.0;
		eStar= -eStar*eyp;
	  
		double alpha=1.0;  // alpha varies between 0.75 and 1.0 (elasto-plasic -- Linear SH)
		double fStarL=Backbone_f(eStar);
		double fStar=fStarL*alpha*(1.1 - 0.016*sqrt(fyp/Esp*2000)*LDratio);
		if (fStar > -0.2*fyp) fStar= -0.2*fyp;

		if ((e< -eyp) && (e>=eStar)) {
			aveStress = Backbone_f(e)*(1.0-(1.0-fStar/fStarL)*(e+eyp)/(eStar+eyp));
		} else if (e<eStar) {
			aveStress = fStar-0.02*Esp*(e-eStar);
			if (aveStress>-0.2*fyp) aveStress=-0.2*fyp;
		}
		aveStress = (fyp-aveStress)/(fyp-Backbone_f(e));
		return pow(aveStress,1.5);
  } else {
	return 1.0;
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
		
	double fs_buck = alpha_*sqrt(32.0/(e_cross-ess))/(9.42477796076938*LDratio);
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
	    double eon = ea*krev+eb*(1-krev);
	    if (eon > Teo_n) Teo_n=eon;
	    Teb=Teo_n+emin;

	    // set stress dependent curve parameters
	    Tfa=CStress;
      Cfa[0]=CStress;
	    TEa=ReturnSlope(Tea-Teo_n-Temin);
  	  
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
  	  
	    T_FDamage[2]=0.0;
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
	    Tfb=Backbone_f(emin);
	    TEb=(1.0/(1.0/Esp+pr*(1.0/Eshp - 1.0/Esp)));
  	  
	    SetTRn();
	    TEsec = (Tfb-Tfa)/(Teb-Tea);
	    if (TEsec<TEb) TEb=TEsec*0.999;
	    if (TEsec>TEa) TEa=TEsec*1.001;
	    res += SetMP();
  	  
	    T_FDamage[2]=0.0;
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
	  if(Temin<0.0) {
	    TFatDamage-=T_FDamage[0];
	    T_FDamage[0]=damage(TStrain-TeAbsMin,TStress-Cfa[1]);
	    TFatDamage+=T_FDamage[0];
	  }
  }
  return res;
}

// Rule 2: Compresion Envelope Branch
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
	  double eop=ea*(1-krev)+eb*krev;
	  if (eop<Teo_p) Teo_p=eop;
	  Teb=Teo_p+emax;
	  
	  // set stress dependent curve parameters 
	  Tfa=CStress;
    Cfa[1]=CStress;
	  TEa=ReturnSlope(Temax + Teo_p -Tea);

	  Tfb= Backbone_f(emax);
	  TEb= Backbone_E(emax);
	  
	  SetTRp();
	  TEsec = (Tfb-Tfa)/(Teb-Tea);
	  res += SetMP();
		
	  T_FDamage[2]=0.0;
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
	  Tfb= Backbone_f(emax);
	  TEb=1.0/(1.0/Esp+pr*(1.0/Eshp - 1.0/Esp));
	  
	  SetTRp();
	  TEsec = (Tfb-Tfa)/(Teb-Tea);
	  if (TEsec<TEb) TEb=TEsec*0.999;
	  if (TEsec>TEa) TEa=TEsec*1.001;
	  res += SetMP();
	  
	  T_FDamage[2]=0.0;
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
	  if(Temax>0.0) {
	    TFatDamage-=T_FDamage[1];
	    T_FDamage[1]=damage(TeAbsMax-TStrain,Cfa[0]-TStress);
	    TFatDamage+=T_FDamage[1];
	  }
  }
  return res;
}

// Rule 3: Unloading Reversal Branch
int
ReinforcingSteel::Rule3(int res)
{
  if (TStrain-CStrain > 0.0) {	
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
	
	Tfb= Backbone_f(Teb-Teo_p);
	TEb= Backbone_E(Teb-Teo_p);

	SetTRp();
	TEsec = (Tfb-Tfa)/(Teb-Tea);
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
	
	T_FDamage[3]=0.0;
	TBranchNum=5;
	Rule5(res);
  } else {
	  if (TStrain - Teb <= ZeroTol) {
	    T_FDamage[1]=T_FDamage[2];
	    TBranchNum=2;
	    Rule2(res);
	  } else {
      TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage -=T_FDamage[2];
	    T_FDamage[2]=damage(TeAbsMax-TStrain,Tfa-TStress);
	    TFatDamage +=T_FDamage[2];
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

	Tfb= Backbone_f(Teb-Teo_n);
	TEb= Backbone_E(Teb-Teo_n);
	
	SetTRn();
	TEsec = (Tfb-Tfa)/(Teb-Tea); 
	if (TEsec<TEb) TEb=TEsec*0.999;
	if (TEsec>TEa) TEa=TEsec*1.001;
	res += SetMP();
	
	T_FDamage[3]=0.0;
	TBranchNum=6;
	Rule6(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
	    T_FDamage[0]=T_FDamage[2];
	    TBranchNum=1;
	    Rule1(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=T_FDamage[2];
	    T_FDamage[2]=damage(TStrain-TeAbsMin,TStress-Tfa);
	    TFatDamage+=T_FDamage[2];
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
	
	Tfa   = Backbone_f(Tea-Teo_p);
	TEa   = CEa[2];
	
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
  
	T_FDamage[4]=0.0;
	TBranchNum=7;
	Rule7(res);	
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
	    TFatDamage-=T_FDamage[3];
	    TFatDamage+=damage(Teb-Tea,Tfb-Tfa);
	    TBranchNum=1;
	    Rule1(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=T_FDamage[3];
	    T_FDamage[3]=damage(TStrain-Tea,TStress-Tfa);
	    TFatDamage+=T_FDamage[3];
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
	
	Tfa   = Backbone_f(Tea-Teo_n);
	TEa   = CEa[2];

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
	
	T_FDamage[4]=0.0;
	TBranchNum=8;
	Rule8(res);
  } else {
	  if (TStrain - Teb <= ZeroTol) {
	    TFatDamage-=T_FDamage[3];
	    TFatDamage+=damage(Tea-Teb,Tfa-Tfb);
	    TBranchNum=2;
	    Rule2(res);
	  } else {
	    TStress  = MP_f(TStrain);
	    TTangent = MP_E(TStrain);
	    TFatDamage-=T_FDamage[3];
	    T_FDamage[3]=damage(Tea-TStrain,Tfa-TStress);
	    TFatDamage+=T_FDamage[3];
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

	T_FDamage[5]=0.0;
	TBranchNum=9;
	Rule9(res);	
  } else {
	  if (TStrain - Teb <= ZeroTol) {
	    TFatDamage-=T_FDamage[4];
	    TFatDamage+=damage(Tea-Teb,Tfa-Tfb);

	    Tea   = Ceb[3]*(Tea-Cea[3])/(Ceb[3]-Cea[3]) + Cea[2]*(Ceb[3]-Tea)/(Ceb[3]-Cea[3]);
	    Teb   = Ceb[2];
	    Tfa   = Backbone_f(Tea-Teo_p);
	    TEa   = CEa[2];
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
	    TFatDamage-=T_FDamage[4];
	    T_FDamage[4]=damage(Tea-TStrain,Tfa-TStress);
	    TFatDamage+=T_FDamage[4];
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

	T_FDamage[5]=0.0;
	TBranchNum=10;
	Rule10(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
	    TFatDamage-=T_FDamage[4];
	    TFatDamage+=damage(Teb-Tea,Tfb-Tfa);
  	  
	    Tea   = Ceb[3]*(Tea-Cea[3])/(Ceb[3]-Cea[3]) + Cea[2]*(Ceb[3]-Tea)/(Ceb[3]-Cea[3]);
	    Teb   = Ceb[2];
	    Tfa   = Backbone_f(Tea-Teo_n);
	    TEa   = CEa[2];
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
	    TFatDamage-=T_FDamage[4];
	    T_FDamage[4]=damage(TStrain-Tea,TStress-Tfa);
	    TFatDamage+=T_FDamage[4];
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
	T_FDamage[TBranchMem]=0.0;
	Rule11(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage -=T_FDamage[TBranchMem];
	    TFatDamage +=damage(Teb-Tea,Tfb-Tfa);
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
	    TFatDamage -=T_FDamage[TBranchMem];
	    T_FDamage[TBranchMem]=damage(TStrain-Tea,TStress-Tfa);
	    TFatDamage +=T_FDamage[TBranchMem];
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
	T_FDamage[TBranchMem]=0.0;
	Rule12(res);
  } else {
	  if (TStrain - Teb <= ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=T_FDamage[TBranchMem];
	    TFatDamage+=damage(Tea-Teb,Tfa-Tfb);

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
	    TFatDamage  -=T_FDamage[TBranchMem];
	    T_FDamage[TBranchMem]=damage(Tea-TStrain,Tfa-TStress);
	    TFatDamage  +=T_FDamage[TBranchMem];
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
	T_FDamage[TBranchMem]=0.0;
	Rule9(res);	
  } else {
	  if (TStrain - Teb <= ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=T_FDamage[TBranchMem-2];
	    TFatDamage+=damage(Tea-Teb,Tfa-Tfb);
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
	    TFatDamage-=T_FDamage[TBranchMem];
	    T_FDamage[TBranchMem]=damage(Tea-TStrain,Tfa-TStress);
	    TFatDamage+=T_FDamage[TBranchMem];
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
	T_FDamage[TBranchMem]=0.0;
	Rule10(res);
  } else {
	  if (TStrain - Teb >= -ZeroTol) {
		  TBranchMem = (TBranchNum+1)/2;
	    TFatDamage-=T_FDamage[TBranchMem-2];
	    TFatDamage+=damage(Teb-Tea,Tfb-Tfa);
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
	    TFatDamage-=T_FDamage[TBranchMem];
	    T_FDamage[TBranchMem]=damage(TStrain-Tea,TStress-Tfa);
	    TFatDamage+=T_FDamage[TBranchMem];
	  }
  }
  return res;
}



/*****************************************************************************************/
/********        Strength and stiffness degradation including buckling         ***********/
/*****************************************************************************************/
double
ReinforcingSteel::damage(double ehalf, double stressAmp){
  double ehalfPlastic = fabs(ehalf)-fabs(stressAmp/Esp);
  if(ehalfPlastic>0.0)
    return pow(ehalfPlastic/Fat1,Fat2);
  else
    return 0.0;
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
  return Esp*(1-2.0*(dea));
}
