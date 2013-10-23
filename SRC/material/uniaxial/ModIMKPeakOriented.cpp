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
                                                                        
// $Revision: 1.7 $
// $Date: 2009/03/05 00:52:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ModIMKPeakOriented.cpp,v $
                                                                        
// Written: Dimitrios G. Lignos, PhD, Assistant Professor, McGill University 
// Created: February, 2011
// Revision: A
//
// Description: This file contains the class implementation for ModIMKPeakOriented Model

#include <math.h>

#include <elementAPI.h>
#include <ModIMKPeakOriented.h>
#include <Vector.h>
#include <Channel.h>
#include <MaterialResponse.h>

#include <OPS_Globals.h>

static int numModIMKPeakOrientedMaterials = 0;

void *
OPS_ModIMKPeakOriented()
{
  if (numModIMKPeakOrientedMaterials == 0) {
    numModIMKPeakOrientedMaterials++;
    opserr << "Modified Ibarra-Medina-Krawinkler Model with Peak-Oriented Hysteretic Response\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[1];
  double dData[23];
  int numData = 1;
  // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  ModIMKPeakOriented tag" << endln;
    return 0;
  }
  
  numData = 23;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial ModIMKPeakOriented tag? Ke?, alfaPos?, alfaNeg?, My_pos?, My_neg?"; 
    opserr << "Ls?, Ld?, La?, Lk?, Cs?, Cd?, Ca?, Ck?, thetaPpos?, thetaPneg?, thetaPCpos?, thetaPCneg? "; 
    opserr << "ResfacPos?, ResfacNeg?, fracDispPos?, fracDispNeg?,DPos?, DNeg?";
    
    return 0;	
  }
  
  // Parsing was successful, allocate the material with zero index
  theMaterial = new ModIMKPeakOriented(iData[0], 
				       dData[0], dData[1], dData[2], dData[3],
				       dData[4], dData[5], dData[6], dData[7],
				       dData[8], dData[9], dData[10], dData[11],
				       dData[12], dData[13], dData[14], dData[15],
				       dData[16], dData[17], dData[18], dData[19],
				       dData[20], dData[21], dData[22]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ModIMKPeakOriented Material\n";
    return 0;
  }
  
  return theMaterial;
}


ModIMKPeakOriented::ModIMKPeakOriented(int tag, double ke, double alfaPos, double alfaNeg, double my_pos, double my_neg, 
				       double ls, double ld, double la, double lk, double cs, double cd, double ca, double ck,
				       double thetaPpos, double thetaPneg, double thetaPCpos, double thetaPCneg,
				       double resfacPos, double resfacNeg, double fracDispPos, double fracDispNeg,
				       double dPos, double dNeg)
 :UniaxialMaterial(tag,MAT_TAG_ModIMKPeakOriented), Ke(ke), AlfaPos(alfaPos), AlfaNeg(alfaNeg), My_pos(my_pos), My_neg(my_neg), 
  Ls(ls), Ld(ld), La(la), Lk(lk), Cs(cs), Cd(cd), Ca(ca), Ck(ck),ThetaPpos(thetaPpos), ThetaPneg(thetaPneg), 
  ThetaPCpos(thetaPCpos), ThetaPCneg(thetaPCneg), ResfacPos(resfacPos), ResfacNeg(resfacNeg), 
  FracDispPos(fracDispPos), FracDispNeg(fracDispNeg),DPos(dPos), DNeg(dNeg)
{
  // Initialize Variables by Calling revertToStart function
  this->revertToStart();
  
}

ModIMKPeakOriented::ModIMKPeakOriented()
:UniaxialMaterial(0,MAT_TAG_ModIMKPeakOriented),
 Ke(0.0), AlfaPos(0.0), AlfaNeg(0.0), My_pos(0.0), My_neg(0.0), Ls(0.0), Ld(0.0), La(0.0), Lk(0.0), 
 Cs(0.0), Cd(0.0), Ca(0.0), Ck(0.0), ThetaPpos(0.0), ThetaPneg(0.0), 
 ThetaPCpos(0.0), ThetaPCneg(0.0), ResfacPos(0.0), ResfacNeg(0.0), 
 FracDispPos(0.0), FracDispNeg(0.0), DPos(0.0), DNeg(0.0)
{
  // Initialize Variables by Calling revertToStart function
  this->revertToStart();
}

ModIMKPeakOriented::~ModIMKPeakOriented()
{
  // does nothing
}

int 
ModIMKPeakOriented::setTrialStrain(double strain, double strainRate)
{
  //all variables to the last commit state
  
  this->revertToLastCommit();
  
  double f,d,deltaD;
  
  d  = strain;
  
  // Determine change in deformation from last converged state
  deltaD = d - dP;
  
  // added by DL 04/22/2013
  
  flagdeg = 0;
  
  if (fabs(deltaD) < 1.0e-18 && strain != 0.0) {
    return 0;
  }
  
  
  // Initialize parameters in the first cycle
  if (kon==0) {
    
    // strain hardening range
    ekhardPos = Ke * AlfaPos;
    ekhardNeg = Ke * AlfaNeg;
    
    // Reference Energy Capacity based on Lignos (2008)
    Enrgts = Ls * My_pos;
    Enrgtk = 2.0 * Lk * My_pos;
    Enrgta = La * My_pos;
    Enrgtd = Ld * My_pos;
    
    Enrgtot = 0.0;
    Enrgc = 0.0;
    
    dmax = My_pos/Ke;
    dmin = My_neg/Ke;
    
    ekP = Ke;
    ekunload = Ke;
    ekexcurs = Ke;
    
    fyPos = My_pos;
    fyNeg = My_neg;
    
    cpPos = My_pos/Ke + ThetaPpos;
    cpNeg = My_neg/Ke + (-ThetaPneg);
    
    // Capping Point
    fPeakPos = My_pos + ekhardPos * (ThetaPpos);
    fPeakNeg = My_neg + ekhardNeg * (-ThetaPneg);
    
    // Capping Slope
    capSlopePos = -fPeakPos/(ThetaPCpos * Ke);
    capSlopeNeg =  fPeakNeg/(ThetaPCneg * Ke);
    
    //Reference Point for Cyclic Deterioration
    
    fCapRefPos = -capSlopePos * Ke * (ThetaPpos + My_pos/Ke) + fPeakPos;
    fCapRefNeg = -capSlopeNeg * Ke * (-ThetaPneg + My_neg/Ke) + fPeakNeg;
    
    fmin = My_neg;
    fmax = My_pos;
    
    flagdeg = 0;
    flagStop = 0;     // if deformation exceeds theta_u for +/- loading direction f = 0;
    
    dlstNeg = My_neg/Ke;
    dlstPos = My_pos/Ke;
    
    flstNeg = My_neg;
    flstPos = My_pos;
    Tangent = Ke;
    
    ek = 0.0;
    sp = 0.0;
    fP = 0.0;
    sn = 0.0;
    f = 0.0;
    
    Unl = 1;
    
    RSE = 0.0;
    
    if(deltaD >= 0.0){
      kon = 1;
    } else{
      kon = 2;
    }
  }
  
  
  // After the first step		
  if(flagStop == 0){
    
    if (deltaD >= 0.0){
      
      if (kon==2){	
	kon = 1;	
	Unl = 0;
	
	RSE = 0.5 * fP * fP / ekunload;
	double a2 = Enrgtk-(Enrgtot-RSE);
	double betak;
	
	if(a2 < 0.0 && Enrgtk != 0.0){
	  flagStop = 1;
	}
	
	if(Lk != 0.0){
	  betak = pow((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE)),Ck);
	  ekunload = ekexcurs*(1.0-betak);
          
	  if(ekunload <= ekhardNeg){
	    flagStop = 1;
	  }           	
	}
	
	//  determination of sn according to the hysteresis status        
	
	if(fP < 0.0){
	  
	  double tst = dP-fP/ekunload;
	  
	  if(fabs(dmax-My_pos/Ke)>=1.e-10 && fabs(tst) <= 1.e-10){
	    
	    sn = 1.e-9;
	    
	  } else {
	    sn = dP-fP/ekunload;
	  }
	}
	if(fabs(dmin - dP) <= 1.e-10){
	  sp = sn + 1.e-10;
	}
      } // if kon = 2
      
      //	LOADING
      
      //     Push envelope
      
      if (d>=dmax){
	envelPosCap2(fyPos,AlfaPos,capSlopePos,cpPos,d, f, ek, Ke, My_pos, ResfacPos,FracDispPos, flagStop);
	dmax = d;
	fmax = f;
	dlstPos = dmax+1.e-10;
	flstPos = f;
	
      } else if (fabs(sn)>1.e-10){
	
	envelPosCap2 (fyPos,AlfaPos,capSlopePos,cpPos,dmax,fmax, ek, Ke, My_pos, ResfacPos,FracDispPos, flagStop);
	
	if ( d <= sn){
	  ek = ekunload;
	  f  = fP+ek*deltaD;
	  
	  if(Unl == 0 && fabs(ek-ekP) > 1.e-10 && dP != dmin){
	    dlstNeg = dP;
	    flstNeg = fP;
	  }
	  //	If d>sn, 	
	} else {
	  ek = fmax/(dmax-sn);
	  double f2 = ek*(d-sn);
	  
	  if(dlstPos > sn && dlstPos < dmax){
	    
	    double ekc = flstPos/(dlstPos-sn);
	    
	    if(ekc>ek && flstPos < fmax){
	      
	      if(d < dlstPos){
		ek = flstPos/(dlstPos-sn);
		f2 = ek*(d-sn);
		
	      } else {
		
		ek = (fmax-flstPos)/(dmax-dlstPos);
		f2 = flstPos+ek*(d-dlstPos);
	      }
	    }
	  }
	  
	  double f1 = fP+ekunload*deltaD;
	  
	  //f = min(f1,f2);
	  if (f1 < f2) {
	    f = f1;
	    
	  } else {
	    
	    f = f2;
	  }
	  
	  if (fabs(f-f1)<1.e-10) {
	    ek=ekunload;
	  }
					
	} 
	
      } else {							
	if (d > 0.0) {
	  envelPosCap2(fyPos,AlfaPos,capSlopePos,cpPos,d, f, ek, Ke, My_pos, ResfacPos,FracDispPos, flagStop);	
	} else {	
	  envelNegCap2(fyNeg,AlfaNeg,capSlopeNeg,cpNeg, d, f, ek, Ke, My_neg, ResfacNeg,-FracDispNeg, flagStop);	
	}
				
      }
      
      // if deltaD < 0
    } else {
      if (kon == 1) {
	kon = 2;
	Unl = 0;
	
	RSE = 0.5 * fP * fP/ekunload;						
	double a2 = Enrgtk-(Enrgtot-RSE);
	
	
	if(Lk!= 0.0 ){
	  double betak = pow((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE)),Ck);
	  ekunload = ekexcurs*(1.0-betak);
	  
	  if(ekunload<= ekhardPos){
	    flagStop=1;
	  }           
	}
	
	//		Determination of sn according to the hysteresis status
	
	if(fP > 0.0) {
	  
	  double tst = dP-fP/ekunload;
	  if(fabs(dmin-My_neg/Ke)>=1.e-10 && fabs(tst) <= 1.e-10){
	    sp = 1.e-9;
	  } else {
	    sp = dP-fP/ekunload;
	  }
	}
	if(fabs(dmax-dP) <=1.e-10){
	  sn = sp-1.e-10;
	}
	
      } // kon = 1
      
      //	UNLOADING
      
      //  Push envelope						
      if (d <= dmin){
	envelNegCap2 (fyNeg,AlfaNeg,capSlopeNeg,cpNeg,d, f, ek, Ke, My_neg, ResfacNeg,-FracDispNeg, flagStop);
	dmin = d;
	fmin = f;
	dlstNeg = dmin - 1.e-10;
	flstNeg = f;
	
      } else if (fabs(sp) > 1.e-10){
	envelNegCap2 (fyNeg,AlfaNeg,capSlopeNeg,cpNeg, dmin, fmin, ek, Ke, My_neg, ResfacNeg,-FracDispNeg, flagStop);
	
	if ( d >= sp ){
	  
	  ek = ekunload;
	  f = fP+ek*deltaD;
	  
	  if(Unl==0 && fabs(ek-ekP) > 1.e-10 && dP != dmax){
	    dlstPos = dP;	
	    flstPos = fP;
	  }
	} else {
	  
	  ek = fmin/(dmin-sp);
	  double f2 = ek * (d-sp);
	  
	  if(dlstNeg < sp && dlstNeg > dmin){
	    
	    double ekc = flstNeg/(dlstNeg-sp);
	    
	    if(ekc > ek && flstNeg > fmin){
	      
	      if(d > dlstNeg) {  
		
		ek = flstNeg/(dlstNeg-sp);
		f2 = ek*(d-sp);
		
	      } else {
		
		ek = (fmin-flstNeg)/(dmin-dlstNeg);								
		f2 = flstNeg+ek*(d-dlstNeg);							
	      }
	    }					
	  }
	  
	  double f1 = fP+ekunload*deltaD;
	  // f=max(f1,f2);
	  if (f1>f2) {
	    f=f1;
	  } else {
	    f=f2;
	  }
	  
	  if (fabs(f-f1) < 1.e-10){ 
	    ek=ekunload;
	  }
	}
	
	
      } else {
	
	if (d > 0.0){
	  envelPosCap2 (fyPos,AlfaPos,capSlopePos,cpPos,d, f, ek, Ke, My_pos, ResfacPos ,FracDispPos, flagStop);
	} else {
	  envelNegCap2 (fyNeg,AlfaNeg,capSlopeNeg,cpNeg,d, f, ek, Ke, My_neg, ResfacNeg ,-FracDispNeg, flagStop);
	}
      }
      
    }
    
    //Energy Computations
    
    double Enrgi = 0.0;
    
    Enrgi = 0.5 * (f+fP) * deltaD;
    Enrgc = Enrgc + Enrgi;
    Enrgtot = Enrgtot + Enrgi;
    RSE = 0.5*f*f/ekunload;
    
    // Added by DL, 04/22/2013
    
    if((f*fP<0.0)) {
      if(((fP>0.0)&&(dmax>(fyPos/Ke)))||((fP<0.0)&&(dmin<(fyNeg/Ke)))) {
	flagdeg = 1;		
	//interup = 1;
      }
    }	
    
    //	UPDATING OF DETERIORATION PARAMETERS ----------------------------
    
    // Check the remaining capacity of the system
    double betas = 0.0;
    double betaa = 0.0;
    double betac = 0.0;
    
    // no cyclic deterioration
    if(flagdeg == 1){
      
      if(Enrgtot >= Enrgts && Enrgts != 0.0){	    
	betas = 1.0;
      } else if (Enrgtot >= Enrgta && Enrgta != 0.0) {  
	betaa = 1.0;
      } else if(Enrgtot >= Enrgtd && Enrgtd == 0.0) {  
	betac = 1.0;
      }
      
      // cyclic deterioration
    } else {
      
      if(Ls != 0.0){
	betas = pow(Enrgc/(Enrgts-Enrgtot),Cs);  	
      } else {
	betas = 0.0;
      }
      
      if(La != 0.0){
	betaa = pow(Enrgc/(Enrgta-Enrgtot),Ca);
      } else {
	betaa = 0.0;
      }
      if(Ld != 0.0){
	betac = pow(Enrgc/(Enrgtd-Enrgtot),Cd);
      } else {
	betac = 0.0;
      }
      if(fabs(betas) >= 1.0) {
	betas = 1.0;
      }
      if(fabs(betaa) >= 1.0){ 
	betaa = 1.0;
      }
      if(fabs(betac) >= 1.0){ 
	betac = 1.0;								
      }
    }
    
    //	Update values for the next half cycle
    ekexcurs = ekunload;
    Enrgc = 0.0;
    //	Deteriorate parameters for the next half cycle
    
    // if(deltaD > 0.0 && f > 0.0 || deltaD < 0.0 && f >=0.0 && d <=sp) {
    if(deltaD < 0.0 ) {			
      fyNeg = fyNeg*(1.0-betas*DNeg);			
      AlfaNeg=AlfaNeg*(1.0-betas*DNeg);			
      fCapRefNeg=fCapRefNeg*(1.0-betac*DNeg);
      dmin = dmin*(1.0+betaa);
      
      double dyNeg = fyNeg/Ke;
      ekhardNeg=AlfaNeg*Ke;
      
      double dCap1Neg=fCapRefNeg/(Ke-capSlopeNeg*Ke);
      double dCap2Neg=(fCapRefNeg+ekhardNeg*dyNeg-fyNeg)/(ekhardNeg-capSlopeNeg*Ke);
      // cpNeg=min(dCap1Neg,dCap2Neg);
      //double fCapNeg = fCapRefNeg + capSlopeNeg*Ke*cpNeg;
      
      if (dCap1Neg < dCap2Neg){
	cpNeg = dCap1Neg;
      } else {
	cpNeg = dCap2Neg;
      }
      
    } else { //else if (deltaD < 0.0 && f < 0.0 || deltaD > 0.0 && f <=0.0 && d >= sn) {
      
      fyPos = fyPos*(1.0-betas*DPos);
      AlfaPos=AlfaPos*(1.0-betas*DPos);
      fCapRefPos=fCapRefPos*(1.0-betac*DPos);
      dmax = dmax*(1.0+betaa);	
      
      double dyPos = fyPos/Ke;
      ekhardPos=AlfaPos*Ke;
      
      double dCap1Pos=fCapRefPos/(Ke-capSlopePos*Ke);
      double dCap2Pos=(fCapRefPos+ekhardPos*dyPos-fyPos)/(ekhardPos-capSlopePos*Ke);
      
      //cpPos=max(dCap1Pos,dCap2Pos);
      
      if (dCap1Pos>dCap2Pos) {
	cpPos = dCap1Pos;
      } else {
	cpPos = dCap2Pos;
      }
      //double fCapPos = fCapRefPos + capSlopePos * Ke * cpPos;
    } 
    
  } else {
    f = 0.0;
    ek = 0.0;
  }
  
  dP = d;
  fP = f;
  ekP = ek;
  Tangent = ek;
  
  return 0;
}

double 
ModIMKPeakOriented::getStress(void)
{

  return  (fP);
}

double 
ModIMKPeakOriented::getTangent(void)
{
  return (Tangent);
}

double 
ModIMKPeakOriented::getInitialTangent(void)
{
  return (Ke);
}

double 
ModIMKPeakOriented::getDampTangent(void)
{
  double DTangent = ekP;
  return DTangent;	
}


double 
ModIMKPeakOriented::getStrain(void)
{
  return dP;
}

double 
ModIMKPeakOriented::getStrainRate(void)
{
  return 0;
}

int 
ModIMKPeakOriented::commitState(void)
{
  Cstrain = dP;
  Cstress = fP;
  Ctangent = Tangent;
  
  CdP = dP;
  CfP = fP;
  Cek = ek;
  
  CUnl = Unl;
  Ckon = kon;
  CflagStop = flagStop;
  
  Cdmax = dmax;
  Cdmin = dmin;
  Cfmax = fmax;
  Cfmin = fmin;
  
  CfyPos = fyPos;
  CfyNeg = fyNeg;
  
  CdlstPos = dlstPos;
  CdlstNeg = dlstNeg;
  
  CflstPos = flstPos;
  CflstNeg = flstNeg;
  
  Csn = sn;
  Csp = sp;
  
  CEnrgc = Enrgc;
  CEnrgtot = Enrgtot;
  
  CEnrgts = Enrgts;
  CEnrgtd = Enrgtd;
  CEnrgtk = Enrgtk;
  CEnrgta = Enrgta;
  
  CfPeakPos = fPeakPos;
  CfPeakNeg = fPeakNeg;
  
  CcapSlopePos = capSlopePos;
  CcapSlopeNeg = capSlopeNeg;
  
  CfCapRefPos = fCapRefPos;
  CfCapRefNeg = fCapRefNeg;
  
  Cekunload = ekunload;
  
  CcpNeg = cpNeg;
  CcpPos = cpPos;
  
  CekhardPos = ekhardPos;
  CekhardNeg = ekhardNeg;
  Cekexcurs = ekexcurs;
  CekP = ekP;
  Cflagdeg = flagdeg;
  CRSE=RSE; 
  
  return 0;
}

int 
ModIMKPeakOriented::revertToLastCommit(void)
{
  dP = Cstrain;
  fP = Cstress;
  Tangent = Ctangent;
  
  dP = CdP;
  fP = CfP;
  ek = Cek;
  
  Unl = CUnl;
  kon = Ckon;
  flagStop = CflagStop;
  
  dmax = Cdmax;
  dmin = Cdmin;
  fmax = Cfmax;
  fmin = Cfmin;
  
  fyPos = CfyPos;
  fyNeg = CfyNeg;
  
  dlstPos = CdlstPos;
  dlstNeg = CdlstNeg;
  
  flstPos = CflstPos;
  flstNeg = CflstNeg;
  
  sn = Csn;
  sp = Csp;
  
  // Energy Variables
  Enrgc = CEnrgc;
  Enrgtot = CEnrgtot;
  
  Enrgts = CEnrgts;
  Enrgtd = CEnrgtd;
  Enrgtk = CEnrgtk;
  Enrgta = CEnrgta;
  
  fPeakPos = CfPeakPos;
  fPeakNeg = CfPeakNeg;
  
  capSlopePos = CcapSlopePos;
  capSlopeNeg = CcapSlopeNeg;
  
  fCapRefPos = CfCapRefPos;
  fCapRefNeg = CfCapRefNeg;
  
  ekunload = Cekunload;
  
  cpNeg = CcpNeg;
  cpPos = CcpPos;
  
  ekhardPos = CekhardPos;
  ekhardNeg = CekhardNeg;
  ekexcurs = Cekexcurs;
  ekP = CekP;
  flagdeg = Cflagdeg;
  RSE=CRSE; 
  
  return 0;
}

int 
ModIMKPeakOriented::revertToStart(void)
{
	
  // Initialize state variables	
  dP=CdP=0.0;
  fP=CfP=0.0;
  ek=Cek=0.0;
  RSE=CRSE=0.0; 
  
  Unl=CUnl=1;
  kon=Ckon=0;
  flagStop=CflagStop=0;
  
  dmax = My_pos/Ke;
  dmin = My_neg/Ke;
  fmax = My_pos;
  fmin = My_neg;
  
  fyPos = My_pos;
  fyNeg = My_neg;
  
  dlstPos = My_pos/Ke;
  dlstNeg = My_neg/Ke;
  
  flstPos = My_pos;
  flstNeg = My_neg;
  
  sn=Csn=0.0;
  sp=Csp=0.0;
  
  // Energy Variables
  Enrgc=CEnrgc=0.0;
  Enrgtot=CEnrgtot=0.0;
  
  Enrgts = Ls * My_pos;
  Enrgtd = Ld * My_pos;
  Enrgtk = 2.0 * Lk * My_pos;
  Enrgta = La * My_pos;
  
  fPeakPos = My_pos + ekhardPos * ThetaPpos;
  fPeakNeg = My_neg + ekhardNeg * (-ThetaPneg);
  
  capSlopePos = -fPeakPos/(ThetaPCpos*Ke);
  capSlopeNeg = fPeakNeg/(ThetaPCneg*Ke);
  
  fCapRefPos = -capSlopePos * Ke * (ThetaPpos + My_pos/Ke) + fPeakPos;
  fCapRefNeg = -capSlopeNeg * Ke * (-ThetaPneg + My_neg/Ke) + fPeakNeg;
  
  ekunload = Ke;
  
  cpPos = My_pos/Ke + ThetaPpos;
  cpNeg = My_neg/Ke + (-ThetaPneg);
  
  ekhardPos = Ke * AlfaPos;
  ekhardNeg = Ke * AlfaNeg;
  ekexcurs = Ke;
  ekP = 0.0;
  
  flagdeg = 0;
  
  Cstrain=0.0;
  Cstress = 0.0;
  Tangent=Ctangent=ek;
  
  Cdmax = My_pos/Ke;
  Cdmin = My_neg/Ke;
  Cfmax = My_pos;
  Cfmin = My_neg;
  
  CfyPos = My_pos;
  CfyNeg = My_neg;
  
  CdlstPos = My_pos/Ke;
  CdlstNeg = My_neg/Ke;
  
  CflstPos = My_pos;
  CflstNeg = My_neg;
  
  CEnrgts = Ls * My_pos;
  CEnrgtd = Ld * My_pos;
  CEnrgtk = 2.0 * Lk * My_pos;
  CEnrgta = La * My_pos;
  
  CfPeakPos = My_pos + ekhardPos * ThetaPpos;
  CfPeakNeg = My_neg + ekhardNeg * (-ThetaPneg);
  
  CcapSlopePos = -fPeakPos/(ThetaPCpos*Ke);
  CcapSlopeNeg = fPeakNeg/(ThetaPCneg*Ke);
  
  CfCapRefPos = -capSlopePos * Ke * (ThetaPpos + My_pos/Ke) + fPeakPos;
  CfCapRefNeg = -capSlopeNeg * Ke * (-ThetaPneg + My_neg/Ke) + fPeakNeg;
  
  Cekunload = Ke;
  
  CcpPos = My_pos/Ke + ThetaPpos;
  CcpNeg = My_neg/Ke + (-ThetaPneg);
  
  CekhardPos = Ke * AlfaPos;
  CekhardNeg = Ke * AlfaNeg;
  Cekexcurs = Ke;
  CekP = 0.0;
  
  Cflagdeg = 0;
  
  return 0;
}

UniaxialMaterial *
ModIMKPeakOriented::getCopy(void)
{
  ModIMKPeakOriented *theCopy = new ModIMKPeakOriented(this->getTag(), Ke, AlfaPos, AlfaNeg, My_pos, My_neg, Ls, Ld, La, Lk, 
						       Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg,ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg,FracDispPos, FracDispNeg,DPos, DNeg);
    // Converged state variables
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Ctangent = Ctangent;
  
  theCopy->CdP = CdP;
  theCopy->CfP = CfP;
  theCopy->Cek = ek;
  
  theCopy->CUnl = CUnl;
  theCopy->Ckon = Ckon;
  theCopy->CflagStop = CflagStop;
  theCopy->Cflagdeg = Cflagdeg;
  
  theCopy->Cdmax = Cdmax;
  theCopy->Cdmin = Cdmin;
  theCopy->Cfmax = Cfmax;
  theCopy->Cfmin = Cfmin;
  
  theCopy->CfyPos = CfyPos;
  theCopy->CfyNeg = CfyNeg;
  
  theCopy->CdlstPos = CdlstPos;
  theCopy->CdlstNeg = CdlstNeg;
  theCopy->CflstPos = CflstPos;
  theCopy->CflstNeg = CflstNeg;
  
  theCopy->Csn= Csn;
  theCopy->Csp= Csp;
  
  // Energy Variables
  theCopy->CEnrgc = CEnrgc;
  theCopy->CEnrgtot = CEnrgtot;
  
  theCopy->CEnrgts = CEnrgts;
  theCopy->CEnrgtd = CEnrgtd;
  theCopy->CEnrgtk = CEnrgtk;
  theCopy->CEnrgta = CEnrgta;
  
  theCopy->CfPeakPos = CfPeakPos;
  theCopy->CfPeakNeg = CfPeakNeg;
  
  theCopy->CcapSlopePos = CcapSlopePos;
  theCopy->CcapSlopeNeg = CcapSlopeNeg;
  
  theCopy->CfCapRefPos = CfCapRefPos;
  theCopy->CfCapRefNeg = CfCapRefNeg;
  
  theCopy->Cekunload = Cekunload;
  
  theCopy->CcpNeg = CcpNeg;
  theCopy->CcpPos = CcpPos;
  
  theCopy->CekhardPos = CekhardPos;
  theCopy->CekhardNeg = CekhardNeg;
  theCopy->Cekexcurs = Cekexcurs;
  theCopy->CekP = CekP;
  
  // Trial state variables
  theCopy->dP = dP;
  theCopy->fP = fP;
  theCopy->ek = ek;
  
  theCopy->Unl = Unl;
  theCopy->kon = kon;
  theCopy->flagStop = flagStop;
  theCopy->Cflagdeg = Cflagdeg;
  
  theCopy->dmax = dmax;
  theCopy->dmin = dmin;
  theCopy->fmax = fmax;
  theCopy->fmin = fmin;
  
  theCopy->fyPos = fyPos;
  theCopy->fyNeg = fyNeg;
  
  theCopy->dlstPos = dlstPos;
  theCopy->dlstNeg = dlstNeg;
  theCopy->flstPos = flstPos;
  theCopy->flstNeg = flstNeg;
  
  theCopy->sn= sn;
  theCopy->sp= sp;
  
  // Energy Variables
  theCopy->Enrgc = Enrgc;
  theCopy->Enrgtot = Enrgtot;
  
  theCopy->Enrgts = Enrgts;
  theCopy->Enrgtd = Enrgtd;
  theCopy->Enrgtk = Enrgtk;
  theCopy->Enrgta = Enrgta;
  
  theCopy->fPeakPos = fPeakPos;
  theCopy->fPeakNeg = fPeakNeg;	
  
  theCopy->capSlopePos = capSlopePos;
  theCopy->capSlopeNeg = capSlopeNeg;
  
  theCopy->fCapRefPos = fCapRefPos;
  theCopy->fCapRefNeg = fCapRefNeg;
  
  theCopy->ekunload = ekunload;
  
  theCopy->cpNeg = cpNeg;
  theCopy->cpPos = cpPos;
  
  theCopy->ekhardPos = ekhardPos;
  theCopy->ekhardNeg = ekhardNeg;
  theCopy->ekexcurs = ekexcurs;
  theCopy->ekP = ekP;
  
  theCopy->RSE=RSE;
  theCopy->CRSE=CRSE; 
  
  return theCopy;
}

int 
ModIMKPeakOriented::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(67);
  data(0) = this->getTag();
  
  // Material properties
  data(1) = Ke;
  data(2) = AlfaPos;
  data(3) = AlfaNeg;
  data(4) = My_pos;
  data(5) = My_neg;
  data(6) = Ls;
  data(7) = Ld;
  data(8) = La;
  data(9) = Lk;
  data(10) = Cs;
  data(11) = Cd;
  data(12) = Ca;
  data(13) = Ck;
  data(14) = ThetaPpos;
  data(15) = ThetaPneg;
  data(16) = ThetaPCpos;
  data(17) = ThetaPCneg;
  data(18) = ResfacPos;
  data(19) = ResfacNeg;
  data(20) = FracDispPos;
  data(21) = FracDispNeg;
  data(22) = DPos;
  data(23) = DNeg;
  
  // State variables from last converged state
  data(24) = Cstrain;
  data(25) = Cstress;
  data(26) = Ctangent;
  
  data(27) = CdP;
  data(28) = CfP;
  data(29) = Cek;
  
  //  Variables 
  data(30) = CUnl;
  data(31) = Ckon;
  data(32) = CflagStop;
  data(33) = Cdmax;
  data(34) = Cdmin;
  data(35) = Cfmax;
  data(36) = Cfmin;
  
  data(37) = CfyPos;
  data(38) = CfyNeg;
  
  data(39) = CdlstPos;
  data(40) = CdlstNeg;
  
  data(41) = Csn;
  data(42) = Csp;
  
  // Energy Variables
  data(43) = CEnrgc;
  data(44) = CEnrgtot;
  
  data(45) = CEnrgts;
  data(46) = CEnrgtd;
  data(47) = CEnrgtk;
  data(48) = CEnrgta;
  
  data(49) = CfPeakPos;
  data(50) = CfPeakNeg;
  
  data(51) = CcapSlopePos;
  data(52) = CcapSlopeNeg;
  
  data(53) = CfCapRefPos;
  data(54) = CfCapRefNeg;
  
  data(55) = Cekunload;
  
  data(56) = CcpNeg;
  data(57) = CcpPos;
  
  data(58) = CekhardPos;
  data(59) = CekhardNeg;
  data(60) = Cekexcurs;
  data(61) = CekP;
  
  data(62) = CflstPos;
  data(63) = CflstNeg;
  data(64) = Cflagdeg;
  
  data(65) = RSE;
  data(66) = CRSE;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ModIMKPeakOriented::sendSelf() - failed to send data\n";

  return res;
}

int 
ModIMKPeakOriented::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(67);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ModIMKPeakOriented::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    // Material properties
    Ke = data(1);
    AlfaPos = data(2);
    AlfaNeg = data(3);
    My_pos = data(4);
    My_neg = data(5);
    Ls = data(6);
    Ld = data(7);
    La = data(8);
    Lk = data(9);
    Cs = data(10);
    Cd = data(11);
    Ca = data(12);
    Ck = data(13);
    ThetaPpos = data(14);
    ThetaPneg = data(15);
    ThetaPCpos = data(16);
    ThetaPCneg = data(17);
    ResfacPos = data(18);
    ResfacNeg = data(19);
    FracDispPos = data(20);
    FracDispNeg = data(21);
    DPos = data(22);
    DNeg = data(23);
    // State variables from last converged state 	
    Cstrain = data(24);
    Cstress = data(25);
    Ctangent = data(26);
    
    CdP = data(27);
    CfP = data(28);
    Cek = data(29);
    
    CUnl = data(30);
    Ckon = data(31);
    CflagStop = data(32);
    
    Cdmax = data(33);
    Cdmin = data(34);
    Cfmax = data(35);
    Cfmin = data(36);
    
    CfyPos = data(37);
    CfyNeg = data(38);
    
    CdlstPos = data(39);
    CdlstNeg = data(40);
    
    Csn = data(41);
    Csp = data(42);
    
    // Energy Variables
    CEnrgc = data(43);
    CEnrgtot = data(44);
    
    CEnrgts = data(45);
    CEnrgtd = data(46);
    CEnrgtk = data(47);
    CEnrgta = data(48);
    
    CfPeakPos = data(49);
    CfPeakNeg = data(50);
    
    CcapSlopePos = data(51);
    CcapSlopeNeg = data(52);
    
    CfCapRefPos = data(53);
    CfCapRefNeg = data(54);
    
    Cekunload = data(55);
    
    CcpNeg = data(56);
    CcpPos = data(57);
    
    CekhardPos = data(58);
    CekhardNeg = data(59);
    Cekexcurs = data(60);
    CekP = data(61);
    
    CflstPos = data(62);
    CflstNeg = data(63);
    Cflagdeg = data(64);
    RSE = data(65);
    CRSE = data(66);
  }
  
  return res;
}

void 
ModIMKPeakOriented::Print(OPS_Stream &s, int flag)
{
    s << "ModIMKPeakOriented tag: " << this->getTag() << endln;
    s << "  Ke: " << Ke << endln;	
    s << "  AlfaPos: " << AlfaPos << endln;
    s << "  AlfaNeg: " << AlfaNeg << endln;
    s << "  My_pos: " << My_pos << endln;
    s << "  My_neg: " << My_neg << endln;
    s << "  Ls: " << Ls << endln;
    s << "  Ld: " << Ld << endln;
    s << "  La: " << La << endln;
    s << "  Lk: " << Lk << endln;
    s << "  Cs: " << Cs << endln;
    s << "  Cd: " << Cd << endln;
    s << "  Ca: " << Ca << endln;
    s << "  Ck: " << Ck << endln;
    s << "  ThetaPpos: " << ThetaPpos << endln;
    s << "  ThetaPneg: " << ThetaPneg << endln;
    s << "  ThetaPCpos: " << ThetaPCpos << endln;
    s << "  ThetaPCneg: " << ThetaPCneg << endln;
    s << "  ResfacPos: " << ResfacPos << endln;
    s << "  ResfacNeg: " << ResfacNeg << endln;
    s << "  FracDispPos: " << FracDispPos << endln;
    s << "  FracDispNeg: " << FracDispNeg << endln;
    s << "  DPos: " << DPos << endln;
    s << "  DNeg: " << DNeg << endln;
}

// Positive and Negative Envelope SubRoutines

void 
ModIMKPeakOriented::envelPosCap2(double fy,double alphaPos,double alphaCap,double cpDsp,double& d,
					double& f, double& ek, double elstk,double fyieldPos,double Resfac, double fracDisp, int& flagStop)
{
		
  double dy = fy/elstk;
  
  double Res,rcap; // dres;
  if(dy<=cpDsp) {
    Res = Resfac*fyieldPos;
    rcap = fy+alphaPos*elstk*(cpDsp-dy);
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);
    
    if (d<0.0){
      f = 0.0;
      ek = 1.0e-7;
    } else if (d<=dy) { 
      ek = elstk;
      f = ek*d;
    } else if(d<=cpDsp) { 
      ek = elstk*alphaPos;
      f = fy+ek*(d-dy);
    } else if(d<=dres) { 
      ek = alphaCap*elstk;
      f = rcap+ek*(d-cpDsp);
    } else {
      ek = 1.0e-7;
      f = Res+d*ek;
    }
    if(d>=fracDisp) {
      ek = 1.0e-7;
      f = 1.0e-10;
      d=fracDisp;
      flagStop=1;
    }
  } else if(dy>cpDsp) { 
    
    rcap = elstk*cpDsp;
    Res = Resfac*rcap;
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);
    
    if (d<0.0) {
      f = 0.0;
      ek = 1.0e-7;
    } else if(d<=cpDsp) {
      ek = elstk;
      f = ek*d;
    } else if(d<=dres) {
      ek = alphaCap*elstk;
      f = rcap+ek*(d-cpDsp);
    } else {
      ek = 1.0e-7;
      f = Res+d*ek;
    }
    if(d>=fracDisp) {
      ek = 1.0e-7;
      f = 1.0e-10;
      d=fracDisp;
      flagStop=1;  
    }
  }
}

void
ModIMKPeakOriented::envelNegCap2(double fy,double alphaNeg,double alphaCap,double cpDsp,double& d,
				double& f, double& ek, double elstk, double fyieldNeg,double Resfac, double fracDisp, int& flagStop)
{

  double dy = fy/elstk;
  double Res,rcap; // dres;
  
  if(dy>=cpDsp){ 
    
    Res = Resfac*fyieldNeg;
    rcap = fy+alphaNeg*elstk*(cpDsp-dy);
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);
    
    if (d>0.0){ 
      f = 0.0;
      ek = 1.0e-7;
    } else if (d>=dy) { 
      ek = elstk;
      f = ek*d;
    } else if (d>=cpDsp) {
      ek = elstk*alphaNeg;
      f = fy+ek*(d-dy);
    } else if (d>=dres) {
      ek = elstk*alphaCap;
      f = rcap+ek*(d-cpDsp);
    } else {
      ek = 1.0e-7;
      f = Res+ek*d;
    }
    if(d<=fracDisp) {
      ek = 1.0e-7;
      f = 1.0e-10;
      d=fracDisp;
      flagStop=1;
    }
    
  } else if(dy<cpDsp) { 
    
    rcap = elstk*cpDsp;
    Res = Resfac*rcap;
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);
    
    if (d>0.0) {
      f = 0.0;
      ek = 1.0e-7;
    } else if (d>=cpDsp) {
      ek = elstk;
      f = ek*d;
    } else if (d>=dres) {
      ek = elstk*alphaCap;
      f = rcap+ek*(d-cpDsp);
    } else {
      ek = 1.0e-7;
      f = Res+ek*d;
    }
    if(d<=fracDisp) {
      ek = 1.0e-7;
      f = 1.0e-10;
      d=fracDisp;
      flagStop=1;
    }
  }
}


Response *
ModIMKPeakOriented::setResponse (const char **argv, int argc,
                                 OPS_Stream &theOutputStream)
{
  //by default, See if the response is one of the defaults
  Response *theResponse = UniaxialMaterial::setResponse(argv, argc, theOutputStream);
  
  if (theResponse != 0)      
    return theResponse;
  
  if ((strcmp(argv[0],"dres") == 0)) {
    theOutputStream.tag("ResponseType", "dres");
    theResponse =  new MaterialResponse(this, 101, dres);
  }
  
  return theResponse;

}

int
ModIMKPeakOriented::getResponse (int responseID, Information &matInformation)
{
  if (responseID == 101) {
    matInformation.setDouble(dres);
  } else
    return UniaxialMaterial::getResponse(responseID, matInformation);
  return 0;
}
