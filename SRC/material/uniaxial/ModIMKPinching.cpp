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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ModIMKPinching.cpp,v $
                                                                        
// Written: Dimitrios G. Lignos, PhD, Assistant Professor, McGill University 
// Created: August, 2012
// Revision: A
//
// Description: This file contains the class implementation for ModIMKPinching Model

#include <math.h>

#include <elementAPI.h>
#include <ModIMKPinching.h>
#include <Vector.h>
#include <Channel.h>

#include <OPS_Globals.h>

static int numModIMKPinchingMaterials = 0;

void *
OPS_ModIMKPinching()
{
  if (numModIMKPinchingMaterials == 0) {
    numModIMKPinchingMaterials++;
    opserr << "Modified Ibarra-Medina-Krawinkler Model with Pinched Hysteretic Response\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[1];
  double dData[26];
  int numData = 1;
  // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  ModIMKPinching tag" << endln;
    return 0;
  }
  
  numData = 26;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial ModIMKPinching tag? Ke?, alfaPos?, alfaNeg?, My_pos?, My_neg?"; 
    opserr << "FprPos?, FprNeg?, A_pinch?, Ls?, Ld?, La?, Lk?, Cs?, Cd?, Ca?, Ck?, thetaPpos?, thetaPneg?"; 
    opserr << "thetaPCpos?, thetaPCneg?, ResfacPos?, ResfacNeg?, fracDispPos?, fracDispNeg?,DPos?, DNeg?";
    
    return 0;	
  }
  
  // Parsing was successful, allocate the material with zero index
  theMaterial = new ModIMKPinching(iData[0], 
				   dData[0], dData[1], dData[2], dData[3],
				   dData[4], dData[5], dData[6], dData[7],
				   dData[8], dData[9], dData[10], dData[11],
				   dData[12], dData[13], dData[14], dData[15],
				   dData[16], dData[17], dData[18], dData[19],
				   dData[20], dData[21], dData[22], dData[23], dData[24], dData[25]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ModIMKPinching Material\n";
    return 0;
  }

  return theMaterial;
}


ModIMKPinching::ModIMKPinching(int tag, double ke, double alfaPos, double alfaNeg, double my_pos, double my_neg, 
			       double fprPos, double fprNeg, double a_pinch, double ls, double ld, double la,
									   double lk, double cs, double cd, double ca, double ck, double thetaPpos, double thetaPneg,
			       double thetaPCpos, double thetaPCneg, double resfacPos, double resfacNeg,
			       double fracDispPos, double fracDispNeg, double dPos, double dNeg)
 :UniaxialMaterial(tag,MAT_TAG_ModIMKPinching), Ke(ke), AlfaPos(alfaPos), AlfaNeg(alfaNeg), My_pos(my_pos), My_neg(my_neg), 
  FprPos(fprPos), FprNeg(fprNeg), A_pinch(a_pinch), Ls(ls), Ld(ld), La(la), Lk(lk), Cs(cs), Cd(cd), Ca(ca), Ck(ck), 
  ThetaPpos(thetaPpos), ThetaPneg(thetaPneg), ThetaPCpos(thetaPCpos), ThetaPCneg(thetaPCneg), ResfacPos(resfacPos), ResfacNeg(resfacNeg), 
  FracDispPos(fracDispPos), FracDispNeg(fracDispNeg),DPos(dPos), DNeg(dNeg)
{
  // Initialize Variables by Calling revertToStart function
  this->revertToStart();
  Tangent = Ke;
}

ModIMKPinching::ModIMKPinching()
 :UniaxialMaterial(0,MAT_TAG_ModIMKPinching),
  Ke(0.0), AlfaPos(0.0), AlfaNeg(0.0), My_pos(0.0), My_neg(0.0), FprPos(0.0), FprNeg(0.0), A_pinch(0.0),
  Ls(0.0), Ld(0.0), La(0.0), Lk(0.0), Cs(0.0), Cd(0.0), Ca(0.0), Ck(0.0), ThetaPpos(0.0), ThetaPneg(0.0), 
  ThetaPCpos(0.0), ThetaPCneg(0.0), ResfacPos(0.0), ResfacNeg(0.0), 
  FracDispPos(0.0), FracDispNeg(0.0), DPos(0.0), DNeg(0.0)
{
  // Initialize Variables by Calling revertToStart function
  //	revertToStart();
}

ModIMKPinching::~ModIMKPinching()
{
  // does nothing
}

int 
ModIMKPinching::setTrialStrain(double strain, double strainRate)
{
  //all variables to the last commit state
  
  this->revertToLastCommit();
  
  double f = 0.0;          // current force
  double d  = strain;      // current deformation
  
  // Determine change in deformation from last converged state
  double deltaD = d - dP;
  
  // added by DL 10/10/2013
  if (fabs(deltaD) < 1.0e-18 && strain != 0.0) {
    return 0;
  }	
  
  // Initialize parameters in the first cycle
  if (kon==0) {
    
    // Post Yield Stiffness
    ekhardPos = Ke * AlfaPos;
    ekhardNeg = Ke * AlfaNeg;
    
    // Reference Energy Capacity as defined Lignos (2008)
    Enrgts = Ls * My_pos;
    Enrgtk = 2.0 * Lk * My_pos;
    Enrgta = La * My_pos;
    Enrgtd = Ld * My_pos;
    
    ekP = 0.0;
    ekunload = Ke;
    ekexcurs = Ke;
    
    Enrgtot = 0.0;
    Enrgc = 0.0;
    
    // updated yield forces per exursion in positive and negative loading direction
    fyPos = My_pos;
    fyNeg = My_neg;
    
    // Capping Deformation
    cpPos = My_pos/Ke + ThetaPpos;
    cpNeg = My_neg/Ke + (-ThetaPneg);
    
    // Capping Slope
    capSlopePos = -(My_pos + ekhardPos * (ThetaPpos))/(ThetaPCpos * Ke);
    capSlopeNeg =  (My_neg + ekhardNeg * (-ThetaPneg))/(ThetaPCneg * Ke);
    
    //Reference Point for Cyclic Deterioration
    
    fCapRefPos = -capSlopePos * Ke * cpPos + (My_pos + ekhardPos * (ThetaPpos));
    fCapRefNeg = -capSlopeNeg * Ke * cpNeg + (My_neg + ekhardNeg * (-ThetaPneg));
    
    // Capping Force (updated per excursion)
    fCapPos = My_pos + ekhardPos * ThetaPpos;
    fCapNeg = My_neg + ekhardNeg * (-ThetaPneg);
    
    dmax = My_pos/Ke;            // max deformation in positive loading direction
    dmin = My_neg/Ke;            // max deformation in negative loading direction
    
    fmax = My_pos;               // maximum force at positive loading direction
    fmin = My_neg;               // maximum force at negative loading direction
    
    fpDegPos = FprPos * fmax;    // Positive pinched force
    fpDegNeg = FprNeg * fmin;    // Negative piched force
    
    sp = 0.0;         // deformation at which unloading range in positive loading direction hits zero
    sn = 0.0;         // deformation at which unloading range in negative loading direction hits zero
    
    ek = 0.0;         // current stiffness
    fP = 0.0;         // previous force
    
    flagdeg = 0;      // Flag that enables cyclic deterioration
    flagStop = 0;     // if deformation exceeds theta_u for +/- loading direction f = 0;
    
    // Previous Deformation
    dP = 0.0;
    
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
	
	double RSE = 0.5 * fP * fP / ekunload;
	double a2 = Enrgtk-(Enrgtot-RSE);
	double betak;
	double ekLim;
	
	if(a2 < 0.0 && Enrgtk != 0.0){
	  flagStop = 1;
	}
	
	if(Lk != 0.0){
	  betak = pow((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE)),Ck);
	  ekunload = ekexcurs*(1.0-betak);
	  ekLim = (fmax - fmin) /(dmax - dmin);
	  
	  if(ekunload <= ekLim){
	    ekunload = ekLim;
	  } 				
	  
	  if(ekunload <= ekhardNeg){
	    flagStop = 1;
	  }           	
	}
	
	//  determination of sn according to the hysteresis status 
	if(ekunload <= 1.e-7){
	  flagStop = 1;
	}           	
	
	if(fP <= 0.0){
	  
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
	
      } else if (fabs(sn)>1.e-10){
	
	envelPosCap2 (fyPos,AlfaPos,capSlopePos,cpPos,dmax,fmax, ek, Ke, My_pos, ResfacPos,FracDispPos, flagStop);
	
	double dch = dmax - fmax/ekunload;
	fpDegPos = fmax * FprPos;
	double ekpinch = fpDegPos / (dmax - sn);
	double fpinch = ekpinch * (A_pinch * dch - sn);
	
	if (sn <= A_pinch * dch) {
	  
	  if ( d <= sn){
	    
	    ek = ekunload;
	    f  = fP+ek*deltaD;
	    
	  } else if (d >= sn && d < A_pinch * dch) {
	    
	    ek = ekpinch;
	    double f2 = ek * (d-sn);
	    double f1 = fP + ekunload * deltaD;
	    //f = min(f1,f2);
	    if (f1 < f2) {
	      f = f1;
	      
	    } else {
	      
	      f = f2;
	    }
	    
	    if (ekunload < ek) {
	      flagStop = 1;
	    }
	    
	    if (fabs(f-f1)<1.e-10) {
	      ek=ekunload;
	    }
	  } else {
	    
	    ek = (fmax - fpinch)/(dmax - A_pinch * dch );
	    
	    double f2 = fpinch + ek * (d-A_pinch*dch);
	    double f1 = fP+ekunload * deltaD;
	    //f = min(f1,f2);
	    if (f1 < f2) {
	      f = f1;
	      
	    } else {
	      
	      f = f2;
	    }
	    
	    if (ekunload < ek) {
	      flagStop = 1;
	    }
	    
	    if (fabs(f-f1)<1.e-10) {
	      ek=ekunload;
	    }
	  }						
	  
	} else if (sn > A_pinch * dch ) {
	  if (d <=sn) {
	    ek = ekunload;
	    f = fP + ek * deltaD;
	  } else {
	    ek = fmax/(dmax-sn);
	    double f2 = ek*(d-sn);
	    double f1 = fP+ekunload * deltaD;
	    //f = min(f1,f2);
	    if (f1 < f2) {
	      f = f1;
	      
	    } else {
	      
	      f = f2;
	    }
	    
	    if (ekunload < ek) {
	      flagStop = 1;
	    }
	    
	    if (fabs(f-f1)<1.e-10) {
	      ek=ekunload;
	    }
	  }
	}
	
      } else {							
	if (d > 0.0) {
	  envelPosCap2(fyPos,AlfaPos,capSlopePos,cpPos,d, f, ek, Ke, My_pos, ResfacPos,FracDispPos, flagStop);	
	} else {	
	  envelNegCap2(fyNeg,AlfaNeg,capSlopeNeg,cpNeg, d, f, ek, Ke, My_neg, ResfacNeg,-FracDispNeg, flagStop);	
	}
	
      }
      
      
    } else {
      if (kon == 1) {
	kon = 2;
	
	double RSE = 0.5 * fP * fP/ekunload;						
	double a2 = Enrgtk-(Enrgtot-RSE);
	double betak;
	double ekLim;
	if (a2 <= 0.0 && Lk !=0.0) {
	  flagStop = 1;
	}
	
	if(Lk!= 0.0 && dmax > My_pos/Ke){
	  betak = pow((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE)),Ck);
	  ekunload = ekexcurs*(1.0-betak);
	  ekLim = (fmax - fmin) / (dmax - dmin);
	  
	  if (ekunload <= ekLim) {
	    ekunload = ekLim;
	  }
	  
	  if(ekunload<= ekhardPos){
	    flagStop=1;
	  }           
	}
	
	//		Determination of sn according to the hysteresis status
	if (ekunload <= 1.e-7) {
	  flagStop = 1;
	}
	
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
	
      } else if (fabs(sp) > 1.e-10){
	envelNegCap2 (fyNeg,AlfaNeg,capSlopeNeg,cpNeg, dmin, fmin, ek, Ke, My_neg, ResfacNeg,-FracDispNeg, flagStop);
	
	double dch = dmin - fmin/ekunload;
	fpDegNeg = fmin * FprNeg;
	double ekpinch = fpDegNeg / (dmin - sp);
	double fpinch = ekpinch * (A_pinch * dch - sp);
	
	if (sp >= A_pinch * dch ) {
	  
	  if ( d >= sp){
	    
	    ek = ekunload;
	    f = fP+ek*deltaD;
	    
	  } else if(d <= sp && d > A_pinch * dch){
	    ek = ekpinch;
	    double f2 = ek * (d-sp);
	    double f1 = fP + ekunload * deltaD;
	    
	    // f=max(f1,f2);
	    if (f1>f2)
	      {
		f=f1;
	      }
	    else
	      {
		f=f2;
	      }
	    if (ekunload < ek) {
	      flagStop = 1;
	    }
	    if (fabs(f-f1) < 1.e-10){ 
	      ek=ekunload;
	    }						
	    
	  } else {
	    
	    ek = (fmin - fpinch)/(dmin-A_pinch * dch);
	    double f2 = fpinch + ek*(d-A_pinch * dch);
	    double f1 = fP + ekunload * deltaD;
	    
	    // f=max(f1,f2);
	    if (f1>f2)
	      {
		f=f1;
	      }
	    else
	      {
		f=f2;
	      }
	    if (ekunload < ek) {
	      flagStop = 1;
	    }
	    if (fabs(f-f1) < 1.e-10){ 
	      ek=ekunload;
	    }						
	  }
	} else if(sp < A_pinch * dch ){
	  if (d >= sp) {
	    ek = ekunload;
	    f = fP + ek * deltaD;
	  } else {
	    ek = fmin/(dmin-sp);
	    double f2 = ek*(d-sp);
	    double f1 = fP+ekunload * deltaD;
	    
	    // f=max(f1,f2);
	    if (f1 > f2) {
	      f = f1;
	    } else {
	      f=f2;
	    }
	    if (ekunload < ek) {
	      flagStop = 1;
	    }
	    if (fabs(f-f1) < 1.e-10){ 
	      ek=ekunload;
	    }
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
    
    //if (Ls!=0.0 && Lk!=0.0 && La!=0.0 && Lc != 0.0) {
    Enrgi = 0.5 * (f+fP) * deltaD;
    Enrgc = Enrgc + Enrgi;
    Enrgtot = Enrgtot + Enrgi;
    //}
    
    //	UPDATING OF DETERIORATION PARAMETERS ----------------------------
    
    // Check the remaining capacity of the system
    double betas = 0.0;
    double betaa = 0.0;
    double betad = 0.0;
    
    if(flagdeg == 1){
      
      if(Enrgtot >= Enrgts && Enrgts != 0.0){	    
	betas = 1.0;
      } else if (Enrgtot >= Enrgta && Enrgta != 0.0) {  
	betaa = 1.0;
      } else if(Enrgtot >= Enrgtd && Enrgtd == 0.0) {  
	betad = 1.0;
      }
      
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
	betad = pow(Enrgc/(Enrgtd-Enrgtot),Cd);
      } else {
	betad = 0.0;
      }
      
      if(fabs(betas) >= 1.0) {
	betas = 1.0;
	flagStop = 1;
      }
      if(fabs(betaa) >= 1.0){ 
	betaa = 1.0;
	flagStop = 1;
      }
      if(fabs(betad) >= 1.0){ 
	betad = 1.0;
	flagStop = 1;
      }
    }
    
    //	Update values for the next half cycle
    ekexcurs = ekunload;
    Enrgc = 0.0;
    //	Deteriorate parameters for the next half cycle
    
    if(deltaD > 0.0 && f > 0.0 || deltaD < 0.0 && f >=0.0 && d <=sp) {
      
      fyNeg = fyNeg * (1.0 - betas * DNeg);	
      fpDegNeg = fpDegNeg * (1.0 - betas * DNeg);
      AlfaNeg=AlfaNeg * (1.0 - betas * DNeg);			
      fCapRefNeg=fCapRefNeg * (1.0-betad * DNeg);
      
      dmin = dmin * (1.0+betaa);
      
      double dyNeg = fyNeg/Ke;
      ekhardNeg = AlfaNeg * Ke;
      
      double dCap1Neg=fCapRefNeg/(Ke-capSlopeNeg*Ke);
      double dCap2Neg=(fCapRefNeg+ekhardNeg*dyNeg-fyNeg)/(ekhardNeg-capSlopeNeg*Ke);
      // cpNeg=min(dCap1Neg,dCap2Neg);
      
      if (dCap1Neg < dCap2Neg){
	cpNeg = dCap1Neg;
      } else {
	cpNeg = dCap2Neg;
      }
      
      fCapNeg = fCapRefNeg + capSlopeNeg * Ke * cpNeg;
      
    } else if (deltaD < 0.0 && f < 0.0 || deltaD > 0.0 && f <=0.0 && d >= sn) {
      
      fyPos = fyPos*(1.0-betas * DPos);
      fpDegPos = fpDegPos * (1.0 - betas * DPos);
      AlfaPos=AlfaPos*(1.0 - betas * DPos);
      fCapRefPos=fCapRefPos*(1.0 - betad * DPos);
      
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
      fCapPos = fCapRefPos + capSlopePos * Ke * cpPos;
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

double ModIMKPinching::getStress(void)
{
  
  return  (fP);
}

double ModIMKPinching::getTangent(void)
{
  return (Tangent);
}

double ModIMKPinching::getInitialTangent(void)
{
  return (Ke);
}

double ModIMKPinching::getDampTangent(void)
{
  double DTangent = ekP;
  return DTangent;	
}


double ModIMKPinching::getStrain(void)
{
  return dP;
}

double ModIMKPinching::getStrainRate(void)
{
    return 0;
}

int ModIMKPinching::commitState(void)
{
  Cstrain = dP;
  Cstress = fP;
  Ctangent = Tangent;
  
  CdP = dP;
  CfP = fP;
  Cek = ek;
  
  Ckon = kon;
  CflagStop = flagStop;
  Cflagdeg = flagdeg;
  
  Cdmax = dmax;
  Cdmin = dmin;
  Cfmax = fmax;
  Cfmin = fmin;
  
  CfyPos = fyPos;
  CfyNeg = fyNeg;
  
  Csn = sn;
  Csp = sp;
  
  CEnrgc = Enrgc;
  CEnrgtot = Enrgtot;
  
  CEnrgts = Enrgts;
  CEnrgtd = Enrgtd;
  CEnrgtk = Enrgtk;
  CEnrgta = Enrgta;
  
  CcapSlopePos = capSlopePos;
  CcapSlopeNeg = capSlopeNeg;
  
  CfCapRefPos = fCapRefPos;
  CfCapRefNeg = fCapRefNeg;
  
  CfCapPos = fCapPos;
  CfCapNeg = fCapNeg;
  
  CfpDegPos = fpDegPos;
  CfpDegNeg = fpDegNeg;
  
  Cekunload = ekunload;
  
  CcpNeg = cpNeg;
  CcpPos = cpPos;
  
  CekhardPos = ekhardPos;
  CekhardNeg = ekhardNeg;
  Cekexcurs = ekexcurs;
  CekP = ekP;
  
  return 0;
}


int ModIMKPinching::revertToLastCommit(void)
{
  dP = Cstrain;
  fP = Cstress;
  Tangent = Ctangent;
  
  dP = CdP;
  fP = CfP;
  ek = Cek;
  
  kon = Ckon;
  flagStop = CflagStop;
  flagdeg = Cflagdeg;
  
  dmax = Cdmax;
  dmin = Cdmin;
  fmax = Cfmax;
  fmin = Cfmin;
  
  fyPos = CfyPos;
  fyNeg = CfyNeg;
  
  sn = Csn;
  sp = Csp;
  
  // Energy Variables
  Enrgc = CEnrgc;
  Enrgtot = CEnrgtot;
  
  Enrgts = CEnrgts;
  Enrgtd = CEnrgtd;
  Enrgtk = CEnrgtk;
  Enrgta = CEnrgta;
  
  capSlopePos = CcapSlopePos;
  capSlopeNeg = CcapSlopeNeg;
  
  fCapRefPos = CfCapRefPos;
  fCapRefNeg = CfCapRefNeg;
  
  fCapPos = CfCapPos;
  fCapNeg = CfCapNeg;
  
  fpDegPos = CfpDegPos;
  fpDegNeg = CfpDegNeg;
  
  ekunload = Cekunload;
  
  cpNeg = CcpNeg;
  cpPos = CcpPos;
  
  ekhardPos = CekhardPos;
  ekhardNeg = CekhardNeg;
  ekexcurs = Cekexcurs;
  ekP = CekP;
  
  return 0;
}

int ModIMKPinching::revertToStart(void)
{
  
  // Initialize state variables	
  dP = 0.0;
  fP = 0.0;
  ek = Ke;
  
  kon = 0;
  flagStop = 0;
  flagdeg = 0;
  
  fyPos = My_pos;
  fyNeg = My_neg;
  
  ekhardPos = Ke * AlfaPos;
  ekhardNeg = Ke * AlfaNeg;
  
  sn = 0.0;
  sp = 0.0;
  
  // Energy Variables
  Enrgc = 0.0;
  Enrgtot = 0.0;
  
  Enrgts = Ls * My_pos;
  Enrgtd = Ld * My_pos;
  Enrgtk = 2.0 * Lk * My_pos;
  Enrgta = La * My_pos;
  
  cpPos = My_pos/Ke + ThetaPpos;
  cpNeg = My_neg/Ke + (-ThetaPneg);
  
  capSlopePos = -(My_pos + ekhardPos * (ThetaPpos))/(ThetaPCpos * Ke);
  capSlopeNeg = (My_neg + ekhardNeg * (-ThetaPneg))/(ThetaPCneg * Ke);
  
  fCapRefPos = -capSlopePos * Ke * (ThetaPpos + My_pos/Ke) +  (My_pos + ekhardPos * (ThetaPpos));
  fCapRefNeg = -capSlopeNeg * Ke * (-ThetaPneg + My_neg/Ke) + (My_neg + ekhardNeg * (-ThetaPneg));
  
  // Capping Force (updated per excursion)
  fCapPos = My_pos + ekhardPos * ThetaPpos;
  fCapNeg = My_neg + ekhardNeg * (-ThetaPneg);	
  
  dmax = My_pos/Ke;
  dmin = My_neg/Ke;
  
  fmax = My_pos;
  fmin = My_neg;
  
  fpDegPos = FprPos * fmax;
  fpDegNeg = FprNeg * fmin;
  
  ekunload = Ke;		
  ekexcurs = Ke;
  ekP = Ke;
  
  Cstrain=0.0;
  Cstress = 0.0;
  Ctangent = Ke;
  Tangent = Ke;
  
  CdP = 0.0;
  CfP = 0.0;
  Cek = ek;
  
  Ckon = 0;
  CflagStop = 0;
  Cflagdeg = 0;
  
  CfyPos = My_pos;
  CfyNeg = My_neg;
  
  CekhardPos = Ke * AlfaPos;
  CekhardNeg = Ke * AlfaNeg;
  
  Csn = 0.0;
  Csp = 0.0;
  
  // Energy Variables
  CEnrgc = 0.0;
  CEnrgtot = 0.0;
  
  CEnrgts = Ls * My_pos;
  CEnrgtd = Ld * My_pos;
  CEnrgtk = 2.0 * Lk * My_pos;
  CEnrgta = La * My_pos;
  
  CcpNeg = My_pos/Ke + ThetaPpos;
  CcpPos = My_neg/Ke + (-ThetaPneg);
  
  CcapSlopePos =  -(My_pos + ekhardPos * (ThetaPpos))/(ThetaPCpos * Ke);
  CcapSlopeNeg = (My_neg + ekhardNeg * (-ThetaPneg))/(ThetaPCneg * Ke);
  
  CfCapRefPos = -capSlopePos * Ke * cpPos + (My_pos + ekhardPos * (ThetaPpos));
  CfCapRefNeg = -capSlopeNeg * Ke * cpNeg + (My_neg + ekhardNeg * (-ThetaPneg));
  
  // Capping Force (updated per excursion)
  CfCapPos = My_pos + ekhardPos * ThetaPpos;
  CfCapNeg = My_neg + ekhardNeg * (-ThetaPneg);	
  
  Cdmax = My_pos/Ke;
  Cdmin = My_neg/Ke;
  
  Cfmax = My_pos;
  Cfmin = My_neg;
  
  CfpDegPos = FprPos * fmax;
  CfpDegNeg = FprNeg * fmin;
  
  Cekunload = Ke;
  Cekexcurs = Ke;
  CekP = Ke;	
  
  return 0;
}

UniaxialMaterial *
ModIMKPinching::getCopy(void)
{
  ModIMKPinching *theCopy = new ModIMKPinching(this->getTag(), Ke, AlfaPos, AlfaNeg, My_pos, My_neg, FprPos, FprNeg, A_pinch, Ls, Ld, La, Lk, 
					       Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg,ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg,FracDispPos, FracDispNeg,DPos, DNeg);
  // Converged state variables
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Ctangent = Ctangent;
  
  theCopy->CdP = CdP;
  theCopy->CfP = CfP;
  theCopy->Cek = ek;
  
  theCopy->Ckon = Ckon;
  theCopy->CflagStop = CflagStop;
  theCopy->Cflagdeg = Cflagdeg;
  
  theCopy->Cdmax = Cdmax;
  theCopy->Cdmin = Cdmin;
  theCopy->Cfmax = Cfmax;
  theCopy->Cfmin = Cfmin;
  
  theCopy->CfyPos = CfyPos;
  theCopy->CfyNeg = CfyNeg;
  
  theCopy->CekhardPos = CekhardPos;
  theCopy->CekhardNeg = CekhardNeg;
  
  theCopy->ekhardPos = ekhardPos;
  theCopy->ekhardNeg = ekhardNeg;
  
  theCopy->Csn= Csn;
  theCopy->Csp= Csp;
  
  // Energy Variables
  theCopy->CEnrgc = CEnrgc;
  theCopy->CEnrgtot = CEnrgtot;
  
  theCopy->CEnrgts = CEnrgts;
  theCopy->CEnrgtd = CEnrgtd;
  theCopy->CEnrgtk = CEnrgtk;
  theCopy->CEnrgta = CEnrgta;
  
  theCopy->CcapSlopePos = CcapSlopePos;
  theCopy->CcapSlopeNeg = CcapSlopeNeg;
  
  theCopy->CfCapRefPos = CfCapRefPos;
  theCopy->CfCapRefNeg = CfCapRefNeg;
  
  theCopy->CfCapPos = CfCapPos;
  theCopy->CfCapNeg = CfCapNeg;
  
  theCopy->CfpDegPos = CfpDegPos;
  theCopy->CfpDegNeg = CfpDegNeg;
  
  theCopy->Cekunload = Cekunload;
  
  theCopy->CcpNeg = CcpNeg;
  theCopy->CcpPos = CcpPos;
  
  theCopy->Cekexcurs = Cekexcurs;
  theCopy->CekP = CekP;
  
  // Trial state variables
  theCopy->dP = dP;
  theCopy->fP = fP;
  theCopy->ek = ek;
  
  theCopy->kon = kon;
  theCopy->flagStop = flagStop;
  theCopy->Cflagdeg = Cflagdeg;
  
  theCopy->dmax = dmax;
  theCopy->dmin = dmin;
  theCopy->fmax = fmax;
  theCopy->fmin = fmin;
  
  theCopy->fyPos = fyPos;
  theCopy->fyNeg = fyNeg;
  
  theCopy->sn= sn;
  theCopy->sp= sp;
  
  // Energy Variables
  theCopy->Enrgc = Enrgc;
  theCopy->Enrgtot = Enrgtot;
  
  theCopy->Enrgts = Enrgts;
  theCopy->Enrgtd = Enrgtd;
  theCopy->Enrgtk = Enrgtk;
  theCopy->Enrgta = Enrgta;
  
  theCopy->capSlopePos = capSlopePos;
  theCopy->capSlopeNeg = capSlopeNeg;
  
  theCopy->fCapRefPos = fCapRefPos;
  theCopy->fCapRefNeg = fCapRefNeg;
  
  theCopy->fpDegPos = fpDegPos;
  theCopy->fpDegNeg = fpDegNeg;
  
  theCopy->fCapPos = fCapPos;
  theCopy->fCapNeg = fCapNeg;
  
  theCopy->ekunload = ekunload;
  
  theCopy->cpNeg = cpNeg;
  theCopy->cpPos = cpPos;
  
  theCopy->ekexcurs = ekexcurs;
  theCopy->ekP = ekP;
  
  return theCopy;
}

int 
ModIMKPinching::sendSelf(int cTag, Channel &theChannel)
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
    data(6) = FprPos;
    data(7) = FprNeg;
    data(8) = A_pinch;
    data(9) = Ls;
    data(10) = Ld;
    data(11) = La;
    data(12) = Lk;
    data(13) = Cs;
    data(14) = Cd;
    data(15) = Ca;
    data(16) = Ck;
    data(17) = ThetaPpos;
    data(18) = ThetaPneg;
    data(19) = ThetaPCpos;
    data(20) = ThetaPCneg;
    data(21) = ResfacPos;
    data(22) = ResfacNeg;
    data(23) = FracDispPos;
    data(24) = FracDispNeg;
    data(25) = DPos;
    data(26) = DNeg;
    
    // State variables from last converged state
    data(27) = Cstrain;
    data(28) = Cstress;
    data(29) = Ctangent;
    
    data(30) = CdP;
    data(31) = CfP;
    data(32) = Cek;
    
    //  Variables 
    data(33) = Ckon;
    data(34) = CflagStop;
    data(35) = Cdmax;
    data(36) = Cdmin;
    data(37) = Cfmax;
    data(38) = Cfmin;
    
    data(39) = CfyPos;
    data(40) = CfyNeg;
    
    data(41) = Csn;
    data(42) = Csp;
    
    // Energy Variables
    data(43) = CEnrgc;
    data(44) = CEnrgtot;
    
    data(45) = CEnrgts;
    data(46) = CEnrgtd;
    data(47) = CEnrgtk;
    data(48) = CEnrgta;
    
    data(49) = CcapSlopePos;
    data(50) = CcapSlopeNeg;
    
    data(51) = CfCapRefPos;
    data(52) = CfCapRefNeg;
    
    data(53) = CfCapPos;
    data(54) = CfCapNeg;
    
    data(55) = CfpDegPos;
    data(56) = CfpDegNeg;
    
    data(57) = Cekunload;
    
    data(58) = CcpNeg;
    data(59) = CcpPos;
    
    data(60) = CekhardPos;
    data(61) = CekhardNeg;
    data(62) = Cekexcurs;
    data(63) = CekP;
    
    data(64) = Cflagdeg;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
      opserr << "ModIMKPinching::sendSelf() - failed to send data\n";
    
    return res;
}

int 
ModIMKPinching::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(67);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
    opserr << "ModIMKPinching::recvSelf() - failed to receive data\n";
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
    FprPos = data(6);
    FprNeg = data(7);
    A_pinch = data(8);
    Ls = data(9);
    Ld = data(10);
    La = data(11);
    Lk = data(12);
    Cs = data(13);
    Cd = data(14);
    Ca = data(15);
    Ck = data(16);
    ThetaPpos = data(17);
    ThetaPneg = data(18);
    ThetaPCpos = data(19);
    ThetaPCneg = data(20);
    ResfacPos = data(21);
    ResfacNeg = data(22);
    FracDispPos = data(23);
    FracDispNeg = data(24);
    DPos = data(25);
    DNeg = data(26);
    
    // State variables from last converged state 	
    Cstrain = data(27);
    Cstress = data(28);
    Ctangent = data(29);
    
    CdP = data(30);
    CfP = data(31);
    Cek = data(32);
    
    Ckon = data(33);
    CflagStop = data(34);
    Cdmax = data(35);
    Cdmin = data(36);
    Cfmax = data(37);
    Cfmin = data(38);
    
    CfyPos = data(39);
    CfyNeg = data(40);
    
    Csn = data(41);
    Csp = data(42);
    
    // Energy Variables
    CEnrgc = data(43);
    CEnrgtot = data(44);
    
    CEnrgts = data(45);
    CEnrgtd = data(46);
    CEnrgtk = data(47);
    CEnrgta = data(48);
    
    CcapSlopePos = data(49);
    CcapSlopeNeg = data(50);
    
    CfCapRefPos = data(51);
    CfCapRefNeg = data(52);
    
    CfCapPos = data(53);
    CfCapNeg = data(54);
    
    CfpDegPos = data(55);
    CfpDegNeg = data(56);
    
    Cekunload = data(57);
    
    CcpNeg = data(58);
    CcpPos = data(59);
    
    CekhardPos = data(60);
    CekhardNeg = data(61);
    Cekexcurs = data(62);
    CekP = data(63);
    
    Cflagdeg = data(64);
  }
  
  return res;
}

void 
ModIMKPinching::Print(OPS_Stream &s, int flag)
{
  s << "ModIMKPinching tag: " << this->getTag() << endln;
  s << "  Ke: " << Ke << endln;	
  s << "  AlfaPos: " << AlfaPos << endln;
  s << "  AlfaNeg: " << AlfaNeg << endln;
  s << "  My_pos: " << My_pos << endln;
  s << "  My_neg: " << My_neg << endln;
  s << "  FprPos: " << FprPos << endln;
  s << "  FprNeg: " << FprNeg << endln;
  s << "  A_Pinch: " << A_pinch << endln;	
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
ModIMKPinching::envelPosCap2(double fy,double alphaPos,double alphaCap,double cpDsp,double& d,
			     double& f, double& ek, double elstk,double fyieldPos,double Resfac, double fracDisp, int& flagStop)
{
  
  double dy = fy/elstk;
  
  double Res,rcap,dres;
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
ModIMKPinching::envelNegCap2(double fy,double alphaNeg,double alphaCap,double cpDsp,double& d,
			     double& f, double& ek, double elstk, double fyieldNeg,double Resfac, double fracDisp, int& flagStop)
{
  
  double dy = fy/elstk;
  double Res,rcap,dres;
  
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
