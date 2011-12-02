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
                                                                        

#include <math.h>
#include <string.h>

#include <Bilin.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>

#include <OPS_Globals.h>
#include <Parameter.h>

static int numBilinMaterials = 0;

void *
OPS_NewBilinMaterial()
{
  if (numBilinMaterials == 0) {
    numBilinMaterials++;
    opserr << "Modified Ibarra-Krawinkler Model with Bilinear Hysteretic Response\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[23];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Bilin tag" << endln;
    return 0;
  }

  numData = 23;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial Bilin tag? Ke? As? AsNeg? My_pos? My_neg? LamdaS? ";
    opserr << "LamdaK?  LamdaA? LamdaD? Cs? Ck? Ca? Cd? Thetap_pos? Thetap_neg? Thetapc_pos? Thetapc_neg?K? ";
    opserr << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg\n";
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new Bilin(iData[0], 
			  dData[0], dData[1], dData[2], dData[3], dData[4], 
			  dData[5], dData[6], dData[7], dData[8], dData[9], 
			  dData[10], dData[11], dData[12], dData[13], dData[14], 
			  dData[15], dData[16], dData[17], dData[18], dData[19], 
			  dData[20], dData[21], dData[22]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Bilin Material\n";
    return 0;
  }

  return theMaterial;
}

Bilin::Bilin(int tag, double p_Ke,double p_As,double p_AsNeg,double p_My_pos,double p_My_neg,double p_LamdaS,
	     double p_LamdaK,double p_LamdaA,double p_LamdaD,double p_Cs,double p_Ck,double p_Ca,double p_Cd,
	     double p_Thetap_pos,double p_Thetap_neg,double p_Thetapc_pos,double p_Thetapc_neg,double p_K,double p_KNeg,
	     double p_Thetau_pos,double p_Thetau_neg,double p_PDPlus,double p_PDNeg) 
:UniaxialMaterial(tag, MAT_TAG_Bilin),Ke(p_Ke), As(p_As), AsNeg(p_AsNeg), My_pos(p_My_pos), My_neg(p_My_neg), 
 LamdaS(p_LamdaS), LamdaK(p_LamdaK),LamdaA(p_LamdaA),LamdaD(p_LamdaD), Cs(p_Cs), Ck(p_Ck), Ca(p_Ca),Cd(p_Cd),
 Thetap_pos(p_Thetap_pos), Thetap_neg(p_Thetap_neg), Thetapc_pos(p_Thetapc_pos),Thetapc_neg(p_Thetapc_neg),
 K(p_K), KNeg(p_KNeg),Thetau_pos(p_Thetau_pos), Thetau_neg(p_Thetau_neg), PDPlus(p_PDPlus), PDNeg(p_PDNeg)
{
  //initialize variables
  this->revertToStart();

  //this->revertToLastCommit();
}

Bilin::Bilin()
:UniaxialMaterial(0, MAT_TAG_Bilin),
 Ke(0), As(0), AsNeg(0), My_pos(0), My_neg(0), 
 LamdaS(0), LamdaK(0),LamdaA(0),LamdaD(0), Cs(0), Ck(0), Ca(0),Cd(0),
 Thetap_pos(0), Thetap_neg(0), Thetapc_pos(0),Thetapc_neg(0),
 K(0), KNeg(0),Thetau_pos(0), Thetau_neg(0), PDPlus(0), PDNeg(0)
{
  this->revertToStart();
}

Bilin::~Bilin()
{
  // does nothing
}

int 
Bilin::setTrialStrain(double strain, double strainRate)
{  
  //all variables to the last commit
  this->revertToLastCommit();

  //Here I declare function variables
  double Ue,ddise,deltaD,d,temp_1_farzin,temp,betas,betak,betad,ekhard,
    dBoundPos,dBoundNeg,f1,f2,fn,e1,e2,a2,fNewLoadPos,
    f,xDevPos1,yDevPos1,xDevPos2,yDevPos2,xDevPos,yDevPos,fNewLoadNeg,xDevNeg1,
    yDevNeg1,xDevNeg2,yDevNeg2,xDevNeg,yDevNeg,Enrgi,v1,d1,ener;
  
  int NoCollapse,kst,jcode =0;

 
 //state determination algorithm: defines the current force and tangent stiffness
 U=strain; //set trial displacement
 
 //********************************matlab code*************************************************************
 //// total displacement
 Ue = U;
 //// incremental displacement
 ddise=Ue-CU;  //the incremental displacement of the element (from the last commited)
 //// update Uprev
 Uprev=Ue; 
 ////
 //// calculation of deltaD
 deltaD=ddise;
 d=Ue; 
 dP = d-deltaD;	   
 if (d>0.0) 
   {
     
     //function call
     interPoint(temp_1_farzin,temp,0.0,fCapRefPos,capSlope*elstk,0.0,Resfac*fyieldPos,0.0); 
     if (d<temp_1_farzin){ 
       iNoFpos = 0;
       LP=0;
     }
   } else {
   interPoint(temp_1_farzin,temp,0.0,fCapRefNeg,capSlopeNeg*elstk,0.0,ResfacNeg*fyieldNeg,0.0); 
   if (d>temp_1_farzin)
     {
       iNoFneg = 0;
       LN=0;
     }
 }
 //// Other variables
 flagdeg = 0;
 betas = 0.0;
 betak = 0.0;  
 betad = 0.0;
 ////  yield siplacements
 dyieldPos = fyieldPos/elstk;
 dyieldNeg = fyieldNeg/elstk;
 //// Initialize parameters in the first cycle
 if (kon==0) {
   //Because I included the statement revertToLastCommit above:
   elstk          = Ke; 
   fyieldPos      = My_pos; 
   fyieldNeg      = My_neg; 
   alpha          = As; 
   alphaN         = AsNeg;  
   ecaps          = LamdaS/(fyieldPos/(elstk)); 
   ecapk          = LamdaK/(fyieldPos/(elstk)); 
   ecapa		    = LamdaA/(fyieldPos/(elstk)); 
   ecapd		    = LamdaD/(fyieldPos/(elstk)); 
   cs		        = Cs; 
   ck		        = Ck; 
   ca		        = Ca; 
   cd		        = Cd; 
   capDispPos     = Thetap_pos+fyieldPos/elstk; 
   capDispNeg     = -Thetap_neg+fyieldNeg/elstk; 
   Resfac         = K; 
   ResfacNeg      = KNeg;                   
   myFcapping     =-(fyieldPos+alpha*elstk*capDispPos); 
   capSlope	    = myFcapping/(Thetapc_pos*elstk); 
   myFcappingNeg  = -(fabs(fyieldNeg)+alpha*elstk*fabs(capDispNeg)); 
   capSlopeNeg    = myFcappingNeg/(Thetapc_neg*elstk); 
   fracDispPos    = Thetau_pos; 
   fracDispNeg    = -Thetau_neg; 
   DPlus          = PDPlus;                   
   DNeg           = PDNeg;                                                                  
   stif = elstk; 
   ekP  = elstk;
   flagControlResponse = 0;    
   Tangent=Ke;     
   //stop
   ekhard = elstk*alpha;
   alphaNeg=alphaN;
   alphaPos=alpha;
   ekhardPos = ekhard;
   ekhardNeg = ekhard;
   dyieldPos = fyieldPos/elstk;
   dyieldNeg = fyieldNeg/elstk;
   Enrgts = fyieldPos*dyieldPos*ecaps;
   Enrgtk = fyieldPos*dyieldPos*ecapk;
   Enrgtd = fyieldPos*dyieldPos*ecapd;
   dmax = dyieldPos;
   dmin = dyieldNeg; 
   ekP = elstk;
   ekunload = elstk;
   ekexcurs = elstk;
   Enrgtot = 0.0;
   Enrgc = 0.0;
   fyPos = fyieldPos;
   fyNeg = fyieldNeg;
   dyNeg=dyieldNeg;
   dyPos=dyieldPos;
   resSn = fyieldPos;
   resSp = fyieldNeg;
   cpPos = capDispPos;
   cpNeg = capDispNeg;
   fPeakPos=fyieldPos+ekhard*(capDispPos-dyieldPos);
   fPeakNeg=fyieldNeg+ekhard*(capDispNeg-dyieldNeg);
   flagControlResponse = flagControlResponse;      // This flag Controls the ultimate rotation 
   DPlus = DPlus;                                  // Slab Effect Positive Direction
   DNeg = DNeg;                                    // Slab Effect Negative Direction
   if (cpPos<dyieldPos) {
     fPeakPos =fyieldPos*cpPos/dyieldPos;
   }
   if (cpNeg>dyieldNeg) {
     fPeakNeg =fyieldNeg*cpNeg/dyieldNeg;
   }
   fCapPos = fPeakPos;
   fCapNeg = fPeakNeg;
   LP = 0;
   LN = 0;
   fLimPos = 0;
   fLimNeg = 0;
   dLimPos = 0;
   dLimNeg = 0;
   dBoundPos = 0; 
   dBoundNeg = 0; 
   iNoFpos = 0;
   iNoFneg = 0;
   interup=0;
   fCapRefPos=-capSlope*elstk*capDispPos+fPeakPos;
   fCapRefNeg=-capSlopeNeg*elstk*capDispNeg+fPeakNeg;
   capSlopeOrig = capSlope;
   capSlopeOrigNeg = capSlopeNeg;
   flagstopdeg	= 0;
   dCap2Pos = 0.0;
   dCap1Pos = 0.0;
   dCap2Neg = 0.0;
   dCap1Neg = 0.0;
   dyPos = 0.0;
   dyNeg = 0.0;
   f1 = 0.0;  
   f2 = 0.0;
   fmin = 0.0;
   fmax = 0.0;
   RSE = 0.0;
   fn = 0.0;
   dn = 0.0;
   e1 = 0.0;
   e2 = 0.0;
   Enrgi = 0.0;
   ekt = 0.0;
   a2 = 0.0;
   NoCollapse = 0;
   iDeg = 0;
   if(deltaD>=0.0){
     kon = 1;
   } else {
     kon = 2;
   }
   kst = 1;
 }
 ////
 // 	******************* S T A R T S   B I G   L O O P  ****************
 // 	IF   D E L T A > 0 - - - - - - - - - - - - - - - - - - - - - - - -  

 if (deltaD>=0.0) {
   if (iNoFpos==1) {
     interPoint(dNewLoadPos,fNewLoadPos,dyNeg,fyNeg,ekhardNeg,0.0,0.0,0.0);
   }
   //If there is a corner changing delta from negative to positive      
   if (kon==2) {
     kon = 1;
     dlstNeg = dP;
     flstNeg = fP;
     
     if (dlstNeg<=dmin){
       fLimNeg = flstNeg;
       dLimNeg = dlstNeg;
     }
     
     if (resSn>0){
       RSE = 0.5*fP*fP/ekunload;
     } else {
       RSE=0.5*(fP+resSn)*(sn-dP)+0.5*resSn*resSn/(elstk*alphaPos);
     }
     
     a2 = Enrgtk-(Enrgtot-RSE);
     if((a2<=0.0) && (Enrgtk!=0.0)) {
       flgstop=1;
     }
     
     if((ecapk!=0.0) && (d<sp) && fabs(capSlope)>=1.0e-3 && (fabs(capSlopeNeg)>=1.0e-3)){
       betak = pow(((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE))),ck);
       if(((Enrgtot-RSE)>=Enrgtk)||(betak>=1.0)) {
	 betak = 1.0;
       }
       if( flagstopdeg !=1 && flagdeg !=1) {	        
	 ekunload = ekexcurs*(1-betak);
       } else {
	 ekunload = ekexcurs;
       }
       if(ekunload<=0.1*elstk) {
	 ekunload = 0.1*elstk;
       }
     }
     
     
     //// sn CALCULATION-------------------------------------------------
     
     if((dmin<dyieldNeg)||(dmax>dyieldPos)) {
       if((dP<sp)||(((dP>sp)&&(ekP==ekunload)))) {
	 
	 // 			call snCalc(sn,resSn,snEnv,resSnEnv,dP,fP,ekunload,alphaPos,dyPos,fyPos,cpPos,fCapPos,capSlope,fCapRefPos,LP,dLimPos,fLimPos,snHor,resSnHor,elstk,Resfac)
	 snCalc();
	 
	 if((fabs(dmax-dyieldPos)>=1.0e-10)&&(fabs(sn-dyieldPos)<=1.0e-10)) {
	   sn=dyieldPos-1.0e-9;
	 }
       }
     }
     
     if (sn>dmax) {
       dmax = sn;
     }
     
     if ((iNoFneg==1)&&(dP<=dNewLoadNeg)) {
       sn = dNewLoadNeg - 1.0e-6;
       resSn = 0;
     }
     
   }
   // 	LOADING ----------------------------------------------------------
   //      Push envelope
   // Compute force and stiffness for this increment
   if ((iNoFneg==1)&&(iNoFpos==1)&&(d<fracDispPos)) {
     f = Resfac*fyieldPos;
     ek = 1.0e-7;
     jcode = 8;
   } else if((iNoFneg==1)&&(iNoFpos==1)&&(d>=fracDispPos)||flagControlResponse==1) {
     f = 1.0e-10;
     ek = 1.0e-7;
     jcode = 9;
     flagstopdeg = 1;
   } else if ((iNoFpos==1)&&(d>sn)&&(d<fracDispPos)) {
     f = Resfac*fyieldPos;
     ek = 1.0e-7;
     jcode = 8;
   } else if ((iNoFpos==1)&&(d>sn)&&(d>=fracDispPos)){
     f = 1.0e-10;
     ek =1.0e-7;
     jcode = 9;
     flagstopdeg	= 1;
   } else if (d>=fracDispPos || flagControlResponse == 1){
     f = 1.0e-10;
     ek =1.0e-7;
     jcode = 9;
     flagstopdeg = 1;
     flagControlResponse = 1;                 // Switches the response to zero strength
   } else if ((iNoFpos==1)&&(d<sn)) {
     ek = ekunload;
     f  = fP+ek*deltaD;
     jcode = 2;	
   } else if ((iNoFneg==1)&&(d<dNewLoadNeg)){
     f = 0;
     ek = 1.0e-7;
     jcode = 8; 
   } else if (d>dmax) {
     //  		call envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,f,ek,elstk,fyieldPos,Resfac)
     f=0.0;
     envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,f,ek,elstk,fyieldPos,Resfac);
     dmax = d;
     fmax = f;
     jcode = 5;
     if(ek==elstk*capSlope) {
       jcode = 7;
     }
     if(ek<=1.0e-7) {
       jcode = 8;
       iCapPos = 1;
     }
     
     // c		COMPUTE MAXIMUM POSSIBLE DISPLACEMENT
     // 		call boundPos (dBoundPos,fCapRefPos,capSlope,fyNeg,dyNeg,alphaNeg,fCapPos,cpPos,elstk)	
     dBoundPos=boundPos();
     if((d>dBoundPos)||(ek==1.0e-7)) {
       iNoFpos = 1;
     }
     
   } else if (fabs(sn)>1.0e-10) {
     
     
     
     if (LP==0) {
       
       
       if (cpPos<=dyPos) {
	 // 				call interPoint(xDevPos1,yDevPos1,sn,resSn,ekhardPos,cpPos,fCapPos,capSlope*elstk)
	 interPoint(xDevPos1,yDevPos1,sn,resSn,ekhardPos,cpPos,fCapPos,capSlope*elstk);
	 // 				call interPoint(xDevPos2,yDevPos2,sn,resSn,ekunload,cpPos,fCapPos,capSlope*elstk)
	 interPoint(xDevPos2,yDevPos2,sn,resSn,ekunload,cpPos,fCapPos,capSlope*elstk);
	 if (xDevPos1>xDevPos2) {
	   xDevPos = xDevPos1;
	   yDevPos = yDevPos1;
	 } else {
	   xDevPos = xDevPos2;
	   yDevPos = yDevPos2;
	 }		
	 
	 if ((d<=sn)&&(d<=xDevPos)) {
	   ek = ekunload;
	   f  = fP+ek*deltaD;
	   jcode = 2;
	 } else if ((d>sn)&&(d<=xDevPos)) {
	   ek = elstk*alphaPos;
	   f2 = resSn+ek*(d-sn);
	   f1 = fP+ekunload*deltaD ; 
	   //f = min(f1,f2);
	   if (f1<f2)
	     {
	       f=f1;
	     }
	   else
	     {
	       f=f2;
	     }
	   
	   jcode = 5;
	   
	 } else if (d>xDevPos) {
	   ek = capSlope*elstk;
	   f2 = yDevPos+ek*(d-xDevPos);
	   f1 = fP+ekunload*deltaD  ;
	   //f = min(f1,f2); 
	   if (f1<f2)
	     {
	       f=f1;
	     }
	   else
	     {
	       f=f2;
	     }
	   if (ek!=ekunload) {
	     //                         call envHitsZero(f,fP,ek)
	     envHitsZero(f);
	   }
	   iCapPos = 1;
	   jcode = 7;
	 } else {
	 }
         
       } else if (cpPos > dyPos) {
	 interPoint(xDevPos1,yDevPos1,sn,resSn,ekhardPos,cpPos,fCapPos,capSlope*elstk);
	 if(d<=sn && d<=xDevPos1) { // Unloading in positive loading direction
	   ek = ekunload;
	   f  = fP+ek*deltaD;
	   jcode = 2;
	 } else if ((d>sn)&&(d<=cpPos)&& (d<=xDevPos1)) {
	   ek = elstk*alphaPos;
	   f2 = resSn+ek*(d-sn);
	   f1 = fP+ekunload*deltaD;
	   //f = min(f1,f2); 
	   if (f1<f2)
	     {
	       f=f1;
	     }
	   else
	     {
	       f=f2;
	     }
	   if (fabs(f-f1)<1.0e-10) {
	     ek=ekunload;
	   }
	   jcode = 5;
	 } else if((d>cpPos)&&(d<fracDispPos)) {
	   ek = capSlope*elstk;
	   f2 = fCapPos+ek*(d-cpPos);
	   f1 = fP+ekunload*deltaD;  
	   //f = min(f1,f2);
	   if (f1<f2)
	     {
	       f=f1;
	     }
	   else
	     {
	       f=f2;
	     }
	   if (fabs(f-f1)<1.0e-10) {
	     ek=ekunload;
	   }
	   if (ek!=ekunload) {
	     //                         call envHitsZero(f,fP,ek)
	     envHitsZero(f);
	   }
	   iCapPos = 1;
	   jcode = 7;
	   // c added by Dimitrios to simulate ductile fracture
	 } else if(d>=fracDispPos){
	   f=1.0e-10;
	   ek = -1.0e-7;
	   jcode = 9;
	   flagstopdeg = 1;
	 }
	 // Added by Dimitrios to stay on the residual path when  dresidual <d< dfracture               
	 if(d < fracDispPos && d > cpPos+(Resfac*fyieldPos-fCapPos)/(capSlope*elstk)) {
	   f = Resfac*fyieldPos;
	   ek = -1.0e-7;
	   jcode = 9;
	 }
	 if( flagControlResponse == 1) { // Response has to be zero since post capping regions hit zero
	   f = 1.0e-10;
	   jcode = 10;                       // this is for DRAIN output path tracking
	 }                
       }
       
       // c		IF LP IS EQUAL TO 1
     } else if (LP==1) {
       if(d<=sn) {
	 ek = ekunload;
	 f  = fP+ek*deltaD;
	 jcode = 2;
       } else if (((d>sn)&&(sn==snEnv)&&(d<=snHor))||((iNoFneg==1)&&(d>sn)&&(d<snHor))) {
	 ek = elstk*alphaPos;
	 f2 = resSn+ek*(d-sn);
	 f1 = fP+ekunload*deltaD;  
	 //f = min(f1,f2); 
	 if (f1<f2)
	   {
	     f=f1;
	   }
	 else
	   {
	     f=f2;
	   }
	 if (fabs(f-f1)<1.0e-10) {
	   ek=ekunload;
	 }
	 jcode = 5;
	 
       } else {
	 ek = 0;
	 f1 = fP+ekunload*deltaD;  
	 f2 = fLimPos;
	 //f = min(f1,f2); 
	 if (f1<f2)
	   {
	     f=f1;
	   }
	 else
	   {
	     f=f2;
	   }
	 if (fabs(f-f1)<1.0e-10) {
	   ek=ekunload;
	 }
	 jcode = 8;
       }
     } else {
     }
     // c	Elastic
   } else {
     if (d>0.0) {
       // 		    call envelPosCap2 (fyPos,alphaPos,capSlope,cpPos,d,f,ek,elstk,fyieldPos,Resfac)
       // [d,f,ek]=envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,0.0,ek,elstk,fyieldPos,Resfac);
       f=0.0;
       envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,f,ek,elstk,fyieldPos,Resfac);
     } else {
       // 		     call envelNegCap2 (fyNeg,alphaNeg,capSlope,cpNeg,d,f,ek,elstk,fyieldNeg,Resfac)
       f=0.0;
       envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,d,f,ek,elstk,fyieldNeg,ResfacNeg);
     }
     jcode = 0;
     if(ek==elstk*alphaPos) {
       jcode = 1;
     }
     if(ek==elstk*capSlope) {
       jcode = 7;
     }
     if(ek<=1.0e-7) {
       jcode = 8;
     }
   }
   
   //// c	IF   D E L T A < 0 - - - - - - - - - - - - - - - - - - - - - - - -  
   
 } else {
   
   if (iNoFneg==1) {
     // 		call interPoint(dNewLoadNeg,fNewLoadNeg,dyPos,fyPos,ekhardPos,0.d0,0.d0,0.d0)
     interPoint(dNewLoadNeg,fNewLoadNeg,dyPos,fyPos,ekhardPos,0.0,0.0,0.0);
   }
   
   
   // c	If there is a corner changing delta from positive to negative ---      
   
   if (kon==1) {
     kon = 2;
     
     dlstPos = dP;
     flstPos = fP;
     
     if (dlstPos>=dmax) {
       fLimPos = flstPos;
       dLimPos = dlstPos;
     }
     
     if (resSp<0) {
       RSE = 0.5*fP*fP/ekunload;
     } else {
       RSE=0.5*(fP+resSn)*(dP-sp)+0.5*resSp*resSp/(elstk*alphaNeg);
     }
     
     a2 = Enrgtk-(Enrgtot-RSE);
     
     
     if((a2<=0.0)&&(Enrgtk!=0.0)) {
       flgstop=1;
     }
     //// Update the Unloading Stiffness deterioration
     if((ecapk!=0.0)&&(d>sn)&&(flagstopdeg!=1)&&(flagdeg!=1)) {
       betak = pow(((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE))),ck);
       if(((Enrgtot-RSE)>=Enrgtk)||(betak>=1.0)) {
                betak = 1.0;
       }
       // If Post caping slopes have not been flat due to residual update stiffness deterioration
       if(flagdeg !=1 || flagstopdeg!=1) {	        
	 ekunload = ekexcurs*(1-betak);
       } else { // Keep the same unloading stiffness
	 ekunload = ekexcurs;
       }            
       if(ekunload<=0.1*elstk) {
	 ekunload = 0.1*elstk;
       }
     }
     
     
     ////     sp CALCULATION----------------------------------------------------
     
     if((dmin<dyieldNeg)||(dmax>dyieldPos)) {
       if((dP>sn)||(((dP<sn)&&(ekP==ekunload)))) {
	 
	 // 			call spCalc(sp,resSp,spEnv,resSpEnv,dP,fP,ekunload,alphaNeg,dyNeg,fyNeg,cpNeg,fCapNeg,capSlope,fCapRefNeg,LN,dLimNeg,fLimNeg,spHor,resSpHor,elstk,Resfac)
	 spCalc();
	 if((fabs(dmin-dyieldNeg)>=1.0e-10)&&(fabs(sp-dyieldNeg)<=1.0e-10)) {
	   sp=dyieldNeg-1.0e-9;
	 }
       }
     }
     
     if (sp<dmin) {
       dmin = sp;
     }
     
     if ((iNoFpos==1)&&(dP>=dNewLoadPos)) {
       sp = dNewLoadPos + 1.0e-6;
       resSp = 0;
     }
   }
   
   ////	UNLOADING
   // c     Push envelope
   
   if ((iNoFneg==1)&&(iNoFpos==1)&&(d>fracDispNeg)) {
     f = ResfacNeg*fyieldNeg;
     ek = 1.0e-7;
     jcode = 8;
   } else if ((iNoFneg==1)&&(iNoFpos==1)&&(d<=fracDispNeg)||flagControlResponse==1){
     f = -1.0e-10;
     ek =1.0e-7;
     jcode = 9;
     flagstopdeg = 1;
   } else if ((iNoFneg==1)&&(d<sp)&&(d>fracDispNeg)) {
     f = ResfacNeg*fyieldNeg;
     ek = 1.0e-7;
     jcode = 8;
   } else if ((iNoFneg==1)&&(d<sp)&&(d<=fracDispNeg)) {
     f = -1.0e-10;
     ek = 1.0e-7;
     jcode = 9;
     flagstopdeg = 1;
   } else if (d<=fracDispNeg || flagControlResponse ==1){
     f = -1.0e-10;
     ek = 1.0e-7;
     jcode = 9;
     flagstopdeg = 1;
     flagControlResponse = 1;               // To control response after I exceed fracture
   } else if ((iNoFneg==1)&&(d>=sp)) {
     ek = ekunload;
     f  = fP+ek*deltaD;
     jcode = 2;
   } else if ((iNoFpos==1)&&(d>dNewLoadPos)) {
     f = 0;
     ek = 1.0e-7;
     jcode = 8;
   } else if (d<dmin) {
     // 		call envelNegCap2(fyNeg,alphaNeg,capSlope,cpNeg,d,f,ek,elstk,fyieldNeg,Resfac)
     f=0.0;
     envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,d,f,ek,elstk,fyieldNeg,ResfacNeg);
     dmin = d;
     fmin = f;
     
     jcode = 5;
     if(ek==elstk*capSlopeNeg) {
       jcode = 7;
     }
     if(ek<=1.0e-7) {
       jcode = 8;
       iCapNeg = 1;
     }
     
     //// c		COMPUTE MINIMUM POSSIBLE DISPLACEMENT
     // 		call boundNeg (dBoundNeg,fCapRefNeg,capSlope,fyPos,dyPos,alphaPos,fCapNeg,cpNeg,elstk)	
     dBoundNeg=boundNeg();
     if((d<dBoundNeg)||(ek==1.0e-7)) {
       iNoFneg = 1;
     }
     
   } else if (fabs(sp)>1.0e-10){
     
     // c		If LN is equal to zero
     if (LN==0) {
       
       if (cpNeg>=dyNeg) {
	 // 				call interPoint(xDevNeg1,yDevNeg1,sp,resSp,elstk*alphaNeg,cpNeg,fCapNeg,capSlope*elstk)
	 interPoint(xDevNeg1,yDevNeg1,sp,resSp,elstk*alphaNeg,cpNeg,fCapNeg,capSlopeNeg*elstk);
	 // call interPoint(xDevNeg2,yDevNeg2,sp,resSp,ekunload,cpNeg,fCapNeg,capSlope*elstk)
	 interPoint(xDevNeg2,yDevNeg2,sp,resSp,ekunload,cpNeg,fCapNeg,capSlopeNeg*elstk);
	 if (xDevNeg1<xDevNeg2) {
	   xDevNeg = xDevNeg1;
	   yDevNeg = yDevNeg1;
	 } else {
	   xDevNeg = xDevNeg2;
	   yDevNeg = yDevNeg2;
	 }	
	 // 			
	 if ((d>=sp)&&(d>=xDevNeg)) {
	   ek = ekunload;
	   f  = fP+ek*deltaD;
	   jcode = 2;
	 } else if ((d<sp)&&(d>=xDevNeg)) {
	   ek = elstk*alphaNeg;
	   f2 = resSp+ek*(d-sp);
	   f1 = fP+ekunload*deltaD;
	   //f = max(f1,f2);
	   if (f1>f2)
	     {
	       f=f1;
	     }
	   else
	     {
	       f=f2;
	     }
	   if (fabs(f-f1)<1.0e-10) {
	     ek=ekunload;
	   }
	   jcode = 5;
				} else if (d<xDevNeg) {
					ek = capSlopeNeg*elstk;
					f2 = yDevNeg+ek*(d-xDevNeg);
					f1 = fP+ekunload*deltaD; 
					//f = max(f1,f2);
                    if (f1>f2)
					{
						f=f1;
					}
					else
					{
						f=f2;
					}
					if (fabs(f-f1)<1.0e-10) {
                        ek=ekunload;
                    }
					if (ek!=ekunload) {
//                         call envHitsZero(f,fP,ek)
                       envHitsZero(f);
                    }
					iCapNeg = 1;
					jcode = 7;
//                } else
                }


       } else if (cpNeg<dyNeg) {
	 interPoint(xDevNeg1,yDevNeg1,sp,resSp,elstk*alphaNeg,cpNeg,fCapNeg,capSlopeNeg*elstk);

	 if ( d>=sp && d>=xDevNeg1) {
	   ek = ekunload;
	   f = fP+ek*deltaD;
	   jcode = 2;
	 } else if((d<sp)&&(d>=cpNeg)&& (d>=xDevNeg1)) {
	   ek = elstk*alphaNeg;
	   f2 = resSp+ek*(d-sp);
	   f1 = fP+ekunload*deltaD;  
	   //f=max(f1,f2);
	   if (f1>f2)
	     {
	       f=f1;
	     }
	   else
					{
						f=f2;
					}
					if (fabs(f-f1)<1.0e-10) {
                        ek=ekunload;
                    }
					jcode = 5;
				} else if((d<cpNeg)&&(d>fracDispNeg)) {
					ek = capSlopeNeg*elstk;
					f2 = fCapNeg+ek*(d-cpNeg);
 					f1 = fP+ekunload*deltaD;  
					//f = max(f1,f2);
					if (f1>f2)
					{
						f=f1;
					}
					else
					{
						f=f2;
					}
					if (fabs(f-f1)<1.0e-10) {
                        ek=ekunload;
                    }
					if (ek!=ekunload) {
//                         call envHitsZero(f,fP,ek)
                       envHitsZero(f);
                    }
					iCapNeg = 1;
					jcode = 7;
// Added by Dimitris to model ductile tearing. Once theta_u(fracDispNeg) is exceeded Strength drops to zero and stays
				} else if(d<=fracDispNeg ||  flagControlResponse == 1){
	 			    ek = -1.0e-7;
	                f = 1.0e-10;
	                jcode = 9;
					flagstopdeg = 1;
                    flagControlResponse = 1;               // To dictate the response after passing fracture             
                }
// Added by Dimitrios to stay on the residual path when  dresidual <d- < dfracture               
				if(d > fracDispNeg && d < cpNeg+(ResfacNeg*fyieldNeg-fCapNeg)/(capSlopeNeg*elstk)) {
                    f = ResfacNeg*fyieldNeg;
                    ek = -1.0e-7;
                    jcode = 9;
                }
				if( flagControlResponse == 1) { // Response has to be zero since post capping regions hit zero
                    f = 1.0e-10;
                    jcode = 10;                       // this is for DRAIN output path tracking
                }
            }

		} else if (LN==1){
			if(d>=sp){
				ek = ekunload;
				f  = fP+ek*deltaD;
				jcode = 2;
			} else if (((d<sp)&&(sp==spEnv)&&(d>spHor))||((iNoFpos==1)&&(d<sp)&&(d>spHor))) {
				ek = elstk*alphaNeg;
				f2 = resSp+ek*(d-sp);
				f1 = fP+ekunload*deltaD;
				//f = max(f1,f2);
				if (f1>f2)
					{
						f=f1;
					}
					else
					{
						f=f2;
					}
				if (fabs(f-f1)<1.0e-10) {
                    ek=ekunload;
                }
				jcode = 5;
			} else {
				ek = 1.0e-7;
				f1 = fP+ekunload*deltaD ; 
				f2 = fLimNeg;
				//f = max(f1,f2);
				if (f1>f2)
					{
						f=f1;
					}
					else
					{
						f=f2;
					}
				if (fabs(f-f1)<1.0e-10) {
                    ek=ekunload;
                }
				jcode = 8;
            }
		} else {
        }
	  } else {
		  if (d>0.0) {
// 			call envelPosCap2 (fyPos,alphaPos,capSlope,cpPos,d,f,ek,elstk,fyieldPos,Resfac)
            //[d,f,ek]=envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,0.0,ek,elstk,fyieldPos,Resfac);
			  f=0.0;
        envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,f,ek,elstk,fyieldPos,Resfac);
		} else {
// 			call envelNegCap2 (fyNeg,alphaNeg,capSlope,cpNeg,d,f,ek,elstk,fyieldNeg,Resfac)
           envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,d,f,ek,elstk,fyieldNeg,ResfacNeg);
        }
          jcode = 0;
		  if(ek==elstk*alphaNeg) {
            jcode = 1;
        }
		  if(ek==elstk*capSlopeNeg) {
            jcode = 7;
        }
		  if(ek<=1.0e-7)  {
            jcode = 8;
        }
      }

}
// c	}S BIG LOOP ****************************************************
//// 	       
//    
    dn = d;
	fn = f;

	Enrgi = 0.5*(f+fP)*deltaD;

	Enrgc = Enrgc + Enrgi;
	Enrgtot = Enrgtot + Enrgi;

      RSE = 0.5*f*f/ekunload;

////	Flag to deteriorate parameters on the opposite side of the loop --	

	  if((f*fP<0.0)&&(interup==0)) {
		  if(((fP>0.0)&&(dmax>dyieldPos))||((fP<0.0)&&(dmin<dyieldNeg))) {
			flagdeg = 1;		
			interup=1;
        }
    }	


////	energy CALCULATIONS ---------------------------------------------

	  if((flagstopdeg==0)&&(flagdeg==1)) {
	    iDeg = iDeg +1;
        
		if((Enrgtot>=Enrgts)&&(Enrgts!=0.0)) {
			betas = 1.0;
		} else if((Enrgtot>=Enrgtd)&&(Enrgtd!=0.0)) {
			betad = 1.0;
		} else {
			if(ecaps!=0.0) {
                betas = pow((Enrgc/(Enrgts-Enrgtot)),cs);
             }
			if(ecapd!=0.0) {
                betad = pow((Enrgc/(Enrgtd-Enrgtot)),cd);
             }
	
			if(fabs(betas)>=1.0) {
                betas = 1.0;
             }
			if(fabs(betad)>=1.0) {
                betad = 1.0;
             }
        } 
////		Initialize energy of the cycle and Kstif for next loop -------

	    Enrgc = 0.0;
		ekexcurs = ekunload;

////		Check the horizontal limit for subsequent cycles -------------
		if ((dmin<cpNeg)&&(LN==1)) {
                iCapNeg = 1;
            }
		if ((dmax>cpPos)&&(LP==1)) {
                iCapPos = 1;
            }


////		Deteriorate parameters for the next half cycle
		if(deltaD<0.0){
		
			if(flagstopdeg == 0){
                fyNeg = fyNeg*(1-betas*DNeg);
                alphaNeg=alphaNeg*(1-betas*DNeg);
                fCapRefNeg=fCapRefNeg*(1-betad*DNeg);
			}else{
                fyNeg = fyNeg;
                alphaNeg=alphaNeg;
                fCapRefNeg=fCapRefNeg;
			}
            // When we reach post capping slope goes to zero due to residual
			if(fyNeg>=ResfacNeg*fyieldNeg) { // If strength drops below residual
				fyNeg = ResfacNeg*fyieldNeg;
				alphaNeg = 10^(-4);
                fCapRefNeg = fyNeg;
                capSlopeNeg = -pow(10.0,-6);
                flagstopdeg = 1;
			} else { //% Keep updating the post capping slope
                capSlopeNeg = capSlopeOrigNeg*(1-fabs((ResfacNeg*fyieldNeg)/fyNeg));
				if(capSlopeNeg >=0){
                    capSlopeNeg = -pow(10.0,-6);
				}
			}

			dyNeg = fyNeg/elstk;
			ekhardNeg=alphaNeg*elstk;

			dCap1Neg=fCapRefNeg/(elstk-capSlopeOrigNeg*elstk);
			dCap2Neg=(fCapRefNeg+ekhardNeg*dyNeg-fyNeg)/(ekhardNeg-capSlopeOrigNeg*elstk);
			//cpNeg=min(dCap1Neg,dCap2Neg);
            if (dCap1Neg<dCap2Neg)
					{
						cpNeg=dCap1Neg;
					}
					else
					{
						cpNeg=dCap2Neg;
					}

			fCapNeg = fCapRefNeg + capSlopeOrigNeg*elstk*cpNeg;
            
// 			call envelNegCap2 (fyNeg,alphaNeg,capSlope,cpNeg,dLimNeg,fLimNeg,ekt,elstk,fyieldNeg,Resfac)
           envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,dLimNeg,fLimNeg,ekt,elstk,fyieldNeg,ResfacNeg);
// 			call spCalc(sp,resSp,spEnv,resSpEnv,dP,fP,ekunload,alphaNeg,dyNeg,fyNeg,cpNeg,fCapNeg,capSlope,fCapRefNeg,LN,dLimNeg,fLimNeg,spHor,resSpHor,elstk,Resfac)
           spCalc();
// c			In case the degradation point moves from the negative to pos. side
			if(resSp>f) {        
                d = sp;
				f = resSp;
				ek = ekhardNeg;   
            }
		} else {
			if(flagstopdeg == 0){
                fyPos = fyPos*(1-betas*DPlus);	
                alphaPos=alphaPos*(1-betas*DPlus);	
                fCapRefPos=fCapRefPos*(1-betad*DPlus);
			} else {
                fyPos = fyPos;	
                alphaPos=alphaPos;	
                fCapRefPos=fCapRefPos;
			}
                
     //   %If post capping slope goes to zero due to residual:
			if(fyPos <= Resfac*fyieldPos) {  //% If yield Strength Pos drops below residual
                fyPos = Resfac*fyieldPos;
				alphaPos = pow(10.0,-4);
                fCapRefPos = fyPos;
                capSlope = -pow(10.0,-6);
                flagstopdeg = 1;              
			}  else { //% keep updating
		    capSlope = capSlopeOrig*(1-fabs((Resfac*fyieldPos)/fyPos));
			if(capSlope >=0) {
                capSlope = -pow(10.0,-6);
			}
			}
			dyPos = fyPos/elstk;
			ekhardPos=alphaPos*elstk;

			dCap1Pos=fCapRefPos/(elstk-capSlopeOrig*elstk);
			dCap2Pos=(fCapRefPos+ekhardPos*dyPos-fyPos)/(ekhardPos-capSlopeOrig*elstk);
			//cpPos=max(dCap1Pos,dCap2Pos);
            if(dCap1Pos>dCap2Pos)
			{
				cpPos=dCap1Pos;
			}
			else
			{
				cpPos=dCap2Pos;
			}
			fCapPos = fCapRefPos + capSlopeOrig*elstk*cpPos;

			//call envelPosCap2 (fyPos,alphaPos,capSlope,cpPos,dLimPos,fLimPos,ekt,elstk,fyieldPos,Resfac)           
           envelPosCap2(fyPos,alphaPos,capSlope,cpPos,dLimPos,fLimPos,ekt,elstk,fyieldPos,Resfac);
        
			//call snCalc(sn,resSn,snEnv,resSnEnv,dP,fP,ekunload,alphaPos,dyPos,fyPos,cpPos,fCapPos,capSlope,fCapRefPos,LP,dLimPos,fLimPos,snHor,resSnHor,elstk,Resfac)
            snCalc();
// c			In case the degradation point moves from the pos. to neg. side
			if(resSn<f) {
				d = sn;
				f = resSn;
				ek = ekhardPos;

            }	
        }
}

// c		Check the horizontal limit in case that dBound is reached after first neg slope
if ((d<0)&&(fabs(ek)<=1.0e-7)) {
                LN = 1;
      }
if ((d>0)&&(fabs(ek)<=1.0e-7)) {
                LP = 1;
      }


// c	Calculation of static resisting force vector for the element -----
   
    double RestoringForce=f;


	


// c	Update envelope values --------------------------------------------	

      f1 = ftot;	
      v1 = vtot;
      d1 = dtot;

	ftot = fn;
	dtot = dn;
	vtot = dn;



  

//// c	Change in Energies -----------------------------------------------

// c	Hysteretic energy
	ener = Enrgi; 	
	
// c	Updating parameters for next cycle ---------------------------------
	stif = ek;
	ekP = ek;
    fP = f;
	dP = d;

	if(jcode!=kcode) {
        kst = 1;
    }
	kcode = jcode;

//priority of logical operators
	if (interup==1&&(ek==elstk*alphaPos) ||(ek==elstk*alphaNeg)||(ek==capSlope*elstk) ||(ek==capSlopeNeg*elstk)) {
		interup = 0;
    }
	


//*******************************} of matlab code*******************************************************

//Define force and tangent
  Force=RestoringForce;
  Tangent=stif;

  return 0;
}
double 
Bilin::getStress(void)
{
   return (Force);
}

double 
Bilin::getTangent(void)
{

    return (Tangent);
}

double 
Bilin::getInitialTangent(void)
{
    return (Ke);
}




double 
Bilin::getStrain(void)
{
    return (U);
}



int 
Bilin::commitState(void)
{
	//commit trial  variables
   //CU=TU; //displacement
   CU=U; 
   CForce=Force; 
   CTangent=Tangent; 

   CdNewLoadPos=dNewLoadPos;
   CdNewLoadNeg=dNewLoadNeg;
   Cflgstop=flgstop; 
   Cflagdeg=flagdeg ; 
   Cflagstopdeg=flagstopdeg; 
   Cekt=ekt;
   Cinterup=interup;  
   Ckcode=kcode; 
   Ckon=kon;
   CiCapNeg=iCapNeg; 
   CiNoFneg=iNoFneg; 
   CiNoFpos=iNoFpos; 
   CiCapPos=iCapPos; 
   CiDeg=iDeg; 
   CLP=CLP;
   CLN=LN; 
   CcapSlope=capSlope; 
   CcapDispPos=capDispPos; 
   CcapDispNeg=capDispNeg; 
   Celstk=elstk; 
   CfyieldPos=fyieldPos; 
   CfyieldNeg=fyieldNeg; 
   Calpha=alpha;  
   Cecaps=ecaps; 
   Cecapk=ecapk; 
   Cecapd=ecapd; 
   Ccs=cs; 
   Cck=ck; 
   Ccd=cd; 
   Cdmax=dmax; 
   Cdmin=dmin; 
   CEnrgtot=Enrgtot; 
   CEnrgc=Enrgc; 
   CfyPos=fyPos; 
   CfLimNeg=fLimNeg; 
   CfyNeg=fyNeg; 
   CekP=ekP; 
   Cekunload=ekunload; 
   Csp=sp; 
   Csn=sn; 
   CdP=dP; 
   CfP=fP; 
   Cek=ek; 
   Cstif=stif; 
   CdLimPos=dLimPos; 
   CdLimNeg=dLimNeg;  
   Cvtot=vtot; 
   Cftot=ftot; 
   Cdtot=dtot; 
   Cdn=dn; 
   CcpPos=cpPos; 
   CcpNeg=cpNeg; 
   CfLimPos=fLimPos; 
   CdlstPos=dlstPos; 
   CflstPos=flstPos; 
   CdlstNeg=dlstNeg; 
   CflstNeg=flstNeg; 
   Cekexcurs=ekexcurs; 
   CRSE=RSE; 
   CfPeakPos=fPeakPos; 
   CfPeakNeg=fPeakNeg; 
   CdCap1Pos=dCap1Pos; 
   CdCap2Pos=dCap2Pos; 
   CdCap1Neg=dCap1Neg; 
   CdCap2Neg=dCap2Neg; 
   CalphaNeg=alphaNeg; 
   CalphaPos=alphaPos; 
   CekhardNeg=ekhardNeg; 
   CekhardPos=ekhardPos; 
   CfCapRefPos=fCapRefPos; 
   CfCapRefNeg=fCapRefNeg; 
   CEnrgts=Enrgts; 
   CEnrgtk=Enrgtk; 
   CEnrgtd=Enrgtd; 
   CdyPos=dyPos; 
   CdyNeg=dyNeg; 
   CdyieldPos=dyieldPos; 
   CdyieldNeg=dyieldNeg; 
   CresSnHor=resSnHor; 
   Cfmax=fmax; 
   Cfmin=fmin; 
   CresSp=resSp; 
   CresSn=resSn; 
   CfCapPos=fCapPos; 
   CfCapNeg=fCapNeg; 
   CsnHor=snHor; 
   CspHor=spHor; 
   CresSpHor=resSpHor;
   CsnEnv=snEnv;
   CresSnEnv=resSnEnv; 
   CspEnv=spEnv; 
   CresSpEnv=resSpEnv; 
   CResfac=Resfac; 
   CcapSlopeOrig=capSlopeOrig;  
   CfracDispPos=fracDispPos; 
   CfracDispNeg=fracDispNeg; 
   CDPlus=DPlus; 
   CDNeg=DNeg; 
   CalphaN=alphaN;
   Cecapa=ecapa;
   Cca=ca;
   CResfacNeg=ResfacNeg;
   CmyFcapping=myFcapping;
   CmyFcappingNeg=myFcappingNeg;
   CcapSlopeNeg=capSlopeNeg;
   CflagControlResponse=flagControlResponse; 
   CcapSlopeOrigNeg=capSlopeOrigNeg;
   CUprev=Uprev;
   return 0;
}

int 
Bilin::revertToLastCommit(void)
{
	//the opposite of commit trial history variables
   U=CU; 
   Force=CForce; 
   Tangent=CTangent; 

   dNewLoadPos=CdNewLoadPos;
   dNewLoadNeg=CdNewLoadNeg;
   flgstop=Cflgstop; 
   flagdeg=Cflagdeg ; 
   flagstopdeg=Cflagstopdeg; 
   ekt=Cekt;
   interup=Cinterup;  
   kcode=Ckcode; 
   kon=Ckon;
   iCapNeg=CiCapNeg; 
   iNoFneg=CiNoFneg; 
   iNoFpos=CiNoFpos; 
   iCapPos=CiCapPos; 
   iDeg=CiDeg; 
   LP=CLP;
   LN=CLN; 
   capSlope=CcapSlope; 
   capDispPos=CcapDispPos; 
   capDispNeg=CcapDispNeg; 
   elstk=Celstk; 
   fyieldPos=CfyieldPos; 
   fyieldNeg=CfyieldNeg; 
   alpha=Calpha;  
   ecaps=Cecaps; 
   ecapk=Cecapk; 
   ecapd=Cecapd; 
   cs=Ccs; 
   ck=Cck; 
   cd=Ccd; 
   dmax=Cdmax; 
   dmin=Cdmin; 
   Enrgtot=CEnrgtot; 
   Enrgc=CEnrgc; 
   fyPos=CfyPos; 
   fLimNeg=CfLimNeg; 
   fyNeg=CfyNeg; 
   ekP=CekP; 
   ekunload=Cekunload; 
   sp=Csp; 
   sn=Csn; 
   dP=CdP; 
   fP=CfP; 
   ek=Cek; 
   stif=Cstif; 
   dLimPos=CdLimPos; 
   dLimNeg=CdLimNeg;  
   vtot=Cvtot; 
   ftot=Cftot; 
   dtot=Cdtot; 
   dn=Cdn; 
   cpPos=CcpPos; 
   cpNeg=CcpNeg; 
   fLimPos=CfLimPos; 
   dlstPos=CdlstPos; 
   flstPos=CflstPos; 
   dlstNeg=CdlstNeg; 
   flstNeg=CflstNeg; 
   ekexcurs=Cekexcurs; 
   RSE=CRSE; 
   fPeakPos=CfPeakPos; 
   fPeakNeg=CfPeakNeg; 
   dCap1Pos=CdCap1Pos; 
   dCap2Pos=CdCap2Pos; 
   dCap1Neg=CdCap1Neg; 
   dCap2Neg=CdCap2Neg; 
   alphaNeg=CalphaNeg; 
   alphaPos=CalphaPos; 
   ekhardNeg=CekhardNeg; 
   ekhardPos=CekhardPos; 
   fCapRefPos=CfCapRefPos; 
   fCapRefNeg=CfCapRefNeg; 
   Enrgts=CEnrgts; 
   Enrgtk=CEnrgtk; 
   Enrgtd=CEnrgtd; 
   dyPos=CdyPos; 
   dyNeg=CdyNeg; 
   dyieldPos=CdyieldPos; 
   dyieldNeg=CdyieldNeg; 
   resSnHor=CresSnHor; 
   fmax=Cfmax; 
   fmin=Cfmin; 
   resSp=CresSp; 
   resSn=CresSn; 
   fCapPos=CfCapPos; 
   fCapNeg=CfCapNeg; 
   snHor=CsnHor; 
   spHor=CspHor; 
   resSpHor=CresSpHor;
   snEnv=CsnEnv;
   resSnEnv=CresSnEnv; 
   spEnv=CspEnv; 
   resSpEnv=CresSpEnv; 
   Resfac=CResfac; 
   capSlopeOrig=CcapSlopeOrig;  
   fracDispPos=CfracDispPos; 
   fracDispNeg=CfracDispNeg; 
   DPlus=CDPlus; 
   DNeg=CDNeg; 
   alphaN=CalphaN;
   ecapa=Cecapa;
   ca=Cca;
   ResfacNeg=CResfacNeg;
   myFcapping=CmyFcapping;
   myFcappingNeg=CmyFcappingNeg;
   capSlopeNeg=CcapSlopeNeg;
   flagControlResponse=CflagControlResponse; 
   capSlopeOrigNeg=CcapSlopeOrigNeg;
   Uprev=CUprev;
    return 0;
}

int 
Bilin::revertToStart(void)
{

  U =0; CU =0; //displacement
  Force =0; CForce =0; //force
  Tangent =0; CTangent =0; //tangent stiffness

  //History variables
  dNewLoadPos =0;CdNewLoadPos =0;
  dNewLoadNeg =0; CdNewLoadNeg =0;
  flgstop =0;Cflgstop =0; 
  flagdeg =0;Cflagdeg  =0; 
  flagstopdeg =0;Cflagstopdeg =0; 
  ekt =0; Cekt =0;
  interup =0; Cinterup =0;  
  kcode =0; Ckcode =0; 
  kon =0; Ckon =0;
  iCapNeg =0; CiCapNeg =0; 
  iNoFneg =0; CiNoFneg =0; 
  iNoFpos =0; CiNoFpos =0; 
  iCapPos =0; CiCapPos =0; 
  iDeg =0; CiDeg =0; 
  LP =0; CLP =0;
  LN =0; CLN =0; 
  capSlope =0; CcapSlope =0; 
  capDispPos =0; CcapDispPos =0; 
  capDispNeg =0; CcapDispNeg =0; 
  elstk =0; Celstk =0; 
  fyieldPos =0; CfyieldPos =0; 
  fyieldNeg =0; CfyieldNeg =0; 
  alpha =0; Calpha =0;  
  ecaps =0; Cecaps =0; 
  ecapk =0; Cecapk =0; 
  ecapd =0; Cecapd =0; 
  cs =0; Ccs =0; 
  ck =0; Cck =0; 
  cd =0; Ccd =0; 
  dmax =0; Cdmax =0; 
  dmin =0; Cdmin =0; 
  Enrgtot =0; CEnrgtot =0; 
  Enrgc =0; CEnrgc =0; 
  fyPos =0; CfyPos =0; 
  fLimNeg =0; CfLimNeg =0; 
  fyNeg =0; CfyNeg =0; 
  ekP =0; CekP =0; 
  ekunload =0; Cekunload =0; 
  sp =0; Csp =0; 
  sn =0; Csn =0; 
  dP =0; CdP =0; 
  fP =0; CfP =0; 
  ek =0; Cek =0; 
  stif =0; Cstif =0; 
  dLimPos =0; CdLimPos =0; 
  dLimNeg =0; CdLimNeg =0;  
  vtot =0; Cvtot =0; 
  ftot =0; Cftot =0; 
  dtot =0; Cdtot =0; 
  dn =0; Cdn =0; 
  cpPos =0; CcpPos =0; 
  cpNeg =0; CcpNeg =0; 
  fLimPos =0; CfLimPos =0; 
  dlstPos =0; CdlstPos =0; 
  flstPos =0; CflstPos =0; 
  dlstNeg =0; CdlstNeg =0; 
  flstNeg =0; CflstNeg =0; 
  ekexcurs =0; Cekexcurs =0; 
  RSE =0; CRSE =0; 
  fPeakPos =0; CfPeakPos =0; 
  fPeakNeg =0; CfPeakNeg =0; 
  dCap1Pos =0; CdCap1Pos =0; 
  dCap2Pos =0; CdCap2Pos =0; 
  dCap1Neg =0; CdCap1Neg =0; 
  dCap2Neg =0; CdCap2Neg =0; 
  alphaNeg =0; CalphaNeg =0; 
  alphaPos =0; CalphaPos =0; 
  ekhardNeg =0; CekhardNeg =0; 
  ekhardPos =0; CekhardPos =0; 
  fCapRefPos =0; CfCapRefPos =0; 
  fCapRefNeg =0; CfCapRefNeg =0; 
  Enrgts =0; CEnrgts =0; 
  Enrgtk =0; CEnrgtk =0; 
  Enrgtd =0; CEnrgtd =0; 
  dyPos =0; CdyPos =0; 
  dyNeg =0; CdyNeg =0; 
  dyieldPos =0; CdyieldPos =0; 
  dyieldNeg =0; CdyieldNeg =0; 
  resSnHor =0; CresSnHor =0; 
  fmax =0; Cfmax =0; 
  fmin =0; Cfmin =0; 
  resSp =0; CresSp =0; 
  resSn =0; CresSn =0; 
  fCapPos =0; CfCapPos =0; 
  fCapNeg =0; CfCapNeg =0; 
  snHor =0; CsnHor =0; 
  spHor =0; CspHor =0; 
  resSpHor =0; CresSpHor =0;
  snEnv =0; CsnEnv =0;
  resSnEnv =0; CresSnEnv =0; 
  spEnv =0; CspEnv =0; 
  resSpEnv =0; CresSpEnv =0; 
  Resfac =0; CResfac =0; 
  capSlopeOrig =0;CcapSlopeOrig =0;  
  fracDispPos =0; CfracDispPos =0; 
  fracDispNeg =0; CfracDispNeg =0; 
  DPlus =0; CDPlus =0; 
  DNeg =0; CDNeg =0; 
  alphaN =0; CalphaN =0;
  ecapa =0; Cecapa =0;
  ca =0; Cca =0;
  ResfacNeg =0; CResfacNeg =0;
  myFcapping =0; CmyFcapping =0;
  myFcappingNeg =0; CmyFcappingNeg =0;
  capSlopeNeg =0; CcapSlopeNeg =0;
  flagControlResponse =0; CflagControlResponse =0; 
  capSlopeOrigNeg =0; CcapSlopeOrigNeg =0;
  Uprev =0;CUprev =0;

  //initially I zero everything
  /*
   U=CU=0.0; 
   Force=CForce=0.0; 
   Tangent=CTangent=0.0; 

   dNewLoadPos=CdNewLoadPos=0.0;
   dNewLoadNeg=CdNewLoadNeg=0.0;
   flgstop=Cflgstop=0; 
   flagdeg=Cflagdeg=0; 
   flagstopdeg=Cflagstopdeg=0; 
   ekt=Cekt=0.0;
   interup=Cinterup=0;  
   kcode=Ckcode=0; 
   kon=Ckon=0;
   iCapNeg=CiCapNeg=0; 
   iNoFneg=CiNoFneg=0; 
   iNoFpos=CiNoFpos=0; 
   iCapPos=CiCapPos=0; 
   iDeg=CiDeg=0; 
   LP=CLP=0;
   LN=CLN=0; 
   capSlope=CcapSlope=0.0; 
   capDispPos=CcapDispPos=0.0; 
   capDispNeg=CcapDispNeg=0.0; 
   elstk=Celstk=0.0; 
   fyieldPos=CfyieldPos=0.0; 
   fyieldNeg=CfyieldNeg=0.0; 
   alpha=Calpha=0.0;  
   ecaps=Cecaps=0.0; 
   ecapk=Cecapk=0.0; 
   ecapd=Cecapd=0.0; 
   cs=Ccs=0.0; 
   ck=Cck=0.0; 
   cd=Ccd=0.0; 
   dmax=Cdmax=0.0; 
   dmin=Cdmin=0.0; 
   Enrgtot=CEnrgtot=0.0; 
   Enrgc=CEnrgc=0.0; 
   fyPos=CfyPos=0.0; 
   fLimNeg=CfLimNeg=0.0; 
   fyNeg=CfyNeg=0.0; 
   ekP=CekP=0.0; 
   ekunload=Cekunload=0.0; 
   sp=Csp=0.0; 
   sn=Csn=0.0; 
   dP=CdP=0.0; 
   fP=CfP=0.0; 
   ek=Cek=0.0; 
   stif=Cstif=0.0; 
   dLimPos=CdLimPos=0.0; 
   dLimNeg=CdLimNeg=0.0;  
   vtot=Cvtot=0.0; 
   ftot=Cftot=0.0; 
   dtot=Cdtot=0.0; 
   dn=Cdn=0.0; 
   cpPos=CcpPos=0.0; 
   cpNeg=CcpNeg=0.0; 
   fLimPos=CfLimPos=0.0; 
   dlstPos=CdlstPos=0.0; 
   flstPos=CflstPos=0.0; 
   dlstNeg=CdlstNeg=0.0; 
   flstNeg=CflstNeg=0.0; 
   ekexcurs=Cekexcurs=0.0; 
   RSE=CRSE=0.0; 
   fPeakPos=CfPeakPos=0.0; 
   fPeakNeg=CfPeakNeg=0.0; 
   dCap1Pos=CdCap1Pos=0.0; 
   dCap2Pos=CdCap2Pos=0.0; 
   dCap1Neg=CdCap1Neg=0.0; 
   dCap2Neg=CdCap2Neg=0.0; 
   alphaNeg=CalphaNeg=0.0; 
   alphaPos=CalphaPos=0.0; 
   ekhardNeg=CekhardNeg=0.0; 
   ekhardPos=CekhardPos=0.0; 
   fCapRefPos=CfCapRefPos=0.0; 
   fCapRefNeg=CfCapRefNeg=0.0; 
   Enrgts=CEnrgts=0.0; 
   Enrgtk=CEnrgtk=0.0; 
   Enrgtd=CEnrgtd=0.0; 
   dyPos=CdyPos=0.0; 
   dyNeg=CdyNeg=0.0; 
   dyieldPos=CdyieldPos=0.0; 
   dyieldNeg=CdyieldNeg=0.0; 
   resSnHor=CresSnHor=0.0; 
   fmax=Cfmax=0.0; 
   fmin=Cfmin=0.0; 
   resSp=CresSp=0.0; 
   resSn=CresSn=0.0; 
   fCapPos=CfCapPos=0.0; 
   fCapNeg=CfCapNeg=0.0; 
   snHor=CsnHor=0.0; 
   spHor=CspHor=0.0; 
   resSpHor=CresSpHor=0.0;
   snEnv=CsnEnv=0.0;
   resSnEnv=CresSnEnv=0.0; 
   spEnv=CspEnv=0.0; 
   resSpEnv=CresSpEnv=0.0; 
   Resfac=CResfac=0.0; 
   capSlopeOrig=CcapSlopeOrig=0.0;  
   fracDispPos=CfracDispPos=0.0; 
   fracDispNeg=CfracDispNeg=0.0; 
   DPlus=CDPlus=0.0; 
   DNeg=CDNeg=0.0; 
   alphaN=CalphaN=0.0;
   ecapa=Cecapa=0.0;
   ca=Cca=0.0;
   ResfacNeg=CResfacNeg=0.0;
   myFcapping=CmyFcapping=0.0;
   myFcappingNeg=CmyFcappingNeg=0.0;
   capSlopeNeg=CcapSlopeNeg=0.0;
   flagControlResponse=CflagControlResponse=0; 
   capSlopeOrigNeg=CcapSlopeOrigNeg=0.0;
   Uprev=CUprev=0.0;
  */

   //then I initialize everything accordingly
   elstk          = Ke; 
   fyieldPos      = My_pos; 
   fyieldNeg      = My_neg; 
   alpha          = As; 
   alphaN         = AsNeg;  
   ecaps          = LamdaS/(fyieldPos/(elstk)); 
   ecapk          = LamdaK/(fyieldPos/(elstk)); 
   ecapa	  = LamdaA/(fyieldPos/(elstk)); 
   ecapd	  = LamdaD/(fyieldPos/(elstk)); 
   cs		  = Cs; 
   ck		  = Ck; 
   ca		  = Ca; 
   cd		  = Cd; 
   capDispPos     = Thetap_pos+fyieldPos/elstk; 
   capDispNeg     = -Thetap_neg+fyieldNeg/elstk; 
   Resfac         = K; 
   ResfacNeg      = KNeg;                   
   myFcapping     =-(fyieldPos+alpha*elstk*capDispPos); 
   capSlope	    = myFcapping/(Thetapc_pos*elstk); 
   myFcappingNeg  = -(fabs(fyieldNeg)+alpha*elstk*fabs(capDispNeg)); 
   capSlopeNeg    = myFcappingNeg/(Thetapc_neg*elstk); 
   fracDispPos    = Thetau_pos; 
   fracDispNeg    = -Thetau_neg; 
   DPlus          = PDPlus;                   
   DNeg           = PDNeg;                                                                  
   stif = elstk; 
   ekP  = elstk;
   flagControlResponse = 0;    
   Tangent=Ke;                             
   return 0;
}

UniaxialMaterial *
Bilin::getCopy(void)
{
    Bilin *theCopy = new Bilin(this->getTag(),Ke,As,AsNeg,My_pos,My_neg,LamdaS,LamdaK,
			       LamdaA,LamdaD,Cs,Ck,Ca,Cd,Thetap_pos,Thetap_neg,Thetapc_pos,Thetapc_neg,K,KNeg,
			       Thetau_pos,Thetau_neg,PDPlus,PDNeg);

    //Fixed Model parameters: need to change according to material properties
   theCopy->U=U;
   theCopy->CU=CU; 
   theCopy->Force=Force;
   theCopy->CForce=CForce; 
   theCopy->Tangent=Tangent;
   theCopy->CTangent=CTangent; 
   theCopy->dNewLoadPos=dNewLoadPos;
   theCopy->CdNewLoadPos=CdNewLoadPos;
   theCopy->dNewLoadNeg=dNewLoadNeg;
   theCopy->CdNewLoadNeg=CdNewLoadNeg;
   theCopy->flgstop=flgstop;
   theCopy->Cflgstop=Cflgstop; 
   theCopy->flagdeg=flagdeg;
   theCopy->Cflagdeg=Cflagdeg; 
   theCopy->flagstopdeg=flagstopdeg;
   theCopy->Cflagstopdeg=Cflagstopdeg; 
   theCopy->ekt=ekt;
   theCopy->Cekt=Cekt;
   theCopy->interup=interup;
   theCopy->Cinterup=Cinterup;  
   theCopy->kcode=kcode;
   theCopy->Ckcode=Ckcode; 
   theCopy->kon=kon;
   theCopy->Ckon=Ckon;
   theCopy->iCapNeg=iCapNeg;
   theCopy->CiCapNeg=CiCapNeg; 
   theCopy->iNoFneg=iNoFneg;
   theCopy->CiNoFneg=CiNoFneg; 
   theCopy->iNoFpos=iNoFpos;
   theCopy->CiNoFpos=CiNoFpos; 
   theCopy->iCapPos=iCapPos;
   theCopy->CiCapPos=CiCapPos; 
   theCopy->iDeg=iDeg;
   theCopy->CiDeg=CiDeg; 
   theCopy->LP=LP;
   theCopy->CLP=CLP;
   theCopy->LN=LN;
   theCopy->CLN=CLN; 
   theCopy->capSlope=capSlope;
   theCopy->CcapSlope=CcapSlope; 
   theCopy->capDispPos=capDispPos;
   theCopy->CcapDispPos=CcapDispPos; 
   theCopy->capDispNeg=capDispNeg;
   theCopy->CcapDispNeg=CcapDispNeg; 
   theCopy->elstk=elstk;
   theCopy->Celstk=Celstk; 
   theCopy->fyieldPos=fyieldPos;
   theCopy->CfyieldPos=CfyieldPos; 
   theCopy->fyieldNeg=fyieldNeg;
   theCopy->CfyieldNeg=CfyieldNeg; 
   theCopy->alpha=alpha;
   theCopy->Calpha=Calpha;  
   theCopy->ecaps=ecaps;
   theCopy->Cecaps=Cecaps; 
   theCopy->ecapk=ecapk;
   theCopy->Cecapk=Cecapk; 
   theCopy->ecapd=ecapd;
   theCopy->Cecapd=Cecapd; 
   theCopy->cs=cs;
   theCopy->Ccs=Ccs; 
   theCopy->ck=ck;
   theCopy->Cck=Cck; 
   theCopy->cd=cd;
   theCopy->Ccd=Ccd; 
   theCopy->dmax=dmax;
   theCopy->Cdmax=Cdmax; 
   theCopy->dmin=dmin;
   theCopy->Cdmin=Cdmin; 
   theCopy->Enrgtot=Enrgtot;
   theCopy->CEnrgtot=CEnrgtot; 
   theCopy->Enrgc=Enrgc;
   theCopy->CEnrgc=CEnrgc; 
   theCopy->fyPos=fyPos;
   theCopy->CfyPos=CfyPos; 
   theCopy->fLimNeg=fLimNeg;
   theCopy->CfLimNeg=CfLimNeg; 
   theCopy->fyNeg=fyNeg;
   theCopy->CfyNeg=CfyNeg; 
   theCopy->ekP=ekP;
   theCopy->CekP=CekP; 
   theCopy->ekunload=ekunload;
   theCopy->Cekunload=Cekunload; 
   theCopy->sp=sp;
   theCopy->Csp=Csp; 
   theCopy->sn=sn;
   theCopy->Csn=Csn; 
   theCopy->dP=dP;
   theCopy->CdP=CdP; 
   theCopy->fP=fP;
   theCopy->CfP=CfP; 
   theCopy->ek=ek;
   theCopy->Cek=Cek; 
   theCopy->stif=stif;
   theCopy->Cstif=Cstif; 
   theCopy->dLimPos=dLimPos;
   theCopy->CdLimPos=CdLimPos; 
   theCopy->dLimNeg=dLimNeg;
   theCopy->CdLimNeg=CdLimNeg;  
   theCopy->vtot=vtot;
   theCopy->Cvtot=Cvtot; 
   theCopy->ftot=ftot;
   theCopy->Cftot=Cftot; 
   theCopy->dtot=dtot;
   theCopy->Cdtot=Cdtot; 
   theCopy->dn=dn;
   theCopy->Cdn=Cdn; 
   theCopy->cpPos=cpPos;
   theCopy->CcpPos=CcpPos; 
   theCopy->cpNeg=cpNeg;
   theCopy->CcpNeg=CcpNeg; 
   theCopy->fLimPos=fLimPos;
   theCopy->CfLimPos=CfLimPos; 
   theCopy->dlstPos=dlstPos;
   theCopy->CdlstPos=CdlstPos; 
   theCopy->flstPos=flstPos;
   theCopy->CflstPos=CflstPos; 
   theCopy->dlstNeg=dlstNeg;
   theCopy->CdlstNeg=CdlstNeg; 
   theCopy->flstNeg=flstNeg;
   theCopy->CflstNeg=CflstNeg; 
   theCopy->ekexcurs=ekexcurs;
   theCopy->Cekexcurs=Cekexcurs; 
   theCopy->RSE=RSE;
   theCopy->CRSE=CRSE; 
   theCopy->fPeakPos=fPeakPos;
   theCopy->CfPeakPos=CfPeakPos; 
   theCopy->fPeakNeg=fPeakNeg;
   theCopy->CfPeakNeg=CfPeakNeg; 
   theCopy->dCap1Pos=dCap1Pos;
   theCopy->CdCap1Pos=CdCap1Pos; 
   theCopy->dCap2Pos=dCap2Pos;
   theCopy->CdCap2Pos=CdCap2Pos; 
   theCopy->dCap1Neg=dCap1Neg;
   theCopy->CdCap1Neg=CdCap1Neg; 
   theCopy->dCap2Neg=dCap2Neg;
   theCopy->CdCap2Neg=CdCap2Neg; 
   theCopy->alphaNeg=alphaNeg;
   theCopy->CalphaNeg=CalphaNeg; 
   theCopy->alphaPos=alphaPos;
   theCopy->CalphaPos=CalphaPos; 
   theCopy->ekhardNeg=ekhardNeg;
   theCopy->CekhardNeg=CekhardNeg; 
   theCopy->ekhardPos=ekhardPos;
   theCopy->CekhardPos=CekhardPos; 
   theCopy->fCapRefPos=fCapRefPos;
   theCopy->CfCapRefPos=CfCapRefPos; 
   theCopy->fCapRefNeg=fCapRefNeg;
   theCopy->CfCapRefNeg=CfCapRefNeg; 
   theCopy->Enrgts=Enrgts;
   theCopy->CEnrgts=CEnrgts; 
   theCopy->Enrgtk=Enrgtk;
   theCopy->CEnrgtk=CEnrgtk; 
   theCopy->Enrgtd=Enrgtd;
   theCopy->CEnrgtd=CEnrgtd; 
   theCopy->dyPos=dyPos;
   theCopy->CdyPos=CdyPos; 
   theCopy->dyNeg=dyNeg;
   theCopy->CdyNeg=CdyNeg; 
   theCopy->dyieldPos=dyieldPos;
   theCopy->CdyieldPos=CdyieldPos; 
   theCopy->dyieldNeg=dyieldNeg;
   theCopy->CdyieldNeg=CdyieldNeg; 
   theCopy->resSnHor=resSnHor;
   theCopy->CresSnHor=CresSnHor; 
   theCopy->fmax=fmax;
   theCopy->Cfmax=Cfmax; 
   theCopy->fmin=fmin;
   theCopy->Cfmin=Cfmin; 
   theCopy->resSp=resSp;
   theCopy->CresSp=CresSp; 
   theCopy->resSn=resSn;
   theCopy->CresSn=CresSn; 
   theCopy->fCapPos=fCapPos;
   theCopy->CfCapPos=CfCapPos; 
   theCopy->fCapNeg=fCapNeg;
   theCopy->CfCapNeg=CfCapNeg; 
   theCopy->snHor=snHor;
   theCopy->CsnHor=CsnHor; 
   theCopy->spHor=spHor;
   theCopy->CspHor=CspHor; 
   theCopy->resSpHor=resSpHor;
   theCopy->CresSpHor=CresSpHor;
   theCopy->snEnv=snEnv;
   theCopy->CsnEnv=CsnEnv;
   theCopy->resSnEnv=resSnEnv;
   theCopy->CresSnEnv=CresSnEnv; 
   theCopy->spEnv=spEnv;
   theCopy->CspEnv=CspEnv; 
   theCopy->resSpEnv=resSpEnv;
   theCopy->CresSpEnv=CresSpEnv; 
   theCopy->Resfac=Resfac;
   theCopy->CResfac=CResfac; 
   theCopy->capSlopeOrig=capSlopeOrig;
   theCopy->CcapSlopeOrig=CcapSlopeOrig;  
   theCopy->fracDispPos=fracDispPos;
   theCopy->CfracDispPos=CfracDispPos; 
   theCopy->fracDispNeg=fracDispNeg;
   theCopy->CfracDispNeg=CfracDispNeg; 
   theCopy->DPlus=DPlus;
   theCopy->CDPlus=CDPlus; 
   theCopy->DNeg=DNeg;
   theCopy->CDNeg=CDNeg; 
   theCopy->alphaN=alphaN;
   theCopy->CalphaN=CalphaN;
   theCopy->ecapa=ecapa; 
   theCopy->Cecapa=Cecapa;
   theCopy->ca=ca;
   theCopy->Cca=Cca;
   theCopy->ResfacNeg=ResfacNeg;
   theCopy->CResfacNeg=CResfacNeg;
   theCopy->myFcapping=myFcapping;
   theCopy->CmyFcapping=CmyFcapping;
   theCopy->myFcappingNeg=myFcappingNeg;
   theCopy->CmyFcappingNeg=CmyFcappingNeg;
   theCopy->capSlopeNeg=capSlopeNeg;
   theCopy->CcapSlopeNeg=CcapSlopeNeg;
   theCopy->flagControlResponse=flagControlResponse;
   theCopy->CflagControlResponse=CflagControlResponse;  

   theCopy->capSlopeOrigNeg=capSlopeOrigNeg;
   theCopy->CcapSlopeOrigNeg=CcapSlopeOrigNeg;
   theCopy->CUprev=CUprev;
   theCopy->Uprev=Uprev;

    return theCopy;
}

int 
Bilin::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(246); ///diorthose to megethos
   data(0) = this->getTag();
   data(1)=Ke;
   data(2)=As;
   data(3)=AsNeg;
   data(4)=My_pos;
   data(5)=My_neg;
   data(6)=LamdaS;
   data(7)=LamdaK;
   data(8)=LamdaA;
   data(9)=LamdaD;
   data(10)=Cs;
   data(11)=Ck;
   data(12)=Ca;
   data(13)=Cd;
   data(14)=Thetap_pos;
   data(15)=Thetap_neg;
   data(16)=Thetapc_pos;
   data(17)=Thetapc_neg;
   data(18)=K;
   data(19)=KNeg;
   data(20)=Thetau_pos;
   data(21)=Thetau_neg;
   data(22)=PDPlus;
   data(23)=PDNeg;
   data(24)=U;
   data(25)=CU; 
   data(26)=Force;
   data(27)=CForce; 
   data(28)=Tangent;
   data(29)=CTangent; 
   data(30)=dNewLoadPos;
   data(31)=CdNewLoadPos;
   data(32)=dNewLoadNeg;
   data(33)=CdNewLoadNeg;
   data(34)=flgstop;
   data(35)=Cflgstop; 
   data(36)=flagdeg;
   data(37)=Cflagdeg; 
   data(38)=flagstopdeg;
   data(39)=Cflagstopdeg; 
   data(40)=ekt;
   data(41)=Cekt;
   data(42)=interup;
   data(43)=Cinterup;  
   data(44)=kcode;
   data(45)=Ckcode; 
   data(46)=kon;
   data(47)=Ckon;
   data(48)=iCapNeg;
   data(49)=CiCapNeg; 
   data(50)=iNoFneg;
   data(51)=CiNoFneg; 
   data(52)=iNoFpos;
   data(53)=CiNoFpos; 
   data(54)=iCapPos;
   data(55)=CiCapPos; 
   data(56)=iDeg;
   data(57)=CiDeg; 
   data(58)=LP;
   data(59)=CLP;
   data(60)=LN;
   data(61)=CLN; 
   data(62)=capSlope;
   data(63)=CcapSlope; 
   data(64)=capDispPos;
   data(65)=CcapDispPos; 
   data(66)=capDispNeg;
   data(67)=CcapDispNeg; 
   data(68)=elstk;
   data(69)=Celstk; 
   data(70)=fyieldPos;
   data(71)=CfyieldPos; 
   data(72)=fyieldNeg;
   data(73)=CfyieldNeg; 
   data(74)=alpha;
   data(75)=Calpha;  
   data(76)=ecaps;
   data(77)=Cecaps; 
   data(78)=ecapk;
   data(79)=Cecapk; 
   data(80)=ecapd;
   data(81)=Cecapd; 
   data(82)=cs;
   data(83)=Ccs; 
   data(84)=ck;
   data(85)=Cck; 
   data(86)=cd;
   data(87)=Ccd; 
   data(88)=dmax;
   data(89)=Cdmax; 
   data(90)=dmin;
   data(91)=Cdmin; 
   data(92)=Enrgtot;
   data(93)=CEnrgtot; 
   data(94)=Enrgc;
   data(95)=CEnrgc; 
   data(96)=fyPos;
   data(97)=CfyPos; 
   data(98)=fLimNeg;
   data(99)=CfLimNeg; 
   data(100)=fyNeg;
   data(101)=CfyNeg; 
   data(102)=ekP;
   data(103)=CekP; 
   data(104)=ekunload;
   data(105)=Cekunload; 
   data(106)=Csp; 
   data(107)=sn;
   data(108)=Csn; 
   data(109)=dP;
   data(110)=CdP; 
   data(111)=fP;
   data(112)=CfP; 
   data(113)=ek;
   data(114)=Cek; 
   data(115)=stif;
   data(116)=Cstif; 
   data(117)=dLimPos;
   data(118)=CdLimPos; 
   data(119)=dLimNeg;
   data(120)=CdLimNeg;  
   data(121)=vtot;
   data(122)=Cvtot; 
   data(123)=ftot;
   data(124)=Cftot; 
   data(125)=dtot;
   data(126)=Cdtot; 
   data(127)=dn;
   data(128)=Cdn; 
   data(129)=cpPos;
   data(130)=CcpPos; 
   data(131)=cpNeg;
   data(132)=CcpNeg; 
   data(133)=fLimPos;
   data(134)=CfLimPos; 
   data(135)=dlstPos;
   data(136)=CdlstPos; 
   data(137)=flstPos;
   data(138)=CflstPos; 
   data(139)=dlstNeg;
   data(140)=CdlstNeg; 
   data(141)=flstNeg;
   data(142)=CflstNeg; 
   data(143)=ekexcurs;
   data(144)=Cekexcurs; 
   data(145)=RSE;
   data(146)=CRSE; 
   data(147)=fPeakPos;
   data(148)=CfPeakPos; 
   data(149)=fPeakNeg;
   data(150)=CfPeakNeg; 
   data(151)=dCap1Pos;
   data(152)=CdCap1Pos; 
   data(153)=dCap2Pos;
   data(154)=CdCap2Pos; 
   data(155)=dCap1Neg;
   data(156)=CdCap1Neg; 
   data(157)=dCap2Neg;
   data(158)=CdCap2Neg; 
   data(159)=alphaNeg;
   data(160)=CalphaNeg; 
   data(161)=alphaPos;
   data(162)=CalphaPos; 
   data(163)=ekhardNeg;
   data(164)=CekhardNeg; 
   data(165)=ekhardPos;
   data(166)=CekhardPos; 
   data(167)=fCapRefPos;
   data(168)=CfCapRefPos; 
   data(169)=fCapRefNeg;
   data(170)=CfCapRefNeg; 
   data(171)=Enrgts;
   data(172)=CEnrgts; 
   data(173)=Enrgtk;
   data(174)=CEnrgtk; 
   data(175)=Enrgtd;
   data(176)=CEnrgtd; 
   data(177)=dyPos;
   data(178)=CdyPos; 
   data(179)=dyNeg;
   data(180)=CdyNeg; 
   data(181)=dyieldPos;
   data(182)=CdyieldPos; 
   data(183)=dyieldNeg;
   data(184)=CdyieldNeg; 
   data(185)=resSnHor;
   data(186)=CresSnHor; 
   data(187)=fmax;
   data(188)=Cfmax; 
   data(189)=fmin;
   data(190)=Cfmin; 
   data(191)=resSp;
   data(192)=CresSp; 
   data(193)=resSn;
   data(194)=CresSn; 
   data(195)=fCapPos;
   data(196)=CfCapPos; 
   data(197)=fCapNeg;
   data(198)=CfCapNeg; 
   data(199)=snHor;
   data(200)=CsnHor; 
   data(201)=spHor;
   data(202)=CspHor; 
   data(203)=resSpHor;
   data(204)=CresSpHor;
   data(205)=snEnv;
   data(206)=CsnEnv;
   data(207)=resSnEnv;
   data(208)=CresSnEnv; 
   data(209)=spEnv;
   data(210)=CspEnv; 
   data(211)=resSpEnv;
   data(212)=CresSpEnv; 
   data(213)=Resfac;
   data(214)=CResfac; 
   data(215)=capSlopeOrig;
   data(216)=CcapSlopeOrig;  
   data(217)=fracDispPos;
   data(218)=CfracDispPos; 
   data(219)=fracDispNeg;
   data(220)=CfracDispNeg; 
   data(221)=DPlus;
   data(222)=CDPlus; 
   data(223)=DNeg;
   data(224)=CDNeg; 
   data(225)=alphaN;
   data(226)=CalphaN;
   data(227)=ecapa; 
   data(228)=Cecapa;
   data(229)=ca;
   data(230)=Cca;
   data(231)=ResfacNeg;
   data(232)=CResfacNeg;
   data(233)=myFcapping;
   data(234)=CmyFcapping;
   data(235)=myFcappingNeg;
   data(236)=CmyFcappingNeg;
   data(237)=capSlopeNeg;
   data(238)=CcapSlopeNeg;
   data(239)=flagControlResponse;
   data(240)=CflagControlResponse; 
   data(241)=sp;
   data(242)=CcapSlopeOrigNeg;
   data(243)=capSlopeOrigNeg;
   data(244)=CUprev;
   data(245)=Uprev;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Bilin::sendSelf() - failed to send data\n";

  return res;
}

int 
Bilin::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(246);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "Bilin::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
   Ke=data(1);
   As=data(2);
   AsNeg=data(3);
   My_pos=data(4);
   My_neg=data(5);
   LamdaS=data(6);
   LamdaK=data(7);
   LamdaA=data(8);
   LamdaD=data(9);
   Cs=data(10);
   Ck=data(11);
   Ca=data(12);
   Cd=data(13);
   Thetap_pos=data(14);
   Thetap_neg=data(15);
   Thetapc_pos=data(16);
   Thetapc_neg=data(17);
   K=data(18);
   KNeg=data(19);
   Thetau_pos=data(20);
   Thetau_neg=data(21);
   PDPlus=data(22);
   PDNeg=data(23);
   U=data(24);
   CU=data(25); 
   Force=data(26);
   CForce=data(27); 
   Tangent=data(28);
   CTangent=data(29); 
   dNewLoadPos=data(30);
   CdNewLoadPos=data(31);
   dNewLoadNeg=data(32);
   CdNewLoadNeg=data(33);
   flgstop=data(34);
   Cflgstop=data(35); 
   flagdeg=data(36);
   Cflagdeg=data(37); 
   flagstopdeg=data(38);
   Cflagstopdeg=data(39); 
   ekt=data(40);
   Cekt=data(41);
   interup=data(42);
   Cinterup=data(43);  
   kcode=data(44);
   Ckcode=data(45); 
   kon=data(46);
   Ckon=data(47);
   iCapNeg=data(48);
   CiCapNeg=data(49); 
   iNoFneg=data(50);
   CiNoFneg=data(51); 
   iNoFpos=data(52);
   CiNoFpos=data(53); 
   iCapPos=data(54);
   CiCapPos=data(55); 
   iDeg=data(56);
   CiDeg=data(57); 
   LP=data(58);
   CLP=data(59);
   LN=data(60);
   CLN=data(61); 
   capSlope=data(62);
   CcapSlope=data(63); 
   capDispPos=data(64);
   CcapDispPos=data(65); 
   capDispNeg=data(66);
   CcapDispNeg=data(67); 
   elstk=data(68);
   Celstk=data(69); 
   fyieldPos=data(70);
   CfyieldPos=data(71); 
   fyieldNeg=data(72);
   CfyieldNeg=data(73); 
   alpha=data(74);
   Calpha=data(75);  
   ecaps=data(76);
   Cecaps=data(77); 
   ecapk=data(78);
   Cecapk=data(79); 
   ecapd=data(80);
   Cecapd=data(81); 
   cs=data(82);
   Ccs=data(83); 
   ck=data(84);
   Cck=data(85); 
   cd=data(86);
   Ccd=data(87); 
   dmax=data(88);
   Cdmax=data(89); 
   dmin=data(90);
   Cdmin=data(91); 
   Enrgtot=data(92);
   CEnrgtot=data(93); 
   Enrgc=data(94);
   CEnrgc=data(95); 
   fyPos=data(96);
   CfyPos=data(97); 
   fLimNeg=data(98);
   CfLimNeg=data(99); 
   fyNeg=data(100);
   CfyNeg=data(101); 
   ekP=data(102);
   CekP=data(103); 
   ekunload=data(104);
   Cekunload=data(105); 
   Csp=data(106); 
   sn=data(107);
   Csn=data(108); 
   dP=data(109);
   CdP=data(110); 
   fP=data(111);
   CfP=data(112); 
   ek=data(113);
   Cek=data(114); 
   stif=data(115);
   Cstif=data(116); 
   dLimPos=data(117);
   CdLimPos=data(118); 
   dLimNeg=data(119);
   CdLimNeg=data(120);  
   vtot=data(121);
   Cvtot=data(122); 
   ftot=data(123);
   Cftot=data(124); 
   dtot=data(125);
   Cdtot=data(126); 
   dn=data(127);
   Cdn=data(128); 
   cpPos=data(129);
   CcpPos=data(130); 
   cpNeg=data(131);
   CcpNeg=data(132); 
   fLimPos=data(133);
   CfLimPos=data(134); 
   dlstPos=data(135);
   CdlstPos=data(136); 
   flstPos=data(137);
   CflstPos=data(138); 
   dlstNeg=data(139);
   CdlstNeg=data(140); 
   flstNeg=data(141);
   CflstNeg=data(142); 
   ekexcurs=data(143);
   Cekexcurs=data(144); 
   RSE=data(145);
   CRSE=data(146); 
   fPeakPos=data(147);
   CfPeakPos=data(148); 
   fPeakNeg=data(149);
   CfPeakNeg=data(150); 
   dCap1Pos=data(151);
   CdCap1Pos=data(152); 
   dCap2Pos=data(153);
   CdCap2Pos=data(154); 
   dCap1Neg=data(155);
   CdCap1Neg=data(156); 
   dCap2Neg=data(157);
   CdCap2Neg=data(158); 
   alphaNeg=data(159);
   CalphaNeg=data(160); 
   alphaPos=data(161);
   CalphaPos=data(162); 
   ekhardNeg=data(163);
   CekhardNeg=data(164); 
   ekhardPos=data(165);
   CekhardPos=data(166); 
   fCapRefPos=data(167);
   CfCapRefPos=data(168); 
   fCapRefNeg=data(169);
   CfCapRefNeg=data(170); 
   Enrgts=data(171);
   CEnrgts=data(172); 
   Enrgtk=data(173);
   CEnrgtk=data(174); 
   Enrgtd=data(175);
   CEnrgtd=data(176); 
   dyPos=data(177);
   CdyPos=data(178); 
   dyNeg=data(179);
   CdyNeg=data(180); 
   dyieldPos=data(181);
   CdyieldPos=data(182); 
   dyieldNeg=data(183);
   CdyieldNeg=data(184); 
   resSnHor=data(185);
   CresSnHor=data(186); 
   fmax=data(187);
   Cfmax=data(188); 
   fmin=data(189);
   Cfmin=data(190); 
   resSp=data(191);
   CresSp=data(192); 
   resSn=data(193);
   CresSn=data(194); 
   fCapPos=data(195);
   CfCapPos=data(196); 
   fCapNeg=data(197);
   CfCapNeg=data(198); 
   snHor=data(199);
   CsnHor=data(200); 
   spHor=data(201);
   CspHor=data(202); 
   resSpHor=data(203);
   CresSpHor=data(204);
   snEnv=data(205);
   CsnEnv=data(206);
   resSnEnv=data(207);
   CresSnEnv=data(208); 
   spEnv=data(209);
   CspEnv=data(210); 
   resSpEnv=data(211);
   CresSpEnv=data(212); 
   Resfac=data(213);
   CResfac=data(214); 
   capSlopeOrig=data(215);
   CcapSlopeOrig=data(216);  
   fracDispPos=data(217);
   CfracDispPos=data(218); 
   fracDispNeg=data(219);
   CfracDispNeg=data(220); 
   DPlus=data(221);
   CDPlus=data(222); 
   DNeg=data(223);
   CDNeg=data(224); 
   alphaN=data(225);
   CalphaN=data(226);
   ecapa=data(227); 
   Cecapa=data(228);
   ca=data(229);
   Cca=data(230);
   ResfacNeg=data(231);
   CResfacNeg=data(232);
   myFcapping=data(233);
   CmyFcapping=data(234);
   myFcappingNeg=data(235);
   CmyFcappingNeg=data(236);
   capSlopeNeg=data(237);
   CcapSlopeNeg=data(238);
   flagControlResponse=data(239);
   CflagControlResponse=data(240); 
   sp=data(241);
   CcapSlopeOrigNeg=data(242);
   capSlopeOrigNeg=data(243);
   CUprev=data(244);
   Uprev=data(245);
  }
    
  return res;
}

void 
Bilin::Print(OPS_Stream &s, int flag)
{
    //s << "Bilin tag: " << this->getTag() << endln;
    //s << "  G: " << G << endln;
    //s << "  t: " << t << endln;
  //s << "  A: " << A << endln;
  s << "sp: " << sp << endln;
  s << "cpNeg: " << cpNeg << endln;
  s << "dyNeg: " << dyNeg << endln;
  
  
}

//my functions
void
Bilin::interPoint(double &xInt, double &yInt, double x1, double y1, double m1, double x2, double y2, double m2)
{
xInt = (-m2*x2+y2+m1*x1-y1) / (m1-m2);
yInt = m1*xInt-m1*x1+y1;
}
void Bilin::snCalc(void)
{
// c	Each time that the hysteretic loop changes from negative to positive
// c	delta displacement, this subroutine calculates the point where ends
// c	ekunload and starts	whether the positive hardening curve or the 
// c	positive cap curve. In cases where "cpPos<fyPos" this point could
// c	not be reached, because	the unloading stiffness could intersect the
// c	positive cap curve. However, this change is reflected in the main 
// c	program.
// c	
// c	Output Variables: sn,resSn,snHard,resSnHard
// c	Input Variables:  dP,fP,ekunload,alphaPos,dyPos,fyPos,cpPos,fCapPos
// c	     			  capStiff,fCapRefPos

	double Resid = Resfac*fyPos;
	double dresid = cpPos+(Resid-fCapPos)/(capSlope*elstk);
	double ekresid = 1.0e-10;
	dyPos = fyPos/elstk;
    double snHard,resSnHard,snLim,resSnLim,snResid,resSnResid;
	if (dyPos<cpPos) {
// 	call interPoint(snHard,resSnHard,dyPos,fyPos,elstk*alphaPos,dP,fP,ekunload)
		
    interPoint(snHard,resSnHard,dyPos,fyPos,elstk*alphaPos,dP,fP,ekunload);
	} else {
// 	call interPoint(snHard,resSnHard,cpPos,fCapPos,elstk*alphaPos,dP,fP,ekunload)
    
    interPoint(snHard,resSnHard,cpPos,fCapPos,elstk*alphaPos,dP,fP,ekunload);
	}

// 	call interPoint(snCap,resSnCap,0.d0,fCapRefPos,capSlope*elstk,dP,fP,ekunload)
	double snCap,resSnCap;
    interPoint(snCap,resSnCap,0.0,fCapRefPos,capSlope*elstk,dP,fP,ekunload);

	//sn = min(snHard,snCap);
   if (snHard<snCap)
   {
	   sn=snHard;
   }
   else
   {
	   sn=snCap;
   }
	//resSn = min(resSnHard,resSnCap);
   if (resSnHard<resSnCap)
   {
	   resSn=resSnHard;
   }
   else
   {
	   resSn=resSnCap;
   }
	snEnv = sn;
	resSnEnv = resSn;
// 	call interPoint(temp_1_farzin,temp,0.d0,fCapRefPos,capSlope*elstk,0.d0,Resid,0.d0)
    //don't need to call that!!!!!!!!!!!!!!!
	if((LP==1)&&(fLimPos==0.0)) { 
		//call interPoint(snLim,resSnLim,dLimPos,fLimPos,0.d0,dP,fP,ekunload)
        interPoint(snLim,resSnLim,dLimPos,fLimPos,0.0,dP,fP,ekunload);
		if (snLim<sn){  
			sn=snLim;
			resSn=resSnLim;
		}

//call interPoint(snHor,resSnHor,dLimPos,fLimPos,0.d0,dyPos,fyPos,elstk*alphaPos)
     interPoint(snHor,resSnHor,dLimPos,fLimPos,0.0,dyPos,fyPos,elstk*alphaPos);
	}

	if (sn>dresid) { 
// 		call interPoint(snResid,resSnResid,dresid,Resid,ekresid,dP,fP,ekunload)
        interPoint(snResid,resSnResid,dresid,Resid,ekresid,dP,fP,ekunload);
		sn = snResid;
		resSn = resSnResid;
	}
}


void 
Bilin::envelPosCap2(double fy,double alphaPos,double alphaCap,double cpDsp,double& d,
				  double& f,double& ek,double elstk,double fyieldPos,double Resfac)
{


    //double fracDispPosV=0.20;
	double dy = fy/elstk;

     double Res,rcap,dres;
	if(dy<=cpDsp) {
		Res = Resfac*fyieldPos;
		rcap = fy+alphaPos*elstk*(cpDsp-dy);
		dres = cpDsp+(Res-rcap)/(alphaCap*elstk);

		if (d<0.0){
			f = 0.0;
			ek = 1.0e-7;
		}else if (d<=dy) { 
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
// c added by Dimitrios to account for fracture	
		if(d>=fracDispPos) {
	        ek = 1.0e-7;
			f = 1.0e-10;
			d=fracDispPos;
            flgstop=1;
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
// c added by Dimitrios to account for fracture	
		if(d>=fracDispPos) {
	        ek = 1.0e-7;
			f = 1.0e-10;
			d=fracDispPos;
            flgstop=1;  
		}
		
	} else {
	}
}
double
Bilin::boundPos(void)
{

    double Resid=0.0;    
     double dBoundPos;
     dyNeg = fyNeg/elstk;
	double dresid = cpPos+(Resid-fCapPos)/(capSlope*elstk);
	double ekresid = 1.0e-10;

// 	call interPoint(d1,f1,E
	double d1,f1;
    interPoint(d1,f1,dyNeg,fyNeg,elstk*alphaNeg,0.0,fCapRefPos,capSlope*elstk);

//	call interPoint(d2,f2,dyNeg,fyNeg,elstk*alphaNeg,dresid,resid,ekresid)
	double d2,f2;
    interPoint(d2,f2,dyNeg,fyNeg,elstk*alphaNeg,dresid,Resid,ekresid);
	//dBoundPos = max(d1,d2);
	if (d1>d2)
	{
		dBoundPos=d1;
	}
	else
	{
		dBoundPos=d2;
	}
	return(dBoundPos);
}
void 
Bilin::envHitsZero(double& f)
{
	if (fP>0) {
		if ((f*fP)<0) {
			f = 0;
			ek = 1.0e-7;
			iNoFpos = 1;
            flagControlResponse = 1;
		}
	}else if (fP<0) {
		if ((f*fP)<0) {
			f = 0;
			ek = 1.0e-7;
			iNoFneg = 1;
            flagControlResponse = 1;
		}
	} else {
	}
}
void
Bilin::envelNegCap2(double fy,double alphaNeg,double alphaCap,double cpDsp,double& d,double& f,double& ek,
							   double elstk,double fyieldNeg,double Resfac)
{

	//double fracDispNegV = -0.20;
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
//% c added by Dimitrios to account for fracture	
		if(d<=fracDispNeg) {
	        ek = 1.0e-7;
			f = 1.0e-10;
			d=fracDispNeg;
            flgstop=1;
		}
	
	} else if(dy<cpDsp) { 
        
		rcap = elstk*cpDsp;
		Res = Resfac*rcap;
		dres = cpDsp+(Res-rcap)/(alphaCap*elstk);

		if (d>0.0) {
			f = 0.0;
			ek = 1.0e-7;
		}	else if (d>=cpDsp) {
			ek = elstk;
			f = ek*d;
		} else if (d>=dres) {
			ek = elstk*alphaCap;
			f = rcap+ek*(d-cpDsp);
		} else {
			ek = 1.0e-7;
			f = Res+ek*d;
		}
// c added by Dimitrios to account for fracture	
       if(d<=fracDispNeg) {
	        ek = 1.0e-7;
			f = 1.0e-10;
			d=fracDispNeg;
            flgstop=1;
	   }

	} else {
	}
}
void
Bilin::spCalc(void)
{
	double Resid = ResfacNeg*fyNeg;
	dyNeg = fyNeg/elstk;
	double dresid = cpNeg+(Resid-fCapNeg)/(capSlopeNeg*elstk);
	double ekresid = 1.0e-10;
     double spHard,resSpHard,spCap,resSpCap,spLim,resSpLim,spResid,resSpResid;
	if (dyNeg>cpNeg) {
// 	call interPoint(spHard,resSpHard,dyNeg,fyNeg,elstk*alphaNeg,dP,fP,ekunload)
    interPoint(spHard,resSpHard,dyNeg,fyNeg,elstk*alphaNeg,dP,fP,ekunload);
	} else {
// 	call interPoint(spHard,resSpHard,cpNeg,fCapNeg,elstk*alphaNeg,dP,fP,ekunload)
    interPoint(spHard,resSpHard,cpNeg,fCapNeg,elstk*alphaNeg,dP,fP,ekunload);
	}

// 	call interpoint(spCap,resSpCap,0.d0,fCapRefNeg,capSlope*elstk,dP,fP,ekunload)
    interPoint(spCap,resSpCap,0.0,fCapRefNeg,capSlopeNeg*elstk,dP,fP,ekunload);
	//sp = max(spHard,spCap);
	if(spHard>spCap)
	{
		sp=spHard;
	}
	else
	{
		sp=spCap;
	}
	//resSp = max(resSpHard,resSpCap);
	if(resSpHard>resSpCap)
	{
		resSp=resSpHard;
	}
	else
	{
		resSp=resSpCap;
	}
	spEnv = sp;
	resSpEnv = resSp;

	if((LN==1)&&(fLimNeg==0.0)) {
// 		call interPoint(spLim,resSpLim,dLimNeg,fLimNeg,0.d0,dP,fP,ekunload)
        interPoint(spLim,resSpLim,dLimNeg,fLimNeg,0.0,dP,fP,ekunload);
		if (spLim>sp) {
			sp=spLim;
			resSp=resSpLim;
		}

// 		call interPoint(spHor,resSpHor,dLimNeg,fLimNeg,0.d0,dyNeg,fyNeg,elstk*alphaNeg)
         interPoint(spHor,resSpHor,dLimNeg,fLimNeg,0.0,dyNeg,fyNeg,elstk*alphaNeg);
	}

	if (sp<dresid) { 
// 		call interPoint(spResid,resSpResid,dresid,Resid,ekresid,dP,fP,ekunload)
        interPoint(spResid,resSpResid,dresid,Resid,ekresid,dP,fP,ekunload);
		sp = spResid;
		resSp = resSpResid;
	}
}
double
Bilin::boundNeg(void)
{

    double Resid = 0;
    double dBoundNeg;
    dyPos = fyPos/elstk;
	double dresid = cpNeg+(Resid-fCapNeg)/(capSlopeNeg*elstk);
	double ekresid = 1.0e-10;
double d1,f1,d2,f2;
// 	call interPoint(d1,f1,dyPos,fyPos,elstk*alphaPos,0.d0,fCapRefNeg,capSlope*elstk)
    interPoint(d1,f1,dyPos,fyPos,elstk*alphaPos,0.0,fCapRefNeg,capSlopeNeg*elstk);
//% 
// 	call interPoint(d2,f2,dyPos,fyPos,elstk*alphaPos,dresid,resid,ekresid)
    interPoint(d2,f2,dyPos,fyPos,elstk*alphaPos,dresid,Resid,ekresid);
	//dBoundNeg = min(d1,d2);
	 if (d1<d2)
	 {
		 dBoundNeg=d1;
	 }
	 else
	 {
		 dBoundNeg=d2;
	 }
	return (dBoundNeg);
}

int
Bilin::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"Ke") == 0)
    return param.addObject(1, this);
  else if (strcmp(argv[0],"As") == 0)
    return param.addObject(2, this);
  else if (strcmp(argv[0],"AsNeg") == 0)
    return param.addObject(3, this);
  else if (strcmp(argv[0],"My_pos") == 0)
    return param.addObject(4, this);
  else if (strcmp(argv[0],"My_neg") == 0)
    return param.addObject(5, this);
  else if (strcmp(argv[0],"LambaS") == 0)
    return param.addObject(6, this);
  else if (strcmp(argv[0],"LamdaK") == 0)
    return param.addObject(7, this);
  else if (strcmp(argv[0],"LamdaA") == 0)
    return param.addObject(8, this);
  else if (strcmp(argv[0],"LambdaD") == 0)
    return param.addObject(9, this);
  else if (strcmp(argv[0],"Cs") == 0)
    return param.addObject(10, this);
  else if (strcmp(argv[0],"Ck") == 0)
    return param.addObject(11, this);
  else if (strcmp(argv[0],"Ca") == 0)
    return param.addObject(12, this);
  else if (strcmp(argv[0],"Cd") == 0)
    return param.addObject(13, this);
  else if (strcmp(argv[0],"Thetap_pos") == 0)
    return param.addObject(14, this);
  else if (strcmp(argv[0],"Thetap_neg") == 0)
    return param.addObject(15, this);
  else if (strcmp(argv[0],"Thetapc_pos") == 0)
    return param.addObject(16, this);
  else if (strcmp(argv[0],"Thetapc_neg") == 0)
    return param.addObject(17, this);
  else if (strcmp(argv[0],"K") == 0)
    return param.addObject(18, this);
  else if (strcmp(argv[0],"KNeg") == 0)
    return param.addObject(19, this);
  else if (strcmp(argv[0],"Thetau_pos") == 0)
    return param.addObject(20, this);
  else if (strcmp(argv[0],"Thetau_neg") == 0)
    return param.addObject(21, this);
  else if (strcmp(argv[0],"PDPlus") == 0)
    return param.addObject(22, this);
  else if (strcmp(argv[0],"PDNeg") == 0)
    return param.addObject(23, this);

  return -1;
}

int 
Bilin::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    Ke = info.theDouble;
    return 0;
  case 2:
    As = info.theDouble;
    return 0;
  case 3:
    AsNeg = info.theDouble;
    return 0;
  case 4:
    My_pos = info.theDouble;
    return 0;
  case 5:
    My_neg = info.theDouble;
    return 0;
  case 6:
    LamdaS = info.theDouble;
    return 0;
  case 7:
    //    LambdaK = info.theDouble;
    return 0;
  case 8:
    // LambdaA = info.theDouble;
    return 0;
  case 9:
    // LambdaD = info.theDouble;
    return 0;
  case 10:
    Cs = info.theDouble;
    return 0;
  case 11:
    Ck = info.theDouble;
    return 0;
  case 12:
    Ca = info.theDouble;
    return 0;
  case 13:
    Cd = info.theDouble;
    return 0;
  case 14:
    Thetap_pos = info.theDouble;
    return 0;
  case 15:
    Thetap_neg = info.theDouble;
    return 0;
  case 16:
    Thetapc_pos = info.theDouble;
    return 0;
  case 17:
    Thetapc_neg = info.theDouble;
    return 0;
  case 18:
    K = info.theDouble;
    return 0;
  case 19:
    KNeg = info.theDouble;
    return 0;
  case 20:
    Thetau_pos = info.theDouble;
    return 0;
  case 21:
    Thetau_neg = info.theDouble;
    return 0;
  case 22:
    PDPlus = info.theDouble;
    return 0;
  case 23:
    PDNeg = info.theDouble;
    return 0;

  default:
    return -1;
  }
}


