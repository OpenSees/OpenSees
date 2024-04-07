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

//Modified Ibarra-Medina-Krawinkler with Bilinear hysteretic response

//**********************************************************************                                                                    
// Code Developed by: Dimitrios G. Lignos
// Assistant Professor, McGill University, Montreal Canada
// Originally Written by: Theodore Karavasilis
// Lecturer
// Department of Engineering Science, University of Oxford, Oxford, U.K.
// Re-Written by: D. G. Lignos, August 29th 2012

//**********************************************************************
// Adapted by: Filipe Ribeiro and Andre Barbosa, Sep 13th 2013
// Oregon State University, OR, USA
//**********************************************************************

#include <math.h>

#include <Bilin.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <cfloat>
#include <MaterialResponse.h>

#include <OPS_Globals.h>


static int numBilinMaterials = 0;

void *
OPS_Bilin()
{
  if (numBilinMaterials == 0) {
    numBilinMaterials++;
    opserr << "WARNING: DO NOT USE THE \"Bilin\" MATERIAL, IT HAS BEEN REPLACED. Use \"IMKBilin\" or \"HystereticSM\" INSTEAD.\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[24];                             // Updated: Filipe Ribeiro and Andre Barbosa
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  Bilin tag" << endln;
    return 0;
  }
 
  //Changed in order to account for the nFactor as an optional input    // Updated: Filipe Ribeiro and Andre Barbosa
  //numData = 24;               // Updated: Filipe Ribeiro and Andre Barbosa
  numData = OPS_GetNumRemainingInputArgs();     // Updated: Filipe Ribeiro and Andre Barbosa 
 
  if (numData != 23 && numData != 24 ) {        // Updated: Filipe Ribeiro and Andre Barbosa
    opserr << "Invalid Args want: uniaxialMaterial Bilin tag? Ke? AsPos? AsNeg? My_pos? My_neg? LamdaS? ";             
    opserr << "LamdaD?  LamdaA? LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? Thetap_neg? Thetapc_pos? Thetapc_neg?KPos? ";              
    opserr << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg?  <nFactor?> \n";           // Updated: Filipe Ribeiro and Andre Barbosa
    return 0;           // Updated: Filipe Ribeiro and Andre Barbosa
  }
 
  //  if (OPS_GetDoubleInput(&numData, dData) != 0) {                   // Updated: Filipe Ribeiro and Andre Barbosa 
  //    opserr << "Invalid Args want: uniaxialMaterial Bilin tag? Ke? nFactor? AsPos? AsNeg? My_pos? My_neg? LamdaS? ";         
  //    opserr << "LamdaD?  LamdaA? LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? Thetap_neg? Thetapc_pos? Thetapc_neg?KPos? "; 
  //    opserr << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg\n";     
  //    return 0;  
  //  }																	// Updated: Filipe Ribeiro and Andre Barbosa 
 
  if (numData == 23) {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {             // Updated: Filipe Ribeiro and Andre Barbosa 
      opserr << "Invalid Args want: uniaxialMaterial Bilin tag? Ke? AsPos? AsNeg? My_pos? My_neg? LamdaS? ";            
      opserr << "LamdaD?  LamdaA? LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? Thetap_neg? Thetapc_pos? Thetapc_neg?KPos? "; 
      opserr << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg? <nFactor?> \n";          // Updated: Filipe Ribeiro and Andre Barbosa
      return 0;  
    }
    // Parsing was successful, allocate the material
    theMaterial = new Bilin(iData[0],
                            dData[0], dData[1], dData[2], dData[3], dData[4],
                            dData[5], dData[6], dData[7], dData[8], dData[9],
                            dData[10], dData[11], dData[12], dData[13], dData[14],
                            dData[15], dData[16], dData[17], dData[18], dData[19],
                            dData[20], dData[21], dData[22]);           // Updated: Filipe Ribeiro and Andre Barbosa
   
  } else if (numData == 24) { // Updated: Filipe Ribeiro and Andre Barbosa 
    if (OPS_GetDoubleInput(&numData, dData) != 0) {             
      opserr << "Invalid Args want: uniaxialMaterial Bilin tag? Ke? AsPos? AsNeg? My_pos? My_neg? LamdaS? ";            
      opserr << "LamdaD?  LamdaA? LamdaK? Cs? Cd? Ca? Ck? Thetap_pos? Thetap_neg? Thetapc_pos? Thetapc_neg?KPos? "; 
      opserr << "KNeg? Thetau_pos? Thetau_neg? PDPlus?  PDNeg? <nFactor?>\n";           // Updated: Filipe Ribeiro and Andre Barbosa 
      return 0;  
    }
    // Parsing was successful, allocate the material
    theMaterial = new Bilin(iData[0],
                            dData[0], dData[1], dData[2], dData[3], dData[4],
                            dData[5], dData[6], dData[7], dData[8], dData[9],
                            dData[10], dData[11], dData[12], dData[13], dData[14],
                            dData[15], dData[16], dData[17], dData[18], dData[19],
                            dData[20], dData[21], dData[22], dData[23]); // Updated: Filipe Ribeiro and Andre Barbosa
   
  }
 
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Bilin Material\n";
    return 0;
  }
 
  return theMaterial;
}

// if nFactor is assigned                                                       // Updated: Filipe Ribeiro and Andre Barbosa 
Bilin::Bilin(int tag, double p_Ke0,double p_AsPos,double p_AsNeg,double p_My_pos,double p_My_neg,double p_LamdaS,              
             double p_LamdaD,double p_LamdaA,double p_LamdaK,double p_Cs,double p_Cd,double p_Ca,double p_Ck,
             double p_Thetap_pos,double p_Thetap_neg,double p_Thetapc_pos,double p_Thetapc_neg,double p_KPos,double p_KNeg,
             double p_Thetau_pos,double p_Thetau_neg,double p_PDPlus,double p_PDNeg,double p_nFactor)  
:UniaxialMaterial(tag, MAT_TAG_Bilin),Ke0(p_Ke0), AsPos(p_AsPos), AsNeg(p_AsNeg), My_pos(p_My_pos), My_neg(p_My_neg),          
 LamdaS(p_LamdaS), LamdaK(p_LamdaK),LamdaA(p_LamdaA),LamdaD(p_LamdaD), Cs(p_Cs), Ck(p_Ck), Ca(p_Ca),Cd(p_Cd),
 Thetap_pos(p_Thetap_pos), Thetap_neg(p_Thetap_neg), Thetapc_pos(p_Thetapc_pos),Thetapc_neg(p_Thetapc_neg),
 KPos(p_KPos), KNeg(p_KNeg),Thetau_pos(p_Thetau_pos), Thetau_neg(p_Thetau_neg), PDPlus(p_PDPlus), PDNeg(p_PDNeg), nFactor(p_nFactor)    // Updated: Filipe Ribeiro and Andre Barbosa 
{
  //initialize variables
  this->revertToStart();
  //this->revertToLastCommit();
}

// if nFactor is NOT assigned																							// Updated: Filipe Ribeiro and Andre Barbosa         
Bilin::Bilin(int tag, double p_Ke0,double p_AsPos,double p_AsNeg,double p_My_pos,double p_My_neg,double p_LamdaS,              
             double p_LamdaD,double p_LamdaA,double p_LamdaK,double p_Cs,double p_Cd,double p_Ca,double p_Ck,
             double p_Thetap_pos,double p_Thetap_neg,double p_Thetapc_pos,double p_Thetapc_neg,double p_KPos,double p_KNeg,
             double p_Thetau_pos,double p_Thetau_neg,double p_PDPlus,double p_PDNeg)   
:UniaxialMaterial(tag, MAT_TAG_Bilin),Ke0(p_Ke0), AsPos(p_AsPos), AsNeg(p_AsNeg), My_pos(p_My_pos), My_neg(p_My_neg),          
 LamdaS(p_LamdaS), LamdaK(p_LamdaK),LamdaA(p_LamdaA),LamdaD(p_LamdaD), Cs(p_Cs), Ck(p_Ck), Ca(p_Ca),Cd(p_Cd),
 Thetap_pos(p_Thetap_pos), Thetap_neg(p_Thetap_neg), Thetapc_pos(p_Thetapc_pos),Thetapc_neg(p_Thetapc_neg),
 KPos(p_KPos), KNeg(p_KNeg),Thetau_pos(p_Thetau_pos), Thetau_neg(p_Thetau_neg), PDPlus(p_PDPlus), PDNeg(p_PDNeg)        // Updated: Filipe Ribeiro and Andre Barbosa 
{
  //initialize variables
  this->revertToStart();
  //this->revertToLastCommit();
  // Default value for no nFactor
  nFactor=0.0;												// Updated: Filipe Ribeiro and Andre Barbosa 
}


Bilin::Bilin()
 :UniaxialMaterial(0, MAT_TAG_Bilin),
  Ke0(0), AsPos(0), AsNeg(0), My_pos(0), My_neg(0),  
  LamdaS(0), LamdaD(0),LamdaA(0),LamdaK(0), Cs(0), Cd(0), Ca(0),Ck(0),
  Thetap_pos(0), Thetap_neg(0), Thetapc_pos(0),Thetapc_neg(0),
  KPos(0), KNeg(0),Thetau_pos(0), Thetau_neg(0), PDPlus(0), PDNeg(0), nFactor(0) // Updated: Filipe Ribeiro and Andre Barbosa
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
  double deltaD,d,temp_1,temp,betas,betak,betad,
    dBoundPos,dBoundNeg,f1,f2,fNewLoadPos,fNewLoadNeg,
    f,xDevPos1,yDevPos1,xDevPos2,yDevPos2,xDevPos,yDevPos,xDevNeg1,
    yDevNeg1,xDevNeg2,yDevNeg2,xDevNeg,yDevNeg;        
 
  //Initialize f (needed in case residual region is attained)
  f=0.0;     //Updated: Filipe Ribeiro and Andre Barbosa
 
  //state determination algorithm: defines the current force and tangent stiffness
  U=strain; //set trial displacement
 
  //********************************matlab code*************************************************************
  // Incremental Displacement
  deltaD = strain - CU;
  d=strain;
  dP = CU;          

  if (fabs(deltaD) < DBL_EPSILON && commitCalledOnce == 1) {
    return 0;
  }


  if (d>0.0) {
   
    interPoint(temp_1,temp,0.0,fCapRefPos,capSlope*Ke,0.0,KPos*My_pos,0.0);
    if (d<temp_1){
      iNoFpos = 0;
      LP=0;
    }
   
  } else {
   
    interPoint(temp_1,temp,0.0,fCapRefNeg,capSlopeNeg*Ke,0.0,KNeg*My_neg,0.0);
   
    if (d>temp_1) {
     
      iNoFneg = 0;
      LN=0;
    }
  }
 
  // Other variables
 
  flagdeg = 0;
  betas = 0.0;
  betak = 0.0;  
  betad = 0.0;
 
  // Initialize parameters in the first cycle
 
  if (kon==0) {
   
    //Compute elastic stiffness, strain hardening and post-capping ratios (Ibarra & Krawinkler, 2005)  
   
    // Amplify the elastic stiffness in case n>0
    Ke = Ke0*(1 + nFactor);  // Updated: Filipe Ribeiro and Andre Barbosa 
   
    //Compute strain hardening ratios (Ibarra & Krawinkler, 2005)
    alphaNeg=AsNeg/(1+nFactor*(1-AsNeg)); // updated per inelastic cycle // Updated: Filipe Ribeiro and Andre Barbosa
    alphaPos=AsPos/(1+nFactor*(1-AsPos)); // updated per inelastic cycle // Updated: Filipe Ribeiro and Andre Barbosa
   
    ekhardPos = Ke*AsPos/(1+nFactor*(1-AsPos));         // Updated: Filipe Ribeiro and Andre Barbosa
    ekhardNeg = Ke*AsNeg/(1+nFactor*(1-AsNeg));         // Updated: Filipe Ribeiro and Andre Barbosa
   
    capSlopeMember = -(My_pos + ekhardPos * Thetap_pos)/(Thetapc_pos*Ke0);      // Updated: Filipe Ribeiro and Andre Barbosa
    capSlope = capSlopeMember/(1+nFactor*(1- capSlopeMember));          // Updated: Filipe Ribeiro and Andre Barbosa
    capSlopeNegMember = -(-My_neg + ekhardNeg * Thetap_neg)/(Thetapc_neg*Ke0); // Updated: Filipe Ribeiro and Andre Barbosa
    capSlopeNeg = capSlopeNegMember/(1+nFactor*(1- capSlopeNegMember));// Updated: Filipe Ribeiro and Andre Barbosa
   
	prodBeta=1.0;				// Updated: Filipe Ribeiro and Andre Barbosa 

    //
   
    ekP  = Ke;
   
    flagControlResponse = 0;    
    Tangent=Ke;    
   
    Enrgts = My_pos*LamdaS;
    Enrgtk = 2.0*My_pos*LamdaK;
    Enrgtd = My_pos*LamdaD;
   
    dmax = My_pos/Ke;
    dmin = (My_neg/Ke);
   
    ekP = Ke;
   
    ekunload = Ke;
    ekexcurs = Ke;
   
    Enrgtot = 0.0;
    Enrgc = 0.0;
   
    fyPos = My_pos; // updated per inelastic cycle
    fyNeg = My_neg; // updated per inelastic cycle
   
    dyPos=My_pos/Ke;
    dyNeg=(My_neg/Ke);
   
    resSn = My_pos;
    resSp = My_neg;
   
    cpPos = Thetap_pos+My_pos/Ke;
    cpNeg = (-Thetap_neg+My_neg/Ke);
   
    fPeakPos=My_pos+ekhardPos*((Thetap_pos+My_pos/Ke)-My_pos/Ke);
    fPeakNeg=My_neg+ekhardNeg*((-Thetap_neg+My_neg/Ke)-(My_neg/Ke));
   
    if (cpPos<My_pos/Ke) {
     
      fPeakPos =My_pos*cpPos/(My_pos/Ke);
    }
    if (cpNeg>(My_neg/Ke)) {
      fPeakNeg =My_neg*cpNeg/(My_neg/Ke);
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
   
    fCapRefPos=-capSlope*Ke*(Thetap_pos+My_pos/Ke)+fPeakPos;
    fCapRefNeg=-capSlopeNeg*Ke*(-Thetap_neg+My_neg/Ke)+fPeakNeg;
         
    capSlopeOrig = capSlope;
    capSlopeOrigNeg = capSlopeNeg;
   
    flagstopdeg     = 0;
   
    f1 = 0.0;  
    f2 = 0.0;
   
    fmin = 0.0;
    fmax = 0.0;
   
    RSE = 0.0;
   
    if(deltaD>=0.0){
     
      kon = 1;
     
    } else {
     
      kon = 2;
     
    }
  }
 
  //      ******************* S T A R T S   B I G   L O O P  ****************
  //     IF   D E L T A > 0 - - - - - - - - - - - - - - - - - - - - - - - -  
  if (deltaD>=0.0) {
    if (iNoFpos==1) {
      interPoint(dNewLoadPos,fNewLoadPos,dyNeg,fyNeg,ekhardNeg,0.0,0.0,0.0);
    }
    //If there is a corner changing delta from negative to positive      
    if (kon==2) {
      kon = 1;
     
      if (dP<=dmin){
        fLimNeg = fP;
        dLimNeg = dP;
      }
     
      if (resSn>0.0){
        RSE = 0.5*fP*fP/ekunload;
      } else {
        RSE=0.5*(fP+resSn)*(sn-dP)+0.5*resSn*resSn/(Ke*alphaPos);
      }
     
      double a2 = Enrgtk-(Enrgtot-RSE);
      if((a2<=0.0) && (Enrgtk!=0.0)) {
        flagControlResponse=1;
      }
     
      if((LamdaK!=0.0) && (d<sp) && fabs(capSlope*(1+nFactor))>=1.0e-3 && (fabs(capSlopeNeg*(1+nFactor))>=1.0e-3)){             // Updated: Filipe Ribeiro and Andre Barbosa
       
        betak = pow(((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE))),Ck);
       
        if(((Enrgtot-RSE)>=Enrgtk)||(betak>=1.0)) {
          betak = 1.0;
        }
       
        if( flagstopdeg !=1 && flagdeg !=1) {        
			
				ekunload = Ke*((prodBeta*(1.0-betak))/(1+nFactor*(1-(prodBeta*(1-betak)))));			// Updated: Filipe Ribeiro and Andre Barbosa

        } else {
          ekunload = ekexcurs; // no update happens
        }
       
        if(ekunload<=(0.1/(1+nFactor*(1-0.1)))*Ke) { // Updated: Filipe Ribeiro and Andre Barbosa
          ekunload = (0.1/(1+nFactor*(1-0.1)))*Ke;              // Updated: Filipe Ribeiro and Andre Barbosa
        }
       
      }
     
     
      //// sn CALCULATION-------------------------------------------------
     
      if((dmin<(My_neg/Ke))||(dmax>(My_pos/Ke))) {
        if((dP<sp)||(((dP>sp)&&(ekP==ekunload)))) {
         
          snCalc();
         
          if((fabs(dmax-(My_pos/Ke))>=1.0e-10)&&(fabs(sn-(My_pos/Ke))<=1.0e-10)) {
            sn=(My_pos/Ke)-1.0e-9;
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
    //   LOADING ----------------------------------------------------------
    //      Push envelope
    // Compute force and stiffness for this increment
    if ((iNoFneg==1)&&(iNoFpos==1)&&(d<Thetau_pos)) {
      f = KPos*My_pos;
      ek = 1.0e-7;
    } else if(((iNoFneg==1)&&(iNoFpos==1)&&(d>=Thetau_pos))||flagControlResponse==1) {  //fmk
      f = 1.0e-10;
      ek = 1.0e-7;
      flagstopdeg = 1;
    } else if ((iNoFpos==1)&&(d>sn)&&(d<Thetau_pos)) {
      f = KPos*My_pos;
      ek = 1.0e-7;
    } else if ((iNoFpos==1)&&(d>sn)&&(d>=Thetau_pos)){
      f = 1.0e-10;
      ek =1.0e-7;
      flagstopdeg        = 1;
    } else if (d>=Thetau_pos || flagControlResponse == 1){
      f = 1.0e-10;
      ek =1.0e-7;
      flagstopdeg = 1;
      flagControlResponse = 1;                 // Switches the response to zero strength
    } else if ((iNoFpos==1)&&(d<sn)) {
      ek = ekunload;
      f  = fP+ek*deltaD;
    } else if ((iNoFneg==1)&&(d<dNewLoadNeg)){
      f = 0;
      ek = 1.0e-7;
    } else if (d>dmax) {
      f=0.0;
      envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,f,ek,Ke,My_pos,KPos);
      dmax = d;
      fmax = f;
     
      // c COMPUTE MAXIMUM POSSIBLE DISPLACEMENT
      dBoundPos=boundPos();
      if((d>dBoundPos)||(ek==1.0e-7)) {
        iNoFpos = 1;
      }
     
    } else if (fabs(sn)>1.0e-10) {
     
     
     
      if (LP==0) {
       
       
       if (cpPos<=dyPos) {
         //                             call interPoint(xDevPos1,yDevPos1,sn,resSn,ekhardPos,cpPos,fCapPos,capSlope*Ke)
         interPoint(xDevPos1,yDevPos1,sn,resSn,ekhardPos,cpPos,fCapPos,capSlope*Ke);
         //                             call interPoint(xDevPos2,yDevPos2,sn,resSn,ekunload,cpPos,fCapPos,capSlope*Ke)
         interPoint(xDevPos2,yDevPos2,sn,resSn,ekunload,cpPos,fCapPos,capSlope*Ke);
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
         } else if ((d>sn)&&(d<=xDevPos)) {
           ek = Ke*alphaPos;
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
           
         } else if (d>xDevPos) {
           ek = capSlope*Ke;
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
           
         }
         
       } else if (cpPos > dyPos) {
         interPoint(xDevPos1,yDevPos1,sn,resSn,ekhardPos,cpPos,fCapPos,capSlope*Ke);
         if(d<=sn && d<=xDevPos1) { // Unloading in positive loading direction
           ek = ekunload;
           f  = fP+ek*deltaD;
           
         } else if ((d>sn)&&(d<=cpPos)&& (d<=xDevPos1)) {
           ek = Ke*alphaPos;
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

         } else if((d>cpPos)&&(d<Thetau_pos)) {
           ek = capSlope*Ke;
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
           // c added by Dimitrios to simulate ductile fracture
         } else if(d>=Thetau_pos){
           f=1.0e-10;
           ek = -1.0e-7;
           flagstopdeg = 1;
         }
         // Added by Dimitrios to stay on the residual path when  dresidual <d< dfracture              
         if(d < Thetau_pos && d > cpPos+(KPos*My_pos-fCapPos)/(capSlope*Ke)) {
           f = KPos*My_pos;
           ek = -1.0e-7;
         }
         if( flagControlResponse == 1) { // Response has to be zero since post capping regions hit zero
           f = 1.0e-10;
         }                
       
       }
       
       // c             IF LP IS EQUAL TO 1
      } else if (LP==1) {
        if(d<=sn) {
          ek = ekunload;
          f  = fP+ek*deltaD;
        } else if (((d>sn)&&(sn==snEnv)&&(d<=snHor))||((iNoFneg==1)&&(d>sn)&&(d<snHor))) {
          ek = Ke*alphaPos;
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
        }
      }
      // c       Elastic
    } else {
      if (d>0.0) {
       
        f=0.0;
        envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,f,ek,Ke,My_pos,KPos);
      } else {
       
        f=0.0;
        envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,d,f,ek,Ke,My_neg,KNeg);
      }
    }
   
    //// c       IF   D E L T A < 0 - - - - - - - - - - - - - - - - - - - - - - - -  
   
  } else {
   
    if (iNoFneg==1) {
      interPoint(dNewLoadNeg,fNewLoadNeg,dyPos,fyPos,ekhardPos,0.0,0.0,0.0);
    }
   
   
    // c If there is a corner changing delta from positive to negative ---      
   
    if (kon==1) {
      kon = 2;
     
      if (dP>=dmax) {
        fLimPos = fP;
        dLimPos = dP;
      }
     
      if (resSp<0.0) {
        RSE = 0.5*fP*fP/ekunload;
      } else {
        RSE=0.5*(fP+resSn)*(dP-sp)+0.5*resSp*resSp/(Ke*alphaNeg);
      }
     
     double a2 = Enrgtk-(Enrgtot-RSE);
     
     
     if((a2<=0.0)&&(Enrgtk!=0.0)) {
       flagControlResponse=1;
     }
     //// Update the Unloading Stiffness deterioration
     if((LamdaK!=0.0)&&(d>sn)&&(flagstopdeg!=1)&&(flagdeg!=1)) {
       
		 betak = pow(((Enrgc-RSE)/(Enrgtk-(Enrgtot-RSE))),Ck); 
       
	   if(((Enrgtot-RSE)>=Enrgtk)||(betak>=1.0)) {     
                betak = 1.0;
       }
       // If Post caping slopes have not been flat due to residual update stiffness deterioration       
         
		   if( flagstopdeg !=1 && flagdeg !=1) {  

				ekunload = Ke*((prodBeta*(1.0-betak))/(1+nFactor*(1-(prodBeta*(1-betak)))));		// Updated: Filipe Ribeiro and Andre Barbosa

	   } else { // Keep the same unloading stiffness
         ekunload = ekexcurs;
       }            
       if(ekunload<=(0.1/(1+nFactor*(1-0.1)))*Ke) {             // Updated: Filipe Ribeiro and Andre Barbosa
               ekunload = (0.1/(1+nFactor*(1-0.1)))*Ke; // Updated: Filipe Ribeiro and Andre Barbosa
       }
     }
     
     
     ////     sp CALCULATION----------------------------------------------------
     
     if((dmin<(My_neg/Ke))||(dmax>(My_pos/Ke))) {
       if((dP>sn)||(((dP<sn)&&(ekP==ekunload)))) {
         
         spCalc();
         if((fabs(dmin-(My_neg/Ke))>=1.0e-10)&&(fabs(sp-(My_neg/Ke))<=1.0e-10)) {
           sp=(My_neg/Ke)-1.0e-9;
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
   
   //// UNLOADING
   // c     Push envelope
   
   if ((iNoFneg==1)&&(iNoFpos==1)&&(d>(-Thetau_neg))) {
     f = KNeg*My_neg;
     ek = 1.0e-7;
   } else if (((iNoFneg==1)&&(iNoFpos==1)&&(d<=(-Thetau_neg)))||flagControlResponse==1){
     f = -1.0e-10;
     ek =1.0e-7;
     flagstopdeg = 1;
   } else if ((iNoFneg==1)&&(d<sp)&&(d>(-Thetau_neg))) {
     f = KNeg*My_neg;
     ek = 1.0e-7;
   } else if ((iNoFneg==1)&&(d<sp)&&(d<=(-Thetau_neg))) {
     f = -1.0e-10;
     ek = 1.0e-7;
     flagstopdeg = 1;
   } else if (d<=(-Thetau_neg) || flagControlResponse ==1){
     f = -1.0e-10;
     ek = 1.0e-7;
     flagstopdeg = 1;
     flagControlResponse = 1;               // To control response after I exceed fracture
   } else if ((iNoFneg==1)&&(d>=sp)) {
     ek = ekunload;
     f  = fP+ek*deltaD;
   } else if ((iNoFpos==1)&&(d>dNewLoadPos)) {
     f = 0;
     ek = 1.0e-7;
   } else if (d<dmin) {
     //                 call envelNegCap2(fyNeg,alphaNeg,capSlope,cpNeg,d,f,ek,Ke,My_neg,Resfac)
     f=0.0;
     envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,d,f,ek,Ke,My_neg,KNeg);
     dmin = d;
     fmin = f;
   
     //// c             COMPUTE MINIMUM POSSIBLE DISPLACEMENT
     //                 call boundNeg (dBoundNeg,fCapRefNeg,capSlope,fyPos,dyPos,alphaPos,fCapNeg,cpNeg,Ke)    
     dBoundNeg=boundNeg();
     if((d<dBoundNeg)||(ek==1.0e-7)) {
       iNoFneg = 1;
     }
     
   } else if (fabs(sp)>1.0e-10){
     
     // c               If LN is equal to zero
     if (LN==0) {
       
       if (cpNeg>=dyNeg) {
         //                             call interPoint(xDevNeg1,yDevNeg1,sp,resSp,Ke*alphaNeg,cpNeg,fCapNeg,capSlope*Ke)
         interPoint(xDevNeg1,yDevNeg1,sp,resSp,Ke*alphaNeg,cpNeg,fCapNeg,capSlopeNeg*Ke);
         // call interPoint(xDevNeg2,yDevNeg2,sp,resSp,ekunload,cpNeg,fCapNeg,capSlope*Ke)
         interPoint(xDevNeg2,yDevNeg2,sp,resSp,ekunload,cpNeg,fCapNeg,capSlopeNeg*Ke);
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
         } else if ((d<sp)&&(d>=xDevNeg)) {
           ek = Ke*alphaNeg;
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
         } else if (d<xDevNeg) {
           ek = capSlopeNeg*Ke;
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
           
         }
         
         
       } else if (cpNeg<dyNeg) {
         interPoint(xDevNeg1,yDevNeg1,sp,resSp,Ke*alphaNeg,cpNeg,fCapNeg,capSlopeNeg*Ke);
         if ( d>=sp && d>=xDevNeg1) {
           ek = ekunload;
           f = fP+ek*deltaD;
         } else if((d<sp)&&(d>=cpNeg)&& (d>=xDevNeg1)) {
           ek = Ke*alphaNeg;
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
         } else if((d<cpNeg)&&(d>(-Thetau_neg))) {
           ek = capSlopeNeg*Ke;
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
           // Added by Dimitris to model ductile tearing. Once theta_u((-Thetau_neg)) is exceeded Strength drops to zero and stays
         } else if(d<=(-Thetau_neg) ||  flagControlResponse == 1){
           ek = -1.0e-7;
           f = 1.0e-10;
           flagstopdeg = 1;
           flagControlResponse = 1;               // To dictate the response after passing fracture            
         }
         // Added by Dimitrios to stay on the residual path when  dresidual <d- < dfracture              
         if(d > (-Thetau_neg) && d < cpNeg+(KNeg*My_neg-fCapNeg)/(capSlopeNeg*Ke)) {
           f = KNeg*My_neg;
           ek = -1.0e-7;
         }
         if( flagControlResponse == 1) { // Response has to be zero since post capping regions hit zero
           f = 1.0e-10;
         }
       }
       
     } else if (LN==1){
       if(d>=sp){
         ek = ekunload;
         f  = fP+ek*deltaD;
       } else if (((d<sp)&&(sp==spEnv)&&(d>spHor))||((iNoFpos==1)&&(d<sp)&&(d>spHor))) {
         ek = Ke*alphaNeg;
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
       }
     }
   } else {
     if (d>0.0) {
       f=0.0;
       envelPosCap2(fyPos,alphaPos,capSlope,cpPos,d,f,ek,Ke,My_pos,KPos);
     } else {
       //                      call envelNegCap2 (fyNeg,alphaNeg,capSlope,cpNeg,d,f,ek,Ke,My_neg,Resfac)
           envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,d,f,ek,Ke,My_neg,KNeg);
     }
     
   }
   
  }
  // BIG LOOP Ends ****************************************************
 
  double Enrgi = 0.0;
 
  Enrgi = 0.5 * (f+fP) * deltaD;
  Enrgc = Enrgc + Enrgi;
  Enrgtot = Enrgtot + Enrgi;
 
  RSE = 0.5*f*f/ekunload;

  prodBeta=prodBeta*(1-betak);			// Updated: Filipe Ribeiro and Andre Barbosa
 
  ////    Flag to deteriorate parameters on the opposite side of the loop --    
 
  if((f*fP<0.0)&&(interup==0)) {
    if(((fP>0.0)&&(dmax>(dyPos)))||((fP<0.0)&&(dmin<(dyNeg)))) {
      flagdeg = 1;          
      interup = 1;
    }
  }    
 
 
  ////    energy CALCULATIONS ---------------------------------------------
 
  if((flagstopdeg==0)&&(flagdeg==1)) {
   
    if((Enrgtot>=Enrgts)&&(Enrgts!=0.0)) {
      betas = 1.0;
    } else if((Enrgtot>=Enrgtd)&&(Enrgtd!=0.0)) {
      betad = 1.0;
    } else {
      if(LamdaS!=0.0) {
        betas = pow((Enrgc/(Enrgts-Enrgtot)),Cs);
      }
      if(LamdaD!=0.0) {
        betad = pow((Enrgc/(Enrgtd-Enrgtot)),Cd);
      }
     
      if(fabs(betas)>=1.0) {
        betas = 1.0;
      }
      if(fabs(betad)>=1.0) {
        betad = 1.0;
      }
    }
    ////            Initialize energy of the cycle and Kstif for next loop -------
   
    Enrgc = 0.0;
    ekexcurs = ekunload;
   
    ////            Deteriorate parameters for the next half cycle
    if(deltaD<0.0){
     
      if(flagstopdeg == 0){
        fyNeg = fyNeg*(1-betas*PDNeg); 

        //change the strain hardening ratio						 // Updated: Filipe Ribeiro and Andre Barbosa
        //1st - recover the strain hardening ratio of the member        
        alphaNeg=alphaNeg*(1+nFactor)/(1+nFactor*alphaNeg);             
        //2nd - apply the redution to the ratio of the member           
        alphaNeg=alphaNeg*(1-betas*PDNeg);      
        //3rd - recompute the strain hardening ratio (updated) 
        alphaNeg=(alphaNeg)/(1+nFactor*(1-alphaNeg));			// Updated: Filipe Ribeiro and Andre Barbosa
        
		fCapRefNeg=fCapRefNeg*(1-betad*PDNeg);
      }else{
        fyNeg = fyNeg;
        alphaNeg=alphaNeg;
        fCapRefNeg=fCapRefNeg;
      }
      // When we reach post capping slope goes to zero due to residual
      if(fyNeg>=KNeg*My_neg) { // If strength drops below residual
        fyNeg = KNeg*My_neg;
        //alphaNeg = 10^(-4); // This evaluates to -10 (bitwise XOR)
	alphaNeg = 1.0e-4;
        fCapRefNeg = fyNeg;
        //capSlopeNeg = -pow(10.0,-6);
	capSlopeNeg = -1.0e-6;
        flagstopdeg = 1;
      } else { //% Keep updating the post capping slope

        //change the post-capping ratio					// Updated: Filipe Ribeiro and Andre Barbosa
        //1st - recover post-capping ratio of the member 
        capSlopeNeg=capSlopeOrigNeg*(1+nFactor)/(1+capSlopeOrigNeg*nFactor);   
        //2nd - apply the redution to the ratio of the member   
        capSlopeNeg=capSlopeNeg*(fabs((KNeg*My_neg-fyNeg)/(KNeg*My_neg-My_neg)));      
        //3rd - recompute the post-capping ratio (updated) 
        capSlopeNeg = capSlopeNeg/(1+nFactor*(1-capSlopeNeg));  // Updated: Filipe Ribeiro and Andre Barbosa
        
		if(capSlopeNeg >=0){
		  //capSlopeNeg = -pow(10.0,-6);
		  capSlopeNeg = -1.0e-6;
        }
      }
     
      dyNeg = fyNeg/Ke;
      ekhardNeg=alphaNeg*Ke;
     
      double dCap1Neg=fCapRefNeg/(Ke-capSlopeOrigNeg*Ke);
      double dCap2Neg=(fCapRefNeg+ekhardNeg*dyNeg-fyNeg)/(ekhardNeg-capSlopeOrigNeg*Ke);
      //cpNeg=min(dCap1Neg,dCap2Neg);
      if (dCap1Neg<dCap2Neg)
        {
          cpNeg=dCap1Neg;
        }
      else
        {
          cpNeg=dCap2Neg;
        }
     
      fCapNeg = fCapRefNeg + capSlopeOrigNeg*Ke*cpNeg;
     
      envelNegCap2(fyNeg,alphaNeg,capSlopeNeg,cpNeg,dLimNeg,fLimNeg,ek,Ke,My_neg,KNeg);
      spCalc();
      // c                    In case the degradation point moves from the negative to pos. side
      if(resSp>f) {        
        d = sp;
        f = resSp;
        ek = ekhardNeg;  
      }
    } else {
      if(flagstopdeg == 0){
        fyPos = fyPos*(1-betas*PDPlus);

        //change the strain hardening ratio						 // Updated: Filipe Ribeiro and Andre Barbosa
        //1st - recover the strain hardening ratio of the member        
        alphaPos=alphaPos*(1+nFactor)/(1+nFactor*alphaPos);             
        //2nd - apply the redution to the ratio of the member           
        alphaPos=alphaPos*(1-betas*PDPlus);     
        //3rd - recompute the strain hardening ratio (updated)          
        alphaPos=(alphaPos)/(1+nFactor*(1-alphaPos));				// Updated: Filipe Ribeiro and Andre Barbosa
        
		fCapRefPos=fCapRefPos*(1-betad*PDPlus);
      } else {
        fyPos = fyPos;
        alphaPos=alphaPos;    
        fCapRefPos=fCapRefPos;
      }
               
      //   %If post capping slope goes to zero due to residual:
      if(fyPos <= KPos*My_pos) {  //% If yield Strength Pos drops below residual
        fyPos = KPos*My_pos;
        //alphaPos = pow(10.0,-4);
        alphaPos = 1.0e-4;
        fCapRefPos = fyPos;
        //capSlope = -pow(10.0,-6);
        capSlope = -1.0e-6;
        flagstopdeg = 1;              
      }  else { //% keep updating

        //change the post-capping ratio					 // Updated: Filipe Ribeiro and Andre Barbosa
        //1st - recover post-capping ratio of the member        
        capSlope=capSlopeOrig*(1+nFactor)/(1+capSlopeOrig*nFactor);     
        //2nd - apply the redution to the ratio of the member
        capSlope=capSlope*(fabs((KPos*My_pos-fyPos)/(KPos*My_pos-My_pos)));    
        //3rd - recompute the post-capping ratio (updated)     
        capSlope = capSlope/(1+nFactor*(1-capSlope));   // Updated: Filipe Ribeiro and Andre Barbosa
        
		if(capSlope >=0) {
		  //capSlope = -pow(10.0,-6);
		  capSlope = -1.0e-6;
        }
      }
      dyPos = fyPos/Ke;
      ekhardPos=alphaPos*Ke;
     
      double dCap1Pos=fCapRefPos/(Ke-capSlopeOrig*Ke);
      double dCap2Pos=(fCapRefPos+ekhardPos*dyPos-fyPos)/(ekhardPos-capSlopeOrig*Ke);
      //cpPos=max(dCap1Pos,dCap2Pos);
      if(dCap1Pos>dCap2Pos)
        {
          cpPos=dCap1Pos;
        }
      else
        {
          cpPos=dCap2Pos;
        }
      fCapPos = fCapRefPos + capSlopeOrig*Ke*cpPos;
     
      envelPosCap2(fyPos,alphaPos,capSlope,cpPos,dLimPos,fLimPos,ek,Ke,My_pos,KPos);
     
      snCalc();
      // c                    In case the degradation point moves from the pos. to neg. side
      if(resSn<f) {
        d = sn;
        f = resSn;
        ek = ekhardPos;
       
      }  
    }
  }
 
  // c            Check the horizontal limit in case that dBound is reached after first neg slope
  if ((d<0)&&(fabs(ek)<=1.0e-7)) {
    LN = 1;
  }
  if ((d>0)&&(fabs(ek)<=1.0e-7)) {
    LP = 1;
  }
 
 
  // c    Update envelope values --------------------------------------------    
 
  f1 = f;  
 
 
  // c    Updating parameters for next cycle ---------------------------------
  ekP = ek;
  fP = f;
  dP = d;
  Tangent=ek;
 
  //priority of logical operators
  if ((interup==1&&(ek==Ke*alphaPos)) ||(ek==Ke*alphaNeg)||(ek==capSlope*Ke) ||(ek==capSlopeNeg*Ke)) { //fmk
    interup = 0;
  }
 
  return 0;
}
double
Bilin::getStress(void)
{
   return (fP);
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
  commitCalledOnce = 1;
  
  //commit trial  variables
   CU=U;
   CTangent=Tangent;

   CdNewLoadPos=dNewLoadPos;
   CdNewLoadNeg=dNewLoadNeg;
   Cflagdeg=flagdeg ;
   Cflagstopdeg=flagstopdeg;
   Cinterup=interup;  
   Ckon=kon;
   CiNoFneg=iNoFneg;
   CiNoFpos=iNoFpos;
   CLP=LP;  //fmk
   CLN=LN;
   CcapSlope=capSlope;
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
   CdLimPos=dLimPos;
   CdLimNeg=dLimNeg;  
   CcpPos=cpPos;
   CcpNeg=cpNeg;
   CfLimPos=fLimPos;
   Cekexcurs=ekexcurs;
   CRSE=RSE;
   CfPeakPos=fPeakPos;
   CfPeakNeg=fPeakNeg;
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
   CcapSlopeOrig=capSlopeOrig;  
   CcapSlopeNeg=capSlopeNeg;
   CflagControlResponse=flagControlResponse;
   CcapSlopeOrigNeg=capSlopeOrigNeg;
   CcapSlopeNegMember=capSlopeNegMember; // Updated: Filipe Ribeiro and Andre Barbosa
   CcapSlopeMember=capSlopeMember;// Updated: Filipe Ribeiro and Andre Barbosa
   CKe=Ke;  // Updated: Filipe Ribeiro and Andre Barbosa
   CprodBeta=prodBeta;	// Updated: Filipe Ribeiro and Andre Barbosa
   return 0;
}

int
Bilin::revertToLastCommit(void)
{
        //the opposite of commit trial history variables
   U=CU;
   Tangent=CTangent;

   dNewLoadPos=CdNewLoadPos;
   dNewLoadNeg=CdNewLoadNeg;
   flagdeg=Cflagdeg ;
   flagstopdeg=Cflagstopdeg;
   interup=Cinterup;  
   kon=Ckon;
   iNoFneg=CiNoFneg;
   iNoFpos=CiNoFpos;
   LP=CLP;
   LN=CLN;
   capSlope=CcapSlope;
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
   dLimPos=CdLimPos;
   dLimNeg=CdLimNeg;  
   cpPos=CcpPos;
   cpNeg=CcpNeg;
   fLimPos=CfLimPos;
   ekexcurs=Cekexcurs;
   RSE=CRSE;
   fPeakPos=CfPeakPos;
   fPeakNeg=CfPeakNeg;
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
   capSlopeOrig=CcapSlopeOrig;  
   capSlopeNeg=CcapSlopeNeg;
   flagControlResponse=CflagControlResponse;
   capSlopeOrigNeg=CcapSlopeOrigNeg;
   capSlopeMember=CcapSlopeMember;			// Updated: Filipe Ribeiro and Andre Barbosa
   capSlopeNegMember=CcapSlopeNegMember;        // Updated: Filipe Ribeiro and Andre Barbosa
   Ke=CKe;					 // Updated: Filipe Ribeiro and Andre Barbosa
   prodBeta=CprodBeta;			// Updated: Filipe Ribeiro and Andre Barbosa
    return 0;
}

int
Bilin::revertToStart(void)
{
//initially I zero everything
   U=CU=0.0;
   Tangent=CTangent=0.0;
   commitCalledOnce = 0;

   dNewLoadPos=CdNewLoadPos=0.0;
   dNewLoadNeg=CdNewLoadNeg=0.0;
   flagdeg=Cflagdeg=0;
   flagstopdeg=Cflagstopdeg=0;
   interup=Cinterup=0;  
   kon=Ckon=0;
   iNoFneg=CiNoFneg=0;
   iNoFpos=CiNoFpos=0;
   LP=CLP=0;
   LN=CLN=0;
   capSlope=CcapSlope=0.0;
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
   dLimPos=CdLimPos=0.0;
   dLimNeg=CdLimNeg=0.0;  
   cpPos=CcpPos=0.0;
   cpNeg=CcpNeg=0.0;
   fLimPos=CfLimPos=0.0;
   ekexcurs=Cekexcurs=0.0;
   RSE=CRSE=0.0;
   fPeakPos=CfPeakPos=0.0;
   fPeakNeg=CfPeakNeg=0.0;
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
   capSlopeOrig=CcapSlopeOrig=0.0;  
   capSlopeNeg=CcapSlopeNeg=0.0;
   flagControlResponse=CflagControlResponse=0;
   capSlopeOrigNeg=CcapSlopeOrigNeg=0.0;
   capSlopeMember=CcapSlopeMember=0.0;  // Updated: Filipe Ribeiro and Andre Barbosa
   capSlopeNegMember=CcapSlopeNegMember=0.0;    // Updated: Filipe Ribeiro and Andre Barbosa
   Ke=CKe=0.0;          // Updated: Filipe Ribeiro and Andre Barbosa
   prodBeta=CprodBeta=1.0;			// Updated: Filipe Ribeiro and Andre Barbosa
                    
    return 0;
}

UniaxialMaterial *
Bilin::getCopy(void)
{
    Bilin *theCopy = new Bilin(this->getTag(),Ke0,AsPos,AsNeg,My_pos,My_neg,LamdaS,LamdaD,             
                LamdaA,LamdaK,Cs,Cd,Ca,Ck,Thetap_pos,Thetap_neg,Thetapc_pos,Thetapc_neg,KPos,KNeg,
                Thetau_pos,Thetau_neg,PDPlus,PDNeg,nFactor);  // Updated: Filipe Ribeiro and Andre Barbosa

    //Fixed Model parameters: need to change according to material properties
   theCopy->U=U;
   theCopy->CU=CU;
   theCopy->Tangent=Tangent;
   theCopy->CTangent=CTangent;
   theCopy->dNewLoadPos=dNewLoadPos;
   theCopy->CdNewLoadPos=CdNewLoadPos;
   theCopy->dNewLoadNeg=dNewLoadNeg;
   theCopy->CdNewLoadNeg=CdNewLoadNeg;
   theCopy->flagdeg=flagdeg;
   theCopy->Cflagdeg=Cflagdeg;
   theCopy->flagstopdeg=flagstopdeg;
   theCopy->Cflagstopdeg=Cflagstopdeg;
   theCopy->interup=interup;
   theCopy->Cinterup=Cinterup;  
   theCopy->kon=kon;
   theCopy->Ckon=Ckon;
   theCopy->iNoFneg=iNoFneg;
   theCopy->CiNoFneg=CiNoFneg;
   theCopy->iNoFpos=iNoFpos;
   theCopy->CiNoFpos=CiNoFpos;
   theCopy->LP=LP;
   theCopy->CLP=CLP;
   theCopy->LN=LN;
   theCopy->CLN=CLN;
   theCopy->capSlope=capSlope;
   theCopy->CcapSlope=CcapSlope;
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
   theCopy->dLimPos=dLimPos;
   theCopy->CdLimPos=CdLimPos;
   theCopy->dLimNeg=dLimNeg;
   theCopy->CdLimNeg=CdLimNeg;  
   theCopy->cpPos=cpPos;
   theCopy->CcpPos=CcpPos;
   theCopy->cpNeg=cpNeg;
   theCopy->CcpNeg=CcpNeg;
   theCopy->fLimPos=fLimPos;
   theCopy->CfLimPos=CfLimPos;
   theCopy->ekexcurs=ekexcurs;
   theCopy->Cekexcurs=Cekexcurs;
   theCopy->RSE=RSE;
   theCopy->CRSE=CRSE;
   theCopy->fPeakPos=fPeakPos;
   theCopy->CfPeakPos=CfPeakPos;
   theCopy->fPeakNeg=fPeakNeg;
   theCopy->CfPeakNeg=CfPeakNeg;
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
   theCopy->capSlopeOrig=capSlopeOrig;
   theCopy->CcapSlopeOrig=CcapSlopeOrig;  
   theCopy->capSlopeNeg=capSlopeNeg;
   theCopy->CcapSlopeNeg=CcapSlopeNeg;
   theCopy->flagControlResponse=flagControlResponse;
   theCopy->CflagControlResponse=CflagControlResponse;  
   theCopy->capSlopeMember=capSlopeMember;      // Updated: Filipe Ribeiro and Andre Barbosa
   theCopy->CcapSlopeMember=CcapSlopeMember;// Updated: Filipe Ribeiro and Andre Barbosa
   theCopy->capSlopeNegMember=capSlopeNegMember;        // Updated: Filipe Ribeiro and Andre Barbosa
   theCopy->CcapSlopeNegMember=CcapSlopeNegMember;      // Updated: Filipe Ribeiro and Andre Barbosa
   theCopy->CKe=CKe;  // Updated: Filipe Ribeiro and Andre Barbosa
   theCopy->capSlopeOrigNeg=capSlopeOrigNeg;
   theCopy->CcapSlopeOrigNeg=CcapSlopeOrigNeg;
   theCopy->CprodBeta=CprodBeta;	// Updated: Filipe Ribeiro and Andre Barbosa

    return theCopy;
}

int
Bilin::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(166);      // Updated: Filipe Ribeiro and Andre Barbosa
   data(0) = this->getTag();
   data(1)=Ke0;         // Updated: Filipe Ribeiro and Andre Barbosa
   data(2)=AsPos;
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
   data(18)=KPos;
   data(19)=KNeg;
   data(20)=Thetau_pos;
   data(21)=Thetau_neg;
   data(22)=PDPlus;
   data(23)=PDNeg;
   data(24)=U;
   data(25)=CU;
   data(26)=Tangent;
   data(27)=CTangent;
   data(28)=dNewLoadPos;
   data(29)=CdNewLoadPos;
   data(30)=dNewLoadNeg;
   data(31)=CdNewLoadNeg;
   data(32)=flagdeg;
   data(33)=Cflagdeg;
   data(34)=flagstopdeg;
   data(35)=Cflagstopdeg;
   data(36)=interup;
   data(37)=Cinterup;  
   data(38)=kon;
   data(39)=Ckon;
   data(40)=iNoFneg;
   data(41)=CiNoFneg;
   data(42)=iNoFpos;
   data(43)=CiNoFpos;
   data(44)=LP;
   data(45)=CLP;
   data(46)=LN;
   data(47)=CLN;
   data(48)=capSlope;
   data(49)=CcapSlope;
   data(50)=dmax;
   data(51)=Cdmax;
   data(52)=dmin;
   data(53)=Cdmin;
   data(54)=Enrgtot;
   data(55)=CEnrgtot;
   data(56)=Enrgc;
   data(57)=CEnrgc;
   data(58)=fyPos;
   data(59)=CfyPos;
   data(60)=fLimNeg;
   data(61)=CfLimNeg;
   data(62)=fyNeg;
   data(63)=CfyNeg;
   data(64)=ekP;
   data(65)=CekP;
   data(66)=ekunload;
   data(67)=Cekunload;
   data(162) = sp; //fmk
   data(68)=Csp;
   data(69)=sn;
   data(70)=Csn;
   data(71)=dP;
   data(72)=CdP;
   data(73)=fP;
   data(74)=CfP;
   data(75)=ek;
   data(76)=Cek;
   data(77)=dLimPos;
   data(78)=CdLimPos;
   data(79)=dLimNeg;
   data(80)=CdLimNeg;  
   data(81)=cpPos;
   data(82)=CcpPos;
   data(83)=cpNeg;
   data(84)=CcpNeg;
   data(85)=fLimPos;
   data(86)=CfLimPos;
   data(87)=ekexcurs;
   data(88)=Cekexcurs;
   data(89)=RSE;
   data(90)=CRSE;
   data(91)=fPeakPos;
   data(92)=CfPeakPos;
   data(93)=fPeakNeg;
   data(94)=CfPeakNeg;
   data(95)=alphaNeg;
   data(96)=CalphaNeg;
   data(97)=alphaPos;
   data(98)=CalphaPos;
   data(99)=ekhardNeg;
   data(100)=CekhardNeg;
   data(101)=ekhardPos;
   data(102)=CekhardPos;
   data(103)=fCapRefPos;
   data(104)=CfCapRefPos;
   data(105)=fCapRefNeg;
   data(106)=CfCapRefNeg;
   data(107)=Enrgts;
   data(108)=CEnrgts;
   data(109)=Enrgtk;
   data(110)=CEnrgtk;
   data(111)=Enrgtd;
   data(112)=CEnrgtd;
   data(113)=dyPos;
   data(114)=CdyPos;
   data(115)=dyNeg;
   data(116)=CdyNeg;
   data(117)=resSnHor;
   data(118)=CresSnHor;
   data(119)=fmax;
   data(120)=Cfmax;
   data(121)=fmin;
   data(122)=Cfmin;
   data(123)=resSp;
   data(124)=CresSp;
   data(125)=resSn;
   data(126)=CresSn;
   data(127)=fCapPos;
   data(128)=CfCapPos;
   data(129)=fCapNeg;
   data(130)=CfCapNeg;
   data(131)=snHor;
   data(132)=CsnHor;
   data(133)=spHor;
   data(134)=CspHor;
   data(135)=resSpHor;
   data(136)=CresSpHor;
   data(137)=snEnv;
   data(138)=CsnEnv;
   data(139)=resSnEnv;
   data(140)=CresSnEnv;
   data(141)=spEnv;
   data(142)=CspEnv;
   data(143)=resSpEnv;
   data(144)=CresSpEnv;

   data(145)=capSlopeOrig;
   data(146)=CcapSlopeOrig;  
   data(147)=capSlopeNeg;
   data(148)=CcapSlopeNeg;

   data(163) = capSlopeOrigNeg; //fmk
   data(164) = CcapSlopeOrigNeg; //fmk
   data(165) = Ke;    //fmk

   data(149)=flagControlResponse;
   data(150)=CflagControlResponse;
   data(151)=sp;
   data(152)=CcapSlopeOrigNeg;
   data(153)=capSlopeOrigNeg;
   data(154)=nFactor;				// Updated: Filipe Ribeiro and Andre Barbosa

   data(155)=capSlopeMember;		// Updated: Filipe Ribeiro and Andre Barbosa
   data(156)=CcapSlopeMember;		// Updated: Filipe Ribeiro and Andre Barbosa
   data(157)=capSlopeNegMember;		// Updated: Filipe Ribeiro and Andre Barbosa
   data(158)=CcapSlopeNegMember;    // Updated: Filipe Ribeiro and Andre Barbosa
   data(159)=CKe;					// Updated: Filipe Ribeiro and Andre Barbosa
   data(160) = prodBeta;
   data(161)=CprodBeta;				// Updated: Filipe Ribeiro and Andre Barbosa
  
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
  static Vector data(166);  // Updated: Filipe Ribeiro and Andre Barbosa
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
 
  if (res < 0) {
      opserr << "Bilin::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
   Ke0=data(1);				// Updated: Filipe Ribeiro and Andre Barbosa
   AsPos=data(2);
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
   KPos=data(18);
   KNeg=data(19);
   Thetau_pos=data(20);
   Thetau_neg=data(21);
   PDPlus=data(22);
   PDNeg=data(23);
   U=data(24);
   CU=data(25);
   Tangent=data(26);
   CTangent=data(27);
   dNewLoadPos=data(28);
   CdNewLoadPos=data(29);
   dNewLoadNeg=data(30);
   CdNewLoadNeg=data(31);
   flagdeg=data(32);
   Cflagdeg=data(33);
   flagstopdeg=data(34);
   Cflagstopdeg=data(35);
   interup=data(36);
   Cinterup=data(37);  
   kon=data(38);
   Ckon=data(39);
   iNoFneg=data(40);
   CiNoFneg=data(41);
   iNoFpos=data(42);
   CiNoFpos=data(43);
   LP=data(44);
   CLP=data(45);
   LN=data(46);
   CLN=data(47);
   capSlope=data(48);
   CcapSlope=data(49);
   dmax=data(50);
   Cdmax=data(51);
   dmin=data(52);
   Cdmin=data(53);
   Enrgtot=data(54);
   CEnrgtot=data(55);
   Enrgc=data(56);
   CEnrgc=data(57);
   fyPos=data(58);
   CfyPos=data(59);
   fLimNeg=data(60);
   CfLimNeg=data(61);
   fyNeg=data(62);
   CfyNeg=data(63);
   ekP=data(64);
   CekP=data(65);
   ekunload=data(66);
   Cekunload=data(67);
   sp = data(162);
   Csp=data(68);
   sn=data(69);
   Csn=data(70);
   dP=data(71);
   CdP=data(72);
   fP=data(73);
   CfP=data(74);
   ek=data(75);
   Cek=data(76);
   dLimPos=data(77);
   CdLimPos=data(78);
   dLimNeg=data(79);
   CdLimNeg=data(80);  
   cpPos=data(81);
   CcpPos=data(82);
   cpNeg=data(83);
   CcpNeg=data(84);
   fLimPos=data(85);
   CfLimPos=data(86);
   ekexcurs=data(87);
   Cekexcurs=data(88);
   RSE=data(89);
   CRSE=data(90);
   fPeakPos=data(91);
   CfPeakPos=data(92);
   fPeakNeg=data(93);
   CfPeakNeg=data(94);
   alphaNeg=data(95);
   CalphaNeg=data(96);
   alphaPos=data(97);
   CalphaPos=data(98);
   ekhardNeg=data(99);
   CekhardNeg=data(100);
   ekhardPos=data(101);
   CekhardPos=data(102);
   fCapRefPos=data(103);
   CfCapRefPos=data(104);
   fCapRefNeg=data(105);
   CfCapRefNeg=data(106);
   Enrgts=data(107);
   CEnrgts=data(108);
   Enrgtk=data(109);
   CEnrgtk=data(110);
   Enrgtd=data(111);
   CEnrgtd=data(112);
   dyPos=data(113);
   CdyPos=data(114);
   dyNeg=data(115);
   CdyNeg=data(116);
   resSnHor=data(117);
   CresSnHor=data(118);
   fmax=data(119);
   Cfmax=data(120);
   fmin=data(121);
   Cfmin=data(122);
   resSp=data(123);
   CresSp=data(124);
   resSn=data(125);
   CresSn=data(126);
   fCapPos=data(127);
   CfCapPos=data(128);
   fCapNeg=data(129);
   CfCapNeg=data(130);
   snHor=data(131);
   CsnHor=data(132);
   spHor=data(133);
   CspHor=data(134);
   resSpHor=data(135);
   CresSpHor=data(136);
   snEnv=data(137);
   CsnEnv=data(138);
   resSnEnv=data(139);
   CresSnEnv=data(140);
   spEnv=data(141);
   CspEnv=data(142);
   resSpEnv=data(143);
   CresSpEnv=data(144);

   capSlopeOrig=data(145);
   CcapSlopeOrig=data(146);  
   capSlopeNeg=data(147);
   CcapSlopeNeg=data(148);

   capSlopeOrigNeg = data(163); //fmk
   CcapSlopeOrigNeg = data(164); //fmk
   Ke = data(165);    //fmk

   flagControlResponse=data(149);
   CflagControlResponse=data(150);
   sp=data(151);
   CcapSlopeOrigNeg=data(152);
   capSlopeOrigNeg=data(153);
   nFactor=data(154);		        // Updated: Filipe Ribeiro and Andre Barbosa
   capSlopeMember=data(155);		// Updated: Filipe Ribeiro and Andre Barbosa
   CcapSlopeMember=data(156);		// Updated: Filipe Ribeiro and Andre Barbosa
   capSlopeNegMember=data(157);		// Updated: Filipe Ribeiro and Andre Barbosa
   CcapSlopeNegMember=data(158);        // Updated: Filipe Ribeiro and Andre Barbosa
   CKe=data(159);			// Updated: Filipe Ribeiro and Andre Barbosa
   prodBeta=data(160);		        // Updated: Filipe Ribeiro and Andre Barbosa
   CprodBeta = data(161);

  }

  return res;
}

void
Bilin::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Bilin tag: " << this->getTag() << endln;
        s << "Ke0: " << Ke0 << ", ";
        s << "AsPos: " << AsPos << ", ";
        s << "AsNeg: " << AsNeg << ", ";
        s << "My_pos: " << My_pos << ", ";
        s << "My_neg: " << My_neg << ", ";
        s << "LamdaS: " << LamdaS << ", ";
        s << "LamdaK: " << LamdaK << ", ";
        s << "LamdaA: " << LamdaA << ", ";
        s << "LamdaD: " << LamdaD << ", ";
        s << "Cs: " << Cs << ", ";
        s << "Ck: " << Ck << ", ";
        s << "Ca: " << Ca << ", ";
        s << "Cd: " << Cd << ", ";
        s << "Thetap_pos: " << Thetap_pos << ", ";
        s << "Thetap_neg: " << Thetap_neg << ", ";
        s << "Thetapc_pos: " << Thetapc_pos << ", ";
        s << "Thetapc_neg: " << Thetapc_neg << ", ";
        s << "KPos: " << KPos << ", ";
        s << "KNeg: " << KNeg << ", ";
        s << "Thetau_pos: " << Thetau_pos << ", ";
        s << "Thetau_neg: " << Thetau_neg << ", ";
        s << "PDPlus: " << PDPlus << ", ";
        s << "PDNeg: " << PDNeg << ", ";
        s << "nFactor: " << nFactor;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"Bilin\", ";
        s << "\"Ke0\": " << Ke0 << ", ";
        s << "\"AsPos\": " << AsPos << ", ";
        s << "\"AsNeg\": " << AsNeg << ", ";
        s << "\"My_pos\": " << My_pos << ", ";
        s << "\"My_neg\": " << My_neg << ", ";
        s << "\"LamdaS\": " << LamdaS << ", ";
        s << "\"LamdaK\": " << LamdaK << ", ";
        s << "\"LamdaA\": " << LamdaA << ", ";
        s << "\"LamdaD\": " << LamdaD << ", ";
        s << "\"Cs\": " << Cs << ", ";
        s << "\"Ck\": " << Ck << ", ";
        s << "\"Ca\": " << Ca << ", ";
        s << "\"Cd\": " << Cd << ", ";
        s << "\"Thetap_pos\": " << Thetap_pos << ", ";
        s << "\"Thetap_neg\": " << Thetap_neg << ", ";
        s << "\"Thetapc_pos\": " << Thetapc_pos << ", ";
        s << "\"Thetapc_neg\": " << Thetapc_neg << ", ";
        s << "\"KPos\": " << KPos << ", ";
        s << "\"KNeg\": " << KNeg << ", ";
        s << "\"Thetau_pos\": " << Thetau_pos << ", ";
        s << "\"Thetau_neg\": " << Thetau_neg << ", ";
        s << "\"PDPlus\": " << PDPlus << ", ";
        s << "\"PDNeg\": " << PDNeg << ", ";
        s << "\"nFactor\": " << nFactor << "}";
    }
}

Response *Bilin::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  // See if the response is one of the defaults
  Response *theResponse = UniaxialMaterial::setResponse(argv, argc, theOutput);

  if (theResponse != 0)
    return theResponse;

  if (strcmp(argv[0], "RSE") == 0)
  {
    theOutput.tag("ResponseType", "RSE");
    theResponse = new MaterialResponse(this, 101, CRSE);
  }

  theOutput.endTag();
  return theResponse;
}

int Bilin::getResponse(int responseID, Information &matInformation)
{
  if (responseID == 101)
  {
    matInformation.setDouble(CRSE);
  }
  else
  {
    return UniaxialMaterial::getResponse(responseID, matInformation);
  }
  return 0;
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
// c    Each time that the hysteretic loop changes from negative to positive
// c    delta displacement, this subroutine calculates the point where ends
// c    ekunload and starts     whether the positive hardening curve or the
// c    positive cap curve. In cases where "cpPos<fyPos" this point could
// c    not be reached, because the unloading stiffness could intersect the
// c    positive cap curve. However, this change is reflected in the main
// c    program.
// c  
// c    Output Variables: sn,resSn,snHard,resSnHard
// c    Input Variables:  dP,fP,ekunload,alphaPos,dyPos,fyPos,cpPos,fCapPos
// c                              capStiff,fCapRefPos

        double Resid = KPos*fyPos;
        double dresid = cpPos+(Resid-fCapPos)/(capSlope*Ke);   
        double ekresid = 1.0e-10;
        dyPos = fyPos/Ke;    
    double snHard,resSnHard,snLim,resSnLim,snResid,resSnResid;
        if (dyPos<cpPos) {
               
    interPoint(snHard,resSnHard,dyPos,fyPos,Ke*alphaPos,dP,fP,ekunload);        
        } else {
   
    interPoint(snHard,resSnHard,cpPos,fCapPos,Ke*alphaPos,dP,fP,ekunload);      
        }

        double snCap,resSnCap;
    interPoint(snCap,resSnCap,0.0,fCapRefPos,capSlope*Ke,dP,fP,ekunload);      

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
        if((LP==1)&&(fLimPos==0.0)) {
                //call interPoint(snLim,resSnLim,dLimPos,fLimPos,0.d0,dP,fP,ekunload)
        interPoint(snLim,resSnLim,dLimPos,fLimPos,0.0,dP,fP,ekunload);
                if (snLim<sn){  
                        sn=snLim;
                        resSn=resSnLim;
                }

     interPoint(snHor,resSnHor,dLimPos,fLimPos,0.0,dyPos,fyPos,Ke*alphaPos);           
        }

        if (sn>dresid) {
//              call interPoint(snResid,resSnResid,dresid,Resid,ekresid,dP,fP,ekunload)
        interPoint(snResid,resSnResid,dresid,Resid,ekresid,dP,fP,ekunload);
                sn = snResid;
                resSn = resSnResid;
        }
}


void
Bilin::envelPosCap2(double fy,double alphaPos,double alphaCap,double cpDsp,double& d,
                                  double& f,double& ek,double elstk,double fyieldPos,double Resfac)
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
                if(d>=Thetau_pos) {
                ek = 1.0e-7;
                        f = 1.0e-10;
                        d=Thetau_pos;
            flagControlResponse=1;
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
                if(d>=Thetau_pos) {
                ek = 1.0e-7;
                        f = 1.0e-10;
                        d=Thetau_pos;
            flagControlResponse=1;  
                }
               
        } else {
        }
}
double
Bilin::boundPos(void)
{

    double Resid=0.0;    
     double dBoundPos;
     dyNeg = fyNeg/Ke;                 
        double dresid = cpPos+(Resid-fCapPos)/(capSlope*Ke);   
        double ekresid = 1.0e-10;

//      call interPoint(d1,f1,E
        double d1,f1;
    interPoint(d1,f1,dyNeg,fyNeg,Ke*alphaNeg,0.0,fCapRefPos,capSlope*Ke);

        double d2,f2;
    interPoint(d2,f2,dyNeg,fyNeg,Ke*alphaNeg,dresid,Resid,ekresid);            
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
                if(d<=(-Thetau_neg)) {
                ek = 1.0e-7;
                        f = 1.0e-10;
                        d=(-Thetau_neg);
            flagControlResponse=1;
                }
       
        } else if(dy<cpDsp) {
       
                rcap = elstk*cpDsp;
                Res = Resfac*rcap;
                dres = cpDsp+(Res-rcap)/(alphaCap*elstk);

                if (d>0.0) {
                        f = 0.0;
                        ek = 1.0e-7;
                }       else if (d>=cpDsp) {
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
       if(d<=(-Thetau_neg)) {
                ek = 1.0e-7;
                        f = 1.0e-10;
                        d=(-Thetau_neg);
            flagControlResponse=1;
           }

        } else {
        }
}
void
Bilin::spCalc(void)
{
        double Resid = KNeg*fyNeg;
        dyNeg = fyNeg/Ke;    
        double dresid = cpNeg+(Resid-fCapNeg)/(capSlopeNeg*Ke);                
        double ekresid = 1.0e-10;
     double spHard,resSpHard,spCap,resSpCap,spLim,resSpLim,spResid,resSpResid;
        if (dyNeg>cpNeg) {
//      call interPoint(spHard,resSpHard,dyNeg,fyNeg,Ke*alphaNeg,dP,fP,ekunload)
    interPoint(spHard,resSpHard,dyNeg,fyNeg,Ke*alphaNeg,dP,fP,ekunload);
        } else {
//      call interPoint(spHard,resSpHard,cpNeg,fCapNeg,Ke*alphaNeg,dP,fP,ekunload)
    interPoint(spHard,resSpHard,cpNeg,fCapNeg,Ke*alphaNeg,dP,fP,ekunload);
        }

//      call interpoint(spCap,resSpCap,0.d0,fCapRefNeg,capSlope*Ke,dP,fP,ekunload)
    interPoint(spCap,resSpCap,0.0,fCapRefNeg,capSlopeNeg*Ke,dP,fP,ekunload);           
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
//              call interPoint(spLim,resSpLim,dLimNeg,fLimNeg,0.d0,dP,fP,ekunload)
        interPoint(spLim,resSpLim,dLimNeg,fLimNeg,0.0,dP,fP,ekunload);
                if (spLim>sp) {
                        sp=spLim;
                        resSp=resSpLim;
                }

//              call interPoint(spHor,resSpHor,dLimNeg,fLimNeg,0.d0,dyNeg,fyNeg,Ke*alphaNeg)
         interPoint(spHor,resSpHor,dLimNeg,fLimNeg,0.0,dyNeg,fyNeg,Ke*alphaNeg);
        }

        if (sp<dresid) {
//              call interPoint(spResid,resSpResid,dresid,Resid,ekresid,dP,fP,ekunload)
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
  dyPos = fyPos/Ke;      
  double dresid = cpNeg+(Resid-fCapNeg)/(capSlopeNeg*Ke);
  double ekresid = 1.0e-10;
  double d1,f1,d2,f2;
  //      call interPoint(d1,f1,dyPos,fyPos,Ke*alphaPos,0.d0,fCapRefNeg,capSlope*Ke)
  interPoint(d1,f1,dyPos,fyPos,Ke*alphaPos,0.0,fCapRefNeg,capSlopeNeg*Ke);     
  //%
  //      call interPoint(d2,f2,dyPos,fyPos,Ke*alphaPos,dresid,resid,ekresid)
  interPoint(d2,f2,dyPos,fyPos,Ke*alphaPos,dresid,Resid,ekresid);  
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
