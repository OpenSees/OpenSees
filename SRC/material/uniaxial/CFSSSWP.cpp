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
// $Date: 12-10-2015 $

// Written by Smail KECHIDI, Ph.D student at University of Blida 1 (s_kechidi@univ-blida.dz), PhD mobility Student at University of Porto FEUP (smail.kechidi@fe.up.pt)
// Created: 12-10-2015 12:24:20 $
//
// Description: This file contains the class implementation for CFSSSWP
// CFSSSWP is based on Pinching4 uniaxialMaterial

#include <elementAPI.h>
#include "CFSSSWP.h"
#include <OPS_Globals.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <OPS_Stream.h>
#include <stdio.h>
#include <string.h>
#include "CubicSpline.h"
#include "TriMatrix.h"
 

static int numCFSSSWP = 0;

#ifdef _WIN32
#define isnan _isnan
#endif

void *
OPS_CFSSSWP(void)
{
  // print out some KUDO's
  if (numCFSSSWP == 0) {
    opserr << "Cold Formed Steel Steel-Sheathed Shear Wall Panel uniaxialMaterial - Written by Smail KECHIDI Ph.D Student at University of Blida 1 - Please when using this make reference as: Smail Kechidi and Nouredine Bourahla (2016), Deteriorating hysteresis model for cold-formed steel shear wall panel based on its physical and mechanical characteristics, Journal of Thin-Walled Structures, DOI: 10.1016/j.tws.2015.09.022\n";
    numCFSSSWP =1;
  }

  opserr << "Due to known issues and unreliable results, this material has been" << endln;
  opserr << "temporarily removed from the compiled versions of OpenSees (Tcl and Py)" << endln;
  opserr << "The material source code remains available. Compile at your own risk." << endln;
  return 0;    
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int    iData[1];
  double dData[16];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial CFSSSWP tag" << endln;
    return 0;
  }

  numData = 15;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Material parameters\n";
    return 0;	
  }

  // 
  // create a new material
  //

  theMaterial = new CFSSSWP(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]
  ,dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type CFSSSWP\n";
    return 0;
  }

  // return the material
  return theMaterial;
}

CFSSSWP::CFSSSWP(int tag, double H, int B, double fuf, double fyf,
                          double tf,double Af,double fus, double fys, double ts,
						  double np, double ds, double Vs,double sc, double A, double L): 
             UniaxialMaterial(tag, MAT_TAG_CFSSSWP), hight(H), width(B), fuf(fuf),
             fyf(fyf), tf(tf),
	         Af(Af), fus(fus), fys(fys), ts(ts),
             np(np), ds(ds), Vs(Vs),
             screw_Spacing(sc), A(A), L(L),
      envlpPosStress(7), envlpPosStrain(7), envlpNegStress(7), envlpNegStrain(7), tagMat(tag),
      gammaDLimit(0.0),
	  gammaFLimit(0.0),
	  gammaE(10.0),
	  TnCycle(0.0), CnCycle(0.0),
	  rDispP(0.488), rForceP(0.183), uForceP(-0.08), rDispN(0.488), rForceN(0.244), uForceN(-0.08),
	  state3Stress(4), state3Strain(4), state4Stress(4), state4Strain(4), 
	  envlpPosDamgdStress(7), envlpNegDamgdStress(7)

    {

	 double ddeg;
	 ddeg = 0.1*((hight/(2*width))*(screw_Spacing/152.0));
	 gammaDLimit = ddeg;
	
	        // set envelope slopes
	        this->lateralShearStrength();
	        this->SetEnvelope();
	        envlpPosDamgdStress = envlpPosStress; envlpNegDamgdStress = envlpNegStress;
	        state3Stress.Zero(); state3Strain.Zero(); state4Stress.Zero(); state4Strain.Zero();

	        // Initialize history variables
	        this->revertToStart();
	        this->revertToLastCommit();
}

void CFSSSWP :: lateralShearStrength(void) 
{
	    double Alpha,Alpha1,Alpha2,Beta,Beta1,Beta2,Beta3,Lambda,Wmax,Pns,Pns1,
		Pns2,Pns3,Pnsed,We,rho,V,V1,V2,Gs,Omega1,Omega2,Omega3,Omega4,Delta1,
        Delta2,Delta3,Delta4,DeltaV,MinPns,MinPns1,MinPns2,N,Pn;
		Pns=0;
		MinPns=0;
        double mu=0.3;
		E=203000.00;
        Dy=0;
		Alpha=hight/width;
		Alpha1=fus/310.27;
		Alpha2=fuf/310.27;
		Beta1=ts/0.4572;
		Beta2=tf/0.4572;
		Beta3=screw_Spacing/152.4;
		Lambda=1.736*(Alpha1*Alpha2)/(Beta1*Beta2*pow(Beta3,2)*Alpha);
		Wmax=width/(hight/(sqrt(pow(hight,2)+(width*width))));
		if (tf/ts<=1.0)
		{
			Pns1=4.2*sqrt(pow(tf,3)*ds)*fuf;
			Pns2=2.7*ts*ds*fus;
			Pns3=2.7*tf*ds*fuf;
			MinPns=Pns1;
			MinPns=(Pns2<MinPns)? Pns2:MinPns;
			MinPns=(Pns3<MinPns)? Pns3:MinPns;
		}
		else if (tf/ts>=2.5)
		{
			Pns1=2.7*ts*ds*fus;
			Pns2=2.7*tf*ds*fuf;
			MinPns=(Pns1<Pns2)? Pns1:Pns2;
		}
		else if ((tf/ts)>1.0 && (tf/ts)<2.5)
		{
			Pns1=4.2*sqrt(pow(tf,3)*ds)*fuf;
			Pns2=2.7*ts*ds*fus;
			Pns3=2.7*tf*ds*fuf;
			MinPns1=Pns1;
		  	MinPns1=(Pns2<MinPns)? Pns2:MinPns1;
 	     	MinPns1=(Pns3<MinPns)? Pns3:MinPns1;
			MinPns2=(Pns1<Pns2)? Pns2:Pns3;
			MinPns=MinPns1+(MinPns2-MinPns1)*((tf/ts)-1)/1.5;
         }
			double dis=3*ds;
			Pnsed=0.5*dis*ts*fus;
			if (Lambda<=0.0819)
			We=Wmax;
		else
		{
			rho=(1-0.05*pow((Lambda-0.08),0.12))/pow(Lambda,0.12);
			We=rho*Wmax;
		}
			Pn=(MinPns<Pnsed)? MinPns:Pnsed;
			V1=(((We/(2*screw_Spacing))*Pn)+((We*width)/(2*screw_Spacing*hight)*Pn)+Vs*(width/(sqrt(pow(hight,2)+(width*width)))));
			V2=(We*ts*fys)*(width/sqrt(pow(hight,2)+(width*width)));
			V=(V1<V2)? V1:V2;
			double r,fo;
	        r=1/(1+A/(hight*(width-L)));
	        fo=r/(3-2*r);
			stress3p=fo*V*np;
			Beta=500*(ts/0.457);
			Gs=E/(2*(1+mu));
    		Omega4=sqrt(227.53/fyf);
			rho=0.075*(ts/0.457);
			Delta1=(2*(stress3p/(width*np))*pow(hight,3)/(3*E*Af*width));
            Omega1=screw_Spacing/152.4;
			Omega2=0.838/tf; 
			Delta2=Omega1*Omega2*((stress3p/(width*np))*hight)/(rho*Gs*ts);
			Omega3=sqrt(hight/(2*width));
			Delta3=pow(Omega1,(5/4))*Omega2*Omega3*Omega4*pow((stress3p/(width*np)/(0.0029*Beta)),2);
			Delta4=2.5*hight/width;
            strain3p=(Delta1+Delta2+Delta3+Delta4)/(1000); 
        	stress4p=0.8*stress3p;
			strain4p=1.4*strain3p;
			stress1p=0.4*stress3p;
			strain1p=strain3p/9.25;
			ke=stress1p/strain1p; 
			stress2p=0.85*stress3p;
			Dy=(stress2p/ke);
			strain2p=(stress2p*(strain3p+Dy-2*strain4p-strain1p)+stress3p*strain4p+stress4p*(strain4p-strain3p))/(0.6*stress3p);
			stress1p=stress1p; stress2p=stress2p; stress3p=stress3p; stress4p=stress4p; 
			strain1p=strain1p; strain2p=strain2p; strain3p=strain3p; strain4p=strain4p; 
			strain1n = -strain1p; stress1n = -stress1p; strain2n = -strain2p; stress2n = -stress2p;
			strain3n = -strain3p; stress3n = -stress3p; strain4n = -strain4p; stress4n = -stress4p;
            envlpPosStress.Zero(); envlpPosStrain.Zero(); envlpNegStress.Zero(); envlpNegStrain.Zero(); 
	        energyCapacity = 0.0; kunload = 0.0; elasticStrainEnergy = 0.0; 
}

 CFSSSWP::CFSSSWP():
   UniaxialMaterial(0, MAT_TAG_CFSSSWP),
   stress1p(0.0), strain1p(0.0), stress2p(0.0), strain2p(0.0),
   stress3p(0.0), strain3p(0.0), stress4p(0.0), strain4p(0.0),
   stress1n(0.0), strain1n(0.0), stress2n(0.0), strain2n(0.0),
   stress3n(0.0), strain3n(0.0), stress4n(0.0), strain4n(0.0),

   gammaDLimit(0.0),
   gammaFLimit(0.0), 
   gammaE(0.0),
   rDispP(0.0), rForceP(0.0), uForceP(0.0), rDispN(0.0), rForceN(0.0), uForceN(0.0)
 {
 
 }
 
 CFSSSWP::~CFSSSWP()
 {
  
 }
 
int CFSSSWP::setTrialStrain(double strain, double CstrainRate)
 {
         Tstate = Cstate;
         Tenergy = Cenergy;
         Tstrain = strain;
         lowTstateStrain = lowCstateStrain;
         hghTstateStrain = hghCstateStrain;
         lowTstateStress = lowCstateStress;
         hghTstateStress = hghCstateStress;
         TminStrainDmnd = CminStrainDmnd;
         TmaxStrainDmnd = CmaxStrainDmnd;
         TgammaF = CgammaF;
		 TgammaFN = CgammaFN;
         TgammaD = CgammaD;
		 TgammaDN = CgammaDN;
 
         dstrain = Tstrain - Cstrain;
         if (dstrain<1e-12 && dstrain>-1e-12){
                 dstrain = 0.0;
         }
 
         // determine new state if there is a change in state
         getstate(Tstrain,dstrain);
 
         switch (Tstate)
         { 
 
         case 0:
                 Ttangent = envlpPosStress(0)/envlpPosStrain(0);
                 Tstress = Ttangent*Tstrain;
                 break;
         case 1:
                 Tstress = posEnvlpStress(strain);
                 Ttangent = posEnvlpTangent(strain);
                 break;
         case 2:
                 Ttangent = negEnvlpTangent(strain);
                 Tstress = negEnvlpStress(strain);
                 break;
         case 3:
                 kunload = (hghTstateStrain<0.0) ? kElasticNeg:kElasticPos;    
                         state3Strain(0) = lowTstateStrain;
                         state3Strain(3) = hghTstateStrain;
                         state3Stress(0) = lowTstateStress;
                         state3Stress(3) = hghTstateStress;
 
                 getState3(state3Strain,state3Stress,kunload); 
				 SetSpline();
                 Ttangent = Envlp3Tangent(state3Strain,state3Stress,strain);   
                 Tstress = Envlp3Stress(state3Strain,state3Stress,strain); 
                 break;
         case 4:
                 kunload = (lowTstateStrain<0.0) ? kElasticNeg:kElasticPos;
                         state4Strain(0) = lowTstateStrain;
                         state4Strain(3) = hghTstateStrain;
                         state4Stress(0) = lowTstateStress;
                         state4Stress(3) = hghTstateStress;

                 getState4(state4Strain,state4Stress,kunload);
				 SetSpline();
                 Ttangent = Envlp4Tangent(state4Strain,state4Stress,strain);
                 Tstress = Envlp4Stress(state4Strain,state4Stress,strain);
                 break;
         }
 
         double denergy = 0.5*(Tstress+Cstress)*dstrain;
         elasticStrainEnergy = (Tstrain>0.0) ? 0.5*Tstress/kElasticPos*Tstress:0.5*Tstress/kElasticNeg*Tstress;
         Tenergy = Cenergy + denergy;
         updateDmg(Tstrain,dstrain);
		 return 0;
 }

static int 
getIndex(Vector v,double value)
{
	for(int i = 0; i < v.Size(); i++)
	{
		if(v[i] > value) return i;
	}
	return -1;
}
static int 
getIndexNeg(Vector v,double value)
{
	for(int i = 0; i < v.Size(); i++)
	{
		if(v[i] < value) return i;
	}
	return -1;
}

 void CFSSSWP::SetSpline(void)
 {
			constexpr int Size = 5;
			double X[Size]; double Y[Size];
			
			int fifth = getIndexNeg(envlpNegStrain,state3Strain(0));
			if(fifth == -1)
			{
				printf("erreur fifth");
				exit(5);
			}
			
			X[0] = state3Strain(0) - 20;
			X[1] = state3Strain(0);
			X[2] = state3Strain(1);
			X[3] = state3Strain(2);
			X[4] = state3Strain(3);

			
			Y[0] = state3Stress(0) - 1;
			Y[1] = state3Stress(0);
			Y[2] = state3Stress(1);
			Y[3] = state3Stress(2);
			Y[4] = state3Stress(3);

			if(X[3] - X[0] < 0)
			{
				printf("erreur1\n");	
			}
			
			float a0,an,b0,bn;

			a0 = GetTangentFromCurve(state3Strain(0));
			an = GetTangentFromCurve(state3Strain(3));
			b0 = state3Stress(0) - a0*state3Strain(0);
			bn = state3Stress(3) - an*state3Strain(3);
			
			Spline3.Fit(X,Size,Y,Size);
		
			fifth = getIndex(envlpPosStrain,state4Strain(3));
			if(fifth == -1)
			{
				printf("erreur fifth1");
				exit(5);
			}
			
			X[0] = state4Strain(0);
			X[1] = state4Strain(1);
			X[2] = state4Strain(2);
			X[3] = state4Strain(3);
			X[4] = state4Strain(3) + 20;

			Y[0] = state4Stress(0);
			Y[1] = state4Stress(1);
			Y[2] = state4Stress(2);
			Y[3] = state4Stress(3);
			Y[4] = state4Stress(3) + 1;
			
			if(X[3] - X[0] < 0)
			{
				printf("erreur2\n");
				//while(1);
			}
			
			a0 = GetTangentFromCurve(state4Strain(0));
			an = GetTangentFromCurve(state4Strain(3));
			b0 = state4Stress(0) - a0 * state4Strain(0);
			bn = state4Stress(3) - an * state4Strain(3);
			
			Spline4.Fit(X,Size,Y,Size);
 }
 
 double CFSSSWP::getStrain(void)
 {
         return Tstrain;
 }
 
 double CFSSSWP::getStress(void)
 {
         return Tstress;
 }
 
 double CFSSSWP::getTangent(void)
 {
         return Ttangent;
 }
 
 double CFSSSWP::getInitialTangent(void)
 {
         return envlpPosStress(0)/envlpPosStrain(0);
 }
 
 int CFSSSWP::commitState(void)                                           
 {
         Cstate = Tstate;
 
         if (dstrain>1e-12||dstrain<-(1e-12)) {
                 CstrainRate = dstrain;}
         else {
                 CstrainRate = TstrainRate;}
		 lowCstateStrain = lowTstateStrain;
         lowCstateStress = lowTstateStress;
         hghCstateStrain = hghTstateStrain;
         hghCstateStress = hghTstateStress;
         CminStrainDmnd = TminStrainDmnd;
         CmaxStrainDmnd = TmaxStrainDmnd;
         Cenergy = Tenergy;
         Cstress = Tstress;
         Cstrain = Tstrain;
         CgammaD = TgammaD;
		 CgammaDN = TgammaDN;
         CgammaF = TgammaF;
		 CgammaFN = TgammaFN;
		 CnCycle = TnCycle;
         
         // define adjusted strength and stiffness parameters

         uMaxDamgd = TmaxStrainDmnd*(1 + CgammaD);   
         uMinDamgd = TminStrainDmnd*(1 + CgammaDN);
 
         envlpPosDamgdStress = envlpPosStress*(1-gammaFUsed);
		 envlpNegDamgdStress = envlpNegStress*(1-gammaFUsed);

         return 0;
 }
 
 int CFSSSWP::revertToLastCommit(void)
 {
         
         Tstate = Cstate;
 
         TstrainRate = CstrainRate;
 
         lowTstateStrain = lowCstateStrain;
         lowTstateStress = lowCstateStress;
         hghTstateStrain = hghCstateStrain;
         hghTstateStress = hghCstateStress;
         TminStrainDmnd = CminStrainDmnd;
         TmaxStrainDmnd = CmaxStrainDmnd;
         Tenergy = Cenergy;
 
         Tstrain = Cstrain; Tstress = Cstress;
 
         TgammaD = CgammaD;
		 TgammaDN = CgammaDN;
         TgammaF = CgammaF;
		 TgammaFN = CgammaFN;
		 TnCycle = CnCycle;
 
         return 0;
 }
 
 int CFSSSWP::revertToStart(void)
 {
         Cstate = 0;
         Cstrain = 0.0;
         Cstress = 0.0;
         CstrainRate = 0.0;
         lowCstateStrain = envlpNegStrain(0);
         lowCstateStress = envlpNegStress(0);
         hghCstateStrain = envlpPosStrain(0);
         hghCstateStress = envlpPosStress(0);
         CminStrainDmnd = envlpNegStrain(1);
         CmaxStrainDmnd = envlpPosStrain(1);
         Cenergy = 0.0;
         CgammaD = 0.0;
		 CgammaDN = 0.0;
         CgammaF = 0.0;
		 CgammaFN = 0.0;
		 TnCycle = 0.0;
		 CnCycle = 0.0;
         Ttangent = envlpPosStress(0)/envlpPosStrain(0);
         dstrain = 0.0;       
         gammaFUsed = 0.0;
         uMaxDamgd = CmaxStrainDmnd;
         uMinDamgd = CminStrainDmnd;
 
         return 0;
 }
 
 UniaxialMaterial* CFSSSWP::getCopy(void)
 {
         CFSSSWP *theCopy = new CFSSSWP (this->getTag(),
                  hight,  width,  fuf,  fyf, tf,  Af,  fus,  
				  fys,  ts, np,  ds,  Vs,  screw_Spacing, A, L);
         
         theCopy->rDispN = rDispN;
         theCopy->rDispP = rDispP;
         theCopy->rForceN = rForceN;
         theCopy->rForceP = rForceP;
         theCopy->uForceN = uForceN;
         theCopy->uForceP = uForceP;
 
         // Trial state variables
         theCopy->Tstress = Tstress;
         theCopy->Tstrain = Tstrain;
         theCopy->Ttangent = Ttangent;
 
         // Coverged material history parameters
         theCopy->Cstate = Cstate;
         theCopy->Cstrain = Cstrain;
         theCopy->Cstress = Cstress;
         theCopy->CstrainRate = CstrainRate;
         theCopy->lowCstateStrain = lowCstateStrain;
         theCopy->lowCstateStress = lowCstateStress;
         theCopy->hghCstateStrain = hghCstateStrain;
         theCopy->hghCstateStress = hghCstateStress;
         theCopy->CminStrainDmnd = CminStrainDmnd;
         theCopy->CmaxStrainDmnd = CmaxStrainDmnd;
         theCopy->Cenergy = Cenergy;
         theCopy->CgammaD = CgammaD;
		 theCopy->CgammaDN = CgammaDN;
         theCopy->CgammaF = CgammaF;
		 theCopy->CgammaFN = CgammaFN;
		 theCopy->CnCycle = CnCycle;
         theCopy->gammaFUsed = gammaFUsed;
 
         // trial material history parameters
         theCopy->Tstate = Tstate;
         theCopy->dstrain = dstrain;
         theCopy->lowTstateStrain = lowTstateStrain;
         theCopy->lowTstateStress = lowTstateStress;
         theCopy->hghTstateStrain = hghTstateStrain;
         theCopy->hghTstateStress = hghTstateStress;
         theCopy->TminStrainDmnd = TminStrainDmnd;
         theCopy->TmaxStrainDmnd = TmaxStrainDmnd;
         theCopy->Tenergy = Tenergy;
         theCopy->TgammaD = TgammaD;
		 theCopy->TgammaDN = TgammaDN;
         theCopy->TgammaF = TgammaF;
		 theCopy->TgammaFN = TgammaFN;
		 theCopy->TnCycle = TnCycle;
 
         // Strength and stiffness parameters
         theCopy->kElasticPos = kElasticPos;
         theCopy->kElasticNeg = kElasticNeg;
         theCopy->uMaxDamgd = uMaxDamgd;
         theCopy->uMinDamgd = uMinDamgd;
 
         for (int i = 0; i<7; i++)
         {
                 theCopy->envlpPosStrain(i) = envlpPosStrain(i);
                 theCopy->envlpPosStress(i) = envlpPosStress(i);
                 theCopy->envlpNegStrain(i) = envlpNegStrain(i);
                 theCopy->envlpNegStress(i) = envlpNegStress(i);
                 theCopy->envlpNegDamgdStress(i) = envlpNegDamgdStress(i);
                 theCopy->envlpPosDamgdStress(i) = envlpPosDamgdStress(i);
         }
 
         for (int j = 0; j<4; j++)
         {
                 theCopy->state3Strain(j) = state3Strain(j);
                 theCopy->state3Stress(j) = state3Stress(j);
                 theCopy->state4Strain(j) = state4Strain(j);
                 theCopy->state4Stress(j) = state4Stress(j);
         }
 
         theCopy->energyCapacity = energyCapacity;
         theCopy->kunload = kunload;
         theCopy->elasticStrainEnergy = elasticStrainEnergy;
 
         return theCopy;
 }
 
 int CFSSSWP::sendSelf(int commitTag, Channel &theChannel)
 {
         return -1;
 }
 
 int CFSSSWP::recvSelf(int commitTag, Channel &theChannel,
                                                            FEM_ObjectBroker & theBroker)
 {
         return -1;
 }
 
 void CFSSSWP::Print(OPS_Stream &s, int flag)
 {
         s << "CFSSSWP, tag: " << this-> getTag() << endln;
         s << "Displacement: " << Tstrain << endln;
         s << "Strength: " << Tstress << endln;
         s << "state: " << Tstate << endln;
 }
 
 void CFSSSWP::SetEnvelope(void)
 { 
         double kPos = stress1p/strain1p;
         double kNeg = stress1n/strain1n;
         double k = (kPos>kNeg) ? kPos:kNeg;
         double u = (strain1p>-strain1n) ? 1e-20*strain1p:-1e-20*strain1n;
     
         envlpPosStrain(0) = u;
         envlpPosStress(0) = u*k;
         envlpNegStrain(0) = -u;
         envlpNegStress(0) = -u*k;
	
         envlpPosStrain(1) = strain1p;
         envlpPosStrain(2) = strain2p;
         envlpPosStrain(3) = strain3p;
         envlpPosStrain(4) = strain4p;
 
         envlpNegStrain(1) = strain1n;
         envlpNegStrain(2) = strain2n;
         envlpNegStrain(3) = strain3n;
         envlpNegStrain(4) = strain4n;
 
         envlpPosStress(1) = stress1p;
         envlpPosStress(2) = stress2p;
         envlpPosStress(3) = stress3p;
         envlpPosStress(4) = stress4p;
 
         envlpNegStress(1) = stress1n;
         envlpNegStress(2) = stress2n;
         envlpNegStress(3) = stress3n;
         envlpNegStress(4) = stress4n;
 
         double k1 = (stress4p - stress3p)/(strain4p - strain3p);
         double k2 = (stress4n - stress3n)/(strain4n - strain3n);
  
        envlpPosStress(5) =0.05*stress3p;
		envlpPosStrain(5) = strain4p + 3.75*(strain4p-strain3p);
	    envlpNegStress(5) = 0.05*stress3n;
		envlpNegStrain(5) = strain4n + 3.75*(strain4n-strain3n);

		envlpPosStrain(6) = 1e+6*envlpPosStress(5);
	    envlpPosStress(6) = (k1>0.0)? envlpPosStress(5)+k1*(envlpPosStrain(6) - envlpPosStrain(5)):envlpPosStress(5)*1.1;
	    envlpNegStrain(6) = 1e+6*strain4n;
	    envlpNegStress(6) = (k2>0.0)? envlpNegStress(5)+k1*(envlpNegStrain(6) - envlpNegStrain(5)):envlpNegStress(5)*1.1;
       
         // define critical material properties
         kElasticPos = envlpPosStress(1)/envlpPosStrain(1);      
         kElasticNeg = envlpNegStress(1)/envlpNegStrain(1);
 
         double energypos = 0.5*envlpPosStrain(0)*envlpPosStress(0);
 
         for (int jt = 0; jt<4; jt++){
                 energypos += 0.5*(envlpPosStress(jt) + envlpPosStress(jt+1))*(envlpPosStrain(jt+1)-envlpPosStrain(jt));
         }
 
         double energyneg = 0.5*envlpNegStrain(0)*envlpNegStress(0);
 
         for (int jy = 0; jy<4; jy++){
                 energyneg += 0.5*(envlpNegStress(jy) + envlpNegStress(jy+1))*(envlpNegStrain(jy+1)-envlpNegStrain(jy));
         }
 
         double max_energy = (energypos>energyneg) ? energypos:energyneg;
 
         energyCapacity = gammaE*max_energy;
 
		 // BSpline Adds

		 constexpr int Size = 9;
		 double X[Size]; double Y[Size];

		 for(int i = 0;i < 2;i++)
		 {
			 X[i] = envlpPosStrain(0);
			 Y[i] = envlpPosStress(0);
			 X[Size - i - 1] = envlpPosStrain(4);
			 Y[Size - i - 1] = envlpPosStress(4);
		 }

		 for(int i = 0;i < Size - 4;i++)
		 {
			 X[i + 2] = envlpPosStrain(i);
			 Y[i + 2] = envlpPosStress(i);
		 }
		 //double *XFit = new double[(Size-3)*Precision+2],*YFit = new double[(Size-3)*Precision+2];
		 double XFit[(Size-3)*Precision+2]; double YFit[(Size-3)*Precision+2]; 
		 double a[4]; double b[4];
			 
		
		double p1X,p1Y,p2X,p2Y,p3X,p3Y,p4X,p4Y;
		
		for(int i = 0;i < Size-3;i++)
	{
          p1X = X[i];
		  p1Y = Y[i];
          p2X = X[i + 1];
		  p2Y = Y[i + 1];
          p3X = X[i + 2];
		  p3Y = Y[i + 2];
          p4X = X[i + 3];
		  p4Y = Y[i + 3];
          a[0] = (-p1X + 3 * p2X - 3 * p3X + p4X) / 6.0f;
          a[1] = (3 * p1X - 6 * p2X + 3 * p3X) / 6.0f;
          a[2] = (-3 * p1X + 3 * p3X) / 6.0f;
          a[3] = (p1X + 4 * p2X + p3X) / 6.0f;
          b[0] = (-p1Y + 3 * p2Y - 3 * p3Y + p4Y) / 6.0f;
          b[1] = (3 * p1Y - 6 * p2Y + 3 * p3Y) / 6.0f;
          b[2] = (-3 * p1Y + 3 * p3Y) / 6.0f;
          b[3] = (p1Y + 4 * p2Y + p3Y) / 6.0f;
          for (int j = 0; j < Precision; j++)
          {
               float t = (float)(j) / (float)(Precision);
               XFit[i*Precision+j] = ((a[2] + t * (a[1] + t * a[0])) * t + a[3]);
               YFit[i*Precision+j] = ((b[2] + t * (b[1] + t * b[0])) * t + b[3]);
          }
	}
	
	double XAdvance = XFit[Precision*(Size-3)-1] - XFit[Precision*(Size-3)-2];
	double YAdvance = YFit[Precision*(Size-3)-1] - YFit[Precision*(Size-3)-2];
	double Tangente = (YFit[Precision*(Size-3)-1] - YFit[Precision*(Size-3)-2])/(XFit[Precision*(Size-3)-1] - XFit[Precision*(Size-3)-2]);
	
	double Epsilon = 0.1f;
		
	YFit[Precision*(Size-3)] = Epsilon;
	XFit[Precision*(Size-3)] = XFit[Precision*(Size-3) - 1] + (Epsilon - YFit[Precision*(Size-3) - 1]) / Tangente;
		
	YFit[Precision*(Size-3)+1] = Epsilon;
	XFit[Precision*(Size-3)+1] = 10000;
	
	BSplineXs = XFit;
	BSplineYs = YFit;

	BSplineXLength = Precision*(Size-3) + 2;
	BSplineYLength = Precision*(Size-3) + 2;
         
 }

 double CFSSSWP::GetTangentFromCurve(double Strain)
 {
		int i = 0;
		int Neg = 0;
		while(i < BSplineXLength && BSplineXs[i] < Strain)
		{
			i++;
		}
		if(i == BSplineXLength && BSplineXs[i-1] < Strain)
		{
			if(Neg == 0) return 1;
			return -1;
			exit(0);
		}
		if(BSplineXs[i] == Strain) 
		{
			return (BSplineYs[i+1] - BSplineYs[i-1]) / (BSplineXs[i+1] - BSplineXs[i-1]);
		}
		else if (i < BSplineXLength - 2 && BSplineXs[i+1] == Strain)
		{
			return (BSplineYs[i+2] - BSplineYs[i]) / (BSplineXs[i+2] - BSplineXs[i]);
		}
		else
		{
			return (BSplineYs[i] - BSplineYs[i-1]) / (BSplineXs[i] - BSplineXs[i-1]);
		}
 }
 
 double CFSSSWP::GetStressFromCurve(double Strain)
 {
		int i = 0;
		int Neg = 0;
		if (Strain < 0)
		{
			Neg = 1;
			Strain = -Strain;
		}
		while(i < BSplineXLength && BSplineXs[i] < Strain)
		{
			i++;
		}
		if(i == BSplineXLength && BSplineXs[i-1] < Strain)
		{
			if(Neg == 0) return -1;
			return 1;
		}
		if(BSplineXs[i] == Strain)
		{
			if(Neg == 1)
				return BSplineYs[i];
				return BSplineYs[i];
		}
		else if (i < BSplineXLength - 1 && BSplineXs[i+1] == Strain)
		{
			return BSplineYs[i+1];
		}
		else
		{
			double Stress = BSplineYs[i-1] + (BSplineYs[i] - BSplineYs[i-1]) / (BSplineXs[i] - BSplineXs[i-1]) * (Strain - BSplineXs[i-1]);
			if(Neg == 1)
			return -Stress;
			return  Stress;
			}
 }
   
 void CFSSSWP::getstate(double u,double du)
 {
         int cid = 0;
         int cis = 0;
         int newState = 0;
         if (du*CstrainRate<=0.0){   
                 cid = 1;
         }
         if (u<lowTstateStrain || u>hghTstateStrain || cid) {                
                 if (Tstate == 0) {                                              
                         if (u>hghTstateStrain) {
                                 cis = 1;
                                 newState = 1;
                                 lowTstateStrain = envlpPosStrain(0);
                                 lowTstateStress = envlpPosStress(0);
                                 hghTstateStrain = envlpPosStrain(5);
                                 hghTstateStress = envlpPosStress(5);
                         }
                         else if (u<lowTstateStrain){
                                 cis = 1;
                                 newState = 2;
                                 lowTstateStrain = envlpNegStrain(5);
                                 lowTstateStress = envlpNegStress(5);
                                 hghTstateStrain = envlpNegStrain(0);
                                 hghTstateStress = envlpNegStress(0);
                         }
                 }
                 else if (Tstate==1 && du<0.0) {
                         cis = 1;
                         if (Cstrain>TmaxStrainDmnd) {
                                 TmaxStrainDmnd = u - du;
                         }
                         if (TmaxStrainDmnd<uMaxDamgd) {
                                 TmaxStrainDmnd = uMaxDamgd;
                         }
                         if (u<uMinDamgd) {
                                 newState = 2;
                                 gammaFUsed = CgammaFN;     
                                 for (int i=0; i<=6; i++) {
                                         envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
                                 }
                                 lowTstateStrain = envlpNegStrain(6);
                                 lowTstateStress = envlpNegStress(6);
                                 hghTstateStrain = envlpNegStrain(0);
                                 hghTstateStress = envlpNegStress(0);
                         }
                         else {
                                 newState = 3;
                                 lowTstateStrain = uMinDamgd;
                                 gammaFUsed = CgammaFN;        
                                 for (int i=0; i<=6; i++) {
                                         envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
                                 }
                                 lowTstateStress = negEnvlpStress(uMinDamgd); 
                                 hghTstateStrain = Cstrain;
                                 hghTstateStress = Cstress;
                         }

                 }
                 else if (Tstate ==2 && du>0.0){
                         cis = 1;
                         if (Cstrain<TminStrainDmnd) {
                                 TminStrainDmnd = Cstrain;
                         }
                         if (TminStrainDmnd>uMinDamgd) {
                                 TminStrainDmnd = uMinDamgd;
                         }
                         if (u>uMaxDamgd) {
                                 newState = 1;
                                 gammaFUsed = CgammaF;      
                                 for (int i=0; i<=6; i++) {
                                         envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
                                 }
                                 lowTstateStrain = envlpPosStrain(0);
                                 lowTstateStress = envlpPosStress(0);
                                 hghTstateStrain = envlpPosStrain(5);
                                 hghTstateStress = envlpPosStress(5);
                         }
                         else {
                                 newState = 4;
                                 lowTstateStrain = Cstrain;
                                 lowTstateStress = Cstress;
                                 hghTstateStrain = uMaxDamgd;
                                 gammaFUsed = CgammaF;         
                                 for (int i=0; i<=6; i++) {
                                         envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
                                 }
                                 hghTstateStress = posEnvlpStress(uMaxDamgd);
                         }
                        
                 }
                         else if (Tstate ==3) {
                                 if (u<lowTstateStrain){
                                         cis = 1;
                                         newState = 2;
                                         lowTstateStrain = envlpNegStrain(5);
                                         hghTstateStrain = envlpNegStrain(0);
                                         lowTstateStress = envlpNegDamgdStress(5);
                                         hghTstateStress = envlpNegDamgdStress(0);
                                 }
                                 else if (u>uMaxDamgd && du>0.0) {
                                         cis = 1;
                                         newState = 1;
                                         lowTstateStrain = envlpPosStrain(0);
                                         lowTstateStress = envlpPosStress(0);
                                         hghTstateStrain = envlpPosStrain(5);
                                         hghTstateStress = envlpPosStress(5);
                                 }
                                 else if (du>0.0) {
                                       cis = 1;
                                         newState = 4;
                                         lowTstateStrain = Cstrain;
                                         lowTstateStress = Cstress;
                                         hghTstateStrain = uMaxDamgd;
                                         gammaFUsed = CgammaF;
                                         for (int i=0; i<=6; i++) {
                                         envlpPosDamgdStress(i) = envlpPosStress(i)*(1.0-gammaFUsed);
                                     }
                                         hghTstateStress = posEnvlpStress(uMaxDamgd);
                                        
                                 }
                         }
                         else if (Tstate == 4){
                                 if (u>hghTstateStrain){
                                         cis = 1; 
                                         newState = 1;
                                         lowTstateStrain = envlpPosStrain(0);
                                         lowTstateStress = envlpPosDamgdStress(0);
                                         hghTstateStrain = envlpPosStrain(5);
                                         hghTstateStress = envlpPosDamgdStress(5);
                                 }
                                 else if (u<uMinDamgd && du <0.0) {
                                         cis = 1;
                                         newState = 2;
                                         lowTstateStrain = envlpNegStrain(5);
                                         lowTstateStress = envlpNegDamgdStress(5);
                                         hghTstateStrain = envlpNegStrain(0);
                                         hghTstateStress = envlpNegDamgdStress(0);
                                 }
                                 else if (du<0.0) { 
                                         cis = 1;
                                         newState = 3;
                                         lowTstateStrain = uMinDamgd;
                                         gammaFUsed = CgammaFN;         
                                         for (int i=0; i<=6; i++) {
                                         envlpNegDamgdStress(i) = envlpNegStress(i)*(1.0-gammaFUsed);
                                      }
                                         lowTstateStress = negEnvlpStress(uMinDamgd);
                                         hghTstateStrain = Cstrain;
                                         hghTstateStress = Cstress;

                                 }
                         }
                         }
                         if (cis) {
                                 Tstate = newState;
                                 
                         }
                 }
 
double CFSSSWP::posEnvlpStress(double u)
                 {
                         double k = 0.0;
                         int i = 0;
                         double f = 0.0;
						 f = GetStressFromCurve(u);
						 f = f*(1-gammaFUsed);
						 return f;
                         while (k==0.0 && i<=5){
                                  
                                  if (u<=envlpPosStrain(i+1)){
                                          k = (envlpPosDamgdStress(i+1)-envlpPosDamgdStress(i))/(envlpPosStrain(i+1)-envlpPosStrain(i));
                                          f = envlpPosDamgdStress(i) + (u-envlpPosStrain(i))*k;
                                  }
                                  i++;
                         }
             
 
                         if (k==0.0){
                                 k = (envlpPosDamgdStress(6) - envlpPosDamgdStress(5))/(envlpPosStrain(6) - envlpPosStrain(5));
                                 f = envlpPosDamgdStress(6) + k*(u-envlpPosStrain(6));
                         }
                         return f;
 
                 }
 
 double CFSSSWP::posEnvlpTangent(double u)
                 {
                         double k = 0.0;
                         int i = 0;
                         while (k==0.0 && i<=5){        
                                  
                                  if (u<=envlpPosStrain(i+1)){
                                          k = (envlpPosDamgdStress(i+1)-envlpPosDamgdStress(i))/(envlpPosStrain(i+1)-envlpPosStrain(i));
                         }
                                  i++;
                         }
 
                         if (k==0.0){
                                 k = (envlpPosDamgdStress(6) - envlpPosDamgdStress(5))/(envlpPosStrain(6) - envlpPosStrain(5));
						 }
						 k = GetTangentFromCurve(u);
 
                         return k;
 
                 }
 
 
 double CFSSSWP::negEnvlpStress(double u)
                 {
                         double k = 0.0;
                         int i = 0;
                         double f = 0.0;
						 f = GetStressFromCurve(u);
						 f = f*(1-gammaFUsed);
						 return f;
                         while (k==0.0 && i<=5){                                  
                                  if (u>=envlpNegStrain(i+1)){
                                          k = (envlpNegDamgdStress(i)-envlpNegDamgdStress(i+1))/(envlpNegStrain(i)-envlpNegStrain(i+1));
                                          f = envlpNegDamgdStress(i+1)+(u-envlpNegStrain(i+1))*k;
                                  }
                                  i++;
                         }
 
                         if (k==0.0){
                                 k = (envlpNegDamgdStress(5) - envlpNegDamgdStress(6))/(envlpNegStrain(5)-envlpNegStrain(6));
                                 f = envlpNegDamgdStress(6) + k*(u-envlpNegStrain(6));
                         }
                         return f;
 
                 }
 
 double CFSSSWP::negEnvlpTangent(double u)
                 {
                         double k = 0.0;
                         int i = 0;
                         while (k==0.0 && i<=5){              
                                  
                                  if (u>=envlpNegStrain(i+1)){
                                          k = (envlpNegDamgdStress(i)-envlpNegDamgdStress(i+1))/(envlpNegStrain(i)-envlpNegStrain(i+1));
                                         }
                                  i++;
                         }
 
                         if (k==0.0){
                                 k = (envlpNegDamgdStress(5) - envlpNegDamgdStress(6))/(envlpNegStrain(5)-envlpNegStrain(6));
                                 }
						 k = GetTangentFromCurve(u);
                         return k;
 
                 }
 
 
 void CFSSSWP::getState3(Vector& state3Strain, Vector& state3Stress, double kunload)
                 {
 
                         double kmax = (kunload>kElasticNeg) ? kunload:kElasticNeg;
 
                         if (state3Strain(0)*state3Strain(3) <0.0){
                                 // trilinear unload reload path expected, first define point for reloading
                                 state3Strain(1) = lowTstateStrain*rDispN;
                                if (rForceN-uForceN > 1e-8) {
                                         state3Stress(1) = lowTstateStress*rForceN;
                                 }
                                 else {
                                         if (TminStrainDmnd < envlpNegStrain(3)) {
                                                 double st1 = lowTstateStress*uForceN*(1.0+1e-6);
                                                 double st2 = envlpNegDamgdStress(4)*(1.0+1e-6);
                                                 state3Stress(1) = (st1<st2) ? st1:st2;
                                         }
                                         else {
                                                 double st1 = envlpNegDamgdStress(3)*uForceN*(1.0+1e-6);
                                                 double st2 = envlpNegDamgdStress(4)*(1.0+1e-6);
                                                 state3Stress(1) = (st1<st2) ? st1:st2;
                                         }
                                 }
                                 // if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
                                if ((state3Stress(1)-state3Stress(0))/(state3Strain(1)-state3Strain(0)) > kElasticNeg) {
                                         state3Strain(1) = lowTstateStrain + (state3Stress(1)-state3Stress(0))/kElasticNeg;
                                 }
                                 // check that reloading point is not behind point 4 
                                if (state3Strain(1)>state3Strain(3)) {
                                         // path taken to be a straight line between points 1 and 4
                                         double du = state3Strain(3) - state3Strain(0);
                                         double df = state3Stress(3) - state3Stress(0);
                                         state3Strain(1) = state3Strain(0) + 0.33*du;
                                         state3Strain(2) = state3Strain(0) + 0.67*du;
                                         state3Stress(1) = state3Stress(0) + 0.33*df;
                                         state3Stress(2) = state3Stress(0) + 0.67*df;
                                 }
                                 else {
                                         if (TminStrainDmnd < envlpNegStrain(3)) {
                                                 state3Stress(2) = uForceN*envlpNegDamgdStress(4);
                                         }
                                         else {
                                                 state3Stress(2) = uForceN*envlpNegDamgdStress(3);
                                         }
                                         state3Strain(2) = hghTstateStrain - (hghTstateStress-state3Stress(2))/kunload;
 
                                         if (state3Strain(2) > state3Strain(3)) {
                                                 // point 3 should be along a line between 2 and 4
                                                 double du = state3Strain(3) - state3Strain(1);
                                                 double df = state3Stress(3) - state3Stress(1);
                                                 state3Strain(2) = state3Strain(1) + 0.5*du;
                                                 state3Stress(2) = state3Stress(1) + 0.5*df;
                                         }
                                         else if ((state3Stress(2) - state3Stress(1))/(state3Strain(2) - state3Strain(1)) > kmax) {
                                                 // linear unload-reload path expected
                                                 double du = state3Strain(3) - state3Strain(0);
                                                 double df = state3Stress(3) - state3Stress(0);
                                                 state3Strain(1) = state3Strain(0) + 0.33*du;
                                                 state3Strain(2) = state3Strain(0) + 0.67*du;
                                                 state3Stress(1) = state3Stress(0) + 0.33*df;
                                                 state3Stress(2) = state3Stress(0) + 0.67*df;
                                         }
                                         else if ((state3Strain(2) < state3Strain(1))||((state3Stress(2)-state3Stress(1))/(state3Strain(2)-state3Strain(1))<0)) {
                                                 if (state3Strain(2)<0.0) {
                                                         // pt 3 should be along a line between 2 and 4
                                                         double du = state3Strain(3)-state3Strain(1);
                                                         double df = state3Stress(3)-state3Stress(1);
                                                         state3Strain(2) = state3Strain(1) + 0.5*du;
                                                         state3Stress(2) = state3Stress(1) + 0.5*df;
                                                 }
                                                 else if (state3Strain(1) > 0.0) {
                                                         // pt 2 should be along a line between 1 and 3
                                                         double du = state3Strain(2)-state3Strain(0);
                                                         double df = state3Stress(2)-state3Stress(0);
                                                         state3Strain(1) = state3Strain(0) + 0.5*du;
                                                         state3Stress(1) = state3Stress(0) + 0.5*df;
                                                 }
                                                 else {
                                                         double avgforce = 0.5*(state3Stress(2) + state3Stress(1));
                                                         double dfr = 0.0;
                                                         if (avgforce < 0.0){
                                                                 dfr = -avgforce/100;
                                                         }
                                                         else {
                                                                 dfr = avgforce/100;
                                                         }
                                                         double slope12 = (state3Stress(1) - state3Stress(0))/(state3Strain(1) - state3Strain(0));
                                                         double slope34 = (state3Stress(3) - state3Stress(2))/(state3Strain(3) - state3Strain(2));
                                                         state3Stress(1) = avgforce - dfr;
                                                         state3Stress(2) = avgforce + dfr;
                                                         state3Strain(1) = state3Strain(0) + (state3Stress(1) - state3Stress(0))/slope12;
                                                         state3Strain(2) = state3Strain(3) - (state3Stress(3) - state3Stress(2))/slope34;
                                                }
                                         }
                                 }
                         }
                                 else {
                                         // linear unload reload path is expected                 
                                         double du = state3Strain(3)-state3Strain(0);
                                         double df = state3Stress(3)-state3Stress(0);
                                         state3Strain(1) = state3Strain(0) + 0.33*du;
                                         state3Strain(2) = state3Strain(0) + 0.67*du;
                                         state3Stress(1) = state3Stress(0) + 0.33*df;
                                         state3Stress(2) = state3Stress(0) + 0.67*df;
                                 }
                         
                                 
                                 double checkSlope = state3Stress(0)/state3Strain(0);
                                 double slope = 0.0;
 
                                 // final check
                                 int i = 0;
                                 while (i<3) {
                                         double du = state3Strain(i+1)-state3Strain(i);
                                         double df = state3Stress(i+1)-state3Stress(i);
                                         if (du<0.0 || df<0.0) {
                                                 double du = state3Strain(3)-state3Strain(0);
                                                 double df = state3Stress(3)-state3Stress(0);
                                                 state3Strain(1) = state3Strain(0) + 0.33*du;
                                                 state3Strain(2) = state3Strain(0) + 0.67*du;
                                                 state3Stress(1) = state3Stress(0) + 0.33*df;
                                                 state3Stress(2) = state3Stress(0) + 0.67*df;
                                                 slope = df/du;
                                                 i = 3;
                                         }
                                         if (slope > 1e-8 && slope < checkSlope) {
                                                 state3Strain(1) = 0.0; state3Stress(1) = 0.0;
                                                 state3Strain(2) = state3Strain(3)/2; state3Stress(2) = state3Stress(3)/2;
                                         } 
                                         i++;
                                 }
 
 
                         }
                                 
 void CFSSSWP::getState4(Vector& state4Strain,Vector& state4Stress, double kunload)
                 {
 
                         double kmax = (kunload>kElasticPos) ? kunload:kElasticPos;
 
                         if (state4Strain(0)*state4Strain(3) <0.0){
                                 // trilinear unload reload path expected
                                 state4Strain(2) = hghTstateStrain*rDispP;
                                 if (uForceP==0.0){
                                         state4Stress(2) = hghTstateStress*rForceP;
                                 }
                                 else if (rForceP-uForceP > 1e-8) {
                                         state4Stress(2) = hghTstateStress*rForceP;
                                 }
                                 else {
                                         if (TmaxStrainDmnd > envlpPosStrain(3)) {
                                                 double st1 = hghTstateStress*uForceP*(1.0+1e-6);
                                                 double st2 = envlpPosDamgdStress(4)*(1.0+1e-6);
                                                 state4Stress(2) = (st1>st2) ? st1:st2;
                                         }
                                         else {
                                                 double st1 = envlpPosDamgdStress(3)*uForceP*(1.0+1e-6);
                                                 double st2 = envlpPosDamgdStress(4)*(1.0+1e-6);
                                                 state4Stress(2) = (st1>st2) ? st1:st2;
                                         }
                                 }
                                 // if reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
                                if ((state4Stress(3)-state4Stress(2))/(state4Strain(3)-state4Strain(2)) > kElasticPos) {
                                        state4Strain(2) = hghTstateStrain - (state4Stress(3)-state4Stress(2))/kElasticPos;
                                 }
                                 // check that reloading point is not behind point 1 
                               if (state4Strain(2)<state4Strain(0)) {
                                         // path taken to be a straight line between points 1 and 4
                                         double du = state4Strain(3) - state4Strain(0);
                                         double df = state4Stress(3) - state4Stress(0);
                                         state4Strain(1) = state4Strain(0) + 0.33*du;
                                         state4Strain(2) = state4Strain(0) + 0.67*du;
                                         state4Stress(1) = state4Stress(0) + 0.33*df;
                                         state4Stress(2) = state4Stress(0) + 0.67*df;
                                }
                                 else {
                                         if (TmaxStrainDmnd > envlpPosStrain(3)) {
                                                 state4Stress(1) = uForceP*envlpPosDamgdStress(4);
                                         }
                                         else {
                                                 state4Stress(1) = uForceP*envlpPosDamgdStress(3);
                                         }
                                         state4Strain(1) = lowTstateStrain + (-lowTstateStress+state4Stress(1))/kunload;
 
                                         if (state4Strain(1) < state4Strain(0)) {
                                                // point 2 should be along a line between 1 and 3
                                                 double du = state4Strain(2) - state4Strain(0);
                                                 double df = state4Stress(2) - state4Stress(0);
                                                 state4Strain(1) = state4Strain(0) + 0.5*du;
                                                 state4Stress(1) = state4Stress(0) + 0.5*df;
                                         }
                                        else if ((state4Stress(2) - state4Stress(1))/(state4Strain(2) - state4Strain(1)) > kmax) {
                                                 // linear unload-reload path expected
                                                 double du = state4Strain(3) - state4Strain(0);
                                                 double df = state4Stress(3) - state4Stress(0);
                                                 state4Strain(1) = state4Strain(0) + 0.33*du;
                                                 state4Strain(2) = state4Strain(0) + 0.67*du;
                                                 state4Stress(1) = state4Stress(0) + 0.33*df;
                                                 state4Stress(2) = state4Stress(0) + 0.67*df;
                                        }
                                       else if ((state4Strain(2) < state4Strain(1))||((state4Stress(2)-state4Stress(1))/(state4Strain(2)-state4Strain(1))<0)) {
                                                 if (state4Strain(1)>0.0) {
                                                         // pt 2 should be along a line between 1 and 3
                                                         double du = state4Strain(2)-state4Strain(0);
                                                         double df = state4Stress(2)-state4Stress(0);
                                                         state4Strain(1) = state4Strain(0) + 0.5*du;
                                                         state4Stress(1) = state4Stress(0) + 0.5*df;
                                                 }
                                                 else if (state4Strain(2) < 0.0) {
                                                         // pt 2 should be along a line between 2 and 4
                                                         double du = state4Strain(3)-state4Strain(1);
                                                         double df = state4Stress(3)-state4Stress(1);
                                                         state4Strain(2) = state4Strain(1) + 0.5*du;
                                                         state4Stress(2) = state4Stress(1) + 0.5*df;
                                                 }
                                                 else {
                                                         double avgforce = 0.5*(state4Stress(2) + state4Stress(1));
                                                         double dfr = 0.0;
                                                         if (avgforce < 0.0){
                                                                 dfr = -avgforce/100;
                                                         }
                                                         else {
                                                                 dfr = avgforce/100;
                                                         }
                                                         double slope12 = (state4Stress(1) - state4Stress(0))/(state4Strain(1) - state4Strain(0));
                                                         double slope34 = (state4Stress(3) - state4Stress(2))/(state4Strain(3) - state4Strain(2));
                                                         state4Stress(1) = avgforce - dfr;
                                                         state4Stress(2) = avgforce + dfr;
                                                         state4Strain(1) = state4Strain(0) + (state4Stress(1) - state4Stress(0))/slope12;
                                                         state4Strain(2) = state4Strain(3) - (state4Stress(3) - state4Stress(2))/slope34;
                                                 }
                                         }
                                 }
                         }
                                else {
                                         // linear unload reload path is expected
                                         double du = state4Strain(3)-state4Strain(0);
                                         double df = state4Stress(3)-state4Stress(0);
                                         state4Strain(1) = state4Strain(0) + 0.33*du;
                                         state4Strain(2) = state4Strain(0) + 0.67*du;
                                         state4Stress(1) = state4Stress(0) + 0.33*df;
                                         state4Stress(2) = state4Stress(0) + 0.67*df;
                                 }
                         
 
                                 
                                 double checkSlope = state4Stress(0)/state4Strain(0);
                                 double slope = 0.0;
 
                                 // final check
                                 int i = 0;
                                 while (i<3) {
                                         double du = state4Strain(i+1)-state4Strain(i);
                                         double df = state4Stress(i+1)-state4Stress(i);
                                         if (du<0.0 || df<0.0) {
                                                 double du = state4Strain(3)-state4Strain(0);
                                                 double df = state4Stress(3)-state4Stress(0);
                                                 state4Strain(1) = state4Strain(0) + 0.33*du;
                                                 state4Strain(2) = state4Strain(0) + 0.67*du;
                                                 state4Stress(1) = state4Stress(0) + 0.33*df;
                                                 state4Stress(2) = state4Stress(0) + 0.67*df;
                                                 slope = df/du;
                                                 i = 3;
                                         }
                                         if (slope > 1e-8 && slope < checkSlope) {
                                                 state4Strain(1) = 0.0; state4Stress(1) = 0.0;
                                                 state4Strain(2) = state4Strain(3)/2; state4Stress(2) = state4Stress(3)/2;
                                         } 
 
                                         i++;
                                 }
                         }
 
 double CFSSSWP::Envlp3Tangent(Vector s3Strain, Vector s3Stress, double u)
                         {
							    double k = 0.0;
								k = Spline3.EvalT(u);
								if(k != 10e8)
								 {
									 return k;
								}
								 int i = 0;
                                 while ((k==0.0||i<=2) && (i<=2)) 
                                 {
                                         if (u>= s3Strain(i)) {
                                                 k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
                                         }
                                         i++;
                                 }
                                 if (k==0.0) {
                                         if (u<s3Strain(0)) {
                                                 i = 0;
                                         }
                                         else {
                                                 i = 2;
                                         }
                                         k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
                                         
                                 }
								printf("Tangente = %f\n",k);
                                 return k;
                         }
 
 double CFSSSWP::Envlp4Tangent(Vector s4Strain, Vector s4Stress, double u)
                         {
                                 double k = 0.0;
                                int i = 0;
								k = Spline4.EvalT(u);
								if(k != 10e8)
								{
								return k;
								}
                                 while ((k==0.0||i<=2) && (i<=2)) 
                                 {
                                         if (u>= s4Strain(i)) {
                                                 k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
                                         }
                                         i++;
                                 }
                                 if (k==0.0) {
                                         if (u<s4Strain(0)) {
                                                 i = 0;
                                         }
                                         else {
                                                 i = 2;
                                         }
                                         k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
                                         
                                 }
								printf("Tangente = %f\n",k);
                                 return k;
                         }
 
 
 double CFSSSWP::Envlp3Stress(Vector s3Strain, Vector s3Stress, double u)
{
  double k = 0.0;
  int i = 0;
  double f = 0.0;
  f = Spline3.Eval(u);
  if(isnan(f))
    {
      printf("erreur3");
      //while(1);
    }
  if(f != 10e8)
    {
      return f;
      return GetStressFromCurve(u);
    }
  while ((k==0.0||i<=2) && (i<=2)) 
    {
      if (u>= s3Strain(i)) {
	k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
	f = s3Stress(i)+(u-s3Strain(i))*k;
      }
      i++;
    }
  if (k==0.0) {
    if (u<s3Strain(0)) {
      i = 0;
    }
    else {
      i = 2;
    }
    k = (s3Stress(i+1)-s3Stress(i))/(s3Strain(i+1)-s3Strain(i));
    f = s3Stress(i)+(u-s3Strain(i))*k;
  }
  printf("Strain = %f	Stress = %f	Min = %f, Max = %f\n",u,f,s3Strain(0),s3Strain(3));
  if(u > s3Strain(3))
    {
      //while(1);
    }
  return f;
}

double CFSSSWP::Envlp4Stress(Vector s4Strain, Vector s4Stress, double u)
{
  double k = 0.0;
  int i = 0;
  double f = 0.0;
  f = Spline4.Eval(u);
  if(isnan(f))
    {
										 printf("erreur4");
										 //while(1);
    }
  if(f != 10e8)
								 {
								   return f;
								   return GetStressFromCurve(u);
								 }
  while ((k==0.0||i<=2) && (i<=2)) 
    {
      if (u>= s4Strain(i)) {
                                                 k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
                                                 f = s4Stress(i)+(u-s4Strain(i))*k;
                                         }
                                         i++;
                                 }
                                 if (k==0.0) {
                                         if (u<s4Strain(0)) {
                                                 i = 0;
                                         }
                                         else {
                                                 i = 2;
                                         }
                                         k = (s4Stress(i+1)-s4Stress(i))/(s4Strain(i+1)-s4Strain(i));
                                         f = s4Stress(i)+(u-s4Strain(i))*k;
                                 }
								 printf("Strain = %f	Stress = %f	Min = %f, Max = %f\n",u,f,s4Strain(0),s4Strain(3));
								 if(u > s4Strain(3))
								 {
									 //while(1);
								 }
                                 return f;
                         }
 
   void CFSSSWP::updateDmg(double strain, double dstrain)
         {
                 double tes = 0.0;
				 double umaxAbs = (TmaxStrainDmnd>-TminStrainDmnd) ? TmaxStrainDmnd:-TminStrainDmnd;
                 double uultAbs = (envlpPosStrain(1)>-envlpNegStrain(1)) ? envlpPosStrain(1):-envlpNegStrain(1);
				 TnCycle = CnCycle;
                 if ((strain<uultAbs && strain>-uultAbs) && Tenergy < elasticStrainEnergy)
                 {
                      
						 TgammaD += TnCycle;
                         TgammaF += TnCycle;

				 }
                         if (Tenergy>elasticStrainEnergy) {
                                 tes = ((Tenergy-elasticStrainEnergy)/energyCapacity);
                                 TgammaD += tes;
                                 TgammaF += tes;
                 }
                       
						 TgammaDN = TgammaD;
                         TgammaD = (TgammaD<gammaDLimit) ? TgammaD:gammaDLimit;
						 TgammaFN = TgammaF;
                         TgammaF = (TgammaF<gammaFLimit) ? TgammaF:gammaFLimit;
         }
