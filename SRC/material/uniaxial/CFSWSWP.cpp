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

// $Revision: 2.0 $
// $Date: 11-10-2015 $

// Written by Smail KECHIDI, Ph.D student at University of Blida 1 (s_kechidi@univ-blida.dz), PhD mobility Student at University of Porto FEUP (smail.kechidi@fe.up.pt)
// Created: 11-10-2015 21:45:20 $
//
// Description: This file contains the class implementation for CFSWSWP
// CFSWSWP is based on Pinching4 uniaxialMaterial

#include <elementAPI.h>
#include "CFSWSWP.h"
#include <OPS_Globals.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <OPS_Stream.h>
#include <stdio.h>
#include <string.h>
#include "CubicSpline.h"
#include "TriMatrix.h"
 
static int numCFSWSWP = 0;

#ifdef _WIN32
#define isnan _isnan
#endif

void *
OPS_CFSWSWP(void)
{
  // print out some KUDO's
  if (numCFSWSWP == 0) {
    opserr << "Cold Formed Steel Wood-Sheathed Shear Wall Panel uniaxialMaterial - Written by Smail KECHIDI Ph.D Student at University of Blida 1 - Please when using this make reference as: Smail Kechidi and Nouredine Bourahla (2016), Deteriorating hysteresis model for cold-formed steel shear wall panel based on its physical and mechanical characteristics, Journal of Thin-Walled Structures, DOI: 10.1016/j.tws.2015.09.022\n";
    numCFSWSWP =1;
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
  double dData[15];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial CFSWSWP tag" << endln;
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

  theMaterial = new CFSWSWP(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]
  , dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type CFSWSWP\n";
    return 0;
  }

  // return the material
  return theMaterial;
}

CFSWSWP::CFSWSWP(int tag,
                  double H, int B, double fuf,
                  double tf, double Ife,
                  double Ifi, double ts,
                  double np, double ds, double Vs,
                  double sc, double nc, double type, double A, double L):
                  UniaxialMaterial(tag, MAT_TAG_CFSWSWP),
                  hight(H), width(B), fuf(fuf),
                  tf(tf), Ife(Ife), Ifi(Ifi), ts(ts), np(np), ds(ds),
                  Vs(Vs),screw_Spacing(sc), nc(nc),type(type), A(A), L(L),
      envlpPosStress(7), envlpPosStrain(7), envlpNegStress(7), envlpNegStrain(7), tagMat(tag),
      gammaDLimit(0.0),
	  gammaFLimit(0.0),
	  gammaE(10.0),
	  rDispP(0.488), rForceP(0.183), uForceP(-0.08), rDispN(0.488), rForceN(0.244), uForceN(-0.08),
	  state3Stress(4), state3Strain(4), state4Stress(4), state4Strain(4), 
	  envlpPosDamgdStress(7), envlpNegDamgdStress(7),
	  TnCycle(0.0), CnCycle(0.0)
	  
{
	double fdeg, ddeg;
	 fdeg = 0.1*((hight/(2*width))*(screw_Spacing/152.0));
	 ddeg = 0.1*((hight/(2*width))*(screw_Spacing/152.0));
	 gammaDLimit = ddeg;
	 gammaFLimit = fdeg;

	
	        // set envelope slopes
	        this->lateralShearStrength();
	        this->SetEnvelope();
	
	        envlpPosDamgdStress = envlpPosStress; envlpNegDamgdStress = envlpNegStress;
	        state3Stress.Zero(); state3Strain.Zero(); state4Stress.Zero(); state4Strain.Zero();
	        // Initialize history variables
	        this->revertToStart();
	        this->revertToLastCommit();
}

void CFSWSWP :: lateralShearStrength(void) {
	E=203000.00;
    int nstud=0;
	int nstude=0;
	double stress;
  
    switch (width)
  	{
		case 610 : nstud=0; nstude=2;  break;
		case 1220: nstud=1; nstude=2; break;
		case 2440: nstud=2; nstude=3; break;
	}
    double Es,Gs,Fu,J,D;
	int inttype=0;
	inttype=floor(type);
	switch (inttype)
	{ 
		case 1 :  Es=10445; Gs=825; Fu=4.5; break;
		case 2 :  Es=9917; Gs=925;Fu=4.2; break;
		case 3 :  Es=7376; Gs=497;Fu=4.5; break;
	}
    J=D=0;
	double l,k;
	double dis = 12.7;
	k=floor ((width/2)/screw_Spacing);
	l=floor ((hight/2)/screw_Spacing);
	double  x[50], y[50], deltay, ey=0,f,a,M,Mp,Cu;
	for (int i=0; i<50; i++)
	{
		x[i]=y[i]=0; 
	}
	for (int i=1; i<=k; i++)
	{
		x[i-1]=i*screw_Spacing;
		y[i-1]=((hight/2)-dis);
		J+=(4*(pow(x[i-1],2)+pow(y[i-1],2)));
	}
    a=(((width/2)/screw_Spacing)-floor((width/2)/screw_Spacing))*screw_Spacing;
	if ((a-dis)>=0)
	J+=((4*pow((width/2)-dis,2)+pow((hight/2)-dis,2)));
    J+=pow((hight/2)-dis,2);
	for (int i=1; i<l; i++)
	{
		x[i-1]=((width/2)-dis);
		y[i-1]=screw_Spacing*i;
		J+=(4*(pow(x[i-1],2)+pow(y[i-1],2)));
	}
    a=(((hight/2)/screw_Spacing)-floor((hight/2)/screw_Spacing))*screw_Spacing;
	if ((a-dis)>=0) 
	 {
	 	a=((width/2)-dis);
	J+=(4*(pow(a,2)+pow((screw_Spacing*l),2))); 
	}
	J+=pow(((width/2)-dis),2);
	deltay=(J/(nc*(hight/2)));
	ey=hight/2+deltay;
	Mp=ey;
    for (int i=0; i<50; i++)
	{
		x[i]=y[i]=0;     
	}
	for (int i=1; i<=k; i++ )
	{
		x[i-1]=i*screw_Spacing;
		y[i-1]=(hight/2-dis)+deltay;
		D+=2*sqrt(pow(x[i-1],2)+pow(y[i-1],2));
	}
    a=(((width/2)/screw_Spacing)-floor((width/2)/screw_Spacing))*screw_Spacing;
	if  ((a-dis)>=0)
	D+=2*sqrt(pow(width/2-dis,2)+pow(hight/2-dis,2));
    D+=hight/2-dis+deltay;
	for (int i=1; i<l; i++)
	{
		x[i-1]=(width/2)-dis;
		y[i-1]=(screw_Spacing*i)+deltay;
		D+=2*sqrt(pow(x[i-1],2)+pow(y[i-1],2));
	}

	D+=2*sqrt(pow(width/2-dis,2)+pow(deltay,2));
    a=(((hight/2)/screw_Spacing)-floor((hight/2)/screw_Spacing))*screw_Spacing;
	if ((a-dis)>=0)
	 {
	 	a=width/2-dis;
     	D+=2*sqrt(pow(a,2)+pow(screw_Spacing*l,2)); 
	 }

	f=floor(deltay/screw_Spacing);
	for (int i=1; i<=f; i++)
	{
		x[i-1]=((width/2)-dis);
		y[i-1]=(screw_Spacing*i);
		D+=(2*sqrt((pow(x[i-1],2)+pow(y[i-1],2))));
	}

    a=((deltay/screw_Spacing)-floor(deltay/screw_Spacing))*screw_Spacing;
	D+=(2*sqrt((pow(a,2)+pow(((width/2)-dis),2))));
	double w;
	w=screw_Spacing-a;
	D+=(2*sqrt((pow(w,2)+pow(((width/2)-dis),2))));

	l=floor(((hight/2)-(deltay+w))/screw_Spacing);
	for (int i=1; i<=k; i++)
	{
		x[i-1]=i*screw_Spacing;
		y[i-1]=(screw_Spacing*l)+w;
		D+=(2*sqrt((pow(x[i-1],2)+pow(y[i-1],2))));
	}
    a=(((width/2)/screw_Spacing)-floor((width/2)/screw_Spacing))*screw_Spacing;
	if ((a-dis)>=0)  
	 {
	    D+=(2*sqrt((pow(((width/2)-dis),2)+pow((screw_Spacing*(l+w)),2))));  
	 }

	for (int i=1; i<l; i++)
	{
		x[i-1]=(width/2)-dis;
		y[i-1]=(screw_Spacing*i)+w;
		D+=(2*sqrt((pow(x[i-1],2)+pow(y[i-1],2))));
	}
	a=(((hight/2)/screw_Spacing)-floor((hight/2)/screw_Spacing))*screw_Spacing;
	if ((a-dis)>=0)
	{
	D+=(2*sqrt((pow(((width/2)-dis),2)+pow((screw_Spacing*l),2)))); 
	}
 	M=0.93*D;
	Cu=M/Mp;
	double n,Bf,Bs,MinPs,Ps,Alphav,Alphab,As,Is,Ks,Kf;
	n=sqrt(8-(hight/width))-1.45;
	Bf=3*tf*ds*fuf;
	Bs=3*ts*ds*Fu;
	MinPs=Bf;
	MinPs=(Bs<MinPs)? Bs:MinPs;
	MinPs=(Vs<MinPs)? Vs:MinPs;
	Ps=Cu*n*MinPs*np;
	Alphav=pow(Cu/(3.3*nc),1.8)*(6/(screw_Spacing/25.4));
	Alphab=pow(6/Cu,2)*pow(6/(screw_Spacing/25.4),(1.3*nc)/Cu);
	if ((Cu>30)&&(Cu<50))
	Alphav=Alphab=0.06;
	As=ts*width;
	double b=width;
	Is=ts*(pow(b,3)/12);
	Ks=((Gs*As)/(1.2*hight))*Alphav+((3*Es*Is)/pow(hight,3))*Alphab;
	Kf=(nstud*3*Ifi*E/pow(hight,3))+((nstude*3*E*Ife)/pow(hight,3));
	double r,fo;
	r=1/(1+A/(hight*(width-L))); 
	fo=r/(3-2*r);
	stress3p=fo*(1+(Kf/Ks))*Ps;  
	strain3p=((stress3p)/(Kf+Ks))/(1000*np);
    stress4p=0.8*stress3p;
	strain4p=1.4*strain3p;
	stress1p=0.4*stress3p;
	strain1p=strain3p/9.25;
	ke=stress1p/strain1p;
	stress2p=0.85*stress3p;
	Dy=(stress2p/ke);
	strain2p=(stress2p*(strain3p+Dy-2*strain4p-strain1p)+stress3p*strain4p+stress4p*(strain4p-strain3p))/(0.6*stress3p);
	strain1n = -strain1p; stress1n = -stress1p; strain2n = -strain2p; stress2n = -stress2p;
	strain3n = -strain3p; stress3n = -stress3p; strain4n = -strain4p; stress4n = -stress4p;
  		
    envlpPosStress.Zero(); envlpPosStrain.Zero(); envlpNegStress.Zero(); envlpNegStrain.Zero(); 
	energyCapacity = 0.0; kunload = 0.0; elasticStrainEnergy = 0.0;
}

 CFSWSWP::CFSWSWP():
   UniaxialMaterial(0, MAT_TAG_CFSWSWP),
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
 
 CFSWSWP::~CFSWSWP()
 {
  
 }
 
int CFSWSWP::setTrialStrain(double strain, double CstrainRate)
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

static int getIndex(Vector v,double value)
{
	for(int i = 0; i < v.Size(); i++)
	{
		if(v[i] > value) return i;
	}
	return -1;
}
static int getIndexNeg(Vector v,double value)
{
	for(int i = 0; i < v.Size(); i++)
	{
		if(v[i] < value) return i;
	}
	return -1;
}


 void CFSWSWP::SetSpline(void)
 {
			
			const int Size = 5;
			double X[Size]; double Y [Size];
			
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
 
 double CFSWSWP::getStrain(void)
 {
         return Tstrain;
 }
 
 double CFSWSWP::getStress(void)
 {
         return Tstress;
 }
 
 double CFSWSWP::getTangent(void)
 {
         return Ttangent;
 }
 
 double CFSWSWP::getInitialTangent(void)
 {
         return envlpPosStress(0)/envlpPosStrain(0);
 }
 
 int CFSWSWP::commitState(void)                                           
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
 
 int CFSWSWP::revertToLastCommit(void)
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
 
 int CFSWSWP::revertToStart(void)
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
 
 UniaxialMaterial* CFSWSWP::getCopy(void)
 {
         CFSWSWP *theCopy = new CFSWSWP (this->getTag(),
                 hight,  width,  fuf,  tf,  Ife, Ifi,  ts, np, ds, Vs, screw_Spacing, nc, type, A, L);
         
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
 
 int CFSWSWP::sendSelf(int commitTag, Channel &theChannel)
 {
         return -1;
 }
 
 int CFSWSWP::recvSelf(int commitTag, Channel &theChannel,
                                                            FEM_ObjectBroker & theBroker)
 {
         return -1;
 }
 
 void CFSWSWP::Print(OPS_Stream &s, int flag)
 {
         s << "CFSWSWP, tag: " << this-> getTag() << endln;
         s << "Displacement: " << Tstrain << endln;
         s << "Strength: " << Tstress << endln;
         s << "state: " << Tstate << endln;
 }
 
 void CFSWSWP::SetEnvelope(void)
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

		 const int Size = 9;
		 double X [Size]; double Y[Size];

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

 double CFSWSWP::GetTangentFromCurve(double Strain)
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
  
  double CFSWSWP::GetStressFromCurve(double Strain)
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
 
 void CFSWSWP::getstate(double u,double du)
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
 
double CFSWSWP::posEnvlpStress(double u)
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
 
 double CFSWSWP::posEnvlpTangent(double u)
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
 
 
  double CFSWSWP::negEnvlpStress(double u)
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
 
 double CFSWSWP::negEnvlpTangent(double u)
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
 
 
 void CFSWSWP::getState3(Vector& state3Strain, Vector& state3Stress, double kunload)
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
                                 
 void CFSWSWP::getState4(Vector& state4Strain,Vector& state4Stress, double kunload)
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
 
 double CFSWSWP::Envlp3Tangent(Vector s3Strain, Vector s3Stress, double u)
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
 
 double CFSWSWP::Envlp4Tangent(Vector s4Strain, Vector s4Stress, double u)
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
 
 
 double CFSWSWP::Envlp3Stress(Vector s3Strain, Vector s3Stress, double u)
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
 
 double CFSWSWP::Envlp4Stress(Vector s4Strain, Vector s4Stress, double u)
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
 
   void CFSWSWP::updateDmg(double strain, double dstrain)
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
