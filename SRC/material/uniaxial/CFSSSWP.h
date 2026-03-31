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
// $Date: 12-10-2015 $

// Written by Smail KECHIDI, Ph.D student at University of Blida 1 (s_kechidi@univ-blida.dz), PhD mobility Student at University of Porto FEUP (smail.kechidi@fe.up.pt)
// Created: 12-10-2015 10:24:20 $
//
// Description: This file contains the class implementation for CFSSSWP
// CFSSSWP is based on Pinching4 uniaxialMaterial
 
#ifndef CFSSSWP_h
#define CFSSSWP_h
 
 #include <UniaxialMaterial.h>
 #include <OPS_Stream.h>
 #include <Vector.h>
 #include"CubicSpline.h"
 
 class CFSSSWP : public UniaxialMaterial
 {
 public :

   CFSSSWP(int tag,
	   double hight, int width, double fuf, double fyf,
	   double tf, double Af, double fus, double fys, double ts,
	   double np, double ds, double Vs, double screw_Spacing, double A, double L);
   
   CFSSSWP();
   ~CFSSSWP();

   const char *getClassType(void) const {return "CFSSSWP";};

   double GetTangentFromCurve(double Strain);
   double GetTangentFromCurve3(double Strain);
   double GetTangentFromCurve4(double Strain);
   double GetStressFromCurve(double Strain);
   double GetStressFromCurve3(double Strain);
   double GetStressFromCurve4(double Strain);
   int setTrialStrain(double strain, double strainRate = 0.0);
   double getStrain(void);
   double getStress(void);
   double getTangent(void);
   double getInitialTangent(void);
   int commitState(void);
   int revertToLastCommit(void);
   int revertToStart(void);
   
   UniaxialMaterial *getCopy(void);
   
   int sendSelf(int commitTag, Channel &theChannel);
   int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);
   void Print(OPS_Stream &s, int flag = 0);
   
 protected:
   
 private:
   
   //  BSpline Adds
   void SetSpline(void);
   
   double* BSplineXs,*BSplineYs,*BSplinePosDamgdYs, *BSplineNegDamgdYs; int BSplineXLength,BSplineYLength;
   CubicSpline Spline3,Spline4;
   
   //  Physical and mechanical characteristics of the panel:
   
   // Shear Wall Panel's Dimensions :
   
   double hight; int width; double A; double L;
   enum {Precision = 100};
   
   // Characteristics and material properties of the steel framing studs :
   
   double fuf; double fyf; double tf; double Ife; double Ifi; double E;
   double Af;
   
   // Characteristics and material properties of sheathing :
   
   double fus; double fys; double ts; double np; double type;
   
   // Characteristics of the screw fasteners :
   
   double ds; double screw_Spacing; double nc; double Vs; 
   
   // Backbone parameters
   
   double stress1p; double strain1p; double stress2p; double strain2p;
   double stress3p; double strain3p; double stress4p; double strain4p;
   double stress1n; double strain1n; double stress2n; double strain2n;
   double stress3n; double strain3n; double stress4n; double strain4n;
   double Dy; double ke;
   
   Vector envlpPosStress; Vector envlpPosStrain; 
   Vector envlpNegStress; Vector envlpNegStrain;
   
   int tagMat;  // material tag
   
   // Damage parameters
   
   double gammaDLimit;
   double gammaFLimit;
   double gammaE;
   double TnCycle, CnCycle;
   
   // unloading-reloading parameters
   double rDispP; double rForceP; double uForceP;
   double rDispN; double rForceN; double uForceN;
   
   Vector state3Stress; Vector state3Strain; Vector state4Stress; Vector state4Strain;
   
   Vector envlpPosDamgdStress; Vector envlpNegDamgdStress;
   
   // Trial State Variables
   double Tstress;
   double Tstrain;
   double Ttangent;
   
   // Converged Material History parameters
   int Cstate;
   double Cstrain;
   double Cstress;
   double CstrainRate;
   double lowCstateStrain;
   double lowCstateStress;
   double hghCstateStrain;
   double hghCstateStress;
   double CminStrainDmnd;
   double CmaxStrainDmnd;
   double Cenergy;
   double CgammaD;
   double CgammaDN;
   double CgammaF;
   double CgammaFN;
   double gammaFUsed;
   
   // Trial Material History Parameters
   int Tstate;
   double dstrain;
   double TstrainRate;
   double lowTstateStrain;
   double lowTstateStress;
   double hghTstateStrain;
   double hghTstateStress;
   double TminStrainDmnd;
   double TmaxStrainDmnd;
   double Tenergy;
   double TgammaD;
   double TgammaDN;
   double TgammaF;
   double TgammaFN;
   
   // strength and stiffness parameters;
   double kElasticPos;
   double kElasticNeg;
   double uMaxDamgd;
   double uMinDamgd;
 
   
   // energy parameters
   double energyCapacity;
   double kunload;
   double elasticStrainEnergy;
   void lateralShearStrength(void);
   void SetEnvelope(void);
   void getstate(double, double);
   double posEnvlpStress(double);
   double posEnvlpTangent(double);
   double negEnvlpStress(double);
   double negEnvlpTangent(double);
   void getState3(Vector& , Vector& , double);
   void getState4(Vector& , Vector& , double);
   double Envlp3Tangent(Vector , Vector , double);
   double Envlp3Stress(Vector , Vector , double);
   double Envlp4Tangent(Vector , Vector , double);
   double Envlp4Stress(Vector , Vector , double);
   void updateDmg(double, double);
 };
 #endif
