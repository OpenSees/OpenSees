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
// Added by Liming Jiang (UoE), Mian Zhou (Brunel)
// Created: 06/13
// --------------------------------------------------------------------
// Description: This file contains the class definition for
// StainlessECThermal. StainlessECThermal is modified on the basis of Steel02Thermal
// and steel01Thermal.StainlessECThermal is developed for modelling steel material
// which strictly satisfies Eurocode regarding the temperature dependent properties.

#include <StainlessECThermal.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include<iostream>
using namespace std;

void *
OPS_StainlessECThermal(void)
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;
	int    iData[2];
	double dData[4];
	
	int numData = 1;
	//string gradeInput;

	if (OPS_GetIntInput(&numData, iData) != 0)
	{
		opserr << "WARNING invalid uniaxialMaterial StainlessECThermal tag?" << endln;
		return 0;
	}

	//
	//numData = OPS_GetNumRemainingInputArgs();
	//if (numData == 4 || numData == 8)
	//{
	//	//if (OPS_GetString(gradeChar, 20) == 0)
	//	//{
	//	//	opserr << "WARNING invalid gradeTag for uniaxialMaterial StainlessECThermal " << iData[1] << endln;
	//	//	//return 0;
	//	// }
	//	
	//	//gradeInput = OPS_GetString(gradeChar, 20);
	const char* gradeChar = OPS_GetString();
	if (strcmp(gradeChar, "Grade14301") == 0)
		{
			iData[1] = 1;
		}
	else if ((strcmp(gradeChar,"Grade14401")==0 || strcmp(gradeChar, "Grade14404") == 0))
		{ 
			iData[1] = 2; 
		}
	else if (strcmp(gradeChar, "Grade14571") ==0) {iData[1]=3;}
	else if (strcmp(gradeChar, "Grade14003") ==0) {iData[1]=4;}
	else if (strcmp(gradeChar, "Grade14462") == 0) {iData[1] = 5;}
	else
		{
			opserr << "WARNING invalid material grade for uniaxialMaterial StainlessECThermal "<<iData[0]<<endln;
			return 0;
		}
	numData = OPS_GetNumRemainingInputArgs();
	if (numData != 3 && numData != 4)
	{
			opserr << "Invalid #args, want: uniaxialMaterial StainlessECThermal " << iData[0] << " fy? E? fu?" << endln;
			return 0;
	}
//

	if (OPS_GetDoubleInput(&numData, dData) != 0)
	{
			opserr << "Invalid #args, want: uniaxialMaterial StainlessECThermal " << iData[0] << " fy? E? fu?" << endln;
			return 0;
	}
//
	if (numData == 3)
	{
		dData[3] = 0.0;  // set Initial stress=0.0
	}

		// Parsing was successful, allocate the material
		theMaterial = new StainlessECThermal(iData[0], iData[1], dData[0], dData[1], dData[2],dData[3]);

		if (theMaterial == 0) {
			opserr << "WARNING could not create uniaxialMaterial of type StainlessECThermal Material\n";
			return 0;
		}

		return theMaterial;
	}

//StainlessECThermal::StainlessECThermal(int tag, int grade, double Fy, double E, double Fu, double A1, double A2, double A3, double A4): UniaxialMaterial(tag,MAT_TAG_StainlessECThermal),gradeTag(grade), fyT(Fy), E0T(E),fuT(Fu), a1(A1), a2(A2), a3(A3), a4(A4)
StainlessECThermal::StainlessECThermal(int tag, int grade, double Fy, double E, double Fu,double sigInit) : UniaxialMaterial(tag, MAT_TAG_StainlessECThermal), gradeTag(grade)
{
   // Sets all history and state variables to initial values
   // History variables
   // converged history variables
	fyT = Fy;
    E0T = E;
	fuT = Fu;
	sigini = sigInit;
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;
   //Trial history variables
   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;
//
   // State variables
   //converged history variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0T;
   Ctemp=0;
   // trial history variables
   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0T;
   Ttemp = 0;
//
// // initialise values
	  ThermalElongation = 0;
	  E0 = E0T;
	  fy = fyT;
	  fu = fuT;
	  
// determine the initial strain
	  if (sigini != 0.0)
	  {
		  epsini = determineYieldSurface(sigini);
		  Cstrain = epsini;
		  Cstress = sigini;
#ifdef _DEBUG
		 // opserr << "  CstrainINI:  " << Cstrain << "  CStressINI: " << Cstress << endln;
#endif
		  // need to initialise Tstrian? Tstress?? or they will be assigned in the setTrial(), setTrialStrain()anyway	  
	  }
	  else
		  epsini = 0.0;


//
//
// Grade 1.4301
	  if (gradeTag == 1)
	  {
		  EpsiUT = 0.4; //  initial  value at temp=20C, grade type & temperature dependent
		  EctT = 0.11*E0; //initial value at temp 20C,grade type & temperature dependent
	  }
	  // Grade 1.4401/1.4404
	  else if (gradeTag == 2)
	  {
		  EpsiUT = 0.4; //  initial  value at temp=20C, grade type & temperature dependent
		  EctT = 0.05*E0; //initial value at temp 20C,grade type & temperature dependent
	  }
	  // Grade 1.4571
	  else if (gradeTag == 3)
	  {
		  EpsiUT = 0.4; //  initial  value at temp=20C, grade type & temperature dependent
		  EctT = 0.06*E0; //initial value at temp 20C,grade type & temperature dependent
	  }
	  // Grade 1.4003
	  else if (gradeTag == 4)
	  {
		  EpsiUT = 0.2; //  initial  value at temp=20C, grade type & temperature dependent
		  EctT = 0.055*E0; //initial value at temp 20C,grade type & temperature dependent
	  }
	  // Grade 1.4462
	  else if (gradeTag == 5)
	  {
		  EpsiUT = 0.2; //  initial  value at temp=20C, grade type & temperature dependent
		  EctT = 0.1*E0; //initial value at temp 20C,grade type & temperature dependent
	  }

	  EpsiU = EpsiUT;
	  Ect = EctT;
//	  

}
// default constructor
StainlessECThermal::StainlessECThermal():UniaxialMaterial(0,MAT_TAG_StainlessECThermal),gradeTag(0),fyT(0.0), E0T(0.0), fuT(0.0),sigini(0.0)
{
	  ThermalElongation = 0; //initialize
	  E0 = E0T;
	  fy = fyT;
	  fu = fuT; // added in ultimate strength, MZ 07/16
	  Ttemp=0;
	  Ctemp = 0;
	  EpsiU = 0;
	  Ect = 0;
}

StainlessECThermal::~StainlessECThermal ()
{
	// do nothing
}
//
double StainlessECThermal::determineYieldSurface(double sigini)  // Added by Mian Zhou 10/16
{

	double fabsSigini = fabs(sigini);
	if (fabsSigini < fyT)
	{
		epsini = sigini / E0T;
	}

	else if (fabsSigini = fyT)
	{
		if (sigini > 0)   // positive/tensile initial stress
		{
			epsini = 0.02;
		}
		else {//negative/compressive initial stress
			epsini = -0.02;
		}
	}
	else
	{
		opserr << "WARNING: Initial Stress Exceeds Plastic Yield strength " << endln;
	}
	return epsini;
}

int StainlessECThermal::setTrialStrain(double strain, double FiberTemperature, double strainRate)
{
  Ttemp = FiberTemperature;
 // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain; 
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
//
// Determine change in strain from last converged state
   double dStrain = strain - Cstrain;
//
//opserr<< "Ttemp: "<<Ttemp<<"  Ctemp: "<< Ctemp<<endln;
   if (fabs(dStrain) > DBL_EPSILON ||Ttemp>Ctemp) {
     // Set trial strain
     Tstrain = strain+ epsini;
     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }
   return 0;
}

int StainlessECThermal::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
 //
   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;
//
   if (fabs(dStrain) > DBL_EPSILON ||Ttemp!= Ctemp) 
   {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);
   }
//
   stress = Tstress;
   tangent = Ttangent;
//
   return 0;
}

void StainlessECThermal::determineTrialState (double dStrain)
{

	if (Tloading ==0)
	{
		if(dStrain >0)
			Tloading =1;
		else
			Tloading =-1;
	}
//Changed on Aug 5 by Liming;
    if(fabs(Ttemp- Ctemp)>1e-5)
	{
		if(Cloading!=0){
			Tloading = Cloading;}
	}
	else
	{
		if(Tstrain>0)
			Tloading =1;
		else if(Tstrain <0)
			Tloading =-1;
		else
			{
			if(Cstrain > 0)
				Tloading = 1;
			else
			    Tloading =-1;
			}
	}
// define stress-strain relationship based on EN1993-1-2 Annex C
// calculating the Positive stress
	/*  E0 = E0T; //trial
	  fy = fyT;
	  fu = fuT;
	  EpsiU = EpsiUT;
	  Ect = EctT;  //trial
	  */
	 double EpsiC = fy/E0+0.002;
	 double EpsiU_C = EpsiU - EpsiC;
	 double fu_y = fu - fy;
	 double ET = fu_y*fu_y / (EpsiU_C*Ect - 2 * fu_y);
	 double DT = pow((ET*EpsiU_C*Ect + ET*ET), 0.5);
	 double CT = pow(EpsiU_C*(EpsiU_C + ET / Ect), 0.5);
	 double BT = (1 - EpsiC*Ect / fy)*E0*EpsiC / ((E0*EpsiC / fy - 1)*fy);
	 double EpsiC_BT = pow(EpsiC, BT);
	 double AT = (E0*EpsiC - fy) / (fy* EpsiC_BT);
//
#ifdef _bDEBUG
	 opserr << "  AT:  " << AT << "  BT: " << BT << "  CT:  " << CT << endln;
	 opserr << "  DT:  " << DT << "  ET: " << ET << endln;
#endif
    double fabsTstrain = fabs(Tstrain);
// opserr<<"fabsTstrain is: "<< fabsTstrain <<" at Temp "<<Ttemp<<endln;
 // opserr<<"EpsiPT: "<<EpsiPT<< " fp: "<< fp <<" E0: "<<E0<<endln;
	 if (fabsTstrain <= EpsiC)
		 {
			 double Epsi_pow_b = pow(fabsTstrain, BT);
			 Tstress = E0*fabsTstrain / (1 + AT*Epsi_pow_b);
			 Ttangent = E0*(1 + AT*Epsi_pow_b - AT*BT*Epsi_pow_b)/((1+AT*Epsi_pow_b)*(1 + AT*Epsi_pow_b));
		  }
	else if (fabsTstrain <= EpsiU)
		  {
			  double EpsiU_Epsi = EpsiU - fabsTstrain;
			  Tstress = fy - ET + (DT / CT)*pow(CT*CT - (EpsiU_Epsi)*(EpsiU_Epsi), 0.5);
			  Ttangent= DT*(EpsiU_Epsi) /(CT * pow(CT*CT - (EpsiU_Epsi)*(EpsiU_Epsi), 0.5));
			
		  }
	else if(fabsTstrain <= EpsiU+0.01){
		
		Tstress = fu*(1 - (fabsTstrain - EpsiU) / 0.01);
		//opserr<<"Error: Stiffness of SteelECthermal is negative"<<endln;
		Ttangent = -fu/ 0.01;
	}
	else{
			  Tstress = 1E-10;
			  Ttangent = 1E-10;
	}
	
//#ifdef 
	 //	opserr<<"Ttangent:"<<Ttangent<<endln;
//#endif
   	if(Tloading == 1)
	{
	  Tstress = Tstress;
	}
	else if(Tloading ==-1) // this needs to get changed, as the stress-strain relation is based on Ulitmate stensile strength Fu
	{
		Tstress = -Tstress;
	}
	else
	{
		Tstress = 0;
	}

	//if (Ttemp> 450)
		//opserr << "Trial strain: " << Tstrain << "   tangent: " << Ttangent << "   Tstress: " << Tstress << "Thelong: " << ThermalElongation << endln;

#ifdef _BDEBUG
	opserr<< "Trial strain: "<<Tstrain<< "   tangent: "<< Ttangent <<"   Tstress: "<<Tstress<<endln;
#endif
	Ttangent = 1e11;
    Ctemp = Ttemp;

}
/*
void StainlessECThermal::detectLoadReversal (double dStrain)
{
   // Determine initial loading condition
   if (Tloading == 0 && dStrain != 0.0)
   {
      if (dStrain > 0.0)
         Tloading = 1;
      else
         Tloading = -1;
   }

   double epsy = fy/E0;

   // Transition from loading to unloading, i.e. positive strain increment
   // to negative strain increment
   if (Tloading == 1 && dStrain < 0.0)
   {
      Tloading = -1;
      if (Cstrain > TmaxStrain)
         TmaxStrain = Cstrain;
      TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
   }

   // Transition from unloading to loading, i.e. negative strain increment
   // to positive strain increment
   if (Tloading == -1 && dStrain > 0.0)
   {
      Tloading = 1;
      if (Cstrain < TminStrain)
         TminStrain = Cstrain;
      TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
   }
}
*/
double StainlessECThermal::getStrain ()
{
   return Tstrain;
}

double StainlessECThermal::getStress ()
{
   return Tstress;
}

double StainlessECThermal::getTangent ()
{
   return Ttangent;
}

double
StainlessECThermal::getThermalElongation(void)
{
  return ThermalElongation;
}

//Updating the reduction factors////////////start
double
StainlessECThermal::getElongTangent(double TempT, double &ET, double &Elong, double TempTmax)
{
	double FyRfactors[12];
	double FuRfactors[12];
	double E0Rfactors[12];
	double EctRfactors[12];
	double EpsiURfactors[12];
	//double FyR[12];
	//double FuR[12];
	//double E0R[12];
	//double EctR[12];
	//double EpsiUR[12];
// Grade 1.4301
	if (gradeTag == 1)
	{
		double FyR[12]= { 0.82,0.68, 0.64, 0.60, 0.54, 0.49, 0.40, 0.27, 0.14 ,0.06, 0.03, 0.0 };
		double FuR[12] = { 0.87,0.77,0.73, 0.72 ,0.67, 0.58, 0.43, 0.27, 0.15, 0.07 ,0.03, 0.0 };
		double E0R[12] = { 0.96, 0.92 ,0.88, 0.84 ,0.80, 0.76,0.71, 0.63, 0.45, 0.20, 0.10 , 0.0 };
		double EctR[12] = { 0.05, 0.02 ,0.02, 0.02 ,0.02, 0.02,0.02, 0.02, 0.02, 0.02, 0.02 , 0.02 };
		double EpsiUR[12] = { 0.4, 0.4,0.4, 0.4 ,0.4, 0.35,0.3, 0.20, 0.2, 0.20, 0.20 , 0.20 };
		for (int i = 0; i < 11; i++)
		{
			FyRfactors[i] = FyR[i];
			FuRfactors[i] = FuR[i];
			E0Rfactors[i] = E0R[i];
			EctRfactors[i] = EctR[i];
			EpsiURfactors[i] = EpsiUR[i];
		}

	}
	// Grade 1.4401/1.4404
	else if(gradeTag==2)
	{
		double FyR[12] = {0.88 ,0.76, 0.71, 0.66, 0.63, 0.61, 0.51, 0.40, 0.19 ,0.10, 0.005, 0.0};
		double FuR[12] = {0.93,0.87,0.84, 0.83,0.79, 0.72, 0.55, 0.34, 0.18, 0.09 ,0.04, 0.0};
		double E0R[12] = {0.96, 0.92 ,0.88, 0.84 ,0.80, 0.76, 0.71,0.63, 0.45, 0.20, 0.10 , 0.0};
		double EctR[12] = { 0.049, 0.047 ,0.045, 0.030 ,0.025, 0.02,0.02, 0.02, 0.02, 0.02, 0.02 , 0.02 };
		double EpsiUR[12] = { 0.4, 0.4,0.4, 0.4 ,0.4, 0.40,0.30, 0.20, 0.2, 0.20, 0.20 , 0.20 };
		for (int i = 0; i < 11; i++)
		{
			FyRfactors[i] = FyR[i];
			FuRfactors[i] = FuR[i];
			E0Rfactors[i] = E0R[i];
			EctRfactors[i] = EctR[i];
			EpsiURfactors[i] = EpsiUR[i];
		}
	}
	// Grade 1.4571
	else if (gradeTag == 3)
	{
		double FyR[12] = { 0.89 ,0.83, 0.77, 0.72, 0.69, 0.66, 0.59, 0.50, 0.28 ,0.15, 0.075, 0.0 };
		double FuR[12] = { 0.88,0.81,0.80, 0.80,0.77, 0.71, 0.57, 0.38, 0.22, 0.11 ,0.055, 0.0 };
		double E0R[12] = { 0.96, 0.92 ,0.88, 0.84 ,0.80, 0.76, 0.71,0.63, 0.45, 0.20, 0.10 , 0.0 };
		double EctR[12] = { 0.06, 0.05 ,0.04, 0.030 ,0.025, 0.02,0.02, 0.02, 0.02, 0.02, 0.02 , 0.02 };
		double EpsiUR[12] = { 0.4, 0.4,0.4, 0.4 ,0.4, 0.35,0.30, 0.20, 0.2, 0.20, 0.20 , 0.20 };
		for (int i = 0; i < 11; i++)
		{
			FyRfactors[i] = FyR[i];
			FuRfactors[i] = FuR[i];
			E0Rfactors[i] = E0R[i];
			EctRfactors[i] = EctR[i];
			EpsiURfactors[i] = EpsiUR[i];
		}
	}
	// grade 1.4003
	else if(gradeTag == 4)
	{
		double FyR[12] = { 1.0 ,1.0, 0.98, 0.91, 0.80, 0.45, 0.19, 0.13, 0.10 ,0.07, 0.035, 0.0 };
		double FuR[12] = { 0.94,0.88,0.86, 0.83 ,0.81, 0.42, 0.21, 0.12, 0.11, 0.09 ,0.045, 0.0 };
		double E0R[12] = { 0.96, 0.92 ,0.88, 0.84 ,0.80, 0.76,0.71, 0.63, 0.45, 0.20, 0.10 , 0.0 };
		double EctR[12] = { 0.03, 0.03 ,0.03, 0.030 ,0.03, 0.03,0.03, 0.03, 0.03, 0.03, 0.03 , 0.03 };
		double EpsiUR[12] = { 0.2, 0.2,0.2, 0.15 ,0.15, 0.15,0.15, 0.15, 0.15, 0.15, 0.15 , 0.15 };
		for (int i = 0; i < 11; i++)
		{
			FyRfactors[i] = FyR[i];
			FuRfactors[i] = FuR[i];
			E0Rfactors[i] = E0R[i];
			EctRfactors[i] = EctR[i];
			EpsiURfactors[i] = EpsiUR[i];
		}
	}
	// grade 1.4462
	else if (gradeTag == 5)
	{
		double FyR[12] = { 0.91 ,0.80, 0.75, 0.72, 0.65, 0.56, 0.37, 0.26, 0.10 ,0.03, 0.015, 0.0 };
		double FuR[12] = { 0.93,0.85,0.83, 0.82,0.71, 0.57, 0.38, 0.29, 0.12, 0.04,0.02, 0.0 };
		double E0R[12] = { 0.96, 0.92 ,0.88, 0.84 ,0.80, 0.76, 0.71,0.63, 0.45, 0.20, 0.10 , 0.0 };
		double EctR[12] = { 0.07, 0.037 ,0.035, 0.030 ,0.03, 0.025,0.025, 0.025, 0.025, 0.025, 0.025, 0.025 };
		double EpsiUR[12] = { 0.2, 0.2,0.2, 0.20 ,0.20, 0.20,0.15, 0.15, 0.15, 0.15, 0.15 , 0.15 };
		for (int i = 0; i < 11; i++)
		{
			FyRfactors[i] = FyR[i];
			FuRfactors[i] = FuR[i];
			E0Rfactors[i] = E0R[i];
			EctRfactors[i] = EctR[i];
			EpsiURfactors[i] = EpsiUR[i];
		}
	}
	else
	{
		opserr << "WARNING StainlessECThermal received an invalid gradeTag: " << gradeTag << endln;
	}
//
   //Now Updating modulus, strengths
	for (int i = 0; i < 12; i++)
		{
			if (TempT <= 80 + 100 * i)
			{
				if (i == 0)
				{
					fy = fyT*(1.0 - TempT*(1.0 - FyRfactors[0]) / 80);
					fu = fuT*(1.0 - TempT*(1.0 - FuRfactors[0]) / 80);
					E0 = E0T*(1.0 - TempT*(1.0 - E0Rfactors[0]) / 80);
					Ect = E0T*(0.11 - TempT*(0.11- EctRfactors[0]) / 80);
					EpsiU = EpsiUT - (EpsiUT - EpsiURfactors[0]) *TempT / 80;
				}
				else if (i == 12) {
					opserr << "Warning:The temperature " << TempT << " for StainlessECThermal is out of range\n";
					return -1;
				}
				else
				{
					fy = fyT*(FyRfactors[i - 1] - (TempT + 20 - 100 * i)*(FyRfactors[i - 1] - FyRfactors[i]) / 100);
					fu = fuT*(FuRfactors[i - 1] - (TempT + 20 - 100 * i)*(FuRfactors[i - 1] - FuRfactors[i]) / 100);
					E0 = E0T*(E0Rfactors[i - 1] - (TempT + 20 - 100 * i)*(E0Rfactors[i - 1] - E0Rfactors[i]) / 100);
					Ect = E0T*(EctRfactors[i - 1] - (TempT + 20 - 100 * i)*(EctRfactors[i - 1] - EctRfactors[i]) / 100);
					EpsiU = EpsiURfactors[i - 1] - (TempT + 20 - 100 * i)*(EpsiURfactors[i - 1] - EpsiURfactors[i]) / 100;
				}
				break;
			}
		}
#ifdef _DEBUG
	//opserr<<", TempT:"<<TempT<< " fy: "<< fy <<" E0T:  "<< E0<<endln;
#endif
// calculation of thermal elongation
///*
  if (TempT <= 1)
	{
		ThermalElongation = 1.61e-10;
	}
  else if (TempT <= 1200)
    {
	  double RealTemp = TempT + 20;
	  //ThermalElongation = 12e-6*RealTemp;
	  //double RealTemp = TempT + 20;
	  ThermalElongation = (16 + 4.79e-3*RealTemp - 1.243e-6*RealTemp*RealTemp)*TempT *1e-6;
    }
  else 
  {
	  opserr << "the temperature is invalid\n";
	  return -1;
  }

  //ThermalElongation = 0 ;   //debug  Liming
  ET = E0;
  Elong = ThermalElongation;
  TemperautreC = TempT;
//
#ifdef _dDEBUG
  opserr <<"  TempT:  "<< TempT<< "  Ttangent: "<< Ttangent <<"  E0:  "<<E0<< endln;
#endif
  return 0;
}

//END  for updating the reduction factors


int StainlessECThermal::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int StainlessECThermal::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int StainlessECThermal::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0)
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* StainlessECThermal::getCopy ()
{
   StainlessECThermal* theCopy = new StainlessECThermal(this->getTag(), gradeTag,fy, E0,fu, sigini);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;

   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->TshiftP = TshiftP;
   theCopy->TshiftN = TshiftN;
   theCopy->Tloading = Tloading;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}

int StainlessECThermal::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(17);
   data(0) = this->getTag();

   // Material properties
   data(1) = gradeTag;
   data(2) = fy;
   data(3) = E0;
   data(4) = fu;

   // History variables from last converged state
   data(5) = CminStrain;
   data(6) = CmaxStrain;
   data(7) = CshiftP;
   data(8) = CshiftN;
   data(9) = Cloading;

   // State variables from last converged state
   data(10) = Cstrain;
   data(11) = Cstress;
   data(12) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0)
      opserr << "StainlessECThermal::sendSelf() - failed to send data\n";

   return res;
}

int StainlessECThermal::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(16);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "StainlessECThermal::recvSelf() - failed to receive data\n";
      this->setTag(0);
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      gradeTag = data (1);
	  fy = data(2);
      E0 = data(3);
	  fu = data(4);
 

      // History variables from last converged state
      CminStrain = data(5);
      CmaxStrain = data(6);
      CshiftP = data(7);
      CshiftN = data(8);
      Cloading = int(data(9));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;

      // State variables from last converged state
      Cstrain = data(10);
      Cstress = data(11);
      Ctangent = data(12);

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

void StainlessECThermal::Print (OPS_Stream& s, int flag)
{
   s << "StainlessECThermal tag: " << this->getTag() << endln;
   s << "gradeTag: " <<gradeTag << " ";
   s << "  fy: " << fy << " ";
   s << "  E0: " << E0 << " ";
   s << "  fu: " << fu << " ";
}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
StainlessECThermal::setParameter(const char **argv, int argc, Parameter &param)
{
/*
  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"b") == 0)
    return param.addObject(3, this);

  if (strcmp(argv[0],"a1") == 0)
    return param.addObject(4, this);

  if (strcmp(argv[0],"a2") == 0)
    return param.addObject(5, this);

  if (strcmp(argv[0],"a3") == 0)
    return param.addObject(6, this);

  if (strcmp(argv[0],"a4") == 0)
    return param.addObject(7, this);
*/
  return -1;
}



int
StainlessECThermal::updateParameter(int parameterID, Information &info)
{
/*	
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->fy = info.theDouble;
		break;
	case 2:
		this->E0 = info.theDouble;
		break;
	case 3:
	this->b = info.theDouble;
		break;
	case 4:
		this->a1 = info.theDouble;
		break;
	case 5:
		this->a2 = info.theDouble;
		break;
	case 6:
		this->a3 = info.theDouble;
		break;
	case 7:
		this->a4 = info.theDouble;
		break;
	default:
		return -1;
	}

	Ttangent = E0;          // Initial stiffness
	*/
	return 0;
}




int
StainlessECThermal::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}



double
StainlessECThermal::getStressSensitivity(int gradIndex, bool conditional)
{
	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;
	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + E0*dStrain;
	double fyOneMinusB = fy * (1.0 - b);
	double Esh = b*E0;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity*b*Tstrain
				 + E0*bSensitivity*Tstrain
				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity
			     + E0Sensitivity*(Tstrain-Cstrain)
				 - E0*CstrainSensitivity;
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity*b*Tstrain
			     + E0*bSensitivity*Tstrain
				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
	}
	
	return gradient;
}




double
StainlessECThermal::getInitialTangentSensitivity(int gradIndex)
{
	// For now, assume that this is only called for initial stiffness
	if (parameterID == 2) {
		return 1.0;
	}
	else {
		return 0.0;
	}
}


int
StainlessECThermal::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2,numGrads);
	}


	// Initialize unconditaional stress sensitivity
	double gradient = 0.0;

	
	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity	 = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;
	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + E0*dStrain;
	double fyOneMinusB = fy * (1.0 - b);
	double Esh = b*E0;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity ('gradient')
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity*b*Tstrain
				 + E0*bSensitivity*Tstrain
				 + E0*b*TstrainSensitivity
				 + TshiftP*(fySensitivity*(1-b)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity
			     + E0Sensitivity*(Tstrain-Cstrain)
				 + E0*(TstrainSensitivity-CstrainSensitivity);
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity*b*Tstrain
			     + E0*bSensitivity*Tstrain
			     + E0*b*TstrainSensitivity
				 - TshiftN*(fySensitivity*(1-b)-fy*bSensitivity);
	}


	// Commit history variables
	(*SHVs)(0,gradIndex) = TstrainSensitivity;
	(*SHVs)(1,gradIndex) = gradient;

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////

//this function is no use, just for the definiation of pure virtual function.
int StainlessECThermal::setTrialStrain (double strain, double strainRate)
{
  opserr << "StainlessECThermal::setTrialStrain (double strain, double strainRate) - should never be called\n";
  return 0;
}


int
StainlessECThermal::getVariable(const char *variable, Information &info)  // is there need to get Ect and EpsiU?? MZ 07/16
{
  if (strcmp(variable,"ThermalElongation") == 0) {
    info.theDouble = ThermalElongation;
    return 0;
  } else if (strcmp(variable,"ElongTangent") == 0) {
    Vector *theVector = info.theVector;
    if (theVector != 0) {
      double tempT, ET, Elong, TempTmax;
      tempT = (*theVector)(0);
	  ET = (*theVector)(1);
	  Elong = (*theVector)(2);
      TempTmax = (*theVector)(3);
      this->getElongTangent(tempT, ET, Elong, TempTmax);
	  (*theVector)(0) = tempT;
      (*theVector)(1) = ET;
      (*theVector)(2) = Elong;
	  (*theVector)(3) = TempTmax;
    }
    return 0;
  }else if (strcmp(variable,"TempAndElong") == 0) {
    Vector *theVector = info.theVector;
	if (theVector!= 0) {
		(*theVector)(0) = Ttemp;
        (*theVector)(1) = ThermalElongation;
	}else{
		opserr<<"null Vector in EC"<<endln;
	}

	return 0;
  }
  return -1;
}


