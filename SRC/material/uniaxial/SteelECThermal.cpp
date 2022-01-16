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
// Added by Liming Jiang (UoE)
// Created: 06/13
// --------------------------------------------------------------------
// Description: This file contains the class definition for 
// SteelECThermal. SteelECThermal is modified on the basis of Steel02Thermal
// and steel01Thermal.SteelECthermal is developed for modelling steel material 
// which strictly satisfies Eurocode regarding the temperature dependent properties.





#include <SteelECThermal.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


void * OPS_ADD_RUNTIME_VPV(OPS_SteelECThermal)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;


  int    iData[2];
  double dData[6];
  const char *typeChar = new char[20];
  int numData = 1;
  
 if (OPS_GetIntInput(&numData, iData) != 0) 
 {
    opserr << "WARNING invalid uniaxialMaterial SteelECThermal tag?" << endln;
    return 0;
  }
  
  if (OPS_GetNumRemainingInputArgs()==3 ||OPS_GetNumRemainingInputArgs()==7)
  {
    typeChar = OPS_GetString();

	  if(strcmp(typeChar,"EC3") == 0 ){
		  iData[1]=3;
		  }
	  else if (strcmp(typeChar, "EC2Nh") ==0 ||strcmp(typeChar, "EC2NH") ==0) {
		  iData[1]=21;
		  }
	  else if (strcmp(typeChar, "EC2NC") ==0 ||strcmp(typeChar, "EC2Nc") ==0) {
		  iData[1]=22;
		  }
	  else if (strcmp(typeChar, "EC2X") ==0 ||strcmp(typeChar, "EC2x") ==0) {
		  iData[1]=23;
		  }
	  else{
		  opserr << "WARNING invalid material type for uniaxialMaterial SteelECThermal "<<iData[0]<<endln;
		  return 0;
		  }
	 
  }
  else if (OPS_GetNumRemainingInputArgs()==2 ||OPS_GetNumRemainingInputArgs()==6) 
  {
    iData[1] = 0;
  }



  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 2 && numData != 6) {
    opserr << "Invalid #args, want: uniaxialMaterial SteelECThermal " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial SteelECThermal " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 2) {
    dData[2] = STEEL_01_DEFAULT_A1;
    dData[3] = STEEL_01_DEFAULT_A2;
    dData[4] = STEEL_01_DEFAULT_A3;
    dData[5] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SteelECThermal(iData[0], iData[1], dData[0], dData[1], 
				   dData[2], dData[3], dData[4], dData[5]);
  
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type SteelECThermal Material\n";
    return 0;
  }

  return theMaterial;
}


SteelECThermal::SteelECThermal
(int tag, int TypeTag , double FY, double E,
 double A1, double A2, double A3, double A4):
   UniaxialMaterial(tag,MAT_TAG_SteelECThermal),typeTag(TypeTag),
   fyT(FY), E0T(E), a1(A1), a2(A2), a3(A3), a4(A4)
{
   // Sets all history and state variables to initial values
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
   Ttemp=0;
   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   //Ctangent = E0;

   Ctangent = E0T;///JZ, 07/10//
   Ctemp=0;
 
   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0T;///JZ, 07/10//

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

	  ThermalElongation = 0; //initialize //JZ, 07/10//
	  E0 = E0T;
	  fy = fyT;
	  fp = fyT;
}

SteelECThermal::SteelECThermal():UniaxialMaterial(0,MAT_TAG_SteelECThermal),
 //fy(0.0), E0(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
 typeTag(0),fyT(0.0), E0T(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0) //JZ, 07/10//
{

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

	  ThermalElongation = 0; //initialize //JZ, 07/10//
	  E0 = E0T;
	  fy = fyT;
	  fp = fyT;
	  Ttemp=0;
	  Ctemp = 0;

}

SteelECThermal::~SteelECThermal ()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int SteelECThermal::setTrialStrain(double strain, double FiberTemperature, double strainRate)
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

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;
   
   //opserr<< "Ttemp: "<<Ttemp<<"  Ctemp: "<< Ctemp<<endln; 
   if (fabs(dStrain) > DBL_EPSILON ||Ttemp>Ctemp) {
     // Set trial strain
     Tstrain = strain;
     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);
     
   }

   return 0;
}

int SteelECThermal::setTrial (double strain, double &stress, double &tangent, double strainRate)
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

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON ||Ttemp!= Ctemp) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void SteelECThermal::determineTrialState (double dStrain)
{
	
	if (Tloading ==0)
	{
		if(dStrain >0)
			Tloading =1;
		else 
			Tloading =-1;
	}
//Changed on Aug 5 by Liming;
    if(Ttemp!= Ctemp)
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

	  double EpsiPT = fp/E0;
	  double EpsiYT = 0.02;
	  double EpsiT = 0.15;
	  double EpsiU = 0.2;
	  double CT = (fy-fp)*(fy-fp)/((EpsiYT-EpsiPT)*E0-2*(fy - fp));
	  double BT = pow(CT*(EpsiYT-EpsiPT)*E0+CT*CT, 0.5);
	  double AT = pow((EpsiYT-EpsiPT)*(EpsiYT-EpsiPT+CT/E0),0.5);
	
	  //Calculating the POSITIVE stress according to EC3
      double fabsTstrain = fabs(Tstrain);
	 // opserr<<"fabsTstrain is: "<< fabsTstrain <<" at Temp "<<Ttemp<<endln;
	 // opserr<<"EpsiPT: "<<EpsiPT<< " fp: "<< fp <<" E0: "<<E0<<endln;
	  if (fabsTstrain <= EpsiPT) 
	  	{
		  Tstress = E0*fabsTstrain;
		  Ttangent = E0;
	  	}   
	  else if (fabsTstrain <= EpsiYT )
	 	{
		  Tstress = fp - CT + (BT/AT)*(pow((AT*AT - (EpsiYT- fabsTstrain)* (EpsiYT-fabsTstrain)),0.5));
		  Ttangent = BT*(EpsiYT - fabsTstrain)/ (AT* (pow(( AT*AT - (EpsiYT - fabsTstrain)* (EpsiYT-fabsTstrain)),0.5)));
	  }
	  else if (fabsTstrain <= EpsiT)
		  {
		  	Tstress = fy+(fabsTstrain-EpsiYT)*(1E-4)*E0;
			Ttangent = (1E-4)*E0;
		  }	  
	  else if (fabsTstrain <= EpsiU)
	 	 {
		  double fy1 = fy+(EpsiU-EpsiYT)*(1E-4)*E0;// modeified to add hardeding and avoid cinvergence problem
		  Tstress = fy1*(1- (fabsTstrain - EpsiT)/(EpsiU -EpsiT));
          //opserr<<"Error: Stiffness of SteelECthermal is negative"<<endln;
		  Ttangent = -fy1/(EpsiU-EpsiT);
	  	}
	  else 
		{ 
			 Tstress = 1E-10;
			 Ttangent = 1E-10;
			 //opserr<<"Error: Trial strain ("<< Tstrain<<") for SteelECthermal is invalid"<<endln;;
		}

      
	if(Tloading == 1)
	{
	  Tstress = Tstress;
	}
	else if(Tloading ==-1)
	{
		Tstress = -Tstress;
	}
	else
	{
		Tstress = 0;
	}
	//opserr<< "Trial strain: "<<Tstrain<< "   dStrain: "<<dStrain<<"   Tstress: "<<Tstress<<endln;
    Ctemp = Ttemp;

}

void SteelECThermal::detectLoadReversal (double dStrain)
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

double SteelECThermal::getStrain ()
{
   return Tstrain;
}

double SteelECThermal::getStress ()
{
   return Tstress;
}

double SteelECThermal::getTangent ()
{
   return Ttangent;
}

double 
SteelECThermal::getThermalElongation(void) 
{
  return ThermalElongation;
}

//Liming for updating the reduction factors////////////start
double 
SteelECThermal::getElongTangent(double TempT, double &ET, double &Elong, double TempTmax) 
{
	double FyRfactors[12];
	double FpRfactors[12];
	double E0Rfactors[12];
	//typeTag:3   EC3 Structural Steel
	//typeTag:21  EC2 Reinforcing Steel EC2 NHotRolled
	//typeTag:22  EC2 Reinforcing Steel EC2 NCold formed
	//typeTag:23  EC2 Reinforcing Steel EC2 X
	if(typeTag == 0||typeTag ==3) {
		double FyRfEC3[12] = {1.0, 1.0 ,1.0, 1.0 ,0.78, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0};
		double FpRfEC3[12] = {1.0, 0.807 ,0.613, 0.420 ,0.36, 0.18, 0.075, 0.050, 0.0375, 0.025 ,0.0125, 0.0};
		double E0RfEC3[12] = {1.0, 0.9, 0.8 ,0.7, 0.6 ,0.31, 0.13, 0.09, 0.0675, 0.045, 0.0225 , 0.0};
		for(int i=0;i<12;i++){
			FyRfactors[i] = FyRfEC3[i];
			FpRfactors[i] = FpRfEC3[i];
			E0Rfactors[i] = E0RfEC3[i];	
		}
		 
		}
	else if(typeTag == 21) {
		double FyRfEC21[12] = {1.0, 1.0 ,1.0, 1.0 ,0.78, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0};
		double FpRfEC21[12] = {1.0, 0.81 ,0.61, 0.42 ,0.36, 0.18, 0.07, 0.05, 0.04, 0.02 ,0.01, 0.0};
		double E0RfEC21[12] = {1.0, 0.9, 0.8 ,0.7, 0.6 ,0.31, 0.13, 0.09, 0.07, 0.04, 0.02 , 0.0};
		for(int i=0;i<12;i++){
			FyRfactors[i] = FyRfEC21[i];
			FpRfactors[i] = FpRfEC21[i];
			E0Rfactors[i] = E0RfEC21[i];	
		}
		
		}
	else if(typeTag == 22) {
		double FyRfEC22[12] = {1.0, 1.0 ,1.0, 0.94 ,0.67, 0.40, 0.12, 0.11, 0.08, 0.05 ,0.03, 0.0};
		double FpRfEC22[12] = {0.96 ,0.92, 0.81 ,0.63, 0.44, 0.26, 0.08, 0.06, 0.05 ,0.03, 0.02, 0.0};
		double E0RfEC22[12] = {1.0, 0.87, 0.72 ,0.56, 0.40 ,0.24, 0.08, 0.06, 0.05, 0.03, 0.02 , 0.0};
		for(int i=0;i<12;i++){
			FyRfactors[i] = FyRfEC22[i];
			FpRfactors[i] = FpRfEC22[i];
			E0Rfactors[i] = E0RfEC22[i];	
		}
		}
	else if(typeTag == 23) {
		double FyRfEC23[12] = {1.0, 1.0 ,1.0, 0.90 ,0.70, 0.47, 0.23, 0.11, 0.06, 0.04 ,0.02, 0.0};
		double FpRfEC23[12] = {1.00 ,0.87, 0.74 ,0.70, 0.51, 0.18, 0.07, 0.05, 0.04 ,0.02, 0.01, 0.0};
		double E0RfEC23[12] = {1.0, 0.95, 0.90 ,0.75, 0.60 ,0.31, 0.13, 0.09, 0.07, 0.04, 0.02 , 0.0};
		for(int i=0;i<12;i++){
			FyRfactors[i] = FyRfEC23[i];
			FpRfactors[i] = FpRfEC23[i];
			E0Rfactors[i] = E0RfEC23[i];	
		}
		}
	else 
		opserr<<"WARNING SteelECThermal received an invalid typeTag: "<<typeTag<<endln;

   //Now Updating modulus, strengths
	for( int i=0; i<13; i++) {
		if (TempT <= 80+100*i) 
		{
			if (i==0) {
				fy = fyT*(1.0 - TempT*(1.0-FyRfactors[0])/80);
				fp = fyT*(1.0 - TempT*(1.0-FpRfactors[0])/80);
				E0 = E0T*(1.0 - TempT*(1.0-E0Rfactors[0])/80);
			}
			else if (i==12) {
				opserr << "Warning:The temperature "<<TempT<<" for SteelECthermal is out of range\n"; 
				return -1;
			}
			else {
				fy = fyT*(FyRfactors[i-1] - (TempT+20-100*i)*(FyRfactors[i-1]-FyRfactors[i])/100);
				fp = fyT*(FpRfactors[i-1] - (TempT+20-100*i)*(FpRfactors[i-1]-FpRfactors[i])/100);
				E0 = E0T*(E0Rfactors[i-1] - (TempT+20-100*i)*(E0Rfactors[i-1]-E0Rfactors[i])/100);
			}
			break;
		}
	}
#ifdef _BDEBUG
	//opserr<<", TempT:"<<TempT<< " fy: "<< fy<< " fp: "<< fp <<" E0T:  "<< E0<<endln;
#endif

  // calculation of thermal elongation of reinforcing steel. JZ
///*	
	if (TempT <= 1) {
		  ThermalElongation = TempT * 1.2164e-5;
	  }
  else if (TempT <= 730) {
      ThermalElongation = -2.416e-4 + 1.2e-5 *(TempT+20) + 0.4e-8 *(TempT+20)*(TempT+20);
  }
  else if (TempT <= 840) {
      ThermalElongation = 11e-3;
  }
  else if (TempT <= 1180) {
      ThermalElongation = -6.2e-3 + 2e-5*(TempT+20);
  }
  else {
	  opserr << " SteelEC Temperature "<< TempT<<" is invalid\n";
	  return -1;
  }

  //ThermalElongation = 0 ;   //debug  Liming
  Elong = ThermalElongation;
  ET = E0;
  TemperautreC = TempT;
  return 0;
}


//END  for updating the reduction factors   by Liming Jiang


int SteelECThermal::commitState ()
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

int SteelECThermal::revertToLastCommit ()
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

int SteelECThermal::revertToStart ()
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

UniaxialMaterial* SteelECThermal::getCopy ()
{
   SteelECThermal* theCopy = new SteelECThermal(this->getTag(), typeTag,fy, E0,
				  a1, a2, a3, a4);

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

int SteelECThermal::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(16);
   data(0) = this->getTag();

   // Material properties
   data(1) = typeTag;
   data(2) = fy;
   data(3) = E0;
   data(4) = a1;
   data(5) = a2;
   data(6) = a3;
   data(7) = a4;

   // History variables from last converged state
   data(8) = CminStrain;
   data(9) = CmaxStrain;
   data(10) = CshiftP;
   data(11) = CshiftN;
   data(12) = Cloading;

   // State variables from last converged state
   data(13) = Cstrain;
   data(14) = Cstress;
   data(15) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "SteelECThermal::sendSelf() - failed to send data\n";

   return res;
}

int SteelECThermal::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(16);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "SteelECThermal::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      typeTag = data (1);
	  fy = data(2);
      E0 = data(3);
      a1 = data(4);
      a2 = data(5);
      a3 = data(6);
      a4 = data(7);

      // History variables from last converged state
      CminStrain = data(8);
      CmaxStrain = data(9);
      CshiftP = data(10);
      CshiftN = data(11);
      Cloading = int(data(12));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;

      // State variables from last converged state
      Cstrain = data(13);
      Cstress = data(14);
      Ctangent = data(15);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}

void SteelECThermal::Print (OPS_Stream& s, int flag)
{
   s << "SteelECThermal tag: " << this->getTag() << endln;
   s << "typeTag: " <<typeTag << " ";
   s << "  fy: " << fy << " ";
   s << "  E0: " << E0 << " ";
   s << "  a1: " << a1 << " ";
   s << "  a2: " << a2 << " ";
   s << "  a3: " << a3 << " ";
   s << "  a4: " << a4 << " ";
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
SteelECThermal::setParameter(const char **argv, int argc, Parameter &param)
{

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

  return -1;
}



int
SteelECThermal::updateParameter(int parameterID, Information &info)
{
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

	return 0;
}




int
SteelECThermal::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}



double
SteelECThermal::getStressSensitivity(int gradIndex, bool conditional)
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
SteelECThermal::getInitialTangentSensitivity(int gradIndex)
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
SteelECThermal::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
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
int SteelECThermal::setTrialStrain (double strain, double strainRate)
{
  opserr << "SteelECThermal::setTrialStrain (double strain, double strainRate) - should never be called\n";
  return 0;
}


int 
SteelECThermal::getVariable(const char *variable, Information &info)
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


