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


// Written: Q. Gu & Y.Peng
// Created: May 2012
//
// Description: This file contains the class definition for steel buckling restrained braces
// refer to: Alessandro Zona, Andrea Dall'Asta, 2012, "Elastoplastic model for 
// steel buckling restrained braces", Journal of Constructional steel research 68(2012), 118-125



#include <SteelBRB.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <string.h>
#include <stdlib.h>

#include <elementAPI.h>


void *OPS_SteelBRB(void) 
{ 
  
  int tag;
  double E, sigmaY0, sigmaY_T, alpha_T, beta_T, delta_T, sigmaY_C, alpha_C, beta_C, delta_C, Tol;
  
  Tol = 1.0e-14; 
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 11 && numArgs != 12) {    // ---- correct!!
    opserr << "Warning Insufficient args: unixialMaterial SteelBRB tag E sigmaY0 sigmaY_T alpha_T beta_T delta_T sigmaY_C alpha_C beta_C delta_C <Tol> \n";
    return 0;
  }

  int iData[1];
  double dData[11];

  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer values: nDMaterial SimplifiedJ2 \n";
    return 0;
  }  
  tag = iData[0]; 

  numData = numArgs -1;  
  if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid double values: nDMaterial SimplifiedJ2 " << tag << endln;
    return 0;
  }
  
  E        = dData[0];
  sigmaY0  = dData[1];
  sigmaY_T = dData[2];
  alpha_T  = dData[3];
  beta_T   = dData[4];
  delta_T  = dData[5];
  sigmaY_C = dData[6];
  alpha_C  = dData[7];
  beta_C   = dData[8];
  delta_C  = dData[9];
  if (numArgs == 12) {
    Tol =   dData[10];
  }

  // Parsing was successful, allocate the material
  UniaxialMaterial *theMaterial = new SteelBRB(tag, E,sigmaY0,sigmaY_T,alpha_T,alpha_C,sigmaY_C,beta_T,beta_C,delta_T,delta_C, Tol);      
  return theMaterial;
}


SteelBRB::SteelBRB(int pTag, double pE,double pSigmaY0, double pSigmaY_T,double pAlpha_T,
	double pAlpha_C,	double pSigmaY_C,double pBeta_T,double pBeta_C,double pDelta_T,	double pDelta_C, double pTol)
  :UniaxialMaterial(pTag,MAT_TAG_SteelBRB),
   strain(0.0), stress(0.0), tangent(0.0),	CStress(0.0), CPlastStrain(0.0), CCumPlastStrain(0.0),
   CYieldStress(0.0), CStrain(0.0)
{
	 E = pE;

	 sigmaY0 = pSigmaY0;  
	 sigmaY_T = pSigmaY_T;
	 alpha_T = pAlpha_T;
	 alpha_C = pAlpha_C;
	 sigmaY_C = pSigmaY_C;
	 beta_T = pBeta_T;
	 beta_C = pBeta_C;
	 delta_T = pDelta_T;
	 delta_C = pDelta_C;
	 Tol = pTol; 
	 
	 CDissipatedEnergy =0.0;
	 dissipatedEnergy =0.0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////




}


SteelBRB::~SteelBRB()
{

//	opserr<<"SteelBRB::~SteelBRB() is called!"<<endln; 

}

int 
SteelBRB::setTrialStrain(double pStrain, double pStrainRate)
{
  strain = pStrain;
  double strainInc = strain -CStrain;
  int N0=20;  // % iteration number
  double plastStrainInc =0.0;
  tangent =E;
 
  if (strainInc==0) {
           plastStrainInc = 0.0;
           plastStrain= CPlastStrain;
           stress = CStress;                       
           cumPlastStrain = CCumPlastStrain;
           yieldStress = CYieldStress;
		   tangent =E;  
		   dissipatedEnergy =CDissipatedEnergy;

  }

  else if (CStress* strainInc>=0) { //%loading 
       
      if (CStress >= 0) { 
           plastStrainInc = Newton_BRB(CStress, beta_T, CPlastStrain, sigmaY_T,CCumPlastStrain, delta_T, alpha_T, strainInc,0, Tol, N0);
           plastStrain= CPlastStrain+ plastStrainInc;
           stress = CStress+E*(strainInc-plastStrainInc);                       
           cumPlastStrain = CCumPlastStrain+fabs(plastStrainInc);
           yieldStress = sigmaY0+(sigmaY_T-sigmaY0)*(1.0-exp(-cumPlastStrain/delta_T));

 // ----------- tangent ---------
          double u = CStress +E*(strainInc-plastStrainInc)-beta_T*E*(CPlastStrain+plastStrainInc);
          double v = sigmaY0+(sigmaY_T-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T));
          double A = u/v;
          double C = alpha_T* pow(fabs(A),(alpha_T-2))*A*strainInc*E/v;
          int sign =-1;
          if (plastStrainInc>=0) 
              sign=1;
          
          double D = alpha_T*pow(fabs(A),(alpha_T-2))*A*strainInc*u/v/v*(sigmaY_T-sigmaY0)/delta_T*sign*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T); 
          double tmp = (C+pow(fabs(A),alpha_T))/(1+C*(1+beta_T)+D);
          tangent = E*(1-tmp);
		  dissipatedEnergy = CDissipatedEnergy+0.5*(stress+CStress-beta_T*E*(CPlastStrain+plastStrain))*plastStrainInc;

		  
	  }
	  else {     // %stress(i-1) < 0
           plastStrainInc = Newton_BRB(CStress, beta_C, CPlastStrain, sigmaY_C,CCumPlastStrain, delta_C, alpha_C, strainInc, 0, Tol, N0);  
           plastStrain= CPlastStrain+ plastStrainInc;
           stress = CStress+E*(strainInc-plastStrainInc);                       
           cumPlastStrain = CCumPlastStrain+fabs(plastStrainInc);
           yieldStress = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-cumPlastStrain/delta_C)); 
           
 // ----------- tangent ---------
          double u = CStress +E*(strainInc-plastStrainInc)-beta_C*E*(CPlastStrain+plastStrainInc);
          double v = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C));
          double A = u/v;
          double C = alpha_C* pow(fabs(A),(alpha_C-2))*A*strainInc*E/v;
          int sign =-1;
          if (plastStrainInc>=0) 
              sign=1;
          
          double D = alpha_C*pow(fabs(A),(alpha_C-2))*A*strainInc*u/v/v*(sigmaY_C-sigmaY0)/delta_C*sign*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C); 
          double tmp = (C+pow(fabs(A),alpha_C))/(1+C*(1+beta_C)+D);
          tangent = E*(1-tmp);
		  dissipatedEnergy = CDissipatedEnergy + 0.5*(stress+CStress-beta_C*E*(CPlastStrain+plastStrain))*plastStrainInc;

	  } // %stress(i-1) < 0
  }
        
  else { //  % (CStress* strainInc(i)<0) assume elastic unloading 
	  if (fabs(strainInc)<=fabs(CStress/E)) {  // % strainInc is not big enough, no special treatment needed.
               plastStrainInc = 0;
               plastStrain= CPlastStrain;
               stress = CStress+E*strainInc;                       
               cumPlastStrain = CCumPlastStrain;
               yieldStress = CYieldStress;
			   tangent =E;
		       dissipatedEnergy =CDissipatedEnergy;
	  }
	  else { // % strainInc is too big enough, need special treatment.
              double strain_unloading = -1.0*CStress/E;
              double strain_loading = strainInc-strain_unloading;
              if (CStress < 0) {
                   plastStrainInc = Newton_BRB(0, beta_T, CPlastStrain, sigmaY_T,CCumPlastStrain, delta_T, alpha_T, strain_loading,0, Tol, N0);
                   plastStrain= CPlastStrain+ plastStrainInc;
                   stress = E*(strain_loading-plastStrainInc);                       
                   cumPlastStrain = CCumPlastStrain+fabs(plastStrainInc);
                   yieldStress = sigmaY0+(sigmaY_T-sigmaY0)*(1-exp(-cumPlastStrain/delta_T));

 // ------------ tangent ---
                  double u = 0 +E*(strain_loading-plastStrainInc)-beta_T*E*(CPlastStrain+plastStrainInc);
                  double v = sigmaY0+(sigmaY_T-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T));
                  double A = u/v;
                  double C = alpha_T* pow(fabs(A),(alpha_T-2))*A*strain_loading*E/v;
                  double sign =-1.0;
                  if (plastStrainInc>=0) 
                      sign=1.0;
                  
                  double D = alpha_T*pow(fabs(A),(alpha_T-2))*A*strain_loading*u/v/v*(sigmaY_T-sigmaY0)/delta_T*sign*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T); 
                  double tmp = (C+pow(fabs(A),alpha_T))/(1+C*(1+beta_T)+D);
                  tangent = E*(1-tmp);
		          dissipatedEnergy = CDissipatedEnergy+ 0.5*(stress+0-beta_T*E*(CPlastStrain+plastStrain))*plastStrainInc;		           
			  }
              else { //%stress(i-1) > 0
                   plastStrainInc = Newton_BRB(0, beta_C, CPlastStrain, sigmaY_C,CCumPlastStrain, delta_C, alpha_C, strain_loading, 0, Tol, N0);  
                   plastStrain= CPlastStrain+ plastStrainInc;
                   stress = E*(strain_loading-plastStrainInc);                       
                   cumPlastStrain = CCumPlastStrain+fabs(plastStrainInc);
                   yieldStress = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-cumPlastStrain/delta_C));  

   // ------------ tangent ---
                  double u = 0 +E*(strain_loading-plastStrainInc)-beta_C*E*(CPlastStrain+plastStrainInc);
                  double v = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C));
                  double A = u/v;
                  double C = alpha_C* pow(fabs(A),(alpha_C-2))*A*strain_loading*E/v;
                  double sign =-1.0;
                  if (plastStrainInc>=0) 
                      sign=1.0;
                  
                  double D = alpha_C*pow(fabs(A),(alpha_C-2))*A*strain_loading*u/v/v*(sigmaY_C-sigmaY0)/delta_C*sign*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C); 
                  double tmp = (C+pow(fabs(A),alpha_C))/(1+C*(1+beta_C)+D);
                  tangent = E*(1-tmp);
		          dissipatedEnergy = CDissipatedEnergy + 0.5*(stress+0-beta_C*E*(CPlastStrain+plastStrain))*plastStrainInc;
              } // %stress(i-1) < 0 
	  } // % if (abs(strainInc(i))<=CStress/E)
              
  }





/*		opserr<<"stress is: " << stress<<", tangent is:"<<tangent<<endln;
 
  opserr<<"The fucked parameters are:   sigmaYC  alphaT  alphaC  deltaT  deltaC  " << sigmaY_C <<", "<< alpha_T <<", "<< alpha_C<<", "<<delta_T<<", "<< delta_C;
  opserr<<"   and stress is "<<stress<<endln; 
*/
  return 0;
}

double 
SteelBRB::getStress(void)
{
  return stress;
}

double 
SteelBRB::getTangent(void)
{
  return tangent;
}

double 
SteelBRB::getInitialTangent(void)
{
  // return the initial tangent
  return E;
}

double 
SteelBRB::getStrain(void)
{
  return strain;
}

int 
SteelBRB::commitState(void)
{
	CStress =stress;
	CPlastStrain=plastStrain;
	CCumPlastStrain=cumPlastStrain;
	CYieldStress=yieldStress;
	CStrain=strain;
	CDissipatedEnergy = dissipatedEnergy; 
  
	return 0;
}

int 
SteelBRB::revertToLastCommit(void)
{
  return 0;
}

int 
SteelBRB::revertToStart(void)
{

	 CStress = 0.0;
	 CPlastStrain= 0.0;
	 CCumPlastStrain= 0.0;
	 CYieldStress= 0.0;
	 CStrain= 0.0;


	 stress= 0.0;
	 plastStrain= 0.0;
	 cumPlastStrain= 0.0;
	 yieldStress= 0.0;
	 strain= 0.0;


	 CDissipatedEnergy= 0.0;
	 dissipatedEnergy= 0.0;

     parameterID= 0.0;

	 if (SHVs !=0) 
	    SHVs->Zero();
	 

     return 0;
}

UniaxialMaterial *
SteelBRB::getCopy(void)
{
  SteelBRB *theCopy = new SteelBRB(this->getTag(),E,sigmaY0, sigmaY_T,alpha_T,
	alpha_C,sigmaY_C,beta_T,beta_C,delta_T,	delta_C, Tol );

  theCopy->strain = 0.0;
  theCopy->stress = 0.0;
  theCopy->tangent = 0.0;
  theCopy->CStress = 0.0;
  theCopy->CPlastStrain = 0.0;
  theCopy->CCumPlastStrain = 0.0;
  theCopy->CYieldStress = 0.0;
  theCopy->CStrain = 0.0;
  

  return theCopy;
}

int 
SteelBRB::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
SteelBRB::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
SteelBRB::Print(OPS_Stream &s, int flag)
{
  return;
}




double SteelBRB::PlastStrainIncRes(double CStress, double beta, double CPlastStrain, double sigmaY,
								   double cumPlastStrain, double delta, double alpha, double strainInc, double plastStrainInc){


	double temp1 = CStress +E*(strainInc-plastStrainInc)-beta*E*(CPlastStrain+plastStrainInc);
	double temp2  = sigmaY0+(sigmaY-sigmaY0)*(1.0-exp(-(cumPlastStrain+fabs(plastStrainInc))/delta));

	double f = plastStrainInc-pow(fabs(temp1/temp2),alpha)*strainInc;

return f ;

};

double SteelBRB::PlastStrainIncResDev(double CStress, double beta, double CPlastStrain, double sigmaY, 
									  double cumPlastStrain, double delta, double alpha, double strainInc, double plastStrainInc ){

	double temp1 = CStress +E*(strainInc-plastStrainInc)-beta*E*(CPlastStrain+plastStrainInc);
	double temp2 = sigmaY0+(sigmaY-sigmaY0)*(1.0-exp(-(cumPlastStrain+fabs(plastStrainInc))/delta));

//	%temp3 = (sigmaY-sigmaY0)*(plastStrainInc/(delta*abs(plastStrainInc)))*exp(-(cumPlastStrain+abs(plastStrainInc))/delta)
	int sign =-1;
	if (plastStrainInc>=0) 
		sign=1;

    
	double temp3 = (sigmaY-sigmaY0)/delta*sign*exp(-(cumPlastStrain+fabs(plastStrainInc))/delta);
	double f = 1.0-strainInc*alpha*pow((fabs(temp1/temp2)),(alpha-2))*(temp1/temp2)*(((-E-beta*E)*temp2-temp1*temp3)/temp2/temp2); 
	return f;

};



double  SteelBRB::Newton_BRB(double CStress, double beta, double CPlastStrain, double sigmaY, double cumPlastStrain,
						  double delta, double alpha, double strainInc, double x0, double Tol, int N0){

	int i =1;
	double lowerbound =0.0;
	double upperbound = 0.0;
 
	if (fabs(strainInc)<1.0e-16){ 
		return 0.0; 
	}
	else if (strainInc>0) {
		lowerbound = 0.0;
		upperbound = strainInc;
	}
	else{
		lowerbound = strainInc;
		upperbound = 0.0;
	}


	double  F_low = PlastStrainIncRes(CStress, beta, CPlastStrain, sigmaY,cumPlastStrain, delta, alpha, strainInc,lowerbound);
	double  F_upp = PlastStrainIncRes(CStress, beta, CPlastStrain, sigmaY,cumPlastStrain, delta, alpha, strainInc,upperbound);
  
	if (F_low*F_upp>0)
		  opserr<< "In SteelBRB::Newton_BRB, lower bound and upper bound have the same sign!\n";
	 
  
	double F = PlastStrainIncRes(CStress, beta, CPlastStrain, sigmaY,cumPlastStrain, delta, alpha, strainInc,x0);
	double F_dev = 0.0;
	double x1=x0;
	
   
	while ((i<=N0)&&(fabs(F)>Tol)) {

		F_dev = PlastStrainIncResDev(CStress, beta, CPlastStrain, sigmaY,cumPlastStrain, delta, alpha, strainInc,x0);
		x1 = x0-F/F_dev;
		if (( x1<lowerbound) ||(x1>upperbound))
			 x1 = (lowerbound+upperbound)/2;
		
      
		F = PlastStrainIncRes(CStress, beta, CPlastStrain, sigmaY,cumPlastStrain, delta, alpha, strainInc,x1);
    
		if (F_low*F<0){
			upperbound = x1;
			F_upp =F;
		}
		else if (F_upp*F<0){
			lowerbound = x1;
			F_low =F;       
		}
		
      
		i=i+1;
    
  
		x0=x1;

	}


	if (fabs(F)>Tol){
		opserr<< "Fatal error: SteelBRB::Newton_BRB does not converge ===============\n";
	    exit(-1); 
	}

	return x0; 
	

};


// AddingSensitivity:BEGIN ///////////////////////////////////


/*
	double E;
	double sigmaY0; 
	
	double sigmaY_T;
	double alpha_T;
	double beta_T;
    double delta_T;
	
	double sigmaY_C;	
	double alpha_C;
	double beta_C;
	double delta_C;
*/
int
SteelBRB::setParameter(const char **argv, int argc, Information &info)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"E") == 0) {
		info.theType = DoubleType;
		return 1;
	}
	if (strcmp(argv[0],"sigmaY0") == 0) {
		info.theType = DoubleType;
		return 2;
	}
	if (strcmp(argv[0],"sigmaY_T") == 0) {
		info.theType = DoubleType;
		return 3;
	}
	if (strcmp(argv[0],"alpha_T") == 0) {
		info.theType = DoubleType;
		return 4;
	}
	if (strcmp(argv[0],"beta_T") == 0) {
		info.theType = DoubleType;
		return 5;
	}
	if (strcmp(argv[0],"delta_T") == 0) {
		info.theType = DoubleType;
		return 6;
	}
	if (strcmp(argv[0],"sigmaY_C") == 0) {
		info.theType = DoubleType;
		return 7;
	}
	if (strcmp(argv[0],"alpha_C") == 0) {
		info.theType = DoubleType;
		return 8;
	}
	if (strcmp(argv[0],"beta_C") == 0) {
		info.theType = DoubleType;
		return 9;
	}
	if (strcmp(argv[0],"delta_C") == 0) {
		info.theType = DoubleType;
		return 10;
	}
	
	else
		opserr << "WARNING: Could not set parameter in SteelBRB. " << endln;
               
	return -1;
}



int
SteelBRB::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->E = info.theDouble;
		break;
	case 2:
		this->sigmaY0 = info.theDouble;
		break;
	case 3:
		this->sigmaY_T = info.theDouble;
		break;
	case 4:
		this->alpha_T = info.theDouble;
		break;
	case 5:
		this->beta_T = info.theDouble;
		break;
	case 6:
		this->delta_T = info.theDouble;
		break;
	case 7:
		this->sigmaY_C = info.theDouble;
		break;
	case 8:
		this->alpha_C = info.theDouble;
		break;
	case 9:
		this->beta_C = info.theDouble;
		break;
	case 10:
		this->delta_C = info.theDouble;
		break;

	default:
		return -1;
	}

//	Ttangent = E0;          // Initial stiffness
	revertToStart();

	return 0;
}




int
SteelBRB::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}


double SteelBRB::getStrainSensitivity(int gradNumber){
	if (SHVs ==0) {
			opserr<<"warning:SteelBRB::getStrainsSensitivity, SHVs =0 "<<endln;	
			return 0.0;
		} 
	double  Tsensitivity =(*SHVs)(0,(gradNumber-1));
	return Tsensitivity;

}

double
SteelBRB::getStressSensitivity(int gradNumber, bool conditional){

//* ----------debug only must change back ------
int iii2 =0;    
    iii2++;

 /*  if (iii2 ==1){
	   strain =1.4696588314348e-006;
       this->setTrialStrain(strain, 0.0);
   }
*/

/*
   if (iii2 ==2) {
        strain = 2.9437234801682e-006;
       this->setTrialStrain(strain, 0.0);
   }

//-----------*/
	   


	// Pick up sensitivity history variables
	double CStrainSensitivity = 0.0;
	double CStressSensitivity = 0.0;
	double CPlastStrainSensitivity = 0.0;
	double CCumPlastStrainSensitivity = 0.0;
	double CYieldStressSensitivity = 0.0;


	double strainSensitivity = 0.0;
	double stressSensitivity = 0.0;
	double plastStrainSensitivity = 0.0;
	double cumPlastStrainSensitivity = 0.0;
	double yieldStressSensitivity = 0.0;
 


	if (SHVs != 0) {
		CStrainSensitivity = (*SHVs)(0,(gradNumber-1));
		CStressSensitivity = (*SHVs)(1,(gradNumber-1));
		CPlastStrainSensitivity = (*SHVs)(2,(gradNumber-1));
		CCumPlastStrainSensitivity = (*SHVs)(3,(gradNumber-1));
		CYieldStressSensitivity = (*SHVs)(4,(gradNumber-1));
 

	}

//	double strainSensitivity = TstrainSensitivity; 
	double strainIncSensitivity = strainSensitivity - CStrainSensitivity; 

	// Assign values to parameter derivatives (depending on what's random)
	double ESensitivity = 0.0;
	double sigmaY0Sensitivity = 0.0;
	double sigmaY_TSensitivity = 0.0;
	double alpha_TSensitivity = 0.0;
	double beta_TSensitivity = 0.0;
    double delta_TSensitivity = 0.0;
	double sigmaY_CSensitivity = 0.0;	
	double alpha_CSensitivity = 0.0;
	double beta_CSensitivity = 0.0;
	double delta_CSensitivity = 0.0;




	if (parameterID == 1) {
		ESensitivity = 1.0;
	}
	else if (parameterID == 2) {
		sigmaY0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		sigmaY_TSensitivity = 1.0;
	}
	else if (parameterID == 4) {
		alpha_TSensitivity = 1.0;
	}
	else if (parameterID == 5) {
		beta_TSensitivity = 1.0;
	}
	else if (parameterID == 6) {
		delta_TSensitivity = 1.0;
	}
	else if (parameterID == 7) {
		sigmaY_CSensitivity = 1.0;
	}
	else if (parameterID == 8) {
		alpha_CSensitivity = 1.0;
	}
	else if (parameterID == 9) {
		beta_CSensitivity = 1.0;
	}
	else if (parameterID == 10) {
		delta_CSensitivity = 1.0;
	}

 double plastStrainInc =0.0;

 double plastStrainIncSensitivity = 0.0;

  
  // determine trial stress and tangent
 double strainInc = strain-CStrain;

  if (strainInc==0) {
           plastStrainInc = 0.0;
           plastStrain= CPlastStrain;
           stress = CStress;                       
           cumPlastStrain = CCumPlastStrain;
           yieldStress = CYieldStress;

		   plastStrainIncSensitivity = 0.0;
		   plastStrainSensitivity = CPlastStrainSensitivity;
//		   stressSensitivity = CStressSensitivity;                      
           stressSensitivity = CStressSensitivity+ESensitivity*strainInc+E*strainIncSensitivity; 
		   cumPlastStrainSensitivity = CCumPlastStrainSensitivity;
		   yieldStressSensitivity = CYieldStressSensitivity;


		 
	//	   tangent =E;

  }

  else if (CStress* strainInc>=0) { //%loading 

  
      if (CStress >= 0) { 
           plastStrainInc = plastStrain - CPlastStrain;  //Newton_BRB(CStress, beta_T, CPlastStrain, sigmaY_T,CCumPlastStrain, delta_T, alpha_T, strainInc,0, Tol, N0);
		   int sign = -1.0;
		   if (plastStrainInc>0)
			   sign =1.0;

           double u = CStress +E*(strainInc-plastStrainInc)-beta_T*E*(CPlastStrain+plastStrainInc);
           double v = sigmaY0+(sigmaY_T-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T));
           double A = u/v;

	       double C = CStressSensitivity+(strainInc-plastStrainInc-beta_T*CPlastStrain-beta_T*plastStrainInc)*ESensitivity +E*strainIncSensitivity
			     -E*(CPlastStrain+plastStrainInc)*beta_TSensitivity-beta_T*E*CPlastStrainSensitivity;
           double D = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T))
			     -(sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*
				 (delta_TSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_T/delta_T-CCumPlastStrainSensitivity/delta_T);

           double F = (sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*sign/delta_T;
			    
		   double G = alpha_T*pow(fabs(A),(alpha_T-2))*A*strainInc*(C/v-u*D/v/v)+strainIncSensitivity*pow(fabs(A),alpha_T)
			             +pow(fabs(A),alpha_T)*log(fabs(A))*strainInc*alpha_TSensitivity;
           double H = alpha_T*pow(fabs(A),(alpha_T-2))*A*strainInc*((E+E*beta_T)/v+u*F/v/v);
			  
           plastStrainIncSensitivity = G/(1+H);


           plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
	       stressSensitivity = CStressSensitivity+ESensitivity*(strainInc-plastStrainInc)+E*(strainIncSensitivity-plastStrainIncSensitivity);                      
	       cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
	       yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_T))
			      +(sigmaY_T-sigmaY0)*exp(-cumPlastStrain/delta_T)*(cumPlastStrainSensitivity*delta_T-cumPlastStrain*delta_TSensitivity)/delta_T/delta_T;



	  }
	  else {     // %stress(i-1) < 0
           plastStrainInc = plastStrain - CPlastStrain;  //Newton_BRB(CStress, beta_C, CPlastStrain, sigmaY_C,CCumPlastStrain, delta_C, alpha_C, strainInc,0, Tol, N0);
		   int sign = -1.0;
		   if (plastStrainInc>0)
			   sign =1.0;
           
           double u = CStress +E*(strainInc-plastStrainInc)-beta_C*E*(CPlastStrain+plastStrainInc);
           double v = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C));
           double A = u/v;

	       double C = CStressSensitivity+(strainInc-plastStrainInc-beta_C*CPlastStrain-beta_C*plastStrainInc)*ESensitivity +E*strainIncSensitivity
			     -E*(CPlastStrain+plastStrainInc)*beta_CSensitivity-beta_C*E*CPlastStrainSensitivity;
           double D = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C))
			     -(sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*
				 (delta_CSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_C/delta_C-CCumPlastStrainSensitivity/delta_C);

           double F = (sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*sign/delta_C;
			    
		   double G = alpha_C*pow(fabs(A),(alpha_C-2))*A*strainInc*(C/v-u*D/v/v)+strainIncSensitivity*pow(fabs(A),alpha_C)
			             +pow(fabs(A),alpha_C)*log(fabs(A))*strainInc*alpha_CSensitivity;
           double H = alpha_C*pow(fabs(A),(alpha_C-2))*A*strainInc*((E+E*beta_C)/v+u*F/v/v);
			  
           plastStrainIncSensitivity = G/(1+H);


           plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
	       stressSensitivity = CStressSensitivity+ESensitivity*(strainInc-plastStrainInc)+E*(strainIncSensitivity-plastStrainIncSensitivity);                      
	       cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
	       yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_C))
			      +(sigmaY_C-sigmaY0)*exp(-cumPlastStrain/delta_C)*(cumPlastStrainSensitivity*delta_C-cumPlastStrain*delta_CSensitivity)/delta_C/delta_C;
	  } // %stress(i-1) < 0
  }
      
  else { //  % (CStress* strainInc(i)<0) assume elastic unloading 
	  if (fabs(strainInc)<=fabs(CStress/E)) {  // % strainInc is not big enough, no special treatment needed.
               plastStrainIncSensitivity = 0;
               plastStrainSensitivity= CPlastStrainSensitivity;
               stressSensitivity = CStressSensitivity+ESensitivity*strainInc+E*strainIncSensitivity; 
          
               cumPlastStrainSensitivity = CCumPlastStrainSensitivity;
               yieldStressSensitivity = CYieldStressSensitivity;
			
	  }
	  else { // % strainInc is too big enough, need special treatment.
		      double strain_unloading = -1.0*CStress/E;
              double strain_loading = strainInc-strain_unloading;

              double strain_unloadingSensitivity =-(E*CStressSensitivity-ESensitivity*CStress)/E/E;
              double strain_loadingSensitivity = strainIncSensitivity-strain_unloadingSensitivity;
             
			  if (CStress < 0) {
                  plastStrainInc = plastStrain - CPlastStrain;//Newton_BRB(CStress, beta_T, CPlastStrain, sigmaY_T,CCumPlastStrain, delta_T, alpha_T, strainInc,0, Tol, N0);

				 int sign = -1.0;
				 if (plastStrainInc>0)
					 sign =1.0;
				 double u = 0.0 +E*(strain_loading-plastStrainInc)-beta_T*E*(CPlastStrain+plastStrainInc);
				 double v = sigmaY0+(sigmaY_T-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T));
				 double A = u/v;

				 double C = 0.0+(strain_loading-plastStrainInc-beta_T*CPlastStrain-beta_T*plastStrainInc)*ESensitivity +E*strain_loadingSensitivity
					  -E*(CPlastStrain+plastStrainInc)*beta_TSensitivity-beta_T*E*CPlastStrainSensitivity;
				 double D = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T))
					  -(sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*
						 (delta_TSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_T/delta_T-CCumPlastStrainSensitivity/delta_T);

				 double F = (sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*sign/delta_T;
			    
			 	 double G = alpha_T*pow(fabs(A),(alpha_T-2))*A*strain_loading*(C/v-u*D/v/v)+strain_loadingSensitivity*pow(fabs(A),alpha_T)
			             +pow(fabs(A),alpha_T)*log(fabs(A))*strain_loading*alpha_TSensitivity;
				 double H = alpha_T*pow(fabs(A),(alpha_T-2))*A*strain_loading*((E+E*beta_T)/v+u*F/v/v);
			  
				 plastStrainIncSensitivity = G/(1+H);

				 plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
				 stressSensitivity = 0.0+ESensitivity*( strain_loading-plastStrainInc)+E*( strain_loadingSensitivity-plastStrainIncSensitivity);                      
				 cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
				 yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_T))
			      +(sigmaY_T-sigmaY0)*exp(-cumPlastStrain/delta_T)*(cumPlastStrainSensitivity*delta_T-cumPlastStrain*delta_TSensitivity)/delta_T/delta_T;
			  }
              else { //%stress(i-1) > 0
                     plastStrainInc = plastStrain - CPlastStrain;//Newton_BRB(CStress, beta_C, CPlastStrain, sigmaY_C,CCumPlastStrain, delta_C, alpha_C, strainInc,0, Tol, N0);
					 int sign = -1.0;
				     if (plastStrainInc>0)
			           sign =1.0;
				     double u = E*(strain_loading-plastStrainInc)-beta_C*E*(CPlastStrain+plastStrainInc);
				     double v = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C));
				     double A = u/v;

					 double C = 0.0+(strain_loading-plastStrainInc-beta_C*CPlastStrain-beta_C*plastStrainInc)*ESensitivity +E*strain_loadingSensitivity
					      -E*(CPlastStrain+plastStrainInc)*beta_CSensitivity-beta_C*E*CPlastStrainSensitivity;
					 double D = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C))
						  -(sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*
							 (delta_CSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_C/delta_C-CCumPlastStrainSensitivity/delta_C);

					 double F = (sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*sign/delta_C;
					
			 		 double G = alpha_C*pow(fabs(A),(alpha_C-2))*A*strain_loading*(C/v-u*D/v/v)+strain_loadingSensitivity*pow(fabs(A),alpha_C)
			             +pow(fabs(A),alpha_C)*log(fabs(A))*strain_loading*alpha_CSensitivity;
					 double H = alpha_C*pow(fabs(A),(alpha_C-2))*A*strain_loading*((E+E*beta_C)/v+u*F/v/v);
				  
					 plastStrainIncSensitivity = G/(1+H);

                     plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
					 stressSensitivity = 0.0+ESensitivity*(strain_loading-plastStrainInc)+E*(strain_loadingSensitivity-plastStrainIncSensitivity);                      
					 cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
	                 yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_C))
			           +(sigmaY_C-sigmaY0)*exp(-cumPlastStrain/delta_C)*(cumPlastStrainSensitivity*delta_C-cumPlastStrain*delta_CSensitivity)/delta_C/delta_C;
              } // %stress(i-1) < 0 
	  } // % if (abs(strainInc(i))<=CStress/E)
              
  }
	
	//*opserr.precision(16);
   /*
   *debug1 << "getStressSensitivity PlastStrainInc:"<< plastStrainInc << ",  plastStrainIncSensitivity: "<<plastStrainIncSensitivity<<"\n"; 
   *debug1 << "getStressSensitivity plastStrain:"<< plastStrain << ",  plastStrainSensitivity: "<<plastStrainSensitivity<<"\n"; 
   *debug1 << "getStressSensitivity stress:"<< stress << ",  stressSensitivity: "<<stressSensitivity<<"\n";    
   *debug1 << "getStressSensitivity cumPlastStrain:"<< cumPlastStrain << ",  cumPlastStrainSensitivity: "<<cumPlastStrainSensitivity<<"\n"; 
   *debug1 << "getStressSensitivity yieldStress:"<< yieldStress << ",  yieldStressSensitivity: "<<yieldStressSensitivity<<"\n"; 
 //*/


//	opserr<<"BRB--getStressSensitivity:" <<stressSensitivity<<endln; 
 int ii = 0;    
	 ii++;

	if (fabs(stressSensitivity) > 1.0e10) {
		opserr<< "error, ii = "<<ii<< endln; 
	}

	return stressSensitivity;

}




int
SteelBRB::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(6,numGrads);
		SHVs->Zero();
	}


	// Pick up sensitivity history variables
	double CStrainSensitivity = 0.0;
	double CStressSensitivity = 0.0;
	double CPlastStrainSensitivity = 0.0;
	double CCumPlastStrainSensitivity = 0.0;
	double CYieldStressSensitivity = 0.0;
	double CDissipatedEnergySensitivity = 0.0;

	double strainSensitivity = TstrainSensitivity;
	double stressSensitivity = 0.0;
	double plastStrainSensitivity = 0.0;
	double cumPlastStrainSensitivity = 0.0;
	double yieldStressSensitivity = 0.0;
	double dissipatedEnergySensitivity = 0.0; 


	if (SHVs != 0) {
		CStrainSensitivity = (*SHVs)(0,(gradNumber-1));
		CStressSensitivity = (*SHVs)(1,(gradNumber-1));
		CPlastStrainSensitivity = (*SHVs)(2,(gradNumber-1));
		CCumPlastStrainSensitivity = (*SHVs)(3,(gradNumber-1));
		CYieldStressSensitivity = (*SHVs)(4,(gradNumber-1));
		CDissipatedEnergySensitivity = (*SHVs)(5,(gradNumber-1));

	}

//	double strainSensitivity = TstrainSensitivity; 
	double strainIncSensitivity = strainSensitivity - CStrainSensitivity; 


	// Assign values to parameter derivatives (depending on what's random)
	double ESensitivity = 0.0;
	double sigmaY0Sensitivity = 0.0;
	double sigmaY_TSensitivity = 0.0;
	double alpha_TSensitivity = 0.0;
	double beta_TSensitivity = 0.0;
    double delta_TSensitivity = 0.0;
	double sigmaY_CSensitivity = 0.0;	
	double alpha_CSensitivity = 0.0;
	double beta_CSensitivity = 0.0;
	double delta_CSensitivity = 0.0;




	if (parameterID == 1) {
		ESensitivity = 1.0;
	}
	else if (parameterID == 2) {
		sigmaY0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		sigmaY_TSensitivity = 1.0;
	}
	else if (parameterID == 4) {
		alpha_TSensitivity = 1.0;
	}
	else if (parameterID == 5) {
		beta_TSensitivity = 1.0;
	}
	else if (parameterID == 6) {
		delta_TSensitivity = 1.0;
	}
	else if (parameterID == 7) {
		sigmaY_CSensitivity = 1.0;
	}
	else if (parameterID == 8) {
		alpha_CSensitivity = 1.0;
	}
	else if (parameterID == 9) {
		beta_CSensitivity = 1.0;
	}
	else if (parameterID == 10) {
		delta_CSensitivity = 1.0;
	}

 double plastStrainInc =0.0;

 double plastStrainIncSensitivity = 0.0;

  
  // determine trial stress and tangent
 double strainInc = strain-CStrain;


 // ======== for debug only =======

//        strainInc 
 // ===============================

  if (strainInc==0) {
           plastStrainInc = 0.0;
           plastStrain= CPlastStrain;
           stress = CStress;                       
           cumPlastStrain = CCumPlastStrain;
           yieldStress = CYieldStress; //sigmaY0+(sigmaY_T-sigmaY0)*(1.0-exp(-cumPlastStrain/delta_T));

		   plastStrainIncSensitivity = 0.0;
		   plastStrainSensitivity = CPlastStrainSensitivity;
		 //stressSensitivity = CStressSensitivity;  
           stressSensitivity = CStressSensitivity+ESensitivity*strainInc+E*strainIncSensitivity; 		   
		   cumPlastStrainSensitivity = CCumPlastStrainSensitivity;
		   yieldStressSensitivity = CYieldStressSensitivity;
           dissipatedEnergySensitivity = CDissipatedEnergySensitivity; 


	//	   tangent =E;

  }

  else if (CStress* strainInc>=0) { //%loading 

  
      if (CStress >= 0) { 
           plastStrainInc = plastStrain - CPlastStrain;  //Newton_BRB(CStress, beta_T, CPlastStrain, sigmaY_T,CCumPlastStrain, delta_T, alpha_T, strainInc,0, Tol, N0);
		   int sign = -1.0;
		   if (plastStrainInc>0)
			   sign =1.0;

           double u = CStress +E*(strainInc-plastStrainInc)-beta_T*E*(CPlastStrain+plastStrainInc);
           double v = sigmaY0+(sigmaY_T-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T));
           double A = u/v;

	       double C = CStressSensitivity+(strainInc-plastStrainInc-beta_T*CPlastStrain-beta_T*plastStrainInc)*ESensitivity +E*strainIncSensitivity
			     -E*(CPlastStrain+plastStrainInc)*beta_TSensitivity-beta_T*E*CPlastStrainSensitivity;
           double D = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T))
			     -(sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*
				 (delta_TSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_T/delta_T-CCumPlastStrainSensitivity/delta_T);

           double F = (sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*sign/delta_T;
			    
		   double G = alpha_T*pow(fabs(A),(alpha_T-2))*A*strainInc*(C/v-u*D/v/v)+strainIncSensitivity*pow(fabs(A),alpha_T)
			             +pow(fabs(A),alpha_T)*log(fabs(A))*strainInc*alpha_TSensitivity;
           double H = alpha_T*pow(fabs(A),(alpha_T-2))*A*strainInc*((E+E*beta_T)/v+u*F/v/v);
			  
           plastStrainIncSensitivity = G/(1+H);


           plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
	       stressSensitivity = CStressSensitivity+ESensitivity*(strainInc-plastStrainInc)+E*(strainIncSensitivity-plastStrainIncSensitivity);                      
	       cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
	       yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_T))
			      +(sigmaY_T-sigmaY0)*exp(-cumPlastStrain/delta_T)*(cumPlastStrainSensitivity*delta_T-cumPlastStrain*delta_TSensitivity)/delta_T/delta_T;

		   dissipatedEnergySensitivity = CDissipatedEnergySensitivity + 0.5*(stress+CStress-beta_T*E*(CPlastStrain+plastStrain))*plastStrainIncSensitivity
			      + 0.5*(stressSensitivity+CStressSensitivity-(beta_TSensitivity*E+beta_T*ESensitivity)*(CPlastStrain+plastStrain)
				  -beta_T*E*(CPlastStrainSensitivity+plastStrainSensitivity))*plastStrainInc; 
//		  dissipatedEnergy += 0.5*(stress+CStress-beta_T*E*(CPlastStrain+plastStrain))*plastStrainInc;

	  }
	  else {     // %stress(i-1) < 0
           plastStrainInc = plastStrain - CPlastStrain;  //Newton_BRB(CStress, beta_C, CPlastStrain, sigmaY_C,CCumPlastStrain, delta_C, alpha_C, strainInc,0, Tol, N0);
		   int sign = -1.0;
		   if (plastStrainInc>0)
			   sign =1.0;
           
           double u = CStress +E*(strainInc-plastStrainInc)-beta_C*E*(CPlastStrain+plastStrainInc);
           double v = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C));
           double A = u/v;

	       double C = CStressSensitivity+(strainInc-plastStrainInc-beta_C*CPlastStrain-beta_C*plastStrainInc)*ESensitivity +E*strainIncSensitivity
			     -E*(CPlastStrain+plastStrainInc)*beta_CSensitivity-beta_C*E*CPlastStrainSensitivity;
           double D = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C))
			     -(sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*
				 (delta_CSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_C/delta_C-CCumPlastStrainSensitivity/delta_C);

           double F = (sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*sign/delta_C;
			    
		   double G = alpha_C*pow(fabs(A),(alpha_C-2))*A*strainInc*(C/v-u*D/v/v)+strainIncSensitivity*pow(fabs(A),alpha_C)
			             +pow(fabs(A),alpha_C)*log(fabs(A))*strainInc*alpha_CSensitivity;

           double H = alpha_C*pow(fabs(A),(alpha_C-2))*A*strainInc*((E+E*beta_C)/v+u*F/v/v);
			  
           plastStrainIncSensitivity = G/(1+H);


           plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
	       stressSensitivity = CStressSensitivity+ESensitivity*(strainInc-plastStrainInc)+E*(strainIncSensitivity-plastStrainIncSensitivity);                      
	       cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
	       yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_C))
			      +(sigmaY_C-sigmaY0)*exp(-cumPlastStrain/delta_C)*(cumPlastStrainSensitivity*delta_C-cumPlastStrain*delta_CSensitivity)/delta_C/delta_C;
		   dissipatedEnergySensitivity = CDissipatedEnergySensitivity + 0.5*(stress+CStress-beta_C*E*(CPlastStrain+plastStrain))*plastStrainIncSensitivity
			      + 0.5*(stressSensitivity+CStressSensitivity-(beta_CSensitivity*E+beta_C*ESensitivity)*(CPlastStrain+plastStrain)
				  -beta_C*E*(CPlastStrainSensitivity+plastStrainSensitivity))*plastStrainInc; 
	  } // %stress(i-1) < 0
  }
      
  else { //  % (CStress* strainInc(i)<0) assume elastic unloading 
	  if (fabs(strainInc)<=fabs(CStress/E)) {  // % strainInc is not big enough, no special treatment needed.
               plastStrainIncSensitivity = 0;
               plastStrainSensitivity= CPlastStrainSensitivity;
               stressSensitivity = CStressSensitivity+ESensitivity*strainInc+E*strainIncSensitivity; 
          
               cumPlastStrainSensitivity = CCumPlastStrainSensitivity;
               yieldStressSensitivity = CYieldStressSensitivity;
			   dissipatedEnergySensitivity = CDissipatedEnergySensitivity; 
			
	  }
	  else { // % strainInc is too big enough, need special treatment.
		      double strain_unloading = -1.0*CStress/E;
              double strain_loading = strainInc-strain_unloading;

              double strain_unloadingSensitivity =-(E*CStressSensitivity-ESensitivity*CStress)/E/E;
              double strain_loadingSensitivity = strainIncSensitivity-strain_unloadingSensitivity;
             
			  if (CStress < 0) {
                  plastStrainInc = plastStrain - CPlastStrain;//Newton_BRB(CStress, beta_T, CPlastStrain, sigmaY_T,CCumPlastStrain, delta_T, alpha_T, strainInc,0, Tol, N0);

				 int sign = -1.0;
				 if (plastStrainInc>0)
					 sign =1.0;
				 double u = 0.0 +E*(strain_loading-plastStrainInc)-beta_T*E*(CPlastStrain+plastStrainInc);
				 double v = sigmaY0+(sigmaY_T-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T));
				 double A = u/v;

				 double C = 0.0+(strain_loading-plastStrainInc-beta_T*CPlastStrain-beta_T*plastStrainInc)*ESensitivity +E*strain_loadingSensitivity
					  -E*(CPlastStrain+plastStrainInc)*beta_TSensitivity-beta_T*E*CPlastStrainSensitivity;
				 double D = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T))
					  -(sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*
						 (delta_TSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_T/delta_T-CCumPlastStrainSensitivity/delta_T);

				 double F = (sigmaY_T-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_T)*sign/delta_T;
			    
			 	 double G = alpha_T*pow(fabs(A),(alpha_T-2))*A*strain_loading*(C/v-u*D/v/v)+strain_loadingSensitivity*pow(fabs(A),alpha_T)
			             +pow(fabs(A),alpha_T)*log(fabs(A))*strain_loading*alpha_TSensitivity;
				 double H = alpha_T*pow(fabs(A),(alpha_T-2))*A*strain_loading*((E+E*beta_T)/v+u*F/v/v);
			  
				 plastStrainIncSensitivity = G/(1+H);

				 plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
				 stressSensitivity = 0.0+ESensitivity*( strain_loading-plastStrainInc)+E*(strain_loadingSensitivity-plastStrainIncSensitivity);                      
				 cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
				 yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_TSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_T))
			      +(sigmaY_T-sigmaY0)*exp(-cumPlastStrain/delta_T)*(cumPlastStrainSensitivity*delta_T-cumPlastStrain*delta_TSensitivity)/delta_T/delta_T;

		         
//		          dissipatedEnergy += 0.5*(stress+0-beta_T*E*(CPlastStrain+plastStrain))*plastStrainInc;
				 dissipatedEnergySensitivity = CDissipatedEnergySensitivity + 0.5*(stress+0-beta_T*E*(CPlastStrain+plastStrain))*plastStrainIncSensitivity
			      + 0.5*(stressSensitivity+0-(beta_TSensitivity*E+beta_T*ESensitivity)*(CPlastStrain+plastStrain)
				  -beta_T*E*(CPlastStrainSensitivity+plastStrainSensitivity))*plastStrainInc; 
			  }
              else { //%stress(i-1) > 0
                     plastStrainInc = plastStrain - CPlastStrain;//Newton_BRB(CStress, beta_C, CPlastStrain, sigmaY_C,CCumPlastStrain, delta_C, alpha_C, strainInc,0, Tol, N0);
					 int sign = -1.0;
				     if (plastStrainInc>0)
			           sign =1.0;
				     double u = E*(strain_loading-plastStrainInc)-beta_C*E*(CPlastStrain+plastStrainInc);
				     double v = sigmaY0+(sigmaY_C-sigmaY0)*(1-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C));
				     double A = u/v;

					 double C = 0.0+(strain_loading-plastStrainInc-beta_C*CPlastStrain-beta_C*plastStrainInc)*ESensitivity+E*strain_loadingSensitivity
					      -E*(CPlastStrain+plastStrainInc)*beta_CSensitivity-beta_C*E*CPlastStrainSensitivity;
					 double D = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C))
						  -(sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*
							 (delta_CSensitivity*(CCumPlastStrain+fabs(plastStrainInc))/delta_C/delta_C-CCumPlastStrainSensitivity/delta_C);

					 double F = (sigmaY_C-sigmaY0)*exp(-(CCumPlastStrain+fabs(plastStrainInc))/delta_C)*sign/delta_C;
					
			 		 double G = alpha_C*pow(fabs(A),(alpha_C-2))*A*strain_loading*(C/v-u*D/v/v)+strain_loadingSensitivity*pow(fabs(A),alpha_C)
			             +pow(fabs(A),alpha_C)*log(fabs(A))*strain_loading*alpha_CSensitivity;
					 double H = alpha_C*pow(fabs(A),(alpha_C-2))*A*strain_loading*((E+E*beta_C)/v+u*F/v/v);
				  
					 plastStrainIncSensitivity = G/(1+H);

                     plastStrainSensitivity = CPlastStrainSensitivity+plastStrainIncSensitivity;
					 stressSensitivity = 0.0+ESensitivity*(strain_loading-plastStrainInc)+E*(strain_loadingSensitivity-plastStrainIncSensitivity);                      
					 cumPlastStrainSensitivity = CCumPlastStrainSensitivity+sign*plastStrainIncSensitivity;
	                 yieldStressSensitivity = sigmaY0Sensitivity+(sigmaY_CSensitivity-sigmaY0Sensitivity)*(1.0-exp(-cumPlastStrain/delta_C))
			           +(sigmaY_C-sigmaY0)*exp(-cumPlastStrain/delta_C)*(cumPlastStrainSensitivity*delta_C-cumPlastStrain*delta_CSensitivity)/delta_C/delta_C;

//      	     dissipatedEnergy += 0.5*(stress+0-beta_C*E*(CPlastStrain+plastStrain))*plastStrainInc;		         
				 dissipatedEnergySensitivity = CDissipatedEnergySensitivity + 0.5*(stress+0-beta_C*E*(CPlastStrain+plastStrain))*plastStrainIncSensitivity
			      + 0.5*(stressSensitivity+0-(beta_CSensitivity*E+beta_C*ESensitivity)*(CPlastStrain+plastStrain)
				  -beta_C*E*(CPlastStrainSensitivity+plastStrainSensitivity))*plastStrainInc; 
              } // %stress(i-1) < 0 
	  } // % if (abs(strainInc(i))<=CStress/E)
              
  }
	
//	return stressSensitivity;

 



	// Commit history variables
		(*SHVs)(0,(gradNumber-1)) = strainSensitivity;
		(*SHVs)(1,(gradNumber-1)) = stressSensitivity;
		
	    (*SHVs)(2,(gradNumber-1)) = plastStrainSensitivity;
		(*SHVs)(3,(gradNumber-1)) = cumPlastStrainSensitivity;
		(*SHVs)(4,(gradNumber-1)) = yieldStressSensitivity;
		(*SHVs)(5,(gradNumber-1)) = dissipatedEnergySensitivity;



   /* debug -----
		opserr.precision(16);
		opserr<<"strain: "<<strain<<",  stress:"<<stress<<",   ";
		opserr<<"strainSensitivity: "<<strainSensitivity<<",  stressSensitivity:"<<stressSensitivity<<endln;
*/

   /*

   *debug1 << "commitSensitivity PlastStrainInc:"<< plastStrainInc << ",  plastStrainIncSensitivity: "<<plastStrainIncSensitivity<<"\n"; 
   *debug1 << "commitSensitivity plastStrain:"<< plastStrain << ",  plastStrainSensitivity: "<<plastStrainSensitivity<<"\n"; 
   *debug1 << "commitSensitivity stress:"<< stress << ",  stressSensitivity: "<<stressSensitivity<<"\n";    
   *debug1 << "commitSensitivity cumPlastStrain:"<< cumPlastStrain << ",  cumPlastStrainSensitivity: "<<cumPlastStrainSensitivity<<"\n"; 
   *debug1 << "commitSensitivity yieldStress:"<< yieldStress << ",  yieldStressSensitivity: "<<yieldStressSensitivity<<"\n";

////////*/

//	opserr<<"BRB--commitStressSensitivity:" <<strainSensitivity<<endln; 
int ii = 0; 
	ii++;

	if (fabs(stressSensitivity) > 1.0e10) {
		opserr<< "error in commitSensitivity, ii = "<<ii<< endln; 
	}


	return 0;
}



double
SteelBRB::getInitialTangentSensitivity(int gradNumber)
{
	if (parameterID == 1) {
		return 1.0;
	}
	else {
		return 0.0;
	}
}


// AddingSensitivity:END /////////////////////////////////////////////*/






Response* 
SteelBRB::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{


  	if (strcmp(argv[0],"plasticStrain") == 0) {
		 
		return new MaterialResponse(this, 11, this->getStrain());  // use only size of matrix

	  }
	else if (strcmp(argv[0],"cumPlasticStrain") == 0) {
	 
		return new MaterialResponse(this, 12, this->getStrain());  // use only size of matrix

	  } 
	else if (strcmp(argv[0],"dissipatedEnergy") == 0) {
	 
		return new MaterialResponse(this, 13, this->getStrain());  // use only size of matrix

	  } 
/*	else if (strcmp(argv[0],"plasticStrainSensitivity") == 0) {
		int gradientNum = atoi(argv[1]);
		return new MaterialResponse(this, gradientNum+100, this->getStrain());  // use only size of matrix

	  }
	else if (strcmp(argv[0],"cumPlasticStrainSensitivity") == 0) {
		int gradientNum = atoi(argv[1]);
		return new MaterialResponse(this, gradientNum+500, this->getStrain());  // limitation: 400 RVs

	  }
	else if (strcmp(argv[0],"stressSensitivity") == 0) {
		int gradientNum = atoi(argv[1]);
		return new MaterialResponse(this, gradientNum+900, this->getStress());  // use only size of matrix

	  }
	  else if (strcmp(argv[0],"strainSensitivity") == 0) {
			int gradientNum = atoi(argv[1]);
			return new MaterialResponse(this, gradientNum+1300, this->getStrain());
	  }
*/

	 else if (strstr(argv[0],"plasticStrainSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradientNum = atoi(token);
  		return new MaterialResponse(this, gradientNum+100, this->getStrain());  // use only size of matrix

	  }
	 else if (strstr(argv[0],"cumPlasticStrainSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradientNum = atoi(token);
  		return new MaterialResponse(this, gradientNum+500, this->getStrain());  // use only size of matrix

	  } 
	 else if (strstr(argv[0],"stressSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradientNum = atoi(token);
  		return new MaterialResponse(this, gradientNum+900, this->getStrain());  // use only size of matrix

	  } 
	 else if (strstr(argv[0],"strainSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradientNum = atoi(token);
  		return new MaterialResponse(this, gradientNum+1300, this->getStrain());  // use only size of matrix

	  } 
	 else if (strstr(argv[0],"dissipatedEnergySensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradientNum = atoi(token);
  		return new MaterialResponse(this, gradientNum+1700, this->getStrain());  // use only size of matrix

	  } 


	  //by default, See if the response is one of the defaults
	  Response *res =  UniaxialMaterial::setResponse(argv, argc, theOutput);

	  if (res != 0)      return res;
	  else { 
		  opserr<<"error in SteelBRB::setResponse"<<endln;
		  return 0;
	  }
}

int 
SteelBRB::getResponse(int responseID, Information &matInfo)
{

	if (responseID==11) {
		return matInfo.setDouble(plastStrain);
	} // if
	else if (responseID==12) {
		return matInfo.setDouble(cumPlastStrain);
	} // if
	else if (responseID==13) {
		return matInfo.setDouble(dissipatedEnergy);
	} // if
	else if (responseID>100 && SHVs ==0) 
		return matInfo.setDouble(0.0);

	else if (responseID>100 && responseID<500) {
		return matInfo.setDouble((*SHVs)(2,(responseID-100-1)));
	} // if

	else if (responseID>500 && responseID<900) {
		return matInfo.setDouble((*SHVs)(3,(responseID-500-1)));

	} // if
	else if (responseID>900 && responseID<1300) {
		return matInfo.setDouble((*SHVs)(1,(responseID-900-1)));

	} // if
	else if (responseID>1300 && responseID<1700) {
		return matInfo.setDouble((*SHVs)(0,(responseID-1300-1)));

	} // if
	else if (responseID>1700) {
		return matInfo.setDouble((*SHVs)(5,(responseID-1700-1)));

	} // if


	else

	  // Just call the base class method ... don't need to define
	  // this function, but keeping it here just for clarity
	  return UniaxialMaterial::getResponse(responseID, matInfo);	
	
}

/*
	    (*SHVs)(2,(gradNumber-1)) = plastStrainSensitivity;
		(*SHVs)(3,(gradNumber-1)) = cumPlastStrainSensitivity;
		(*SHVs)(0,(gradNumber-1)) = strainSensitivity;
		(*SHVs)(1,(gradNumber-1)) = stressSensitivity;
*/
