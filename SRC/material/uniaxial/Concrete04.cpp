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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-07-07 22:58:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete04.cpp,v $
                                                                        
// Written: N.Mitra (nmitra@u.washington.edu) 
// Created: 09/04
// Revision: A
//
// Description: This file contains the class implementation for 
// Concrete04 based on Popovics pre and post compression curve 
// for concrete. 
//
// What: "@(#) Concrete04.C, revA"
// Revision 1. Adding in exponential tensile envelope for concrete
// Dt. 05-16-05


#include <Concrete04.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

void* OPS_Concrete04()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? Ec0?";
	opserr << " <ft? etu? <beta?> >\n";
	return 0;
    }

    int type = 1;

    int tag;
    numdata = 1;
    if(OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[4];
    numdata = 4;
    if(OPS_GetDoubleInput(&numdata,data) < 0) {
	opserr << "WARNING invalid double data\n";
	return 0;
    }

    numdata = OPS_GetNumRemainingInputArgs();
    double data2[2];
    if (numdata > 1) {
	numdata = 2;
	if (OPS_GetDoubleInput(&numdata,data2)<0) {
	    opserr << "WARNING invalid double data\n";
	    return 0;
	}
	type = 2;
    }

    numdata = OPS_GetNumRemainingInputArgs();
    double beta;
    if (numdata > 0) {
	numdata = 1;
	if (OPS_GetDoubleInput(&numdata,&beta)) {
	    opserr << "WARNING invalid double data\n";
	    return 0;
	}
	type = 3;
    }

    UniaxialMaterial* mat = 0;
    if (type == 1) {
	mat = new Concrete04(tag,data[0],data[1],data[2],data[3]);
    } else if (type == 2) {
	mat = new Concrete04(tag,data[0],data[1],data[2],data[3],data2[0],data2[1]);
    } else if (type == 3) {
	mat = new Concrete04(tag,data[0],data[1],data[2],data[3],data2[0],data2[1],beta);
    }

    if (mat == 0) {
	opserr << "WARNING: failed to create Concrete04 material\n";
	return 0;
    }

    return mat;
}

Concrete04::Concrete04
(int tag, double FPC, double EPSC0, double EPSCU, double EC0, double FCT, double ETU)
  :UniaxialMaterial(tag, MAT_TAG_Concrete04),
   fpc(FPC), epsc0(EPSC0), epscu(EPSCU), Ec0(EC0), fct(FCT), etu(ETU), beta(0.1),
   CminStrain(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(FCT),
   Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0) 
{
  // Make all concrete parameters negative
  if (fpc > 0.0 || epsc0 > 0.0 || epscu > 0.0) {
    opserr << "error: negative values required for concrete stress-strain model" << endln;
  }
  
  if (fct < 0.0) {
    fct = 0.0;
    opserr << "warning: fct less than 0.0 so the tensile response part is being set to 0" << endln;
  }
  
  Ctangent = Ec0;
  CunloadSlope = Ec0;
  CUtenSlope = Ec0;
  
  // Set trial values
  this->revertToLastCommit();
}

Concrete04::Concrete04
(int tag, double FPC, double EPSC0, double EPSCU, double EC0, double FCT, double ETU, double BETA)
  :UniaxialMaterial(tag, MAT_TAG_Concrete04),
   fpc(FPC), epsc0(EPSC0), epscu(EPSCU), Ec0(EC0), fct(FCT), etu(ETU), beta(BETA),
   CminStrain(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(FCT),
   Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0) 
{
  // Make all concrete parameters negative
  if (fpc > 0.0 || epsc0 > 0.0 || epscu > 0.0) {
    opserr << "error: negative values required for concrete stress-strain model" << endln;
  }
  
  if (fct < 0.0) {
    fct = 0.0;
    opserr << "warning: fct less than 0.0 so the tensile response part is being set to 0" << endln;
  }
  
  Ctangent = Ec0;
  CunloadSlope = Ec0;
  CUtenSlope = Ec0;
  
  // Set trial values
  this->revertToLastCommit();
}


Concrete04::Concrete04
(int tag, double FPC, double EPSC0, double EPSCU, double EC0)
  :UniaxialMaterial(tag, MAT_TAG_Concrete04),
   fpc(FPC), epsc0(EPSC0), epscu(EPSCU), Ec0(EC0), fct(0.0), etu(0.0), beta(0.0),
   CminStrain(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(0.0),
   Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0) 
{
  // Make all concrete parameters negative
  if (fpc > 0.0 || epsc0 > 0.0 || epscu > 0.0) {
    opserr << "error: negative values required for concrete stress-strain model" << endln;
  }
  
  Ctangent = Ec0;
  CunloadSlope = Ec0;
  CUtenSlope = 0.0;
  
  // Set trial values
  this->revertToLastCommit();
  
}

Concrete04::Concrete04():UniaxialMaterial(0, MAT_TAG_Concrete04),
			 fpc(0.0), epsc0(0.0), epscu(0.0), Ec0(0.0), fct(0.0), etu(0.0), beta(0.0),
			 CminStrain(0.0), CunloadSlope(0.0), CendStrain(0.0), CcompStrain(0.0), CUtenStress(0.0),
			 CUtenSlope(0.0), Cstrain(0.0), Cstress(0.0), CmaxStrain(0.0)
{
  // Set trial values
  this->revertToLastCommit();
  
}

Concrete04::~Concrete04()
{
  // Does nothing
}

int Concrete04::setTrialStrain (double strain, double strainRate)
{

  /*// Reset trial history variables to last committed state*/
  TminStrain = CminStrain;   
  TmaxStrain = CmaxStrain;   
  TendStrain = CendStrain;   
  TunloadSlope = CunloadSlope;   
  TUtenSlope = CUtenSlope;
  Tstrain = Cstrain;   
  Tstress = Cstress;   
  Ttangent = Ctangent;

  /* // Set trial strain*/  
  if (fct == 0.0 && strain > 0.0) {    
    Tstrain = strain;
    Tstress = 0.0;    
    Ttangent = 0.0;    
    TUtenSlope = 0.0;    
    return 0;  
  }

  /*// Determine change in strain from last converged state*/  
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)       
    return 0;

  Tstrain = strain;  

  /*// Calculate the trial state given the change in strain  // determineTrialState (dStrain);*/
  TunloadSlope = CunloadSlope;  
  TUtenSlope = CUtenSlope;
  if (dStrain <= 0.0) {	  /*// Material can be either in Compression-Reloading	  // or Tension-Unloading state.*/
    if (Tstrain > 0.0) {         
      /*// Material is in Tension-Unloading State*/		  
      Ttangent = TUtenSlope;		  
      Tstress = Tstrain * TUtenSlope; 	  
    } else {
      /*// Material is in Compression-Reloading State*/      
      TminStrain = CminStrain;      
      TendStrain = CendStrain;      
      TunloadSlope = CunloadSlope;      
      CompReload();
    }
  } else {
    /*// Material can be either in Compression-Unloading	  // or Tension-Reloading State.*/
    if (Tstrain >= 0.0) {    /*// Material is in Tension-Reloading State*/      
      TmaxStrain = CmaxStrain;                  
      if (Tstrain < TmaxStrain) {        
	Tstress = Tstrain * CUtenSlope;        
	Ttangent = CUtenSlope;        
	TUtenSlope = CUtenSlope;      
      } else {        
	TmaxStrain = Tstrain;        
	TensEnvelope();        
	setTenUnload();      
      }        
    } else {
      if (Tstrain <= TendStrain) {          
	Ttangent = TunloadSlope;          
	Tstress = Ttangent * (Tstrain - TendStrain);        
      } else {          
	Tstress = 0.0;          
	Ttangent = 0.0;        
      }        
    }  
  }    
  return 0;
}

void Concrete04::CompReload()
{
  if (Tstrain <= TminStrain) {
    
    TminStrain = Tstrain;
    
    /*// Determine point on envelope*/
    CompEnvelope ();
    setCompUnloadEnv ();
    
  }
  else if (Tstrain < TendStrain) {
    Ttangent = TunloadSlope;
    Tstress = Ttangent*(Tstrain-TendStrain);
  }
  else if (Tstrain <= 0.0) {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
}

void Concrete04::CompEnvelope()
{
  if (Tstrain >= epscu) {
    double Esec = fpc/epsc0;
    double r = 0.0;
    if (Esec >= Ec0) {
      r = 400.0;
    } else {
      r = Ec0/(Ec0-Esec);
    }
    double eta = Tstrain/epsc0;
    Tstress = fpc*eta*r/(r-1+pow(eta,r));
    Ttangent = fpc*r*(r-1)*(1-pow(eta,r))/(pow((r-1+pow(eta,r)),2)*epsc0);
  } else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  
}

void Concrete04::setCompUnloadEnv()	{
  double tempStrain = TminStrain;
  
  if (tempStrain < epscu)
    tempStrain = epscu;
  
  double eta = tempStrain/epsc0;
  
  double ratio = 0.707*(eta-2.0) + 0.834; // unloading parameter as per Karsan-Jirsa
  
  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;
  
  TendStrain = ratio*epsc0;
  
  double temp1 = TminStrain - TendStrain;
  
  double temp2 = Tstress/Ec0;
  
  if (temp1 > -DBL_EPSILON) {	// temp1 should always be negative
    TunloadSlope = Ec0;
  }
  else if (temp1 <= temp2) {
    TendStrain = TminStrain - temp1;
    TunloadSlope = Tstress/temp1;
  }
  else {
    TendStrain = TminStrain - temp2;
    TunloadSlope = Ec0;
  }
  
  
  if (Tstrain >= 0.0) {
    /*opserr << "actually made it in here" << endln;*/
    /*TunloadSlope = Ec0;*/
  }
  
}

void Concrete04::TensReload()
{  TensEnvelope();  setTenUnload();}

void Concrete04::TensEnvelope()
{  double ect = fct / Ec0;    
  if (Tstrain <= ect) {    
    Tstress = Tstrain * Ec0;    
    Ttangent = Ec0;
  } else if (Tstrain > etu) {    
    Tstress = 0.0;    
    Ttangent = 0.0;  
  } else {    
    Tstress = fct * pow(beta, (Tstrain - ect) / (etu - ect));    
    Ttangent = fct * pow(beta, (Tstrain - ect) / (etu - ect)) * log(beta) / (etu - ect);  
  }
}

void Concrete04::setTenUnload(){
  TUtenStress = Tstress;
  TUtenSlope = Tstress / Tstrain;
}
double Concrete04::getStress ()
{     return Tstress;}

double Concrete04::getStrain (){   return Tstrain;}

double Concrete04::getTangent ()
{   return Ttangent;}

int Concrete04::commitState ()
{
  /*// History variables*/   
  CminStrain = TminStrain;   
  CmaxStrain = TmaxStrain;   
  CunloadSlope = TunloadSlope;   
  CendStrain = TendStrain;   
  CUtenSlope = TUtenSlope;
  CcompStrain = TcompStrain;
  CUtenStress = TUtenStress;
  CUtenSlope = TUtenSlope;
  
  /*// State variables*/   
  Cstrain = Tstrain;   
  Cstress = Tstress;   
  Ctangent = Ttangent;      
  return 0;
}

int Concrete04::revertToLastCommit ()
{
  /*// Reset trial history variables to last committed state*/
  TminStrain = CminStrain;   
  TmaxStrain = CmaxStrain;   
  TunloadSlope = CunloadSlope;   
  TendStrain = CendStrain;   
  TUtenSlope = CUtenSlope;
  TcompStrain = CcompStrain;
  TUtenStress = CUtenStress;
  TUtenSlope = CUtenSlope;
  
  /*// Recompute trial stress and tangent*/   
  Tstrain = Cstrain;   
  Tstress = Cstress;   
  Ttangent = Ctangent;
  return 0;
}

int Concrete04::revertToStart ()
{
  
  /*// History variables*/
  CminStrain = 0.0;   CmaxStrain = 0.0;   CunloadSlope = Ec0;   CendStrain = 0.0;   CUtenSlope = Ec0;
  /*// State variables*/
  Cstrain = 0.0;   Cstress = 0.0;   Ctangent = Ec0;
  /*// Reset trial variables and state*/   this->revertToLastCommit();      return 0;
}

UniaxialMaterial* Concrete04::getCopy ()
{
  Concrete04* theCopy = new Concrete04(this->getTag(),
				       fpc, epsc0, epscu, Ec0, fct, etu, beta);
  
  /*// Converged history variables*/
  theCopy->CminStrain = CminStrain;   theCopy->CmaxStrain = CmaxStrain;   theCopy->CunloadSlope = CunloadSlope;   theCopy->CendStrain = CendStrain;   theCopy->CUtenSlope = CUtenSlope;
  
  /*// Converged state variables*/   theCopy->Cstrain = Cstrain;   theCopy->Cstress = Cstress;   theCopy->Ctangent = Ctangent;
  
  return theCopy;
}

int Concrete04::sendSelf (int commitTag, Channel& theChannel)
{   
  int res = 0;   

  static Vector data(16);   
  data(0) = this->getTag();
  
  /* Material properties*/   
  data(1) = fpc;   
  data(2) = epsc0;   
  data(3) = epscu;   
  data(4) = Ec0;   
  data(5) = fct;

  /*// History variables from last converged state*/
  data(6) = CminStrain;   
  data(7) = CmaxStrain;
  data(8) = CunloadSlope;   
  data(9) = CendStrain;   
  data(10) = CcompStrain;
  data(11) = CUtenStress;   
  data(12) = CUtenSlope;   

  
  /*// State variables from last converged state*/   
  data(13) = Cstrain;   
  data(14) = Cstress;   
  data(15) = Ctangent;   
  /*// Data is only sent after convergence, so no trial variables   // need to be sent through data vector*/
  
  res = theChannel.sendVector(this->getDbTag(), commitTag, data);   
  if (res < 0)       
    opserr << "Concrete04::sendSelf() - failed to send data\n";
  return res;
}

int Concrete04::recvSelf (int commitTag, Channel& theChannel,
			  FEM_ObjectBroker& theBroker)
{
  int res = 0;
  static Vector data(16);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "Concrete04::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    
    /*// Material properties */
    fpc = data(1);
    epsc0 = data(2);
    epscu = data(3);
    Ec0 = data(4);
    fct = data(5);
    
    /*// History variables from last converged state*/
    CminStrain = data(6);      
    CmaxStrain = data(7); 
    CunloadSlope = data(8);      
    CendStrain = data(9);      
    CcompStrain = data(10);      
    CUtenStress = data(11);      
    CUtenSlope = data(12);      

    
    /*// State variables from last converged state*/      
    Cstrain = data(13);      
    Cstress = data(14);      
    Ctangent = data(15);
    
    /*// Set trial state variables*/      
    this->revertToLastCommit();
  }
  
  return res;
}

void Concrete04::Print (OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "Concrete04, tag: " << this->getTag() << endln;
		s << "  fpc: " << fpc << endln;
		s << "  epsc0: " << epsc0 << endln;
		s << "  fct: " << fct << endln;
		s << "  epscu: " << epscu << endln;
		s << "  Ec0:  " << Ec0 << endln;
		s << "  etu:  " << etu << endln;
		s << "  beta: " << beta << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"Concrete04\", ";
		s << "\"Ec\": " << Ec0 << ", ";
		s << "\"fc\": " << fpc << ", ";
		s << "\"epsc\": " << epsc0 << ", ";
		s << "\"ft\": " << fct << ", ";
		s << "\"epstu\": " << etu << ", ";
		s << "\"epscu\": " << epscu << ", ";
		s << "\"beta\": " << beta << "}";
	}
}

/*// LOWES: add functions for variable hinge-length model*/
int
Concrete04::getMaterialType()
{
	return 0;
}
/*// LOWES: end*/
