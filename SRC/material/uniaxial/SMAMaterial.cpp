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

// Description: This file contains the class implementation for SMAMaterial.

#include <SMAMaterial.h> 
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_SMAMaterial)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 7) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial SMA matTag? E? eps_L? sig_AM_s? sig_AM_f? sig_MA_s? sig_MA_f?" << endln;
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING: failed to read tag\n";
	return 0;
    }

    double data[6];
    numdata = 6;
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING: failed to read data\n";
	return 0;
    }

    UniaxialMaterial* mat = new SMAMaterial(tag,data[0],data[1],data[2],data[3],data[4],data[5]);
    if (mat == 0) {
	opserr << "WARNING: failed to create SMAMaterial\n";
	return 0;
    }

    return mat;
}

SMAMaterial::SMAMaterial(int tag, double e, double el, double s_as_s, double s_as_f, 
			 double s_sa_s, double s_sa_f)
  :UniaxialMaterial(tag,MAT_TAG_SMA), 
   E(e), eps_L(el), sig_AS_s(s_as_s), sig_AS_f(s_as_f), sig_SA_s(s_sa_s), sig_SA_f(s_sa_f), 
   Tcsi(0.0), Tstrain(0.0), Tstress(0.0), Ttangent(0.0)
{
  // Initialize varibles
  this->revertToStart();
}

SMAMaterial::SMAMaterial()
  :UniaxialMaterial(0,MAT_TAG_SMA), 
   E(0.0), eps_L(0.0), sig_AS_s(0.0), sig_AS_f(0.0), sig_SA_s(0.0), sig_SA_f(0.0), 
   Tcsi(0.0), Tstrain(0.0), Tstress(0.0), Ttangent(0.0)
{
  // Initialize variables
  this->revertToStart();
}

SMAMaterial::~SMAMaterial()
{
   // does nothing
}

int 
SMAMaterial::setTrialStrain(double strain, double strainRate)
{
  Tstrain = strain;
  
  //-------------------------------------------------------------------------------------
  // A ---> S transformation (loading), sigma > 0 and epsilon > 0
  //-------------------------------------------------------------------------------------
  
  if ((Tstrain - Cstrain) > 0 && Tstrain > 0) {  

    double eps_AS_s = sig_AS_s/E + Ccsi*eps_L;
    double eps_AS_f = sig_AS_f/E +      eps_L;
    
    if (Tstrain <= eps_AS_s) {
      
      Tcsi = Ccsi;
      Tstress = E*(Tstrain - Tcsi*eps_L);
      Ttangent = E;
      
    }
    
    else if (Tstrain > eps_AS_s && Tstrain < eps_AS_f) { 
      
      if (Cstress <= sig_AS_s) {
	
	Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_AS_f - E*Tstrain + sig_AS_s)/(-sig_AS_f + Ccsi*E*eps_L - E*eps_L + sig_AS_s);
	Tstress = E*(Tstrain - Tcsi*eps_L);
	double H = -(1 - Ccsi)*E/((1 - Ccsi)*(-E*eps_L) + sig_AS_s - sig_AS_f);
	Ttangent = E*(1 - eps_L*H);
	
      }
      
      else {
	
	Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_AS_f - E*Tstrain + Cstress)/(-sig_AS_f + Ccsi*E*eps_L - E*eps_L + Cstress);
	Tstress = E*(Tstrain - Tcsi*eps_L);
	double H = -(1 - Ccsi)*E/((1 - Ccsi)*(-E*eps_L) + Cstress - sig_AS_f);
	Ttangent = E*(1 - eps_L*H);
        
      }
      
    }
    
    else {
      
      Tcsi = 1;
      Tstress = E*(Tstrain - Tcsi*eps_L);
      Ttangent = E;
      
    }

  }

  //-------------------------------------------------------------------------------------
  // S ---> A  transformation (unloading), sigma > 0 and epsilon > 0
  //-------------------------------------------------------------------------------------
  
  if ((Tstrain - Cstrain) < 0 && Tstrain > 0) { 
    
    double eps_SA_s = sig_SA_s/E + Ccsi*eps_L;
    double eps_SA_f = sig_SA_f/E;
    
    if (Tstrain >= eps_SA_s) {
      
      Tcsi = Ccsi;
      Tstress = E*(Tstrain - Tcsi*eps_L);
      Ttangent = E;
      
    }
    
    else if (Tstrain < eps_SA_s && Tstrain >= eps_SA_f) { 
      
      if (Cstress > sig_SA_s) {
	
	Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_SA_f)/(-sig_SA_f + Ccsi*E*eps_L + sig_SA_s);
	Tstress = E*(Tstrain - Tcsi*eps_L);
	double H = (Ccsi*E)/(- Ccsi*(-E*eps_L) + sig_SA_s - sig_SA_f);
	Ttangent = E*(1 - eps_L*H);
	
      }
      
      else {
	
	Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_SA_f)/(-sig_SA_f + Ccsi*E*eps_L + Cstress);
	Tstress = E*(Tstrain - Tcsi*eps_L);
	double H = (Ccsi*E)/(-Ccsi*(-E*eps_L) + Cstress - sig_SA_f);
	Ttangent = E*(1 - eps_L*H);
        
      }
      
    }
    
    else {
      
      Tcsi = 0;
      Tstress = E*(Tstrain);
      Ttangent = E;
      
    }
    
  }

  //-------------------------------------------------------------------------------------
  // S ---> A  transformation (loading), sigma < 0 and epsilon < 0
  //-------------------------------------------------------------------------------------
  
  if ((Tstrain - Cstrain) < 0 && Tstrain < 0.0) { 
    
    double eps_AS_s = -sig_AS_s/E - Ccsi*eps_L;
    double eps_AS_f = -sig_AS_f/E - eps_L;
    
    if (Tstrain >= eps_AS_s) {
      
      Tcsi = Ccsi;
      Tstress = E*(Tstrain + Tcsi*eps_L);
      Ttangent = E;
      
    }
    
    else if (Tstrain < eps_AS_s && Tstrain > eps_AS_f) { 
      
      if (Cstress >= -sig_AS_s) {
	
	Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_AS_f - E*Tstrain - sig_AS_s)/(sig_AS_f - Ccsi*E*eps_L + E*eps_L - sig_AS_s);
	Tstress = E*(Tstrain + Tcsi*eps_L);
	double H = (1 - Ccsi)*E/((1 - Ccsi)*(E*eps_L) - sig_AS_s + sig_AS_f);
	Ttangent = E*(1 - eps_L*H);
	
      }
      
      else {
	
	Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_AS_f - E*Tstrain + Cstress)/(sig_AS_f - Ccsi*E*eps_L + E*eps_L + Cstress);
	Tstress = E*(Tstrain + Tcsi*eps_L);
	double H = (1 - Ccsi)*E/((1 - Ccsi)*(E*eps_L) + Cstress + sig_AS_f);
	Ttangent = E*(1 - eps_L*H);
        
      }
      
    }
    
    else {
      
      Tcsi = 1;
      Tstress = E*(Tstrain + Tcsi*eps_L);
      Ttangent = E;
      
    }
    
  }
  
  //-------------------------------------------------------------------------------------
  // S ---> A  transformation (unloading), sigma < 0 and epsilon < 0
  //-------------------------------------------------------------------------------------
  
  if  ((Tstrain - Cstrain) > 0 && Tstrain < 0.0) { 
    
    double eps_SA_s = -sig_SA_s/E - Ccsi*eps_L;
    double eps_SA_f = -sig_SA_f/E;
    
    if (Tstrain <= eps_SA_s) {
      
      Tcsi = Ccsi;
      Tstress = E*(Tstrain + Tcsi*eps_L);
      Ttangent = E;
      
    }
    
    else if (Tstrain > eps_SA_s && Tstrain <= eps_SA_f) { 
      
      if (Cstress < -sig_SA_s) {
	
	Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_SA_f)/(sig_SA_f - Ccsi*E*eps_L - sig_SA_s);
	Tstress = E*(Tstrain + Tcsi*eps_L);
	double H = -(Ccsi*E)/(-Ccsi*(E*eps_L) - sig_SA_s + sig_SA_f);
	Ttangent = E*(1 - eps_L*H);
	
      }
      
      else {
	
	Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_SA_f)/(sig_SA_f - Ccsi*E*eps_L + Cstress);
	Tstress = E*(Tstrain + Tcsi*eps_L);
	double H = -(Ccsi*E)/(-Ccsi*(E*eps_L) + Cstress + sig_SA_f);
	Ttangent = E*(1 - eps_L*H);
        
      }
      
    }
    
    else {
      
      Tcsi = 0;
      Tstress = E*(Tstrain); 
      Ttangent = E;
      
    }
    
  }

  return 0;
  
}

double 
SMAMaterial::getStress(void)
{
  return Tstress;
}

double 
SMAMaterial::getTangent(void)
{
  return Ttangent;
}

double 
SMAMaterial::getInitialTangent(void)
{
  return E;
}

double 
SMAMaterial::getStrain(void)
{
  return Tstrain;
}

int 
SMAMaterial::commitState(void)
{

  // Commit trial state variables

  Cstrain = Tstrain;
  Cstress = Tstress;
  Ccsi = Tcsi;

  return 0;
}

int 
SMAMaterial::revertToLastCommit(void)
{
  return 0;
}

int 
SMAMaterial::revertToStart(void)
{
  
  // Reset committed history variables

  Cstrain = 0.0;
  Cstress = 0.0;
  Ccsi = 0.0;

  // Reset trial history variables

  // Initialize state variables

  Tcsi = 0.0;
  Tstrain = 0.0;
  Tstress = 0.0;
  Ttangent = E;

  return 0;
}

UniaxialMaterial *
SMAMaterial::getCopy(void)
{
  SMAMaterial *theCopy = new SMAMaterial(this->getTag(), E, eps_L, sig_AS_s, sig_AS_f, sig_SA_s, sig_SA_f);

  // Copy committed history variables
 
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Ccsi = Ccsi;

  // Copy trial state variables
 
  theCopy->Tcsi = Tcsi;
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;

  return theCopy;
}

int 
SMAMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(11);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = eps_L;
  data(3) = sig_AS_s;
  data(4) = sig_AS_f;
  data(5) = sig_SA_s;
  data(6) = sig_SA_f;
  data(7) = Cstrain;
  data(8) = Cstress;
  data(9) = Ccsi;
  data(10) = Ttangent;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "SMAMaterial::sendSelf() - failed to send data\n";

  return res;  
}

int 
SMAMaterial::recvSelf(int cTag, Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(11);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
    opserr << "SMAMaterial::recvSelf() - failed to receive data\n";
    E = 0; 
    this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    eps_L = data(2);
    sig_AS_s = data(3);
    sig_AS_f = data(4);
    sig_SA_s = data(5);
    sig_SA_f = data(6);
    Cstrain = data(7);
    Cstress = data(8);
    Ccsi = data(9);
    Ttangent = data(10);
	  
    Tstrain = Cstrain;
    Tstress = Cstress;
    Tcsi = Ccsi;
  }
    
  return res;
}

void 
SMAMaterial::Print(OPS_Stream &s, int flag)
{
  s << "SMAMaterial, tag: " << this->getTag() << endln;
  s << " E:         " << E << endln;
  s << " eps_L:     " << eps_L << endln;
  s << " sig_AS_s:  " << sig_AS_s << endln;
  s << " sig_AS_f:  " << sig_AS_f << endln;
  s << " sig_SA_s:  " << sig_SA_s << endln;
  s << " sig_SA_f:  " << sig_SA_f << endln;
  return;
}


