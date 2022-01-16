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
                                                                        
// $Revision: 1.3 $
// $Date: 2007-06-08 00:38:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete02IS.cpp,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of Concrete02IS. 
// This Concrete02IS is based on an f2c of the FEDEAS material
// Concr2.f which is:
//-----------------------------------------------------------------------
// concrete model with damage modulus    
//       by MOHD YASSIN (1993)
// adapted to FEDEAS material library
// by D. Sze and Filip C. Filippou in 1994
// Modified by Nasser A. Marafi (2018) so that the concrete intiial stiffness
// can be specified. 
//-----------------------------------------------------------------------


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Concrete02IS.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

void * OPS_ADD_RUNTIME_VPV(OPS_Concrete02IS)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[8];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Concrete02IS tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 5 && numData != 8) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete02IS " << iData[0] << " E0? fpc? epsc0? fpcu? epscu? <rat? ft? Ets?>\n";
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete02IS " << iData[0] << " E0? fpc? epsc0? fpcu? epscu? <rat? ft? Ets?>\n";
    return 0;
  }


  // Parsing was successful, allocate the material
  if (numData == 8)
    theMaterial = new Concrete02IS(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7]);
  else
    theMaterial = new Concrete02IS(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Concrete02IS Material\n";
    return 0;
  }

  return theMaterial;
}

Concrete02IS::Concrete02IS(int tag, double _E0, double _fc, double _epsc0, double _fcu,
		       double _epscu, double _rat, double _ft, double _Ets):
  UniaxialMaterial(tag, MAT_TAG_Concrete02IS),
  fc(_fc), epsc0(_epsc0), fcu(_fcu), epscu(_epscu), rat(_rat), ft(_ft), Ets(_Ets), E0(_E0)
{
  ecminP = 0.0;
  deptP = 0.0;

  eP = E0;//2.0*fc / epsc0; // Marafi Change 2018/01/31
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = E0;//2.0*fc/epsc0; // Marafi Change 2018/01/31
  //E0 = E0;
}

Concrete02IS::Concrete02IS(int tag, double _E0, double _fc, double _epsc0, double _fcu,
		       double _epscu):
  UniaxialMaterial(tag, MAT_TAG_Concrete02IS),
  fc(_fc), epsc0(_epsc0), fcu(_fcu), epscu(_epscu), E0(_E0)
{
  ecminP = 0.0;
  deptP = 0.0;

  eP = E0;//2.0*fc / epsc0; // Marafi Change 2018/01/31
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = E0;//2.0*fc/epsc0; // Marafi Change 2018/01/31
  //E0 = E0;

  rat = 0.1;
  ft = 0.1*fc;
  if (ft < 0.0)
    ft = -ft;
  Ets = 0.1*fc/epsc0;  
}

Concrete02IS::Concrete02IS(void):
  UniaxialMaterial(0, MAT_TAG_Concrete02IS)
{
 
}

Concrete02IS::~Concrete02IS(void)
{
  // Does nothing
}

UniaxialMaterial*
Concrete02IS::getCopy(void)
{
  Concrete02IS *theCopy = new Concrete02IS(this->getTag(), E0, fc, epsc0, fcu, epscu, rat, ft, Ets);
  
  return theCopy;
}

double
Concrete02IS::getInitialTangent(void)
{
  //return 2.0*fc/epsc0; // Marafi Change 2018/01/31
  return E0;
}

int
Concrete02IS::setTrialStrain(double trialStrain, double strainRate)
{
	double  ec0 = E0; // fc * 2. / epsc0; // Marafi Change 2018/01/31

  // retrieve concrete hitory variables

  ecmin = ecminP;
  dept = deptP;

  // calculate current strain

  eps = trialStrain;
  double deps = eps - epsP;

  if (fabs(deps) < DBL_EPSILON)
    return 0;

  // if the current strain is less than the smallest previous strain 
  // call the monotonic envelope in compression and reset minimum strain 

  if (eps < ecmin) {
    this->Compr_Envlp(eps, sig, e);
    ecmin = eps;
  } else {;

    // else, if the current strain is between the minimum strain and ept 
    // (which corresponds to zero stress) the material is in the unloading- 
    // reloading branch and the stress remains between sigmin and sigmax 
    
    // calculate strain-stress coordinates of point R that determines 
    // the reloading slope according to Fig.2.11 in EERC Report 
    // (corresponding equations are 2.31 and 2.32 
    // the strain of point R is epsR and the stress is sigmR 
    
    double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
    double sigmr = ec0 * epsr;
    
    // calculate the previous minimum stress sigmm from the minimum 
    // previous strain ecmin and the monotonic envelope in compression 
    
    double sigmm;
    double dumy;
    this->Compr_Envlp(ecmin, sigmm, dumy);
    
    // calculate current reloading slope Er (Eq. 2.35 in EERC Report) 
    // calculate the intersection of the current reloading slope Er 
    // with the zero stress axis (variable ept) (Eq. 2.36 in EERC Report) 
    
    double er = (sigmm - sigmr) / (ecmin - epsr);
    double ept = ecmin - sigmm / er;
    
    if (eps <= ept) {
      double sigmin = sigmm + er * (eps - ecmin);
      double sigmax = er * .5f * (eps - ept);
      sig = sigP + ec0 * deps;
      e = ec0;
      if (sig <= sigmin) {
		sig = sigmin;
		e = er;
      }
      if (sig >= sigmax) {
		sig = sigmax;
		e = 0.5 * er;
      }
    } else {
      
      // else, if the current strain is between ept and epn 
      // (which corresponds to maximum remaining tensile strength) 
      // the response corresponds to the reloading branch in tension 
      // Since it is not saved, calculate the maximum remaining tensile 
      // strength sicn (Eq. 2.43 in EERC Report) 
      
      // calculate first the strain at the peak of the tensile stress-strain 
      // relation epn (Eq. 2.42 in EERC Report) 
      
      double epn = ept + dept;
      double sicn;
      if (eps <= epn) {
		this->Tens_Envlp(dept, sicn, e);
		if (dept != 0.0) {
		  e = sicn / dept;
		} else {
		  e = ec0;
		}
		sig = e * (eps - ept);
		  } else {
	
	// else, if the current strain is larger than epn the response 
	// corresponds to the tensile envelope curve shifted by ept 
	
		double epstmp = eps - ept;
		this->Tens_Envlp(epstmp, sig, e);
		dept = eps - ept;
		  }
    }
  }

  return 0;
}



double 
Concrete02IS::getStrain(void)
{
  return eps;
}

double 
Concrete02IS::getStress(void)
{
  return sig;
}

double 
Concrete02IS::getTangent(void)
{
  return e;
}

int 
Concrete02IS::commitState(void)
{
  ecminP = ecmin;
  deptP = dept;
  
  eP = e;
  sigP = sig;
  epsP = eps;
  return 0;
}

int 
Concrete02IS::revertToLastCommit(void)
{
  ecmin = ecminP;;
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;
  return 0;
}

int 
Concrete02IS::revertToStart(void)
{
  ecminP = 0.0;
  deptP = 0.0;

  eP = E0; // 2.0*fc / epsc0; // Marafi Change 2018/01/31
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = E0;//2.0*fc/epsc0; // Marafi Change 2018/01/31

  return 0;
}

int 
Concrete02IS::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(14);
  data(0) =E0;    
  data(1) =fc;    
  data(2) =epsc0; 
  data(3) =fcu;   
  data(4) =epscu; 
  data(5) =rat;   
  data(6) =ft;    
  data(7) =Ets;   
  data(8) =ecminP;
  data(9) =deptP; 
  data(10) =epsP;  
  data(11) =sigP; 
  data(12) =eP;   
  data(13) = this->getTag();

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Concrete02IS::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
Concrete02IS::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{

  static Vector data(14);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Concrete02IS::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  E0 = data(0);
  fc = data(1);
  epsc0 = data(2);
  fcu = data(3);
  epscu = data(4);
  rat = data(5);
  ft = data(6);
  Ets = data(7);
  ecminP = data(8);
  deptP = data(9);
  epsP = data(10);
  sigP = data(11);
  eP = data(12);
  this->setTag(data(13));

  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}

void 
Concrete02IS::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {      
    s << "Concrete02IS:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"Concrete02IS\", ";
	s << "\"E0\": " << E0 << ", ";
	s << "\"fc\": " << fc << ", ";
    s << "\"epsc\": " << epsc0 << ", ";
    s << "\"fcu\": " << fcu << ", ";
    s << "\"epscu\": " << epscu << ", ";
    s << "\"ratio\": " << rat << ", ";
    s << "\"ft\": " << ft << ", ";
    s << "\"Ets\": " << Ets << "}";
  }
}


void
Concrete02IS::Tens_Envlp (double epsc, double &sigc, double &Ect)
{
/*-----------------------------------------------------------------------
! monotonic envelope of concrete in tension (positive envelope)
!
!   ft    = concrete tensile strength
!   Ec0   = initial tangent modulus of concrete 
!   Ets   = tension softening modulus
!   eps   = strain
!
!   returned variables
!    sigc  = stress corresponding to eps
!    Ect  = tangent concrete modulus
!-----------------------------------------------------------------------*/
  
	double Ec0 = E0; // 2.0*fc / epsc0; // Marafi Change 2018/01/31

  double eps0 = ft/Ec0;
  double epsu = ft*(1.0/Ets+1.0/Ec0);
  if (epsc<=eps0) {
    sigc = epsc*Ec0;
    Ect  = Ec0;
  } else {
    if (epsc<=epsu) {
      Ect  = -Ets;
      sigc = ft-Ets*(epsc-eps0);
    } else {
      //      Ect  = 0.0
      Ect  = 1.0e-10;
	  sigc = 1.0e-10;//0.0; // Marafi
    }
  }
  return;
}

  
void
Concrete02IS::Compr_Envlp (double epsc, double &sigc, double &Ect) 
{
/*-----------------------------------------------------------------------
! monotonic envelope of concrete in compression (negative envelope)
!
!   fc    = concrete compressive strength
!   epsc0 = strain at concrete compressive strength
!   fcu   = stress at ultimate (crushing) strain 
!   epscu = ultimate (crushing) strain
!   Ec0   = initial concrete tangent modulus
!   epsc  = strain
!
!   returned variables
!   sigc  = current stress
!   Ect   = tangent concrete modulus
-----------------------------------------------------------------------*/

	double Ec0 = E0; //2.0*fc / epsc0; // Marafi Change 2018/01/31

  //double ratLocal = epsc/epsc0;
  if (epsc>=epsc0) {
	// Modify this to use Popovics Prepeak Curve
	  double Esec = fc / epsc0;
	  double r = 0.0;
	  if (Esec >= Ec0) {
		  r = 400.0;
	  }
	  else {
		  r = Ec0 / (Ec0 - Esec);
	  }
	  double eta = epsc / epsc0;
	  double fpc = fc;
      sigc = fpc * eta*r / (r - 1 + pow(eta, r));
	  Ect = fpc * r*(r - 1)*(1 - pow(eta, r)) / (pow((r - 1 + pow(eta, r)), 2)*epsc0);
      //sigc = fc*ratLocal*(2.0-ratLocal); // Marafi Change 2018/01/31
      //Ect  = Ec0*(1.0-ratLocal); // Marafi Change 2018/01/31
  } else {
    
    //   linear descending branch between epsc0 and epscu
    if (epsc>epscu) {
      sigc = (fcu-fc)*(epsc-epsc0)/(epscu-epsc0)+fc;
      Ect  = (fcu-fc)/(epscu-epsc0);
    } else {
	   
      // flat friction branch for strains larger than epscu
      
      sigc = fcu;
      Ect  = 1.0e-10;
      //       Ect  = 0.0
    }
  }
  return;
}

int
Concrete02IS::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}
