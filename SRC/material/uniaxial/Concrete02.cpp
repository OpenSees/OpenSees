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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete02.cpp,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of Concrete02. 
// This Concrete02 is based on an f2c of the FEDEAS material
// Concr2.f which is:
//-----------------------------------------------------------------------
// concrete model with damage modulus    
//       by MOHD YASSIN (1993)
// adapted to FEDEAS material library
// by D. Sze and Filip C. Filippou in 1994
//-----------------------------------------------------------------------


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Concrete02.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <Matrix.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

void *
OPS_Concrete02()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Concrete02 tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 4 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete02 " << iData[0] << " fpc? epsc0? fpcu? epscu? <rat? ft? Ets?>\n";
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial Concrete02 " << iData[0] << " fpc? epsc0? fpcu? epscu? <rat? ft? Ets?>\n";
    return 0;
  }


  // Parsing was successful, allocate the material
  if (numData == 7)
    theMaterial = new Concrete02(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);
  else
    theMaterial = new Concrete02(iData[0], dData[0], dData[1], dData[2], dData[3]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Concrete02 Material\n";
    return 0;
  }

  return theMaterial;
}

Concrete02::Concrete02(int tag, double _fc, double _epsc0, double _fcu,
		       double _epscu, double _rat, double _ft, double _Ets):
  UniaxialMaterial(tag, MAT_TAG_Concrete02),
  fc(_fc), epsc0(_epsc0), fcu(_fcu), epscu(_epscu), rat(_rat), ft(_ft), Ets(_Ets),
  parameterID(0), SHVs(0)
{
  ecminP = 0.0;
  deptP = 0.0;

  if (fc > 0) fc = -fc;
  if (epsc0 > 0) epsc0 = -epsc0;
  if (fcu > 0) fcu = -fcu;
  if (epscu > 0) epscu = -epscu;

  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;
}

Concrete02::Concrete02(int tag, double _fc, double _epsc0, double _fcu,
		       double _epscu):
  UniaxialMaterial(tag, MAT_TAG_Concrete02),
  fc(_fc), epsc0(_epsc0), fcu(_fcu), epscu(_epscu),
  parameterID(0), SHVs(0)
{
  ecminP = 0.0;
  deptP = 0.0;

  if (fc > 0) fc = -fc;
  if (epsc0 > 0) epsc0 = -epsc0;
  if (fcu > 0) fcu = -fcu;
  if (epscu > 0) epscu = -epscu;
	  
  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;

  rat = 0.1;
  ft = 0.1*fc;
  if (ft < 0.0)
    ft = -ft;
  Ets = 0.1*fc/epsc0;
}

Concrete02::Concrete02(void):
  UniaxialMaterial(0, MAT_TAG_Concrete02),
  parameterID(0), SHVs(0)
{
 
}

Concrete02::~Concrete02(void)
{
  if (SHVs != 0)
    delete SHVs;
}

UniaxialMaterial*
Concrete02::getCopy(void)
{
  Concrete02 *theCopy = new Concrete02(this->getTag(), fc, epsc0, fcu, epscu, rat, ft, Ets);
  
  return theCopy;
}

double
Concrete02::getInitialTangent(void)
{
  return 2.0*fc/epsc0;
}

int
Concrete02::setTrialStrain(double trialStrain, double strainRate)
{
  double  ec0 = fc * 2. / epsc0;

  // retrieve concrete history variables

  ecmin = ecminP;
  dept = deptP;
  TEnergy = CEnergy;

  // calculate current strain

  eps = trialStrain;
  double deps = eps - epsP;

  if (fabs(deps) < DBL_EPSILON)
    return this->revertToLastCommit();

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

  TEnergy += 0.5 * (sigP + sig) * (eps - epsP);
  //opserr << "CE: " << CEnergy << " -> TE: " << TEnergy << " | s (" << sigP << ", " << sig << ", M: " << 0.5 * (sigP + sig) << "); dE: " << (eps - epsP) << "\n";

  return 0;
}



double 
Concrete02::getStrain(void)
{
  return eps;
}

double 
Concrete02::getStress(void)
{
  return sig;
}

double 
Concrete02::getTangent(void)
{
  return e;
}

int 
Concrete02::commitState(void)
{
  ecminP = ecmin;
  deptP = dept;
  
  eP = e;
  sigP = sig;
  epsP = eps;

  CEnergy = TEnergy;

  return 0;
}

int 
Concrete02::revertToLastCommit(void)
{
  ecmin = ecminP;;
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;

  TEnergy = CEnergy;

  return 0;
}

int 
Concrete02::revertToStart(void)
{
  ecminP = 0.0;
  deptP = 0.0;

  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;

  TEnergy = CEnergy = 0.0;

  if (SHVs != 0)
    SHVs->Zero();

  return 0;
}

int 
Concrete02::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(15);
  data(0) =fc;    
  data(1) =epsc0; 
  data(2) =fcu;   
  data(3) =epscu; 
  data(4) =rat;   
  data(5) =ft;    
  data(6) =Ets;   
  data(7) =ecminP;
  data(8) =deptP; 
  data(9) =epsP;  
  data(10) =sigP; 
  data(11) =eP;   
  data(12) = this->getTag();
  data(13) = parameterID;
  data(14) = -1;
  if (SHVs != 0)
    data(14) = SHVs->noCols();

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Concrete02::sendSelf() - failed to sendSelf\n";
    return -1;
  }

  if (SHVs != 0) {
    if (theChannel.sendMatrix(this->getDbTag(), commitTag, *SHVs) < 0) {
      opserr << "Concrete02::sendSelf() - failed to send SHVs matrix\n";
      return -2;
    }
  }

  return 0;
}

int 
Concrete02::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{

  static Vector data(15);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "Concrete02::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  fc = data(0);
  epsc0 = data(1);
  fcu = data(2);
  epscu = data(3);
  rat = data(4);
  ft = data(5);
  Ets = data(6);
  ecminP = data(7);
  deptP = data(8);
  epsP = data(9);
  sigP = data(10);
  eP = data(11);
  this->setTag(data(12));

  e = eP;
  sig = sigP;
  eps = epsP;

  parameterID = (int)data(13);

  // Receive sensitivity history variables (SHVs)
  int noCols = (int)data(14);

  if (SHVs != 0) {
    delete SHVs;
    SHVs = 0;
  }

  if (noCols > 0) {
    SHVs = new Matrix(4, noCols);
    if (theChannel.recvMatrix(this->getDbTag(), commitTag, *SHVs) < 0) {
      delete SHVs;
      SHVs = 0;
      opserr << "Concrete02::recvSelf() - failed to receive SHVs matrix\n";
      return -2;
    }
  }
  
  return 0;
}

void 
Concrete02::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {      
    s << "Concrete02:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"Concrete02\", ";
	s << "\"Ec\": " << 2.0*fc/epsc0 << ", ";
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
Concrete02::Tens_Envlp (double epsc, double &sigc, double &Ect)
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
  
  double Ec0  = 2.0*fc/epsc0;

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
      sigc = 0.0;
    }
  }
  return;
}

  
void
Concrete02::Compr_Envlp (double epsc, double &sigc, double &Ect) 
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

  double Ec0  = 2.0*fc/epsc0;

  double ratLocal = epsc/epsc0;
  if (epsc>=epsc0) {
    sigc = fc*ratLocal*(2.0-ratLocal);
    Ect  = Ec0*(1.0-ratLocal);
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
Concrete02::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
Concrete02::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // Compressive strength
  if (strcmp(argv[0],"fc") == 0 || strcmp(argv[0],"fpc") == 0) {
    param.setValue(fc);
    return param.addObject(1, this);
  }
  // Strain at compressive strength
  if (strcmp(argv[0],"epsc0") == 0 || strcmp(argv[0],"epsco") == 0) {
    param.setValue(epsc0);
    return param.addObject(2, this);
  }
  // Crushing strength
  if (strcmp(argv[0],"fcu") == 0 || strcmp(argv[0],"fpcu") == 0) {
    param.setValue(fcu);
    return param.addObject(3, this);
  }
  // Strain at crushing strength
  if (strcmp(argv[0],"epscu") == 0 || strcmp(argv[0],"epsU") == 0 ||
      strcmp(argv[0],"epsu") == 0) {
    param.setValue(epscu);
    return param.addObject(4, this);
  }
  // Ratio between unloading slope at epscu and initial slope
  if (strcmp(argv[0],"rat") == 0 || strcmp(argv[0],"lambda") == 0) {
    param.setValue(rat);
    return param.addObject(5, this);
  }
  // Tensile strength
  if (strcmp(argv[0],"ft") == 0) {
    param.setValue(ft);
    return param.addObject(6, this);
  }
  // Tension softening stiffness
  if (strcmp(argv[0],"Ets") == 0) {
    param.setValue(Ets);
    return param.addObject(7, this);
  }

  return -1;
}

int
Concrete02::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    this->fc = info.theDouble;
    break;
  case 2:
    this->epsc0 = info.theDouble;
    break;
  case 3:
    this->fcu = info.theDouble;
    break;
  case 4:
    this->epscu = info.theDouble;
    break;
  case 5:
    this->rat = info.theDouble;
    break;
  case 6:
    this->ft = info.theDouble;
    break;
  case 7:
    this->Ets = info.theDouble;
    break;
  default:
    return -1;
  }

  // Compression parameters are stored as negative values,
  // same convention as the constructors
  if (fc > 0.0) fc = -fc;
  if (epsc0 > 0.0) epsc0 = -epsc0;
  if (fcu > 0.0) fcu = -fcu;
  if (epscu > 0.0) epscu = -epscu;

  // Refresh the initial tangent if the material has not been loaded yet
  if (epsP == 0.0 && sigP == 0.0 && ecminP == 0.0 && deptP == 0.0) {
    eP = 2.0*fc/epsc0;
    e = eP;
  }

  return 0;
}

int
Concrete02::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;
}

void
Concrete02::Compr_EnvlpSensitivity(double epsc, double depsc,
				   double dfc, double depsc0,
				   double dfcu, double depscu,
				   double &sigc, double &dsigc)
{
  // Same branches as Compr_Envlp
  if (epsc >= epsc0) {
    // parabolic ascending branch
    double ratLocal = epsc/epsc0;
    double dratLocal = (depsc - ratLocal*depsc0)/epsc0;
    sigc = fc*ratLocal*(2.0-ratLocal);
    dsigc = dfc*ratLocal*(2.0-ratLocal) + fc*(2.0-2.0*ratLocal)*dratLocal;
  } else {
    if (epsc > epscu) {
      // linear descending branch between epsc0 and epscu
      double slope = (fcu-fc)/(epscu-epsc0);
      double dslope = ((dfcu-dfc) - slope*(depscu-depsc0))/(epscu-epsc0);
      sigc = fc + slope*(epsc-epsc0);
      dsigc = dfc + dslope*(epsc-epsc0) + slope*(depsc-depsc0);
    } else {
      // flat friction branch for strains larger than epscu
      sigc = fcu;
      dsigc = dfcu;
    }
  }
}

void
Concrete02::Tens_EnvlpSensitivity(double epsc, double depsc,
				  double dfc, double depsc0,
				  double dft, double dEts,
				  double &sigc, double &dsigc)
{
  // Same branches as Tens_Envlp
  double Ec0 = 2.0*fc/epsc0;
  double dEc0 = 2.0*dfc/epsc0 - Ec0*depsc0/epsc0;
  double eps0 = ft/Ec0;
  double deps0 = (dft - eps0*dEc0)/Ec0;
  double epsu = ft*(1.0/Ets+1.0/Ec0);

  if (epsc <= eps0) {
    sigc = epsc*Ec0;
    dsigc = depsc*Ec0 + epsc*dEc0;
  } else {
    if (epsc <= epsu) {
      sigc = ft - Ets*(epsc-eps0);
      dsigc = dft - dEts*(epsc-eps0) - Ets*(depsc-deps0);
    } else {
      sigc = 0.0;
      dsigc = 0.0;
    }
  }
}

// Concrete02 is piecewise smooth: at the corners of the envelopes and where
// the active branch changes, the stress has only one-sided derivatives.
// The derivative of the branch selected by setTrialStrain is returned.
double
Concrete02::computeStressSensitivity(double depsdh, int gradIndex,
				     double &decmindh, double &ddeptdh)
{
  // Committed sensitivity history variables
  double depsPdh = 0.0;
  double dsigPdh = 0.0;
  double decminPdh = 0.0;
  double ddeptPdh = 0.0;
  if (SHVs != 0 && gradIndex >= 0 && gradIndex < SHVs->noCols()) {
    depsPdh   = (*SHVs)(0,gradIndex);
    dsigPdh   = (*SHVs)(1,gradIndex);
    decminPdh = (*SHVs)(2,gradIndex);
    ddeptPdh  = (*SHVs)(3,gradIndex);
  }

  // Direct derivatives of the material parameters
  double dfc    = (parameterID == 1) ? 1.0 : 0.0;
  double depsc0 = (parameterID == 2) ? 1.0 : 0.0;
  double dfcu   = (parameterID == 3) ? 1.0 : 0.0;
  double depscu = (parameterID == 4) ? 1.0 : 0.0;
  double drat   = (parameterID == 5) ? 1.0 : 0.0;
  double dft    = (parameterID == 6) ? 1.0 : 0.0;
  double dEts   = (parameterID == 7) ? 1.0 : 0.0;

  // History sensitivities are unchanged unless a branch below updates them
  decmindh = decminPdh;
  ddeptdh = ddeptPdh;

  double ec0 = fc * 2. / epsc0;
  double dec0 = 2.0*dfc/epsc0 - ec0*depsc0/epsc0;

  double deps = eps - epsP;

  // No strain increment: setTrialStrain keeps the committed stress and
  // tangent, so a strain perturbation moves the stress along that tangent
  if (fabs(deps) < DBL_EPSILON)
    return dsigPdh + e*(depsdh - depsPdh);

  double sigc, dsigdh;

  if (eps < ecminP) {

    // monotonic envelope in compression; the minimum strain is reset
    this->Compr_EnvlpSensitivity(eps, depsdh, dfc, depsc0, dfcu, depscu,
				 sigc, dsigdh);
    decmindh = depsdh;

  } else {

    // point R that determines the reloading slope
    double numer = fcu - rat * ec0 * epscu;
    double dnumer = dfcu - drat * ec0 * epscu - rat * dec0 * epscu
      - rat * ec0 * depscu;
    double denom = ec0 * (1.0 - rat);
    double ddenom = dec0 * (1.0 - rat) - ec0 * drat;
    double epsr = numer / denom;
    double depsr = (dnumer - epsr * ddenom) / denom;
    double sigmr = ec0 * epsr;
    double dsigmr = dec0 * epsr + ec0 * depsr;

    // previous minimum stress from the compression envelope
    double sigmm, dsigmm;
    this->Compr_EnvlpSensitivity(ecminP, decminPdh, dfc, depsc0, dfcu, depscu,
				 sigmm, dsigmm);

    // reloading slope and its intersection with the zero stress axis
    double er = (sigmm - sigmr) / (ecminP - epsr);
    double der = ((dsigmm - dsigmr) - er * (decminPdh - depsr))
      / (ecminP - epsr);
    double ept = ecminP - sigmm / er;
    double dept_dh = decminPdh - dsigmm / er + sigmm * der / (er * er);

    if (eps <= ept) {

      // unloading-reloading branch in compression
      double sigmin = sigmm + er * (eps - ecminP);
      double dsigmin = dsigmm + der * (eps - ecminP)
	+ er * (depsdh - decminPdh);
      double sigmax = er * .5f * (eps - ept);
      double dsigmax = .5f * (der * (eps - ept) + er * (depsdh - dept_dh));

      sigc = sigP + ec0 * deps;
      dsigdh = dsigPdh + dec0 * deps + ec0 * (depsdh - depsPdh);
      if (sigc <= sigmin) {
	sigc = sigmin;
	dsigdh = dsigmin;
      }
      if (sigc >= sigmax) {
	sigc = sigmax;
	dsigdh = dsigmax;
      }

    } else {

      double epn = ept + deptP;
      if (eps <= epn) {

	// reloading branch in tension along the secant to the
	// maximum remaining tensile strength
	double sicn, dsicn;
	this->Tens_EnvlpSensitivity(deptP, ddeptPdh, dfc, depsc0, dft, dEts,
				    sicn, dsicn);
	double esec, desec;
	if (deptP != 0.0) {
	  esec = sicn / deptP;
	  desec = (dsicn - esec * ddeptPdh) / deptP;
	} else {
	  esec = ec0;
	  desec = dec0;
	}
	sigc = esec * (eps - ept);
	dsigdh = desec * (eps - ept) + esec * (depsdh - dept_dh);

      } else {

	// tensile envelope shifted by ept; the tensile offset is reset
	double epstmp = eps - ept;
	double depstmp = depsdh - dept_dh;
	this->Tens_EnvlpSensitivity(epstmp, depstmp, dfc, depsc0, dft, dEts,
				    sigc, dsigdh);
	ddeptdh = depstmp;
      }
    }
  }

  return dsigdh;
}

double
Concrete02::getStressSensitivity(int gradIndex, bool conditional)
{
  // Conditional derivative: strain held fixed
  double decmindh, ddeptdh;

  return this->computeStressSensitivity(0.0, gradIndex, decmindh, ddeptdh);
}

double
Concrete02::getInitialTangentSensitivity(int gradIndex)
{
  double dfc    = (parameterID == 1) ? 1.0 : 0.0;
  double depsc0 = (parameterID == 2) ? 1.0 : 0.0;
  double ec0 = 2.0*fc/epsc0;

  return 2.0*dfc/epsc0 - ec0*depsc0/epsc0;
}

int
Concrete02::commitSensitivity(double TstrainSensitivity, int gradIndex,
			      int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(4,numGrads);
    SHVs->Zero();
  }

  if (gradIndex < 0 || gradIndex >= SHVs->noCols())
    return -1;

  double decmindh, ddeptdh;
  double dsigdh = this->computeStressSensitivity(TstrainSensitivity, gradIndex,
						 decmindh, ddeptdh);

  // Commit sensitivity history variables
  (*SHVs)(0,gradIndex) = TstrainSensitivity;
  (*SHVs)(1,gradIndex) = dsigdh;
  (*SHVs)(2,gradIndex) = decmindh;
  (*SHVs)(3,gradIndex) = ddeptdh;

  return 0;
}
// AddingSensitivity:END /////////////////////////////////////////////
