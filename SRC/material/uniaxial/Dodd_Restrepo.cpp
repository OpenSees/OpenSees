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
                                                                        
 // Written: fmk
 // Created: 10/2011
 //
 // Description: This file contains the class definition for 
 // Dodd_Restrepo. Dodd_Restrepo is a c++ wrapper to the steel material
 // found in DoddRestrepo.f that was written by L.L. Dodd and J. Restrepo  1991

#include <elementAPI.h>

#include <OPS_Globals.h>
#include "Dodd_Restrepo.h"
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

static int numDoddRestrepo = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_Dodd_Restrepo)
{
  if (numDoddRestrepo == 0) {
    numDoddRestrepo++;
    opserr << "Dodd_Restrepo unaxial material - Written by L.L. Dodd & J. Restepo\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 8 || numArgs > 10) {
    opserr << "WARNING wrong # args: uniaxialMaterial $tag $Fy $Fsu $ESH $ESU $Youngs $ESHI $FSHI <$OmegaFac>" << endln;
    return 0;
  }

  int    iData[1];
  double dData[9];
  int numData;

  dData[7]=1.0; // omegaFac
  dData[8]=1.0; // Conv

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endln;
    return 0;
  }

  numData = numArgs-1;  
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E & ep\n";
    return 0;	
  }

  theMaterial = new Dodd_Restrepo(iData[0], dData[0], dData[1], dData[2],
				  dData[3], dData[4], dData[5], dData[6],
				  dData[7], dData[8]);       
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticPPCpp\n";
    return 0;
  }
  
  return theMaterial;
}


Dodd_Restrepo::Dodd_Restrepo(int tag, 
			     double fy, 
			     double fsu, 
			     double eSH, 
			     double eSU,
			     double youngs, 
			     double eSHI, 
			     double fSHI, 
			     double omegaFac,
			     double conv)
  :UniaxialMaterial(tag,0),
   Fy(fy), Fsu(fsu), ESH(eSH), ESU(eSU), Youngs(youngs), 
   ESHI(eSHI), FSHI(fSHI), Conv(conv), OmegaFac(omegaFac)
{
  numDoddRestrepo++;
  myTag = numDoddRestrepo;


  if (OmegaFac < 0.65) OmegaFac = 0.65;
  if (OmegaFac > 1.15) OmegaFac = 1.15;
  
  double C1       ; // Temporary constant 
  double EpSHI    ; // Intermediate strain hardening curve natural strain
  double FpSH     ; // True stress at initiation of strain hardening curve
  double FpSHI    ; // Intermediate strain hardening curve true stress

  Epy       = Fy/Youngs;
  EpSH      = log(1+ESH/Conv);
  Epsu      = log(1+ESU/Conv);
  Fpsu      = Fsu*(1+ESU/Conv);

  EpsuSh[0] =  Epsu;
  EpsuSh[1] = -Epsu;
  YoungsUn  = Youngs;

  LMR       = 0;
  BFlag[0]  = 0;
  BFlag[1]  = 0;

  Epa[0]    = 0.0;
  Epa[1]    = 0.0;
  EpaM[0]   = 0.0;
  EpaM[1]   = 0.0;
  Epo[0]   = 0.0;
  Epo[1]    = 0.0;
  EpoMax    = 0.0;
  Epr[0]    = 0.0;
  Epr[1]    = 0.0;
  EprM[0]   = 0.0;
  EprM[1]   = 0.0;

  Fpr[0]   =0.0;
  Fpa[0]   =0.0;
  Power[0] =0.0;
  FprM[0]  =0.0;
  FpaM[0]  =0.0;
  YpTanM[0]=0.0;
  PowerM[0]=0.0; 
  Fpr[1]   =0.0;
  Fpa[1]   =0.0;
  Power[1] =0.0;
  FprM[1]  =0.0;
  FpaM[1]  =0.0;
  YpTanM[1]=0.0;
  PowerM[1]=0.0; 

  //
  // Calculate the power term for the strain hardening branch
  //

  EpSHI   = log(1+ESHI/Conv);
  FpSH    = Fy*(1+ESH/Conv);
  FpSHI   = FSHI*(1+ESHI/Conv);
  C1      = FpSH-Fpsu+Fpsu*(Epsu-EpSH);
  SHPower = log((FpSHI+Fpsu*(Epsu-EpSHI)-Fpsu)/C1) / log((Epsu-EpSHI)/(Epsu-EpSH));
  
  tStrain = 0.0;
  tTangent = Youngs;
  tStress = 0.0;

  Eps = 0.0;
  EpsOld = 0.0;
  EpsLast = 0.0;
  Fps = 0.0;
  FpsLast = 0.0;
  YpTan = Youngs;
  YpTanLast = Youngs;
    
  this->commitState();
}

UniaxialMaterial *
Dodd_Restrepo::getCopy(void)
{
  return new Dodd_Restrepo(this->getTag(), Fy, Fsu, ESH, ESU, Youngs, ESHI, FSHI, OmegaFac, Conv);
}

#ifdef _WIN32

#define steel_ STEEL

extern "C" int  steel_ (double *Es, double *EpsLast, double *FpsLast, double *YpTanLast, 
		       double *EpsOld, double *Fy, double *Epy, double * EpSH, double *Epsu, 
		       double *Fpsu, double *Youngs, double *SHPower, double *Epr, double *Fpr, 
		       double *Epa, double *Fpa, double *Epo, double *EpoMax, double *EpsuSh, 
		       double *YoungsUn, double *Power, int *BFlag, int *LMR, double *EprM, 
		       double *FprM, double *EpaM, double *FpaM, double *YpTanM, double *PowerM, 
		       double * Eps, double *Fps, double *Fs, double *YpTan, double *YTan, double *OmegFac);

// Add more declarations as needed


#else

extern "C" int  steel_ (double *Es, double *EpsLast, double *FpsLast, double *YpTanLast, 
			   double *EpsOld, double *Fy, double *Epy, double * EpSH, double *Epsu, 
			   double *Fpsu, double *Youngs, double *SHPower, double *Epr, double *Fpr, 
			   double *Epa, double *Fpa, double *Epo, double *EpoMax, double *EpsuSh, 
			   double *YoungsUn, double *Power, int *BFlag, int *LMR, double *EprM, 
			   double *FprM, double *EpaM, double *FpaM, double *YpTanM, double *PowerM, 
			   double * Eps, double *Fps, double *Fs, double *YpTan, double *YTan, double *OmegFac);

#endif

int
Dodd_Restrepo::setTrialStrain(double strain, double strainRate)
{
  if (fabs(strain-tStrain) > DBL_EPSILON) {
    tStrain = strain;

    Epr[1]= cEpr[1];     Epr[0]= cEpr[0];   
    Fpr[1]= cFpr[1];     Fpr[0]= cFpr[0];   
    Epa[1]= cEpa[1];     Epa[0]= cEpa[0];   
    Fpa[1]= cFpa[1];     Fpa[0]= cFpa[0];   
    Epo[1]= cEpo[1];     Epo[0]= cEpo[0];   
    EpoMax  = cEpoMax; 
    EpsuSh[1]= cEpsuSh[1];  EpsuSh[0]= cEpsuSh[0];
    YoungsUn = cYoungsUn;
    Power[1]= cPower[1];  Power[0]= cPower[0]; 
    BFlag[1]= cBFlag[1];   BFlag[0]= cBFlag[0];
    LMR = cLMR;
    EprM[1]= cEprM[1];    EprM[0]= cEprM[0];  
    FprM[1]= cFprM[1];    FprM[0]= cFprM[0];  
    EpaM[1]= cEpaM[1];    EpaM[0]= cEpaM[0];  
    FpaM[1]= cFpaM[1];    FpaM[0]= cFpaM[0];  
    YpTanM[1]= cYpTanM[1]; YpTanM[0]= cYpTanM[0];
    PowerM[1]= cPowerM[1]; PowerM[0]= cPowerM[0];
    
    steel_ (&tStrain, &EpsLast, &FpsLast, &YpTanLast, 
	     &EpsOld, &Fy, &Epy,  &EpSH, &Epsu, 
	     &Fpsu, &Youngs, &SHPower, Epr, Fpr, 
	     Epa, Fpa, Epo, &EpoMax, EpsuSh,&YoungsUn, Power, BFlag, &LMR, EprM,
	     FprM, EpaM, FpaM, YpTanM, PowerM,
	     &Eps, &Fps, &Fs, &YpTan, &YTan, &OmegaFac);
    
    tStress = Fs;
    tTangent = YTan;
  }

  return 0;
}

int
Dodd_Restrepo::setTrial(double strain, double &stress, double &stiff, double strainRate)
{
  if (fabs(strain-tStrain) > DBL_EPSILON) {

    // Store the strain
    tStrain = strain;

    Epr[1]= cEpr[1];     Epr[0]= cEpr[0];   
    Fpr[1]= cFpr[1];     Fpr[0]= cFpr[0];   
    Epa[1]= cEpa[1];     Epa[0]= cEpa[0];   
    Fpa[1]= cFpa[1];     Fpa[0]= cFpa[0];   
    Epo[1]= cEpo[1];     Epo[0]= cEpo[0];   
    EpoMax  = cEpoMax; 
    EpsuSh[1]= cEpsuSh[1];  EpsuSh[0]= cEpsuSh[0];
    YoungsUn = cYoungsUn;
    Power[1]= cPower[1];  Power[0]= cPower[0]; 
    BFlag[1]= cBFlag[1];   BFlag[0]= cBFlag[0];
    LMR = cLMR;
    EprM[1]= cEprM[1];    EprM[0]= cEprM[0];  
    FprM[1]= cFprM[1];    FprM[0]= cFprM[0];  
    EpaM[1]= cEpaM[1];    EpaM[0]= cEpaM[0];  
    FpaM[1]= cFpaM[1];    FpaM[0]= cFpaM[0];  
    YpTanM[1]= cYpTanM[1]; YpTanM[0]= cYpTanM[0];
    PowerM[1]= cPowerM[1]; PowerM[0]= cPowerM[0];

    steel_ (&tStrain, &EpsLast, &FpsLast, &YpTanLast, 
	    &EpsOld, &Fy, &Epy,  &EpSH, &Epsu, 
	    &Fpsu, &Youngs, &SHPower, Epr, Fpr, 
	    Epa, Fpa, Epo, &EpoMax, EpsuSh, 
	    &YoungsUn, Power, BFlag, &LMR, EprM, 
	    FprM, EpaM, FpaM, YpTanM, PowerM, 
	    &Eps, &Fps, &Fs, &YpTan, &YTan, &OmegaFac);
    
    tStress = Fs;
    tTangent = YTan;
  }
  
  /*
  if (myTag == 41 || myTag == 59) {
    this->Print(opserr);
    opserr << tStrain << endln;
    opserr << EpsLast << endln;
    opserr << FpsLast << endln;
    opserr << "YpTanLast: " << YpTanLast << endln;
    opserr << EpsOld << endln;
    opserr << Fy << endln;
    opserr << Epy<< endln;
    opserr << EpSH << endln;
    opserr << Epsu<< endln;
    opserr << Fpsu << endln;
    opserr << Youngs << endln;
    opserr << SHPower << endln;
    opserr << Epr[1]<< " " << cEpr[1] << " " <<     
    Epr[0]<< " " << cEpr[0] << endln;   
    opserr << Fpr[1]<< " " << cFpr[1] << " " <<     
    Fpr[0]<< " " << cFpr[0] << endln;
    opserr << Epa[1]<< " " << cEpa[1] << " " <<  
    Epa[0]<< " " << cEpa[0] << endln;
    opserr << Fpa[1]<< " " << cFpa[1] << " " <<  
    Fpa[0]<< " " << cFpa[0] << endln;
    opserr << Epo[1]<< " " << cEpo[1] << " " <<
    Epo[0]<< " " << cEpo[0];   
    opserr << EpoMax  << " " << cEpoMax << endln;
    opserr << EpsuSh[1]<< " " << cEpsuSh[1] << " " <<
    EpsuSh[0]<< " " << cEpsuSh[0] << endln;
    opserr << YoungsUn << " " << cYoungsUn;
    opserr << Power[1]<< " " << cPower[1] << " " << 
    Power[0]<< " " << cPower[0] << endln;
    opserr << BFlag[1]<< " " << cBFlag[1] << " " << 
    BFlag[0]<< " " << cBFlag[0] << endln;
    opserr << LMR << " " << cLMR;
    opserr << EprM[1]<< " " << cEprM[1] << " " <<
    EprM[0]<< " " << cEprM[0] << endln;
    opserr << FprM[1]<< " " << cFprM[1] << " " <<  
    FprM[0]<< " " << cFprM[0] << endln;
    opserr << EpaM[1]<< " " << cEpaM[1] << " " <<  
    EpaM[0]<< " " << cEpaM[0] << endln;
    opserr << FpaM[1]<< " " << cFpaM[1] << " " << FpaM[0]<< " " << cFpaM[0] << endln;
    opserr << YpTanM[1]<< " " << cYpTanM[1] << " " << YpTanM[0]<< " " << cYpTanM[0] << endln;
    opserr << PowerM[1]<< " " << cPowerM[1] <<  " " << PowerM[0] << endln;
    opserr << "Fs: " << Fs << " Ytan: " << YTan << endln;
  }
  */
    
  stress = tStress;
  stiff = tTangent;


  
  return 0;
}

double
Dodd_Restrepo::getStrain(void)
{
  return tStrain;
}

double
Dodd_Restrepo::getStress(void)
{
  return tStress;
}

double
Dodd_Restrepo::getTangent(void)
{
  return tTangent;
}

double
Dodd_Restrepo::getInitialTangent(void)
{
  return Youngs;
}

int
Dodd_Restrepo::commitState(void)
{
  if (cStrain != tStrain) {
    EpsOld    = EpsLast;
    EpsLast   = Eps;
    FpsLast   = Fps;
    YpTanLast = YpTan;
  }
  cStrain = tStrain;
  cStress = tStress;
  cTangent = tTangent;
  
  cEpr[1]=Epr[1];     cEpr[0]=Epr[0];   
  cFpr[1]=Fpr[1];     cFpr[0]=Fpr[0];   
  cEpa[1]= Epa[1];     cEpa[0]= Epa[0];   
  cFpa[1]= Fpa[1];     cFpa[0]= Fpa[0];   
  cEpo[1]= Epo[1];     cEpo[0]= Epo[0];   
  cEpoMax  = EpoMax; 
  cEpsuSh[1]= EpsuSh[1];  cEpsuSh[0]= EpsuSh[0];
  cYoungsUn = YoungsUn;
  cPower[1]= Power[1];   cPower[0]= Power[0]; 
  cBFlag[1]= BFlag[1];  cBFlag[0]= BFlag[0];
  cLMR = LMR;
  cEprM[1]= EprM[1];    cEprM[0]= EprM[0];  
  cFprM[1]= FprM[1];    cFprM[0]= FprM[0];  
  cEpaM[1]= EpaM[1];    cEpaM[0]= EpaM[0];  
  cFpaM[1]= FpaM[1];    cFpaM[0]= FpaM[0];  
  cYpTanM[1]= YpTanM[1];  cYpTanM[0]= YpTanM[0];
  cPowerM[1]= PowerM[1];  cPowerM[0]= PowerM[0];

  
  return 0;
}

int
Dodd_Restrepo::revertToLastCommit(void)
{
  tStrain = cStrain;
  tStress = cStress;
  tTangent = cTangent;

  return 0;
}

int
Dodd_Restrepo::revertToStart(void)
{
  double C1       ; // Temporary constant
  double EpSHI    ; // Intermediate strain hardening curve natural strain
  double ESH      ; // Engineering coordinate strain hardening strain
  double ESHI     ; // Intermediate strain hardening curve engineering strain
  double FpSH     ; // True stress at initiation of strain hardening curve
  double FpSHI    ; // Intermediate strain hardening curve true stress

  Epy       = Fy/Youngs;
  EpSH      = log(1+ESH/Conv);
  Epsu      = log(1+ESU/Conv);
  Fpsu      = Fsu*(1+ESU/Conv);
  EpsuSh[0] =  Epsu;
  EpsuSh[1] = -Epsu;
  YoungsUn  = Youngs;

  LMR       = 0;
  BFlag[0]  = 0;
  BFlag[1]  = 0;

  Epa[0]    = 0.0;
  Epa[1]    = 0.0;
  EpaM[0]   = 0.0;
  EpaM[1]   = 0.0;
  Epo[0]   = 0.0;
  Epo[1]    = 0.0;
  EpoMax    = 0.0;
  Epr[0]    = 0.0;
  Epr[1]    = 0.0;
  EprM[0]   = 0.0;
  EprM[1]   = 0.0;

  Fpr[0]   =0.0;
  Fpa[0]   =0.0;
  Power[0] =0.0;
  FprM[0]  =0.0;
  FpaM[0]  =0.0;
  YpTanM[0]=0.0;
  PowerM[0]=0.0; 
  Fpr[1]   =0.0;
  Fpa[1]   =0.0;
  Power[1] =0.0;
  FprM[1]  =0.0;
  FpaM[1]  =0.0;
  YpTanM[1]=0.0;
  PowerM[1]=0.0; 

  //
  // Calculate the power term for the strain hardening branch
  //

  EpSHI   = log(1+ESHI/Conv);
  FpSH    = Fy  *(1+ESH/Conv);
  FpSHI   = FSHI*(1+ESHI/Conv);
  C1      = FpSH-Fpsu+Fpsu*(Epsu-EpSH);
  SHPower = log((FpSHI+Fpsu*(Epsu-EpSHI)-Fpsu)/C1) / log((Epsu-EpSHI)/(Epsu-EpSH));

  tStrain = 0.0;
  tTangent = Youngs;
  tStress = 0.0;
    
  this->commitState();

  return 0;
}

int 
Dodd_Restrepo::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  return res;
}

int
Dodd_Restrepo::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  return res;
}

void
Dodd_Restrepo::Print(OPS_Stream &s, int flag)
{
  s << "Dodd_Restrepo: " << this->getTag() << endln;
  s << "tStrain: " << tStrain << "  tStress: " << tStress << " tTangent: " << tTangent << endln;
}






