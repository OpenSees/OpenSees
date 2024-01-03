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
                                                                        
 //----------------------------------------------------------------------------------------------------------------------------
 // Developed by:
 // Michael H. Scott
 //
 // Based on TDConcrete implementations by:
 // Adam M. Knaack (adam.knaack@schaefer-inc.com) 
 // Schaefer-Inc, Cincinnati, Ohio, USA
 // Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
 // Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
 // Yahya C. Kurama (ykurama@nd.edu)
 // Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Description: This file contains the source code of CreepMaterial. 
 // CreepMaterial is a wrapper that imposes creep and shrinkage evoluation equations
 // to any uniaxialMaterial.
 //----------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "CreepMaterial.h"
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>
#include <Domain.h>
#include <MaterialResponse.h>
#include <Vector.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>

#include <Concrete02IS.h>
#include <ElasticMaterial.h>

static int numCreepMaterial = 0;

void *
OPS_CreepMaterial() {
  // Print description of material model:
  if (numCreepMaterial == 0) {
    //opserr << "Time-Dependent Concrete Material Model - Written by Adam Knaack, University of Notre Dame, 2012 \n";
    numCreepMaterial = 1;
  }
  
  // Pointer to a uniaxial material that will be returned:
  UniaxialMaterial *theMaterial = 0;
  
  // Parse the input line for the material parameters:
  int iData;
  int numData;
  int numArgs;
  
  numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs == 15) {
    //CreepMaterial(int tag, double _fc, double _epsc0, double _fcu,
    //double _epscu, double _tcr, double _ft, double _Ets, double _Ec, double _age, double _epsshu)
    double dData[14];
    
    //Collect material tag:
    numData = 1;
    if (OPS_GetIntInput(&numData, &iData) != 0) {
      opserr << "WARNING: invalid uniaxialMaterial CreepMaterial tag\n";
      return 0;
    }
    
    //Collect input parameters:
    numData = 14;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING: invalid material property definition\n";
      return 0;
    }
    
    //Create a new materiadouble
    theMaterial = new CreepMaterial(iData,dData[0],dData[1],dData[2],dData[3],dData[4],dData[5],dData[6],dData[7],dData[8],dData[9],dData[10],dData[11],dData[12],dData[13]);
    if (theMaterial == 0) {
      opserr << "WARNING: could not create uniaxialMaterial of type CreepMaterial \n";
      return 0;
    }
    
    //Return new material:
    return theMaterial;
  }
  if (numArgs == 10) {
    //CreepMaterial(int tag, double _fc, double _epsc0, double _fcu,
    //double _epscu, double _tcr, double _ft, double _Ets, double _Ec, double _age, double _epsshu)
    double dData[14];
    
    //Collect material tag:
    numData = 1;
    if (OPS_GetIntInput(&numData, &iData) != 0) {
      opserr << "WARNING: invalid uniaxialMaterial CreepMaterial tag\n";
      return 0;
    }

    int wrappedMatl;
    if (OPS_GetIntInput(&numData, &wrappedMatl) != 0) {
      opserr << "WARNING: invalid uniaxialMaterial CreepMaterial wrapped material tag\n";
      return 0;
    }

    UniaxialMaterial *matl = OPS_getUniaxialMaterial(wrappedMatl);
    if (matl == 0) {
      opserr << "WARNING: CreepMaterial - unable to find material with tag " << wrappedMatl << endln;
      return 0;
    }
    
    //Collect input parameters:
    numData = 8;
    if (OPS_GetDoubleInput(&numData, &dData[6]) != 0) {
      opserr << "WARNING: invalid material property definition\n";
      return 0;
    }
    
    //Create a new materiadouble
    theMaterial = new CreepMaterial(iData,*matl,dData[6],dData[7],dData[8],dData[9],dData[10],dData[11],dData[12],dData[13]);
    if (theMaterial == 0) {
      opserr << "WARNING: could not create uniaxialMaterial of type CreepMaterial \n";
      return 0;
    }
    
    //Return new material:
    return theMaterial;
  }  
  
  return 0;
}

//-----------------------------------------------------------------------


CreepMaterial::CreepMaterial(int tag, double _fc, double _fcu, double _epscu, double _ft, double _Ec, double _beta, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast): 
  UniaxialMaterial(tag, MAT_TAG_CreepMaterial), wrappedMaterial(0),
  fc(_fc), fcu(_fcu), epscu(_epscu), ft(_ft), Ec(_Ec), beta(_beta), age(_age), epsshu(_epsshu), epssha(_epssha), tcr(_tcr), epscru(_epscru), epscra(_epscra), epscrd(_epscrd), tcast(_tcast), maxSize(startSize),
  PHI_i(0), E_i(0), DSIG_i(0), TIME_i(0), DTIME_i(0)
{
  wrappedMaterial = new Concrete02IS(0,Ec,fc,2*fc/Ec,fcu,epscu);
  //wrappedMaterial = new ElasticMaterial(0,Ec);

  ecminP = 0.0;
  deptP = 0.0;
  
  //sigCr = fabs(sigCr);
  eP = Ec; //Added by AMK
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = Ec; //Added by AMK
  Et = Ec;
  count = 0; //Added by AMK
  epsInit = 0.0; //Added by AMK
  sigInit = 0.0; //Added by AMK
  eps_total = 0.0; //Added by AMK
  epsP_total = 0.0; //Added by AMK
  
  eps_m = 0.0; //Added by AMK
  eps_cr = 0.0; //Added by AMK
  eps_sh = 0.0;
  epsP_cr = 0.0; //Added by AMK
  epsP_sh = 0.0; 
  epsP_m = 0.0; //Added by AMK
  
  t_load = -1.0; //Added by AMK
  crack_flag = 0;
  iter = 0;
  
  
  
  //Change inputs into the proper sign convention:
  fc = -fabs(fc); 
  epsshu = -fabs(epsshu);
  epscru = fabs(epscru);

  this->expandArrays();
}

void
CreepMaterial::expandArrays()
{
  if (PHI_i == 0)
    PHI_i = new float[maxSize];

  if (E_i == 0)
    E_i = new float[maxSize];

  if (DSIG_i == 0)
    DSIG_i = new float[maxSize];

  if (TIME_i == 0)
    TIME_i = new float[maxSize];

  if (DTIME_i == 0)
    DTIME_i = new float[maxSize];

  if (count+1 >= maxSize) {
    maxSize += growSize;

    float *a = new float[maxSize];
    float *b = new float[maxSize];
    float *c = new float[maxSize];
    float *e = new float[maxSize];
    float *f = new float[maxSize];

    for (int i = 0; i <= count; i++) {
      a[i] = PHI_i[i];
      b[i] = E_i[i];
      c[i] = DSIG_i[i];
      e[i] = TIME_i[i];
      f[i] = DTIME_i[i];
    }

    if (PHI_i != 0)
      delete [] PHI_i;
    if (E_i != 0)
      delete [] E_i;
    if (DSIG_i != 0)
      delete [] DSIG_i;
    if (TIME_i != 0)
      delete [] TIME_i;
    if (DTIME_i != 0)
      delete [] DTIME_i;

    PHI_i = a;
    E_i = b;
    DSIG_i = c;
    TIME_i = e;
    DTIME_i = f;    
  }
  
}

CreepMaterial::CreepMaterial(int tag, UniaxialMaterial &matl, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast): 
  UniaxialMaterial(tag, MAT_TAG_CreepMaterial), wrappedMaterial(0),
  age(_age), epsshu(_epsshu), epssha(_epssha), tcr(_tcr), epscru(_epscru), epscra(_epscra), epscrd(_epscrd), tcast(_tcast), maxSize(startSize),
  PHI_i(0), E_i(0), DSIG_i(0), TIME_i(0), DTIME_i(0)  
{
  wrappedMaterial = matl.getCopy();
  if (wrappedMaterial == 0) {
    opserr << "CreepMaterial::CreepMaterial - failed to get copy of material" << endln;
    exit(-1);
  }  

  // Set initial tangent
  Ec = wrappedMaterial->getInitialTangent();
  
  ecminP = 0.0;
  deptP = 0.0;
  
  //sigCr = fabs(sigCr);
  eP = Ec; //Added by AMK
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = Ec; //Added by AMK
  Et = Ec;
  count = 0; //Added by AMK
  epsInit = 0.0; //Added by AMK
  sigInit = 0.0; //Added by AMK
  eps_total = 0.0; //Added by AMK
  epsP_total = 0.0; //Added by AMK
  
  eps_m = 0.0; //Added by AMK
  eps_cr = 0.0; //Added by AMK
  eps_sh = 0.0;
  epsP_cr = 0.0; //Added by AMK
  epsP_sh = 0.0; 
  epsP_m = 0.0; //Added by AMK
  
  t_load = -1.0; //Added by AMK
  crack_flag = 0;
  iter = 0;
  
  //Change inputs into the proper sign convention:
  epsshu = -fabs(epsshu);
  epscru = fabs(epscru);

  this->expandArrays();  
}

CreepMaterial::CreepMaterial(void):
  UniaxialMaterial(0, MAT_TAG_CreepMaterial), wrappedMaterial(0), maxSize(startSize),
  PHI_i(0), E_i(0), DSIG_i(0), TIME_i(0), DTIME_i(0)
{
 
}

CreepMaterial::~CreepMaterial(void)
{
  if (wrappedMaterial != 0)
    delete wrappedMaterial;

  if (PHI_i != 0)
    delete [] PHI_i;
  if (E_i != 0)
    delete [] E_i;
  if (DSIG_i != 0)
    delete [] DSIG_i;
  if (TIME_i != 0)
    delete [] TIME_i;
  if (DTIME_i != 0)
    delete [] DTIME_i;
}

UniaxialMaterial*
CreepMaterial::getCopy(void)
{
  CreepMaterial *theCopy = new CreepMaterial(this->getTag(), *wrappedMaterial, age, epsshu, epssha, tcr, epscru, epscra, epscrd, tcast); 

  theCopy->maxSize = maxSize;
  theCopy->count = count;
  
  return theCopy;
}

double
CreepMaterial::getInitialTangent(void)
{
  return Ec; //Added by AMK
}

double
CreepMaterial::getCurrentTime(void)
{
  double currentTime;
  Domain * theDomain = ops_TheActiveDomain;
  
  if (theDomain != 0) {
    currentTime = theDomain->getCurrentTime();
  }
  
  return currentTime;
}	

double
CreepMaterial::setCreepStrain(double time, double stress)
{
  double creep;
  double runSum = 0.0;
  
  DTIME_i[count] = ops_Dt;
  
  for (int i = 1; i<=count; i++) {
    PHI_i[i] = setPhi(time,TIME_i[i]); //Determine PHI
    runSum += PHI_i[i]*DSIG_i[i]/Ec; //CONSTANT STRESS within Time interval
  }
  
  phi_i = PHI_i[count];
  creep = runSum;
  return creep;
}

double 
CreepMaterial::setPhi(double time, double tp)
{	
  // ACI Equation:
  double tmtp = time-tp;
  double f1 = pow((4+0.85*tp)/tp,0.5);
  double f2 = pow(tmtp,epscra)/(epscrd+pow(tmtp,epscra))*epscru;
  double f3 = (1.25*pow((tp-tcast),-0.118))/(1.25*pow(tcr,-0.118));
  double phi = f2*f3;
  return phi;
}

double 
CreepMaterial::setShrink(double time)
{
  double tD = age; //Age at initiation of drying
  double shrink = 0.0;
  if (time-(tD) < 0) {
    shrink = 0.0;
  } else {
    shrink = (time-(tD)) / (epssha + (time - (tD))) * epsshu;
  }
  return shrink;
}

int
CreepMaterial::setTrialStrain(double trialStrain, double strainRate)
{
  double t = getCurrentTime();
  double tol = 1.e-4; // 9/13
  double test = 10.0; // 9/13
  double sigI = 0.0;  // 9/13
  int niter = 500;  // 9/13
  
  // Check casting age:
  if (t-tcast<(2.0-0.0001)) { //Assumed that concrete can only carry load once hardened at 2 days following casting
    eps_cr = 0.0;
    eps_sh = 0.0;
    eps_m = 0.0;
    eps_total = trialStrain;
    sig = 0.0;
  } else { // Concrete has hardened and is ready to accept load
    // Initialize total strain:
    eps_total = trialStrain;
    
    // Calculate shrinkage Strain:
    if (iter < 1) {
      eps_sh = setShrink(t);
    }
    
    // Calculate creep and mechanical strain, assuming stress remains constant in a time step:
    if (ops_Creep == 1) {
      if (fabs(t-TIME_i[count]) <= 0.0001) { //If t = t(i-1), use creep/shrinkage from last calculated time step
	eps_cr = epsP_cr;
	eps_sh = epsP_sh;
	eps_m = eps_total - eps_cr - eps_sh;
	//sig = setStress(eps_m, e);
	wrappedMaterial->setTrialStrain(eps_m, strainRate);
	sig = wrappedMaterial->getStress();
	e = wrappedMaterial->getTangent();
        
      } else { // if the current calculation is a new time step
	if (iter < 1) {
	  eps_cr = setCreepStrain(t,sig); 
	}
	eps_m = eps_total - eps_cr - eps_sh;
	//sig = setStress(eps_m, e);
	wrappedMaterial->setTrialStrain(eps_m, strainRate);
	sig = wrappedMaterial->getStress();
	e = wrappedMaterial->getTangent();	
      }
    } else { //Static Analysis using previously converged time-dependent strains
      eps_cr = epsP_cr;
      eps_sh = epsP_sh;
      eps_m = eps_total-eps_cr-eps_sh;
      //sig = setStress(eps_m, e);
      wrappedMaterial->setTrialStrain(eps_m, strainRate);
      sig = wrappedMaterial->getStress();
      e = wrappedMaterial->getTangent();      
    }
  }
  iter ++;
  return 0;
}

double
CreepMaterial::getStrain(void)
{
  return eps_total; //Added by AMK
  //return eps;
}

double
CreepMaterial::getPHI_i(void)
{
  return phi_i;
}

double 
CreepMaterial::getStress(void)
{
  return sig;
}

double 
CreepMaterial::getTangent(void)
{
  return e;
}

double
CreepMaterial::getCreep(void)
{
  return eps_cr;
}

double
CreepMaterial::getShrink(void)
{
  return eps_sh;
}

double
CreepMaterial::getMech(void)
{
  return eps_m;
}

int 
CreepMaterial::commitState(void)
{
  iter = 0;
  ecminP = ecmin;
  ecmaxP = ecmax;
  deptP = dept;

  // Make sure enough room to write into count+1 -- MHS
  this->expandArrays();
    
  //dsig_i[count]=sig-sigP; // Unused -- MHS
  /* 5/8/2013: commented the following lines so that the DSIG_i[count+1]=sig-sigP;*/
  //if (crack_flag == 1) {// DSIG_i will be different depending on how the fiber is cracked
  //	if (sig < 0 && sigP > 0) { //if current step puts concrete from tension to compression, DSIG_i will be only the comp. stress
  //		DSIG_i[count+1] = sig;
  //	}
  //	if (sig > 0) {// Concrete should not creep when crack is opened
  //		DSIG_i[count+1] = 0.0;
  //	}
  //	if (sig > 0 && sigP < 0) {//if current step goes from compression to tension, DSIG_i will be the stress difference
  //		DSIG_i[count+1] = sig-sigP;
  //	}
  //} else { //concrete is uncracked, DSIG = sig - sigP
  //	DSIG_i[count+1] = sig-sigP;
  //}

  DSIG_i[count+1] = sig-sigP;
  
  //Secant Stiffness for determination of creep strain:
  if (fabs(eps_m/sig)>Ec) {
    E_i[count+1] = Ec;
  } else {
    E_i[count+1] = fabs(sig/eps_m); //ADDED 7/22
  }
  
  if (isnan(E_i[count+1])) {
    E_i[count+1] = Ec;
  }
  
  
  TIME_i[count+1] = getCurrentTime();
  
  eP = e;
  sigP = sig;
  epsP = eps;
	
 //Added by AMK:
  epsP_total = eps_total; //Added by AMK;
  epsP_sh = eps_sh;
  epsP_cr = eps_cr;
  epsP_m = eps_m;
  if (eps_m < 0 && fabs(eps_m)>0.50*fabs(fc/Ec)) {
    double s = fabs(eps_m/fc)*Ec;
    s = 0.5*fabs(fc/Ec);
    //opserr << "Strain Compression Limit Exceeded: " << eps_m << ' ' << -s << endln;
  }
  
  //Cracking flags:
  crackP_flag = crack_flag;
  
  //cracked reloading/unloading stiffness:
  if (crackP_flag==1) {
    if (sig/eps_m<Et) {
      Et = sig/eps_m;
    }
  }
  
  if (count==0) {
    epsInit = epsP_total;
    sigInit = sigP;
  }
  
  if (sigInit<0.0 && t_load<0.0) {
    t_load = getCurrentTime();
    sigInit = sigP;
    epsInit = epsP_m;
  } else if (sigInit>0.0 && sigP<0.0 && t_load<0.0) {
    t_load = getCurrentTime();
    sigInit = sigP;
    epsInit = epsP_m;
  }

  wrappedMaterial->commitState();
  
  //if (ops_Creep==1) {
  //	count++;
  //}
  count++;
  
  return 0;
}

int 
CreepMaterial::revertToLastCommit(void)
{
  eps_total = epsP_total; //Added by AMK;
  eps_sh = epsP_sh;
  eps_cr = epsP_cr;
  eps_m = epsP_m;  
  
  ecmin = ecminP;;
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;

  wrappedMaterial->revertToLastCommit();
  
  return 0;
}

int 
CreepMaterial::revertToStart(void)
{
  ecminP = 0.0;
  deptP = 0.0;
  
  eP = Ec;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = Ec;
  
  if (ops_Creep==0) {
    count = 0;
  } else {
    count = 1;
  }

  wrappedMaterial->revertToStart();
  
  return 0;
}

int 
CreepMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  int dbTag = this->getDbTag();

  static ID classTags(3);

  classTags(0) = wrappedMaterial->getClassTag();

  int matDbTag = wrappedMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      wrappedMaterial->setDbTag(matDbTag);
  }
  classTags(1) = matDbTag;
  classTags(2) = this->getTag();

  res = theChannel.sendID(dbTag, commitTag, classTags);
  if (res < 0) {
    opserr << "CreepMaterial::sendSelf -- could not send ID" << endln;
    return res;
  }

  static Vector data(11);
  data(0) =ft;    
  data(1) =Ec; 
  data(2) =beta;   
  data(3) =age; 
  data(4) =epsshu;   
  data(5) =epssha;    
  data(6) =tcr;   
  data(7) =epscru;
  data(8) =epscra; 
  data(9) =epscrd;     
  data(10) = this->getTag();

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CreepMaterial::sendSelf - failed to send Vector" << endln;
    return res;
  }
  
  res = wrappedMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "CreepMaterial::sendSelf -- could not send UniaxialMaterial" << endln;
    return res;
  }
	
  return res;
}

int 
CreepMaterial::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  int dbTag = this->getDbTag();

  res = theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "CreepMaterial::recvSelf() - failed to receive data\n";
    return res;
  }

  static Vector vdata(11);

  if (theChannel.recvVector(this->getDbTag(), commitTag, vdata) < 0) {
    opserr << "CreepMaterial::recvSelf() - failed to recvSelf\n";
    return -1;
  }
  
  ft = vdata(0);   
  Ec = vdata(1);
  beta = vdata(2);   
  age = vdata(3); 
  epsshu = vdata(4);   
  epssha = vdata(5);    
  tcr = vdata(6);   
  epscru = vdata(7);
  epscra = vdata(8); 
  epscrd = vdata(9);   
  this->setTag(vdata(10));

  e = eP;
  sig = sigP;
  eps = epsP;
  
  
  this->setTag(int(data(2)));
  int matClassTag = data(0);
  
  // Check if material is null
  if (wrappedMaterial == 0) {
    wrappedMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (wrappedMaterial == 0) {
      opserr << "CreepMaterial::recvSelf -- could not get a UniaxialMaterial" << endln;
      return -1;
    }
  }

  dbTag = data(1);
  if (wrappedMaterial->getClassTag() != matClassTag) {
    delete wrappedMaterial;
    wrappedMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (wrappedMaterial == 0) {
      opserr << "CreepMaterial::recvSelf -- could not get a UniaxialMaterial" << endln;
      return -1;
    }
  }
  
  wrappedMaterial->setDbTag(dbTag);
  res = wrappedMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "CreepMaterial::recvSelf -- count not receive Uniaxialmaterial" << endln;
    return res;
  }
  
  return res;
}

void 
CreepMaterial::Print(OPS_Stream &s, int flag)
{
  s << "CreepMaterial:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}


int
CreepMaterial::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}

/* Methods added by AMK: */

Response* 
CreepMaterial::setResponse(const char **argv, int argc,
							  OPS_Stream &theOutput)
{	
  Response *theResponse = 0;
  
  theOutput.tag("UniaxialMaterialOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());
  
  // stress
  if (strcmp(argv[0],"stress") == 0) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse =  new MaterialResponse(this, 1, this->getStress());
  }  
  // tangent
  else if (strcmp(argv[0],"tangent") == 0) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 2, this->getTangent());
  }
  
  // strain
  else if (strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 3, this->getStrain());
  }
  
  // strain
  else if ((strcmp(argv[0],"stressStrain") == 0) || 
	   (strcmp(argv[0],"stressANDstrain") == 0) ||
	   (strcmp(argv[0],"stressAndStrain") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 4, Vector(2));
  }
  
  else if (strcmp(argv[0],"CreepStressStrainTangent")==0) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theOutput.tag("ResponseType", "C11");
    theOutput.tag("ResponseType", "CreepStrain");
    theOutput.tag("ResponseType", "MechStrain");
    theOutput.tag("ResponseType", "ShrinkStrain");
    theOutput.tag("ResponseType", "t_load");
    theResponse = new MaterialResponse(this, 6, Vector(6));
  }
  
  else if ((strcmp(argv[0],"stressStrainTangent") == 0) || 
	   (strcmp(argv[0],"stressANDstrainANDtangent") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 5, Vector(3));
  }
  
  // stress sensitivity for local sensitivity recorder purpose.  Quan 2009
  // limit:  no more than 10000 random variables/sensitivity parameters
  else if (strstr(argv[0],"stressSensitivity") != 0) {
    char *token = strtok((char *) argv[0], " ");
    if (token != NULL) token = strtok(NULL, " ");
    int gradient = atoi(token);
    theOutput.tag("ResponseType", "sigsens11");
    theResponse =  new MaterialResponse(this, gradient+10000, this->getStress());
  }
  // strain sensivitiy
  else if (strstr(argv[0],"strainSensitivity") != 0) {
    char *token = strtok((char *) argv[0], " ");
    if (token != NULL) token = strtok(NULL, " ");
    int gradient = atoi(token);
    theOutput.tag("ResponseType", "epssens11");
    theResponse =  new MaterialResponse(this, gradient+20000, this->getStrain());
  }
  
  
  theOutput.endTag();
  return theResponse;
}

int 
CreepMaterial::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);
  static Vector stressStrainTangent(3);
  static Vector CreepStressStrainTangent(6); //Added by AMK
  // each subclass must implement its own stuff   
  
  // added for sensitivity recorder. Quan 2009
  if ((responseID>10000)&&(responseID<20000)){
    matInfo.setDouble(this->getStressSensitivity(responseID-10000,false));
    return 0;
  }
  else if (responseID>20000){
    matInfo.setDouble(this->getStrainSensitivity(responseID-20000));
    return 0;
  }
  
  switch (responseID) {
  case 1:
    matInfo.setDouble(this->getStress());
    return 0;
    
  case 2:
    matInfo.setDouble(this->getTangent());
    return 0;      
    
  case 3:
    matInfo.setDouble(this->getStrain());
    return 0;      
    
  case 4:
    stressStrain(0) = this->getStress();
    stressStrain(1) = this->getStrain();
    matInfo.setVector(stressStrain);
    return 0;
    
  case 5:
    stressStrainTangent(0) = this->getStress();
    stressStrainTangent(1) = this->getStrain();
    stressStrainTangent(2) = this->getTangent();
    matInfo.setVector(stressStrainTangent);
    return 0;
    
  case 6:
    CreepStressStrainTangent(0) = this->getStress();
    CreepStressStrainTangent(1) = this->getStrain();
    CreepStressStrainTangent(2) = this->getTangent();
    CreepStressStrainTangent(3) = this->getCreep();
    CreepStressStrainTangent(4) = this->getMech();
    CreepStressStrainTangent(5) = this->getShrink();
    matInfo.setVector(CreepStressStrainTangent);
    return 0;
    
  default:      
    return -1;
  }
}
